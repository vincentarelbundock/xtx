/*
 *  xtx: column-structured algebra on a plain design.
 *
 *  Row i of term t of the dense matrix is the Kronecker product of its
 *  components' rows.  Split the components into factors and numerics:
 *  the factors give one combined level F_t(i) (0 .. L_t - 1, the first
 *  factor varying fastest) and the numerics a vector N_t(i) of ns_t
 *  products (the first numeric varying fastest).  Then
 *
 *      X_t[i, colmap_t[u, k]] = N_t(i)[u] * C_t[F_t(i), k]
 *
 *  where C_t is the Kronecker product of the factors' contrast matrices.
 *  Products with the design reduce to sums over rows of N's within
 *  combined levels (moments), followed by small dense products with the
 *  C_t, done in R.  The moments are computed here, in blocks of rows so
 *  that each block's levels and products stay in cache.
 *
 *  Only for designs whose inputs are all finite, with no missing factor
 *  codes, and whose products cannot overflow (checked by xtx_struct_ok):
 *  then the results agree with the dense computation up to rounding.
 *
 *  Copyright (C) 2026 Vincent Arel-Bundock.  GPL (>= 2).
 */

#define USE_FC_LEN_T
#include <string.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include "xtx.h"
#include "structured.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef FCONE
# define FCONE
#endif

#define SAFE (DBL_MAX / 16)

/* Terms of a plain design (the intercept is a term with no component).
   Leaves one object PROTECTed (the coerced variables). */
void build_sctx(sctx *S, SEXP setup, SEXP vars)
{
    SEXP factors = VECTOR_ELT(setup, XTX_SETUP_FACTORS),
	contr1 = VECTOR_ELT(setup, XTX_SETUP_CONTR1),
	contr2 = VECTOR_ELT(setup, XTX_SETUP_CONTR2);
    const int *nlevs = INTEGER(VECTOR_ELT(setup, XTX_SETUP_NLEVS)),
	*columns = INTEGER(VECTOR_ELT(setup, XTX_SETUP_COLUMNS)),
	*count = INTEGER(VECTOR_ELT(setup, XTX_SETUP_COUNT));
    int rhs_response = asInteger(VECTOR_ELT(setup, XTX_SETUP_RHS_RESPONSE)),
	intercept = asInteger(VECTOR_ELT(setup, XTX_SETUP_INTERCEPT)),
	nVar = asInteger(VECTOR_ELT(setup, XTX_SETUP_NVAR)),
	nterms = asInteger(VECTOR_ELT(setup, XTX_SETUP_NTERMS));
    const int *fac = nterms > 0 ? INTEGER(factors) : NULL;
    SEXP variable = PROTECT(xtx_vars(setup, vars));
    S->n = nrows(VECTOR_ELT(vars, 0));
    S->p = LENGTH(VECTOR_ELT(setup, XTX_SETUP_XNAMES));
    S->t = (sterm *) R_alloc(nterms + 1, sizeof(sterm));
    S->nt = 0;
    int off = 0;
    if (intercept) {
	sterm *T = &S->t[S->nt++];
	memset(T, 0, sizeof(sterm));
	T->L = 1; T->ns = 1; T->kc = 1; T->count = 1; T->off = 0;
	off = 1;
    }
    for (int k = 0; k < nterms; k++) {
	if (k == rhs_response || count[k] <= 0) continue;
	sterm *T = &S->t[S->nt++];
	memset(T, 0, sizeof(sterm));
	T->off = off;
	T->count = count[k];
	off += count[k];
	int stride = 1, L = 1, ns = 1, kc = 1;
	/* the same counts in double, to detect int overflow */
	double Ld = 1.0, kcd = 1.0, nsd = 1.0;
	for (int i = 0; i < nVar; i++) {
	    int fik = fac[i + k * nVar];
	    if (!fik || columns[i] == 0) continue;
	    if (T->nf + T->nn >= SMAXC) error("internal error: term with too many variables");
	    SEXP var_i = VECTOR_ELT(variable, i);
	    if (nlevs[i] > 0) {
		SEXP cm = VECTOR_ELT(fik == 1 ? contr1 : contr2, i);
		int c = T->nf++;
		T->code[c] = INTEGER_RO(var_i);
		T->fvar[c] = i;
		T->fcm[c] = cm;
		T->adj[c] = isLogical(var_i) ? 1 : 0;
		T->cm[c] = REAL_RO(cm);
		T->nrc[c] = nrows(cm);
		T->fnc[c] = ncols(cm);
		T->lstride[c] = L;
		T->fstride[c] = stride;
		Ld *= nrows(cm);
		kcd *= ncols(cm);
		L = Ld > INT_MAX / 2 ? 1 : (int) Ld;
		kc = kcd > INT_MAX / 2 ? 1 : (int) kcd;
		stride *= ncols(cm);
	    } else {
		int c = T->nn++;
		T->num[c] = REAL_RO(var_i);
		T->nnc[c] = ncols(var_i);
		T->nstride[c] = stride;
		nsd *= ncols(var_i);
		ns = nsd > INT_MAX / 2 ? 1 : (int) nsd;
		stride *= ncols(var_i);
	    }
	}
	T->L = L; T->ns = ns; T->kc = kc;
	if (Ld > INT_MAX / 2 || kcd > INT_MAX / 2 || nsd > INT_MAX / 2) {
	    /* too large for the structured kernels (xtx_struct_ok declines) */
	    T->L = T->kc = T->ns = -1;
	}
    }
}

/* Offset of numeric column u (u < T->ns) of term T among the term's
   columns: u's digits in the radices nnc (first fastest) times the
   numeric components' strides. */
static inline int num_offset(const sterm *T, int u)
{
    int part = 0;
    for (int q = 0; q < T->nn; q++) {
	part += (u % T->nnc[q]) * T->nstride[q];
	u /= T->nnc[q];
    }
    return part;
}

/* Combined levels and numeric products of rows i0 .. i0 + nb - 1:
   N[u * BLOCK + r] (column-major, leading dimension BLOCK), so that every
   loop over rows reads contiguous memory and the dense terms' columns can
   be passed to BLAS in place. */
void term_block(const sterm *T, R_xlen_t n, R_xlen_t i0, int nb, int *F, double *N)
{
    for (int r = 0; r < nb; r++) F[r] = 0;
    for (int c = 0; c < T->nf; c++) {
	const int *code = T->code[c] + i0;
	int a = T->adj[c] - 1, ls = T->lstride[c];
	for (int r = 0; r < nb; r++) F[r] += (code[r] + a) * ls;
    }
    if (T->nn == 0) {
	for (int r = 0; r < nb; r++) N[r] = 1.0;
	return;
    }
    /* built component by component, first fastest; columns of later
       components are written from the highest down, so that each product
       reads columns not yet overwritten */
    int len = 1;
    for (int c = 0; c < T->nn; c++) {
	int nc = T->nnc[c];
	const double *v = T->num[c] + i0;
	if (c == 0) {
	    for (int j = 0; j < nc; j++)
		memcpy(N + (size_t) j * BLOCK, v + (R_xlen_t) j * n, sizeof(double) * nb);
	} else {
	    for (int j = nc - 1; j >= 0; j--) {
		const double *x = v + (R_xlen_t) j * n;
		for (int u = len - 1; u >= 0; u--) {
		    double *dst = N + (size_t) (u + len * j) * BLOCK;
		    const double *src = N + (size_t) u * BLOCK;
		    for (int r = 0; r < nb; r++) dst[r] = src[r] * x[r];
		}
	    }
	}
	len *= nc;
    }
}

typedef struct {
    int *F;
    double *N;
} tblock;

/* The dense part of a design: the terms with a single level (the
   intercept and products of numeric variables). Their block columns are
   stored side by side in one BLOCK x q buffer (leading dimension BLOCK),
   so that their products with each other, with Y and with B use BLAS
   (cache blocking and vector instructions) in place. */
typedef struct {
    int q;          /* columns of the dense part */
    int *off;       /* per term: first column in the buffer, or -1 */
    double *buf;
} dpart;

/* Number of dense columns (terms with one level), as dense_part() counts
   them, without allocating its buffer. */
static int dense_count(const sctx *S)
{
    int q = 0;
    for (int t = 0; t < S->nt; t++)
	if (S->t[t].L == 1) q += S->t[t].ns;
    return q;
}

static dpart dense_part(const sctx *S)
{
    dpart P;
    P.off = (int *) R_alloc(S->nt > 0 ? S->nt : 1, sizeof(int));
    P.q = 0;
    for (int t = 0; t < S->nt; t++) {
	if (S->t[t].L == 1) { P.off[t] = P.q; P.q += S->t[t].ns; }
	else P.off[t] = -1;
    }
    P.buf = P.q > 0 ? (double *) R_alloc((size_t) BLOCK * P.q, sizeof(double)) : NULL;
    return P;
}

static tblock *alloc_blocks(const sctx *S, const dpart *P)
{
    tblock *B = (tblock *) R_alloc(S->nt > 0 ? S->nt : 1, sizeof(tblock));
    for (int t = 0; t < S->nt; t++) {
	B[t].F = (int *) R_alloc(BLOCK, sizeof(int));
	B[t].N = P->off[t] >= 0 ? P->buf + (size_t) BLOCK * P->off[t]
	    : (double *) R_alloc((size_t) BLOCK * S->t[t].ns, sizeof(double));
    }
    return B;
}

/* Can the structured kernels handle this design?  All factor codes
   present, all inputs finite, moments of a size that fits (`cap` doubles
   in all), and magnitudes bounded so that regrouping products cannot
   overflow or underflow: in every term, the non-zero products of its
   numeric inputs and of its contrast entries lie within [1e-50, 1e50]
   (bounds from each component's smallest and largest non-zero absolute
   value).  An entry of X'X is then a sum of products of four such factors,
   within [1e-200, 1e200], and the callers bound row scales and the other
   operand of a product so that every intermediate stays within 1e+-300. */
#define MAG_LO 1e-50
#define MAG_HI 1e50
static void range_of(const double *v, R_xlen_t len, double *lo, double *hi, int *finite)
{
    double mn = R_PosInf, mx = 0.0;
    for (R_xlen_t e = 0; e < len; e++) {
	double x = v[e];
	if (!R_FINITE(x)) { *finite = 0; return; }
	double a = fabs(x);
	if (a > mx) mx = a;
	if (a > 0.0 && a < mn) mn = a;
    }
    *lo = R_FINITE(mn) ? mn : 1.0;
    *hi = mx > 0.0 ? mx : 1.0;
}

SEXP xtx_struct_ok(SEXP setup, SEXP vars, SEXP cap_)
{
    sctx S;
    build_sctx(&S, setup, vars);
    double cap = asReal(cap_), total = 0.0;
    int ok = 1;
    for (int t = 0; t < S.nt; t++)
	if (S.t[t].L < 0) ok = 0;
    for (int t = 0; t < S.nt && ok; t++) {
	const sterm *T = &S.t[t];
	double clo = 1.0, chi = 1.0, nlo = 1.0, nhi = 1.0;
	for (int c = 0; c < T->nf && ok; c++) {
	    const int *code = T->code[c];
	    int lo = 1 - T->adj[c], hi = T->nrc[c] - T->adj[c];
	    for (R_xlen_t i = 0; i < S.n; i++)
		if (code[i] == NA_INTEGER || code[i] < lo || code[i] > hi) { ok = 0; break; }
	    double a, b;
	    range_of(T->cm[c], (R_xlen_t) T->nrc[c] * T->fnc[c], &a, &b, &ok);
	    clo *= a; chi *= b;
	}
	for (int c = 0; c < T->nn && ok; c++) {
	    double a, b;
	    range_of(T->num[c], S.n * T->nnc[c], &a, &b, &ok);
	    nlo *= a; nhi *= b;
	}
	if (!(clo >= MAG_LO && chi <= MAG_HI && nlo >= MAG_LO && nhi <= MAG_HI)) ok = 0;
	for (int s = 0; s <= t; s++) {
	    const sterm *U = &S.t[s];
	    total += (s == t) ? (double) T->L * T->ns * T->ns
		: (double) T->L * U->L * T->ns * U->ns;
	}
	total += (double) T->L * T->kc;
    }
    if (total > cap) ok = 0;
    UNPROTECT(1);
    return ScalarLogical(ok);
}

/* Term descriptions: list per term of list(off, colmap, fvar, fcm), where
   colmap is the ns x kc matrix of 0-based column indices within the term,
   fvar the 0-based variable index of each factor (first fastest in the
   combined level) and fcm the list of their contrast matrices: the term's
   combined contrast matrix C_t is the Kronecker product of the fcm, the
   first factor varying fastest. */
SEXP xtx_struct_info(SEXP setup, SEXP vars)
{
    sctx S;
    build_sctx(&S, setup, vars);
    SEXP out = PROTECT(allocVector(VECSXP, S.nt));
    for (int t = 0; t < S.nt; t++) {
	const sterm *T = &S.t[t];
	SEXP el = SET_VECTOR_ELT(out, t, allocVector(VECSXP, 4));
	SET_VECTOR_ELT(el, 0, ScalarInteger(T->off));
	SEXP M = SET_VECTOR_ELT(el, 1, allocMatrix(INTSXP, T->ns, T->kc));
	SEXP fv = SET_VECTOR_ELT(el, 2, allocVector(INTSXP, T->nf));
	SEXP fc = SET_VECTOR_ELT(el, 3, allocVector(VECSXP, T->nf));
	for (int f = 0; f < T->nf; f++) {
	    INTEGER(fv)[f] = T->fvar[f];
	    SET_VECTOR_ELT(fc, f, T->fcm[f]);
	}
	int *m = INTEGER(M);
	for (int k = 0; k < T->kc; k++) {
	    /* factor columns of combined column k (first factor fastest) */
	    int rem = k, fpart = 0;
	    for (int f = 0; f < T->nf; f++) {
		fpart += (rem % T->fnc[f]) * T->fstride[f];
		rem /= T->fnc[f];
	    }
	    for (int u = 0; u < T->ns; u++)
		m[u + k * T->ns] = num_offset(T, u) + fpart;
	}
    }
    UNPROTECT(2);
    return out;
}

/* The factors of term B that term A does not have: their number, their
   positions in B and the strides of their combined level (first fastest);
   returns the number of combined levels (Lx). */
static int extra_factors(const sterm *A, const sterm *B, int *nx, int *pos, int *stride)
{
    int L = 1;
    *nx = 0;
    for (int g = 0; g < B->nf; g++) {
	int shared = 0;
	for (int f = 0; f < A->nf; f++)
	    if (A->fvar[f] == B->fvar[g]) { shared = 1; break; }
	if (shared) continue;
	pos[*nx] = g;
	stride[*nx] = L;
	(*nx)++;
	L *= B->nrc[g];
    }
    return L;
}

/* One mode of a product with a Kronecker product of contrast matrices.
   X holds an l x R matrix (the first mode first); A is l x ka and B, if
   not NULL, l x kb (B NULL: kb = 1, B = 1).  Writes to o the R x (ka kb)
   matrix with entries sum_a X[a, r] A[a, i] B[a, j] (column i + ka j), so
   that the next mode comes first.  Only the non-zero entries of A and B
   are visited: a treatment contrast costs R per level, not R per entry.
   x (R), ia (ka) and ib (kb) are work space. */
static void kron_mode(const double *X, R_xlen_t R, int l, const double *A, int ka,
		      const double *B, int kb, double *o, double *x, int *ia, int *ib)
{
    if (!B) kb = 1;
    memset(o, 0, sizeof(double) * (size_t) R * ka * kb);
    for (int a = 0; a < l; a++) {
	int na = 0, nb = 0, any = 0;
	for (int i = 0; i < ka; i++) if (A[a + (R_xlen_t) l * i] != 0.0) ia[na++] = i;
	if (B) { for (int j = 0; j < kb; j++) if (B[a + (R_xlen_t) l * j] != 0.0) ib[nb++] = j; }
	else { ib[0] = 0; nb = 1; }
	if (na == 0 || nb == 0) continue;
	for (R_xlen_t r = 0; r < R; r++) {
	    x[r] = X[a + (R_xlen_t) l * r];
	    if (x[r] != 0.0) any = 1;
	}
	if (!any) continue;
	for (int jj = 0; jj < nb; jj++) {
	    int j = ib[jj];
	    double bj = B ? B[a + (R_xlen_t) l * j] : 1.0;
	    for (int ii = 0; ii < na; ii++) {
		int i = ia[ii];
		double c = A[a + (R_xlen_t) l * i] * bj;
		double *oc = o + R * (i + (R_xlen_t) ka * j);
		for (R_xlen_t r = 0; r < R; r++) oc[r] += x[r] * c;
	    }
	}
    }
}

/* kron_mode() from R: X an array with l = nrow(A) rows first, A, B (or
   NULL) contrast matrices; an R x (ka kb) matrix. */
SEXP xtx_kron_mode(SEXP X_, SEXP A_, SEXP B_)
{
    int l = nrows(A_), ka = ncols(A_), kb = isNull(B_) ? 1 : ncols(B_);
    if (!isNull(B_) && nrows(B_) != l) error("internal error: contrasts of different levels");
    R_xlen_t len = XLENGTH(X_);
    if (l == 0 || len % l) error("internal error: moments do not match the contrasts");
    R_xlen_t R = len / l;
    if (R > INT_MAX) error("internal error: too many moments");
    SEXP out = PROTECT(allocMatrix(REALSXP, (int) R, ka * kb));
    double *x = (double *) R_alloc(R > 0 ? R : 1, sizeof(double));
    int *ia = (int *) R_alloc(ka > 0 ? ka : 1, sizeof(int)), *ib = (int *) R_alloc(kb > 0 ? kb : 1, sizeof(int));
    kron_mode(REAL(X_), R, l, REAL(A_), ka, isNull(B_) ? NULL : REAL(B_), kb, REAL(out), x, ia, ib);
    UNPROTECT(1);
    return out;
}

/* Does every row of an l x k contrast matrix have at most one non-zero
   entry (treatment, SAS or indicator contrasts)?  Then its column for
   each level and that entry. */
static int selection(const double *A, int l, int k, int *col, double *val)
{
    for (int a = 0; a < l; a++) {
	col[a] = -1; val[a] = 0.0;
	for (int i = 0; i < k; i++) {
	    double v = A[a + (R_xlen_t) l * i];
	    if (v == 0.0) continue;
	    if (col[a] >= 0) return 0;
	    col[a] = i; val[a] = v;
	}
    }
    return 1;
}

/* The block of X'WX for terms s <= t from their moments M (layout of
   xtx_struct_gram()), added to G (p x p, zero on entry).  The factors of
   the union of the two terms are taken in turn (the moments' dims: first
   those of s, then those of t that s lacks), each with the contrast of
   s, of t, or both for a shared factor.  A factor whose contrasts are
   selections (selection()) needs no product: its levels stay a dim of
   the moments, scaled by the selected entries, and map to columns when
   the block is written.  Otherwise kron_mode() contracts it, visiting only
   non-zero entries; the last such contraction adds straight into G.  For
   s < t the block is added both above and below the diagonal, by the
   same operations in the same order (so G is exactly symmetric). */
static void gram_block(const sterm *A, const sterm *B, int same, const double *M,
		       int nx, const int *xpos, double *G, int p)
{
    const void *vmax = vmaxget();
    int ns = A->ns, nt = B->ns;
    /* modes: their contrasts, and where their columns go: the first
       contrast's to side ma (0: the rows of the block, s; 1: its columns,
       t) with stride sa, the second's (shared factors) to t with stride
       sb */
    int nm = 0;
    const double *mA[2 * SMAXC], *mB[2 * SMAXC];
    int ml[2 * SMAXC], mka[2 * SMAXC], mkb[2 * SMAXC], ma[2 * SMAXC], sa[2 * SMAXC], sb[2 * SMAXC];
    for (int f = 0; f < A->nf; f++) {
	int g = -1;
	if (same) g = f;
	else for (int h = 0; h < B->nf; h++) if (B->fvar[h] == A->fvar[f]) { g = h; break; }
	mA[nm] = A->cm[f]; ml[nm] = A->nrc[f]; mka[nm] = A->fnc[f];
	ma[nm] = 0; sa[nm] = A->fstride[f];
	mB[nm] = g >= 0 ? B->cm[g] : NULL; mkb[nm] = g >= 0 ? B->fnc[g] : 1;
	sb[nm] = g >= 0 ? B->fstride[g] : 0;
	nm++;
    }
    for (int c = 0; c < nx; c++) {
	int g = xpos[c];
	mA[nm] = B->cm[g]; ml[nm] = B->nrc[g]; mka[nm] = B->fnc[g];
	ma[nm] = 1; sa[nm] = B->fstride[g];
	mB[nm] = NULL; mkb[nm] = 1; sb[nm] = 0;
	nm++;
    }
    /* selections: per level, the scale and the row and column offsets */
    int sel[2 * SMAXC];
    double *scl[2 * SMAXC];
    R_xlen_t *srow[2 * SMAXC], *scol[2 * SMAXC];
    int last = -1; /* the last contracted mode */
    int lmax = 1, lsum = 0;
    for (int m = 0; m < nm; m++) { lsum += ml[m]; if (ml[m] > lmax) lmax = ml[m]; }
    int *ca = (int *) R_alloc(2 * (size_t) lmax, sizeof(int)), *cb = ca + lmax;
    double *va = (double *) R_alloc(2 * (size_t) lmax + lsum, sizeof(double)), *vb = va + lmax,
	*sclw = vb + lmax;
    R_xlen_t *sw = (R_xlen_t *) R_alloc(2 * (size_t) (lsum > 0 ? lsum : 1), sizeof(R_xlen_t));
    for (int m = 0, used = 0; m < nm; m++) {
	int l = ml[m];
	sel[m] = selection(mA[m], l, mka[m], ca, va) &&
	    (!mB[m] || selection(mB[m], l, mkb[m], cb, vb));
	if (!sel[m]) { last = m; continue; }
	scl[m] = sclw + used;
	srow[m] = sw + 2 * used;
	scol[m] = sw + 2 * used + l;
	used += l;
	for (int a = 0; a < l; a++) {
	    int ok = ca[a] >= 0 && (!mB[m] || cb[a] >= 0);
	    scl[m][a] = ok ? va[a] * (mB[m] ? vb[a] : 1.0) : 0.0;
	    R_xlen_t oa = ok ? (R_xlen_t) sa[m] * ca[a] : 0;
	    srow[m][a] = ma[m] ? 0 : oa;
	    scol[m][a] = (ma[m] ? oa : 0) + (ok && mB[m] ? (R_xlen_t) sb[m] * cb[a] : 0);
	}
    }
    /* only selections: no product, each moment goes to G in place */
    if (last < 0 && nm > 0) {
	int idx[2 * SMAXC + 2], dsz[2 * SMAXC + 2], nd = 0;
	for (int m = 1; m < nm; m++) dsz[nd++] = ml[m];
	dsz[nd++] = ns;
	dsz[nd++] = nt;
	for (int d = 0; d < nd; d++) idx[d] = 0;
	int *nps = (int *) R_alloc(ns + nt, sizeof(int)), *npt = nps + ns;
	for (int side = 0; side < 2; side++) {
	    const sterm *T = side ? B : A;
	    int *np = side ? npt : nps;
	    for (int u = 0; u < T->ns; u++) np[u] = num_offset(T, u);
	}
	int l0 = ml[0];
	const double *sc0 = scl[0];
	const R_xlen_t *sr0 = srow[0], *sc0c = scol[0];
	for (const double *Me = M; ; Me += l0) {
	    double so = 1.0;
	    R_xlen_t rb = A->off + nps[idx[nd - 2]], cbs = B->off + npt[idx[nd - 1]];
	    for (int m = 1; m < nm; m++) {
		int a = idx[m - 1];
		so *= scl[m][a];
		rb += srow[m][a];
		cbs += scol[m][a];
	    }
	    if (so != 0.0) {
		for (int a = 0; a < l0; a++) {
		    double c = sc0[a] * so, v = Me[a];
		    if (c == 0.0 || v == 0.0) continue;
		    v *= c;
		    R_xlen_t row = rb + sr0[a], col = cbs + sc0c[a];
		    G[row + p * col] += v;
		    if (!same) G[col + p * row] += v;
		}
	    }
	    int d = 0;
	    for (; d < nd; d++) {
		if (++idx[d] < dsz[d]) break;
		idx[d] = 0;
	    }
	    if (d == nd) break;
	}
	vmaxset(vmax);
	return;
    }
    /* take the modes in turn; the last contracted one is left for G */
    R_xlen_t len = (R_xlen_t) ns * nt;
    for (int m = 0; m < nm; m++) len *= ml[m];
    R_xlen_t maxlen = len;
    {
	R_xlen_t cur = len;
	for (int m = 0; m < nm; m++) {
	    if (!sel[m]) cur = cur / ml[m] * mka[m] * mkb[m];
	    if (cur > maxlen) maxlen = cur;
	}
    }
    int kmax = 1;
    for (int m = 0; m < nm; m++) { if (mka[m] > kmax) kmax = mka[m]; if (mkb[m] > kmax) kmax = mkb[m]; }
    int *ia = (int *) R_alloc(kmax, sizeof(int)), *ib = (int *) R_alloc(kmax, sizeof(int));
    double *x = (double *) R_alloc((size_t) maxlen, sizeof(double));
    double *buf[2];
    buf[0] = (double *) R_alloc((size_t) maxlen, sizeof(double));
    buf[1] = (double *) R_alloc((size_t) maxlen, sizeof(double));
    const double *X = M;
    int nb = 0;
    int nproc = last >= 0 ? last : nm; /* modes before the fused one */
    for (int m = 0; m < nm; m++) {
	if (m == last) break; /* written to G below */
	int l = ml[m];
	R_xlen_t R = len / l;
	double *o = buf[nb++ % 2];
	if (sel[m]) {
	    /* move the dim to the end, scaled */
	    const double *sc = scl[m];
	    for (int a = 0; a < l; a++) {
		double c = sc[a];
		double *oa = o + R * a;
		if (c == 0.0) { memset(oa, 0, sizeof(double) * (size_t) R); continue; }
		for (R_xlen_t r = 0; r < R; r++) oa[r] = X[a + (R_xlen_t) l * r] * c;
	    }
	} else {
	    kron_mode(X, R, l, mA[m], mka[m], mB[m], mkb[m], o, x, ia, ib);
	    len = R * mka[m] * mkb[m];
	}
	X = o;
    }
    /* the modes after the fused one are all selections: move each to the
       end too, so that the fused mode comes first */
    for (int m = last + 1; last >= 0 && m < nm; m++) {
	/* X is [l_last, l_{last+1}, ..., rest]: swap the second dim to the
	   end */
	int l0 = ml[last], l = ml[m];
	R_xlen_t R = len / ((R_xlen_t) l0 * l);
	double *o = buf[nb++ % 2];
	const double *sc = scl[m];
	for (int a = 0; a < l; a++) {
	    double c = sc[a];
	    double *oa = o + (R_xlen_t) l0 * R * a;
	    if (c == 0.0) { memset(oa, 0, sizeof(double) * (size_t) l0 * R); continue; }
	    for (R_xlen_t r = 0; r < R; r++) {
		const double *xr = X + (R_xlen_t) l0 * (a + (R_xlen_t) l * r);
		double *orr = oa + (R_xlen_t) l0 * r;
		for (int b = 0; b < l0; b++) orr[b] = xr[b] * c;
	    }
	}
	X = o;
    }
    /* the order of the dims after the first (l, 1 without a fused mode):
       u, v, then the modes before the fused one, then those after it */
    int order[2 * SMAXC], no = 0;
    for (int m = 0; m < nproc; m++) order[no++] = m;
    if (last >= 0) for (int m = last + 1; m < nm; m++) order[no++] = m;
    int l = last >= 0 ? ml[last] : 1;
    R_xlen_t R = len / l;
    /* row and column offsets of each element of the rest */
    int *nps = (int *) R_alloc(ns, sizeof(int)), *npt = (int *) R_alloc(nt, sizeof(int));
    for (int side = 0; side < 2; side++) {
	const sterm *T = side ? B : A;
	int *np = side ? npt : nps;
	for (int u = 0; u < T->ns; u++) np[u] = num_offset(T, u);
    }
    R_xlen_t *ro = (R_xlen_t *) R_alloc((size_t) R, sizeof(R_xlen_t)),
	*co = (R_xlen_t *) R_alloc((size_t) R, sizeof(R_xlen_t));
    {
	/* dims: per mode, a selection's levels or a contraction's first
	   and (shared) second column; offsets as tables */
	int dsize[4 * SMAXC + 2], nd = 0;
	const R_xlen_t *drow[4 * SMAXC + 2], *dcol[4 * SMAXC + 2];
	for (int q = 0; q < no; q++) {
	    int m = order[q];
	    if (sel[m]) {
		dsize[nd] = ml[m]; drow[nd] = srow[m]; dcol[nd] = scol[m]; nd++;
		continue;
	    }
	    for (int h = 0; h < (mB[m] ? 2 : 1); h++) {
		int k = h ? mkb[m] : mka[m], side = h ? 1 : ma[m];
		R_xlen_t st = h ? sb[m] : sa[m];
		R_xlen_t *zr = (R_xlen_t *) R_alloc(k, sizeof(R_xlen_t)),
		    *zc = (R_xlen_t *) R_alloc(k, sizeof(R_xlen_t));
		for (int i = 0; i < k; i++) {
		    zr[i] = side ? 0 : st * i;
		    zc[i] = side ? st * i : 0;
		}
		dsize[nd] = k; drow[nd] = zr; dcol[nd] = zc; nd++;
	    }
	}
	int idx[4 * SMAXC + 2];
	for (int d = 0; d < nd; d++) idx[d] = 0;
	R_xlen_t fr = 0, fc = 0, r = 0;
	for (int d = 0; d < nd; d++) { fr += drow[d][0]; fc += dcol[d][0]; }
	for (;;) {
	    for (int v = 0; v < nt; v++)
		for (int u = 0; u < ns; u++, r++) {
		    ro[r] = A->off + nps[u] + fr;
		    co[r] = B->off + npt[v] + fc;
		}
	    int d = 0;
	    for (; d < nd; d++) {
		int old = idx[d], nw = old + 1 < dsize[d] ? old + 1 : 0;
		idx[d] = nw;
		fr += drow[d][nw] - drow[d][old];
		fc += dcol[d][nw] - dcol[d][old];
		if (nw) break;
	    }
	    if (d == nd) break;
	}
    }
    if (last < 0) {
	for (R_xlen_t r = 0; r < R; r++) {
	    if (X[r] == 0.0) continue;
	    G[ro[r] + p * co[r]] += X[r];
	    if (!same) G[co[r] + p * ro[r]] += X[r];
	}
	vmaxset(vmax);
	return;
    }
    /* the fused mode, into G */
    int m = last, ka = mka[m], kb = mB[m] ? mkb[m] : 1;
    const double *Am = mA[m], *Bm = mB[m];
    for (int a = 0; a < l; a++) {
	int na = 0, nbb = 0, any = 0;
	for (int i = 0; i < ka; i++) if (Am[a + (R_xlen_t) l * i] != 0.0) ia[na++] = i;
	if (Bm) { for (int j = 0; j < kb; j++) if (Bm[a + (R_xlen_t) l * j] != 0.0) ib[nbb++] = j; }
	else { ib[0] = 0; nbb = 1; }
	if (na == 0 || nbb == 0) continue;
	for (R_xlen_t r = 0; r < R; r++) {
	    x[r] = X[a + (R_xlen_t) l * r];
	    if (x[r] != 0.0) any = 1;
	}
	if (!any) continue;
	for (int jj = 0; jj < nbb; jj++) {
	    int j = ib[jj];
	    double bj = Bm ? Bm[a + (R_xlen_t) l * j] : 1.0;
	    R_xlen_t cj = Bm ? (R_xlen_t) sb[m] * j : 0;
	    for (int ii = 0; ii < na; ii++) {
		int i = ia[ii];
		double c = Am[a + (R_xlen_t) l * i] * bj;
		R_xlen_t ri = ma[m] ? 0 : (R_xlen_t) sa[m] * i,
		    ci = cj + (ma[m] ? (R_xlen_t) sa[m] * i : 0);
		double *g = G + ri + p * ci, *gt = G + ci + p * ri;
		if (same) {
		    for (R_xlen_t r = 0; r < R; r++) g[ro[r] + p * co[r]] += x[r] * c;
		} else {
		    for (R_xlen_t r = 0; r < R; r++) {
			double e = x[r] * c;
			g[ro[r] + p * co[r]] += e;
			gt[co[r] + p * ro[r]] += e;
		    }
		}
	    }
	}
    }
    vmaxset(vmax);
}

/* Per-thread work space of the moment kernels: the blocks of each term
   (with the dense part), the accumulators (moments by pair of terms or by
   term, and the dense part's BLAS product) and scratch. */
typedef struct {
    dpart P;
    tblock *Bk;
    double **acc;
    double *Gd, *WB, *wN;
    int *Fx;
} mws;

static void alloc_mws(mws *ws, const sctx *S, int nacc)
{
    ws->P = dense_part(S);
    ws->Bk = alloc_blocks(S, &ws->P);
    ws->acc = (double **) R_alloc(nacc > 0 ? nacc : 1, sizeof(double *));
    ws->Gd = ws->WB = ws->wN = NULL;
    ws->Fx = NULL;
}

/* The pairs of terms of X'WX (s <= t, in order t = 0.., s = 0..t). */
typedef struct {
    int np;
    int *nx, *Lx, **xpos, **xstr;
    R_xlen_t *len;
} gpairs;

/* Rows i0 .. i1 - 1 (i0 a multiple of BLOCK) of the moments of X'WX, added
   to ws->acc and ws->Gd. */
static void sgram_rows(const sctx *S, const gpairs *G, const double *w,
		       R_xlen_t i0s, R_xlen_t i1, mws *ws)
{
    const dpart *P = &ws->P;
    tblock *Bk = ws->Bk;
    double **acc = ws->acc, *Gd = ws->Gd, *WB = ws->WB, *wN = ws->wN;
    int *Fx = ws->Fx;
    const int *nx = G->nx, *Lx = G->Lx;
    int *const *xpos = G->xpos, *const *xstr = G->xstr;
    int qd = P->q;
    for (R_xlen_t i0 = i0s; i0 < i1; i0 += BLOCK) {
	int nb = (int) (i1 - i0 < BLOCK ? i1 - i0 : BLOCK);
	for (int t = 0; t < S->nt; t++) term_block(&S->t[t], S->n, i0, nb, Bk[t].F, Bk[t].N);
	if (qd > 0) {
	    double one = 1.0;
	    int ld = BLOCK;
	    if (w) {
		for (int c = 0; c < qd; c++)
		    for (int r = 0; r < nb; r++)
			WB[r + (size_t) BLOCK * c] = w[i0 + r] * P->buf[r + (size_t) BLOCK * c];
		F77_CALL(dgemm)("T", "N", &qd, &qd, &nb, &one, P->buf, &ld, WB, &ld,
				&one, Gd, &qd FCONE FCONE);
	    } else {
		F77_CALL(dsyrk)("U", "T", &qd, &nb, &one, P->buf, &ld, &one, Gd, &qd FCONE FCONE);
	    }
	}
	int q = 0;
	for (int t = 0; t < S->nt; t++) {
	    const sterm *B = &S->t[t];
	    const int *Ft = Bk[t].F;
	    const double *Nt = Bk[t].N;
	    int nt = B->ns;
	    for (int s = 0; s <= t; s++, q++) {
		if (P->off[s] >= 0 && P->off[t] >= 0) continue; /* dense part: BLAS */
		const sterm *A = &S->t[s];
		const int *Fs = Bk[s].F;
		const double *Ns = Bk[s].N;
		int ns = A->ns, Ls = A->L;
		double *M = acc[q];
		if (s == t) {
		    for (int u = 0; u < ns; u++) {
			if (w) for (int r = 0; r < nb; r++) wN[r] = w[i0 + r] * Ns[(size_t) u * BLOCK + r];
			else for (int r = 0; r < nb; r++) wN[r] = Ns[(size_t) u * BLOCK + r];
			for (int v = 0; v < ns; v++) {
			    double *Muv = M + (R_xlen_t) Ls * (u + ns * v);
			    if (Ls == 1) {
				double sum = 0.0;
				for (int r = 0; r < nb; r++) sum += wN[r] * Ns[(size_t) v * BLOCK + r];
				Muv[0] += sum;
			    } else {
				for (int r = 0; r < nb; r++) Muv[Fs[r]] += wN[r] * Ns[(size_t) v * BLOCK + r];
			    }
			}
		    }
		} else {
		    /* the cell of the second term: the combined level of
		       its factors that the first term does not have */
		    int Lt = Lx[q];
		    const int *Fc = Ft;
		    if (nx[q] < B->nf) {
			Fc = Fx;
			for (int r = 0; r < nb; r++) Fx[r] = 0;
			for (int c = 0; c < nx[q]; c++) {
			    int g = xpos[q][c], a = B->adj[g] - 1, ls = xstr[q][c];
			    const int *code = B->code[g] + i0;
			    for (int r = 0; r < nb; r++) Fx[r] += (code[r] + a) * ls;
			}
		    }
		    R_xlen_t LL = (R_xlen_t) Ls * Lt;
		    for (int u = 0; u < ns; u++) {
			if (w) for (int r = 0; r < nb; r++) wN[r] = w[i0 + r] * Ns[(size_t) u * BLOCK + r];
			else for (int r = 0; r < nb; r++) wN[r] = Ns[(size_t) u * BLOCK + r];
			for (int v = 0; v < nt; v++) {
			    double *Muv = M + LL * (u + ns * v);
			    if (LL == 1) {
				double sum = 0.0;
				for (int r = 0; r < nb; r++) sum += wN[r] * Nt[(size_t) v * BLOCK + r];
				Muv[0] += sum;
			    } else if (Lt == 1) {
				for (int r = 0; r < nb; r++) Muv[Fs[r]] += wN[r] * Nt[(size_t) v * BLOCK + r];
			    } else if (Ls == 1) {
				for (int r = 0; r < nb; r++) Muv[Fc[r]] += wN[r] * Nt[(size_t) v * BLOCK + r];
			    } else {
				for (int r = 0; r < nb; r++)
				    Muv[Fs[r] + (R_xlen_t) Ls * Fc[r]] += wN[r] * Nt[(size_t) v * BLOCK + r];
			    }
			}
		    }
		}
	    }
	}
    }
}

/* Chunk accumulators: ws[t].acc[a] (a < nacc, `len[a]` doubles) and ws[t].Gd
   (lg doubles) for t >= 1, carved out of one buffer per thread, which the
   chunk zeroes on its thread. */
static double **moment_buffers(mws *ws, int T, int nacc, const R_xlen_t *len, R_xlen_t lg)
{
    R_xlen_t tot = lg;
    for (int a = 0; a < nacc; a++) tot += len[a];
    double **buf = (double **) R_alloc(T, sizeof(double *));
    buf[0] = NULL;
    for (int t = 1; t < T; t++) {
	double *b = (double *) R_alloc(tot > 0 ? tot : 1, sizeof(double));
	buf[t] = b;
	for (int a = 0; a < nacc; a++) { ws[t].acc[a] = b; b += len[a]; }
	ws[t].Gd = lg > 0 ? b : NULL;
    }
    return buf;
}

static R_xlen_t moment_total(int nacc, const R_xlen_t *len, R_xlen_t lg)
{
    R_xlen_t tot = lg;
    for (int a = 0; a < nacc; a++) tot += len[a];
    return tot;
}

/* ws[0].acc[a] += ws[t].acc[a], t = 1 .. T - 1 in order; the same for Gd. */
static void moment_combine(mws *ws, int T, int nacc, const R_xlen_t *len, R_xlen_t lg)
{
    double **p = (double **) R_alloc(T, sizeof(double *));
    for (int a = 0; a < nacc; a++) {
	for (int t = 0; t < T; t++) p[t] = ws[t].acc[a];
	xtx_combine(p, T, len[a]);
    }
    if (lg > 0) {
	for (int t = 0; t < T; t++) p[t] = ws[t].Gd;
	xtx_combine(p, T, lg);
    }
}

/* X'WX (p x p), from moments by pair of terms (s <= t, in order t = 0..,
   s = 0..t), assembled by gram_block(); the moments are arrays: s < t: [Ls, Lx, ns, nt] sums of w N_s[u] N_t[v] by cell of
   the union of the two terms' factors, with Lx the combined levels of the
   factors of t that s does not have (extra_factors(); the levels of the
   factors they share agree within a row, so only these cells can be
   non-zero); s == t: [Ls, ns, ns].  On several threads, fixed chunks of
   rows (whole blocks) with their own moments, added in order (see
   xtx_reduce_threads()). */
SEXP xtx_struct_gram(SEXP setup, SEXP vars, SEXP w_, SEXP nthreads_)
{
    sctx S;
    build_sctx(&S, setup, vars);
    const double *w = isNull(w_) ? NULL : REAL_RO(w_);
    if (w && XLENGTH(w_) != S.n) error("internal error: one weight per row");
    gpairs GP;
    int np = S.nt * (S.nt + 1) / 2;
    GP.np = np;
    SEXP out = PROTECT(allocVector(VECSXP, np));
    int npa = np > 0 ? np : 1;
    GP.nx = (int *) R_alloc(npa, sizeof(int));
    GP.Lx = (int *) R_alloc(npa, sizeof(int));
    GP.xpos = (int **) R_alloc(npa, sizeof(int *));
    GP.xstr = (int **) R_alloc(npa, sizeof(int *));
    GP.len = (R_xlen_t *) R_alloc(npa, sizeof(R_xlen_t));
    double **acc0 = (double **) R_alloc(npa, sizeof(double *));
    double work = 0.0;  /* multiply-adds per row */
    int q = 0;
    for (int t = 0; t < S.nt; t++)
	for (int s = 0; s <= t; s++, q++) {
	    const sterm *A = &S.t[s], *B = &S.t[t];
	    if (s < t) {
		GP.xpos[q] = (int *) R_alloc(SMAXC, sizeof(int));
		GP.xstr[q] = (int *) R_alloc(SMAXC, sizeof(int));
		GP.Lx[q] = extra_factors(A, B, &GP.nx[q], GP.xpos[q], GP.xstr[q]);
	    } else {
		GP.nx[q] = 0; GP.Lx[q] = 1; GP.xpos[q] = NULL; GP.xstr[q] = NULL;
	    }
	    R_xlen_t len = (s == t) ? (R_xlen_t) A->L * A->ns * A->ns
		: (R_xlen_t) A->L * GP.Lx[q] * A->ns * B->ns;
	    GP.len[q] = len;
	    SEXP m = SET_VECTOR_ELT(out, q, allocVector(REALSXP, len));
	    acc0[q] = REAL(m);
	    memset(acc0[q], 0, sizeof(double) * (size_t) len);
	    work += (double) A->ns * B->ns + A->ns;
	}
    int qd = dense_count(&S);
    R_xlen_t lg = (R_xlen_t) qd * qd;
    int T = 1, req = asInteger(nthreads_);
    if (req > 1)
	T = xtx_reduce_threads(req, S.n, (double) S.n * (work + 4.0 * S.nt + 1.0),
			      (double) moment_total(np, GP.len, lg));
    mws *ws = (mws *) R_alloc(T, sizeof(mws));
    for (int t = 0; t < T; t++) {
	alloc_mws(&ws[t], &S, np);
	ws[t].wN = (double *) R_alloc(BLOCK, sizeof(double));
	ws[t].Fx = (int *) R_alloc(BLOCK, sizeof(int));
	if (qd > 0 && w) ws[t].WB = (double *) R_alloc((size_t) BLOCK * qd, sizeof(double));
    }
    for (int a = 0; a < np; a++) ws[0].acc[a] = acc0[a];
    if (qd > 0) {
	ws[0].Gd = (double *) R_alloc((size_t) lg, sizeof(double));
	memset(ws[0].Gd, 0, sizeof(double) * (size_t) lg);
    }
    if (T == 1) {
	sgram_rows(&S, &GP, w, 0, S.n, &ws[0]);
    } else {
	double **buf = moment_buffers(ws, T, np, GP.len, lg);
	R_xlen_t tot = moment_total(np, GP.len, lg);
#ifdef _OPENMP
#pragma omp parallel for num_threads(T) schedule(static)
#endif
	for (int t = 0; t < T; t++) {
	    if (t > 0) xtx_zero(buf[t], tot);
	    sgram_rows(&S, &GP, w, xtx_chunk_start(S.n, T, t, BLOCK),
		       xtx_chunk_start(S.n, T, t + 1, BLOCK), &ws[t]);
	}
	moment_combine(ws, T, np, GP.len, lg);
    }
    double **acc = ws[0].acc, *Gd = ws[0].Gd;
    const dpart *P = &ws[0].P;
    /* moments of pairs of dense terms, from the BLAS accumulation (upper
       triangle without weights) */
    if (qd > 0) {
	q = 0;
	for (int t = 0; t < S.nt; t++)
	    for (int s = 0; s <= t; s++, q++) {
		if (P->off[s] < 0 || P->off[t] < 0) continue;
		int ns = S.t[s].ns, nt = S.t[t].ns;
		for (int u = 0; u < ns; u++)
		    for (int v = 0; v < nt; v++) {
			int i = P->off[s] + u, j = P->off[t] + v;
			double g = (w || i <= j) ? Gd[i + (size_t) qd * j] : Gd[j + (size_t) qd * i];
			acc[q][u + ns * v] = g;
		    }
	    }
    }
    SEXP G = PROTECT(allocMatrix(REALSXP, S.p, S.p));
    memset(REAL(G), 0, sizeof(double) * (size_t) S.p * S.p);
    q = 0;
    for (int t = 0; t < S.nt; t++)
	for (int s = 0; s <= t; s++, q++)
	    gram_block(&S.t[s], &S.t[t], s == t, acc[q], GP.nx[q], GP.xpos[q], REAL(G), S.p);
    UNPROTECT(3);
    return G;
}

/* Rows i0 .. i1 - 1 (i0 a multiple of BLOCK) of the moments of X'Y, added
   to ws->acc and ws->Gd (the dense part, qd x K). */
static void sxty_rows(const sctx *S, const double *Y, int K, R_xlen_t i0s,
		      R_xlen_t i1, mws *ws)
{
    const dpart *P = &ws->P;
    tblock *Bk = ws->Bk;
    int qd = P->q;
    for (R_xlen_t i0 = i0s; i0 < i1; i0 += BLOCK) {
	int nb = (int) (i1 - i0 < BLOCK ? i1 - i0 : BLOCK);
	for (int t = 0; t < S->nt; t++) term_block(&S->t[t], S->n, i0, nb, Bk[t].F, Bk[t].N);
	if (qd > 0) {
	    double one = 1.0;
	    int ldy = (int) S->n, ld = BLOCK;
	    F77_CALL(dgemm)("T", "N", &qd, &K, &nb, &one, P->buf, &ld, Y + i0, &ldy,
			    &one, ws->Gd, &qd FCONE FCONE);
	}
	for (int t = 0; t < S->nt; t++) {
	    if (P->off[t] >= 0) continue; /* dense part: BLAS */
	    const sterm *T = &S->t[t];
	    const int *F = Bk[t].F;
	    const double *N = Bk[t].N;
	    int ns = T->ns, L = T->L;
	    for (int k = 0; k < K; k++) {
		const double *y = Y + i0 + (R_xlen_t) k * S->n;
		for (int u = 0; u < ns; u++) {
		    double *M = ws->acc[t] + (R_xlen_t) L * (u + (R_xlen_t) ns * k);
		    if (L == 1) {
			double sum = 0.0;
			for (int r = 0; r < nb; r++) sum += N[(size_t) u * BLOCK + r] * y[r];
			M[0] += sum;
		    } else {
			for (int r = 0; r < nb; r++) M[F[r]] += N[(size_t) u * BLOCK + r] * y[r];
		    }
		}
	    }
	}
    }
}

/* Moments for X'Y (Y n x K): list per term of [L, ns, K] sums of
   N[u] Y[, k] by level.  On several threads, as xtx_struct_gram(). */
SEXP xtx_struct_xty(SEXP setup, SEXP vars, SEXP Y_, SEXP nthreads_)
{
    sctx S;
    build_sctx(&S, setup, vars);
    if (nrows(Y_) != S.n) error("non-conformable arguments");
    int K = ncols(Y_);
    const double *Y = REAL_RO(Y_);
    SEXP out = PROTECT(allocVector(VECSXP, S.nt));
    int nta = S.nt > 0 ? S.nt : 1;
    R_xlen_t *len = (R_xlen_t *) R_alloc(nta, sizeof(R_xlen_t));
    double **acc0 = (double **) R_alloc(nta, sizeof(double *));
    double work = 0.0;
    for (int t = 0; t < S.nt; t++) {
	len[t] = (R_xlen_t) S.t[t].L * S.t[t].ns * K;
	acc0[t] = REAL(SET_VECTOR_ELT(out, t, allocVector(REALSXP, len[t])));
	memset(acc0[t], 0, sizeof(double) * (size_t) len[t]);
	work += (double) S.t[t].ns * (K + 1);
    }
    int qd = dense_count(&S);
    R_xlen_t lg = (R_xlen_t) qd * K;
    int T = 1, req = asInteger(nthreads_);
    if (req > 1)
	T = xtx_reduce_threads(req, S.n, (double) S.n * (work + 4.0 * S.nt + 1.0),
			      (double) moment_total(S.nt, len, lg));
    mws *ws = (mws *) R_alloc(T, sizeof(mws));
    for (int t = 0; t < T; t++) alloc_mws(&ws[t], &S, S.nt);
    for (int a = 0; a < S.nt; a++) ws[0].acc[a] = acc0[a];
    if (qd > 0) {
	ws[0].Gd = (double *) R_alloc((size_t) lg, sizeof(double));
	memset(ws[0].Gd, 0, sizeof(double) * (size_t) lg);
    }
    if (T == 1) {
	sxty_rows(&S, Y, K, 0, S.n, &ws[0]);
    } else {
	double **buf = moment_buffers(ws, T, S.nt, len, lg);
	R_xlen_t tot = moment_total(S.nt, len, lg);
#ifdef _OPENMP
#pragma omp parallel for num_threads(T) schedule(static)
#endif
	for (int t = 0; t < T; t++) {
	    if (t > 0) xtx_zero(buf[t], tot);
	    sxty_rows(&S, Y, K, xtx_chunk_start(S.n, T, t, BLOCK),
		      xtx_chunk_start(S.n, T, t + 1, BLOCK), &ws[t]);
	}
	moment_combine(ws, T, S.nt, len, lg);
    }
    double **acc = ws[0].acc, *Rd = ws[0].Gd;
    const dpart *P = &ws[0].P;
    if (qd > 0) {
	for (int t = 0; t < S.nt; t++) {
	    if (P->off[t] < 0) continue;
	    int ns = S.t[t].ns;
	    for (int k = 0; k < K; k++)
		for (int u = 0; u < ns; u++)
		    acc[t][u + (R_xlen_t) ns * k] = Rd[P->off[t] + u + (size_t) qd * k];
	}
    }
    UNPROTECT(2);
    return out;
}

/* X B from per-term tables P_t [L, ns, K] (P_t = C_t B_t, built in R):
   row i, column k is sum over t, u of N_t(i)[u] P_t[F_t(i), u, k].  Rows
   are computed in blocks of BLOCK rows, each the same on any number of
   threads. */
SEXP xtx_struct_mult(SEXP setup, SEXP vars, SEXP P, SEXP K_, SEXP nthreads_)
{
    sctx S;
    build_sctx(&S, setup, vars);
    int K = asInteger(K_);
    if (LENGTH(P) != S.nt) error("internal error: one table per term");
    SEXP out = PROTECT(allocMatrix(REALSXP, S.n, K));
    double *o = REAL(out);
    memset(o, 0, sizeof(double) * (size_t) S.n * K);
    R_xlen_t nblocks = xtx_block_count(S.n, BLOCK);
    int nthreads = xtx_worker_count(asInteger(nthreads_), nblocks);
    dpart *DP = (dpart *) R_alloc(nthreads, sizeof(dpart));
    tblock **Bk = (tblock **) R_alloc(nthreads, sizeof(tblock *));
    for (int t = 0; t < nthreads; t++) {
	DP[t] = dense_part(&S);
	Bk[t] = alloc_blocks(&S, &DP[t]);
    }
    int qd = DP[0].q;
    double *Pd = NULL;
    if (qd > 0) {
	/* the dense terms' tables, stacked: qd x K */
	Pd = (double *) R_alloc((size_t) qd * K, sizeof(double));
	for (int t = 0; t < S.nt; t++) {
	    if (DP[0].off[t] < 0) continue;
	    int ns = S.t[t].ns;
	    const double *p = REAL(VECTOR_ELT(P, t));
	    for (int k = 0; k < K; k++)
		for (int u = 0; u < ns; u++)
		    Pd[DP[0].off[t] + u + (size_t) qd * k] = p[u + (R_xlen_t) ns * k];
	}
    }
    const double **ptab = (const double **) R_alloc(S.nt > 0 ? S.nt : 1, sizeof(double *));
    for (int t = 0; t < S.nt; t++) ptab[t] = REAL(VECTOR_ELT(P, t));
#ifdef _OPENMP
#pragma omp parallel for if(nthreads > 1) num_threads(nthreads) schedule(static)
#endif
    for (R_xlen_t b = 0; b < nblocks; b++) {
	int th = 0;
#ifdef _OPENMP
	th = omp_get_thread_num();
#endif
	tblock *Bt = Bk[th];
	const dpart *Dt = &DP[th];
	R_xlen_t i0 = b * BLOCK;
	int nb = (int) (S.n - i0 < BLOCK ? S.n - i0 : BLOCK);
	for (int t = 0; t < S.nt; t++) term_block(&S.t[t], S.n, i0, nb, Bt[t].F, Bt[t].N);
	if (qd > 0) {
	    double one = 1.0;
	    int ldo = (int) S.n, ld = BLOCK;
	    F77_CALL(dgemm)("N", "N", &nb, &K, &qd, &one, Dt->buf, &ld, Pd, &qd,
			    &one, o + i0, &ldo FCONE FCONE);
	}
	for (int t = 0; t < S.nt; t++) {
	    if (Dt->off[t] >= 0) continue; /* dense part: BLAS */
	    const sterm *T = &S.t[t];
	    const int *F = Bt[t].F;
	    const double *N = Bt[t].N;
	    const double *p = ptab[t];
	    int ns = T->ns, L = T->L;
	    for (int k = 0; k < K; k++) {
		double *ok = o + i0 + (R_xlen_t) k * S.n;
		for (int u = 0; u < ns; u++) {
		    const double *pu = p + (R_xlen_t) L * (u + (R_xlen_t) ns * k);
		    if (L == 1) {
			double c = pu[0];
			for (int r = 0; r < nb; r++) ok[r] += N[(size_t) u * BLOCK + r] * c;
		    } else {
			for (int r = 0; r < nb; r++) ok[r] += N[(size_t) u * BLOCK + r] * pu[F[r]];
		    }
		}
	    }
	}
    }
    UNPROTECT(2);
    return out;
}

/* `len` doubles at an address that is a multiple of 64 bytes.  A BLAS
   routine can round differently according to the alignment of its
   arguments (it handles the unaligned head of a vector apart); with
   operands always aligned the same way, the same row gives the same bits
   in every call. */
static double *alloc_aligned(size_t len)
{
    char *p = R_alloc(len * sizeof(double) + 64, 1);
    return (double *) (p + ((64 - ((uintptr_t) p & 63)) & 63));
}

/* diag(X S X'), S symmetric: row i is x_i' S x_i.  The columns split into
   the dense part z (terms with one level: the intercept and products of
   numeric variables, non-zero in every row) and the terms with factors
   x_t = N_t(i) (x) C_t[F_t(i), ], so that

     x' S x = z' S_dd z + 2 sum_t z' S_dt x_t + sum_{s, t} x_s' S_st x_t.

   The first part is a matrix-vector product and a
   dot product for each row (BLAS); the other
   two read tables built in R from S and the contrasts:
     A[[t]]     [qd, L_t, ns_t]: S_dt contracted with C_t,
     Mself[[t]] [L_t, ns_t, ns_t]: C_t[F, ] S_tt C_t[F, ]' for each level F,
     Mpair[[j]] [L_s, L_t, ns_s, ns_t] for the pair pairs[, j] = (s, t),
                s < t (0-based terms): C_s S_st C_t'.
   Entries of the lists for terms of the dense part are ignored.  Rows are
   computed in blocks of BLOCK rows, each the same on any number of
   threads.  Sdd is qd x qd, in the order of dense_part(). */
SEXP xtx_struct_quad(SEXP setup, SEXP vars, SEXP Sdd_, SEXP A_, SEXP Mself_,
		     SEXP pairs_, SEXP Mpair_, SEXP nthreads_)
{
    sctx S;
    build_sctx(&S, setup, vars);
    if (LENGTH(A_) != S.nt || LENGTH(Mself_) != S.nt)
	error("internal error: one table per term");
    int np = LENGTH(Mpair_);
    if (np > 0 && (nrows(pairs_) != 2 || ncols(pairs_) != np))
	error("internal error: one table per pair of terms");
    const int *pr = np > 0 ? INTEGER(pairs_) : NULL;
    SEXP out = PROTECT(allocVector(REALSXP, S.n));
    double *o = REAL(out);
    R_xlen_t nblocks = xtx_block_count(S.n, BLOCK);
    int nthreads = xtx_worker_count(asInteger(nthreads_), nblocks);
    dpart *DP = (dpart *) R_alloc(nthreads, sizeof(dpart));
    tblock **Bk = (tblock **) R_alloc(nthreads, sizeof(tblock *));
    for (int t = 0; t < nthreads; t++) {
	DP[t] = dense_part(&S);
	Bk[t] = alloc_blocks(&S, &DP[t]);
    }
    int qd = DP[0].q, nfac = 0;
    for (int t = 0; t < S.nt; t++) nfac += DP[0].off[t] < 0;
    if (qd > 0 && (nrows(Sdd_) != qd || ncols(Sdd_) != qd))
	error("internal error: the dense block of 'S' must be %d x %d", qd, qd);
    double *Sdd = NULL;
    if (qd > 0) {
	Sdd = alloc_aligned((size_t) qd * qd);
	memcpy(Sdd, REAL(Sdd_), sizeof(double) * (size_t) qd * qd);
    }
    /* per thread: Z by row (qd x BLOCK), one row z, and S_dd z */
    double **Zt = (double **) R_alloc(nthreads, sizeof(double *));
    double **Zr = (double **) R_alloc(nthreads, sizeof(double *));
    double **Zs = (double **) R_alloc(nthreads, sizeof(double *));
    for (int t = 0; t < nthreads; t++) {
	Zt[t] = qd > 0 ? (double *) R_alloc((size_t) BLOCK * qd, sizeof(double)) : NULL;
	Zr[t] = qd > 0 ? alloc_aligned((size_t) qd) : NULL;
	Zs[t] = qd > 0 ? alloc_aligned((size_t) qd) : NULL;
    }
    const double **atab = (const double **) R_alloc(S.nt > 0 ? S.nt : 1, sizeof(double *));
    const double **stab = (const double **) R_alloc(S.nt > 0 ? S.nt : 1, sizeof(double *));
    for (int t = 0; t < S.nt; t++) {
	int fac = DP[0].off[t] < 0;
	atab[t] = fac && qd > 0 ? REAL(VECTOR_ELT(A_, t)) : NULL;
	stab[t] = fac ? REAL(VECTOR_ELT(Mself_, t)) : NULL;
    }
    const double **ptab = (const double **) R_alloc(np > 0 ? np : 1, sizeof(double *));
    for (int j = 0; j < np; j++) ptab[j] = REAL(VECTOR_ELT(Mpair_, j));
#ifdef _OPENMP
#pragma omp parallel for if(nthreads > 1) num_threads(nthreads) schedule(static)
#endif
    for (R_xlen_t b = 0; b < nblocks; b++) {
	int th = 0;
#ifdef _OPENMP
	th = omp_get_thread_num();
#endif
	tblock *Bt = Bk[th];
	const dpart *Dt = &DP[th];
	R_xlen_t i0 = b * BLOCK;
	int nb = (int) (S.n - i0 < BLOCK ? S.n - i0 : BLOCK);
	double *ob = o + i0;
	for (int r = 0; r < nb; r++) ob[r] = 0.0;
	for (int t = 0; t < S.nt; t++) term_block(&S.t[t], S.n, i0, nb, Bt[t].F, Bt[t].N);
	const double *zt = Zt[th];
	if (qd > 0) {
	    double *w = Zt[th];
	    for (int j = 0; j < qd; j++) {
		const double *z = Dt->buf + (size_t) j * BLOCK;
		for (int r = 0; r < nb; r++) w[j + (size_t) qd * r] = z[r];
	    }
	    /* z' S_dd z, one row at a time on a copy of the row at a fixed
	       address: the result for a row does not depend on its position
	       or on the other rows, so that it is the same bit for bit when
	       only the distinct rows of a design are computed. */
	    double one = 1.0, zero = 0.0, *zc = Zr[th], *sz = Zs[th];
	    int inc = 1;
	    for (int r = 0; r < nb; r++) {
		memcpy(zc, w + (size_t) qd * r, sizeof(double) * qd);
		F77_CALL(dgemv)("N", &qd, &qd, &one, Sdd, &qd, zc, &inc, &zero, sz, &inc FCONE);
		ob[r] += F77_CALL(ddot)(&qd, zc, &inc, sz, &inc);
	    }
	}
	for (int t = 0; t < S.nt; t++) {
	    if (Dt->off[t] >= 0) continue;
	    const sterm *T = &S.t[t];
	    const int *F = Bt[t].F;
	    const double *N = Bt[t].N;
	    int ns = T->ns, L = T->L;
	    /* 2 z' S_dt x_t */
	    if (qd > 0) {
		const double *A = atab[t];
		for (int u = 0; u < ns; u++) {
		    const double *Nu = N + (size_t) u * BLOCK;
		    for (int r = 0; r < nb; r++) {
			const double *a = A + (size_t) qd * (F[r] + (size_t) L * u);
			const double *z = zt + (size_t) qd * r;
			double acc = 0.0;
			for (int j = 0; j < qd; j++) acc += z[j] * a[j];
			ob[r] += 2.0 * Nu[r] * acc;
		    }
		}
	    }
	    /* x_t' S_tt x_t */
	    const double *M = stab[t];
	    for (int u = 0; u < ns; u++) {
		const double *Nu = N + (size_t) u * BLOCK;
		for (int v = 0; v < ns; v++) {
		    const double *Nv = N + (size_t) v * BLOCK;
		    const double *m = M + (size_t) L * (u + (size_t) ns * v);
		    for (int r = 0; r < nb; r++) ob[r] += Nu[r] * Nv[r] * m[F[r]];
		}
	    }
	}
	/* 2 x_s' S_st x_t, s < t */
	for (int j = 0; j < np; j++) {
	    int s = pr[2 * j], t = pr[2 * j + 1];
	    const sterm *Ts = &S.t[s], *Tt = &S.t[t];
	    const int *Fs = Bt[s].F, *Ft = Bt[t].F;
	    const double *Ns = Bt[s].N, *Nt = Bt[t].N;
	    size_t Ls = (size_t) Ts->L, Lt = (size_t) Tt->L;
	    int nss = Ts->ns, nst = Tt->ns;
	    const double *M = ptab[j];
	    for (int u = 0; u < nss; u++) {
		const double *Nu = Ns + (size_t) u * BLOCK;
		for (int v = 0; v < nst; v++) {
		    const double *Nv = Nt + (size_t) v * BLOCK;
		    const double *m = M + Ls * Lt * (u + (size_t) nss * v);
		    for (int r = 0; r < nb; r++)
			ob[r] += 2.0 * Nu[r] * Nv[r] * m[Fs[r] + Ls * Ft[r]];
		}
	    }
	}
    }
    UNPROTECT(2);
    return out;
}

/* A lower bound on the absolute value of every non-zero entry of the
   dense matrix of a plain design: the minimum over terms of the product
   of its components' smallest non-zero absolute values (input values or
   contrast entries).  The intercept counts as 1.  Non-finite inputs are
   ignored (the caller handles them). */
SEXP xtx_min_nonzero(SEXP setup, SEXP vars)
{
    sctx S;
    build_sctx(&S, setup, vars);
    double out = 1.0;
    for (int t = 0; t < S.nt; t++) {
	const sterm *T = &S.t[t];
	double prod = 1.0;
	for (int c = 0; c < T->nf; c++) {
	    double mn = R_PosInf;
	    R_xlen_t len = (R_xlen_t) T->nrc[c] * T->fnc[c];
	    for (R_xlen_t e = 0; e < len; e++) {
		double x = fabs(T->cm[c][e]);
		if (x > 0.0 && x < mn) mn = x;
	    }
	    if (R_FINITE(mn)) prod *= mn;
	}
	for (int c = 0; c < T->nn; c++) {
	    double mn = R_PosInf;
	    const double *v = T->num[c];
	    for (R_xlen_t e = 0; e < S.n * T->nnc[c]; e++) {
		double x = fabs(v[e]);
		if (x > 0.0 && x < mn) mn = x;
	    }
	    if (R_FINITE(mn)) prod *= mn;
	}
	if (prod < out) out = prod;
    }
    UNPROTECT(1);
    return ScalarReal(out);
}

/* Are all values of a double vector or matrix finite?  (No allocation.)
   x * 0 is 0 for finite x and NaN otherwise: sums of those in four
   independent chains, over chunks, are much faster than a test per
   value. */
SEXP xtx_all_finite(SEXP x)
{
    const double *v = REAL(x);
    R_xlen_t n = XLENGTH(x);
    for (R_xlen_t i0 = 0; i0 < n; i0 += 4096) {
	R_xlen_t i1 = n - i0 < 4096 ? n : i0 + 4096, i = i0;
	double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
	for (; i + 4 <= i1; i += 4) {
	    a0 += v[i] * 0.0;
	    a1 += v[i + 1] * 0.0;
	    a2 += v[i + 2] * 0.0;
	    a3 += v[i + 3] * 0.0;
	}
	for (; i < i1; i++) a0 += v[i] * 0.0;
	if (a0 + a1 + a2 + a3 != 0.0) return ScalarLogical(0);
    }
    return ScalarLogical(1);
}

/* Euclidean norms of the columns of a finite double matrix, scaled so
   that they neither overflow nor underflow (as LAPACK's dnrm2). */
SEXP xtx_col_norms(SEXP x)
{
    int nr = nrows(x), nc = ncols(x);
    const double *v = REAL(x);
    SEXP out = PROTECT(allocVector(REALSXP, nc));
    double *o = REAL(out);
    for (int j = 0; j < nc; j++) {
	const double *c = v + (R_xlen_t) j * nr;
	double scale = 0.0, ssq = 1.0;
	for (int i = 0; i < nr; i++) {
	    double a = fabs(c[i]);
	    if (a == 0.0) continue;
	    if (scale < a) {
		ssq = 1.0 + ssq * (scale / a) * (scale / a);
		scale = a;
	    } else {
		ssq += (a / scale) * (a / scale);
	    }
	}
	o[j] = scale * sqrt(ssq);
    }
    UNPROTECT(1);
    return out;
}

/* Are all values finite, with non-zero absolute values within [lo, hi]?
   (No allocation.) */
SEXP xtx_in_range(SEXP x, SEXP lo_, SEXP hi_)
{
    const double *v = REAL(x);
    double lo = asReal(lo_), hi = asReal(hi_);
    R_xlen_t n = XLENGTH(x);
    for (R_xlen_t i = 0; i < n; i++) {
	double a = fabs(v[i]);
	if (!R_FINITE(v[i]) || a > hi || (a != 0.0 && a < lo)) return ScalarLogical(0);
    }
    return ScalarLogical(1);
}
