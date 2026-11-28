#ifndef XTX_H
#define XTX_H
#include <Rinternals.h>

/* Ceiling division for nonnegative row counts and positive block sizes. */
static inline R_xlen_t xtx_block_count(R_xlen_t rows, R_xlen_t size)
{
    return rows == 0 ? 0 : 1 + (rows - 1) / size;
}

static inline int xtx_worker_count(int requested, R_xlen_t units)
{
#ifdef _OPENMP
    if (requested > 1 && units > 1)
        return requested < units ? requested : (int) units;
#endif
    return 1;
}

/* Reductions (sums over rows) on several threads: the rows are split into
   T chunks whose boundaries depend only on the number of rows and T, each
   chunk is summed into its own accumulator, and the accumulators are added
   in chunk order afterwards.  Results are identical for a given T (which
   depends only on the requested threads and the problem's dimensions);
   with T = 1 the code path is the single-threaded one.

   xtx_reduce_threads(): T for a reduction over `rows` rows costing `work`
   multiply-adds in all, with accumulators of `acc` doubles.  One thread for
   small problems (fewer than XTX_RED_MIN_ROWS rows or XTX_RED_MIN_WORK
   multiply-adds per thread), and fewer threads when the extra accumulators
   would take more than XTX_RED_MAX_BYTES, or adding them up would cost more
   than a quarter of the work. */
#define XTX_RED_MIN_ROWS 8192
#define XTX_RED_MIN_WORK 2.5e5
#define XTX_RED_MAX_BYTES 268435456.0

static inline int xtx_reduce_threads(int requested, R_xlen_t rows, double work, double acc)
{
    int T = xtx_worker_count(requested, rows / XTX_RED_MIN_ROWS);
    while (T > 1) {
	double extra = (double) (T - 1) * acc;
	if (work >= XTX_RED_MIN_WORK * T && extra * sizeof(double) <= XTX_RED_MAX_BYTES &&
	    4.0 * extra <= work)
	    break;
	T--;
    }
    return T;
}

/* First row of chunk t of T (t = T gives `rows`): chunks are whole blocks
   of `unit` rows, as even as possible. */
static inline R_xlen_t xtx_chunk_start(R_xlen_t rows, int T, int t, R_xlen_t unit)
{
    R_xlen_t nb = xtx_block_count(rows, unit);
    R_xlen_t i = (nb * t / T) * unit;
    return i < rows ? i : rows;
}

/* acc[0][j] += acc[1][j] + ... + acc[T - 1][j], in that order, j < len. */
static inline void xtx_combine(double **acc, int T, R_xlen_t len)
{
    if (T <= 1) return;
#ifdef _OPENMP
#pragma omp parallel for num_threads(T) schedule(static) if(len >= 65536)
#endif
    for (R_xlen_t j = 0; j < len; j++) {
	double s = acc[0][j];
	for (int t = 1; t < T; t++) s += acc[t][j];
	acc[0][j] = s;
    }
}

/* Accumulators for chunks 1 .. T - 1 (chunk 0 uses `first`): acc[t]
   points to `len` doubles, which the chunk zeroes (xtx_zero()) on its
   thread. */
static inline double **xtx_chunk_acc(double *first, int T, R_xlen_t len)
{
    double **acc = (double **) R_alloc(T, sizeof(double *));
    acc[0] = first;
    for (int t = 1; t < T; t++)
	acc[t] = (double *) R_alloc(len > 0 ? len : 1, sizeof(double));
    return acc;
}

static inline void xtx_zero(double *x, R_xlen_t len)
{
    for (R_xlen_t j = 0; j < len; j++) x[j] = 0.0;
}

/* Slots of the list returned by xtx_setup(), in order.  Designs saved
   before the contrast cache was introduced stop at XTX_SETUP_NTERMS: check
   XLENGTH(setup) > XTX_SETUP_CONTRAST_ROWS before reading that slot. */
enum {
    XTX_SETUP_FACTORS = 0,
    XTX_SETUP_CONTR1,
    XTX_SETUP_CONTR2,
    XTX_SETUP_NLEVS,
    XTX_SETUP_COLUMNS,
    XTX_SETUP_COUNT,
    XTX_SETUP_RHS_RESPONSE,
    XTX_SETUP_INTERCEPT,
    XTX_SETUP_XNAMES,
    XTX_SETUP_ASSIGN,
    XTX_SETUP_NVAR,
    XTX_SETUP_NTERMS,
    XTX_SETUP_CONTRAST_ROWS,
    XTX_SETUP_LEN
};

SEXP xtx_setup(SEXP terms, SEXP vars, SEXP rho);
SEXP xtx_vars(SEXP setup, SEXP vars);
SEXP xtx_fill(SEXP setup, SEXP vars, SEXP nthreads);
SEXP xtx_jac_mult(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_jac_crossprod(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_jac_quad(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_jac_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_jac_dense(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_row_groups(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_jac_entries_mean(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_jac_finite(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_jac_gram(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_chol_pivot(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_struct_ok(SEXP, SEXP, SEXP);
SEXP xtx_struct_info(SEXP, SEXP);
SEXP xtx_struct_gram(SEXP, SEXP, SEXP, SEXP);
SEXP xtx_struct_xty(SEXP, SEXP, SEXP, SEXP);
SEXP xtx_struct_mult(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_struct_quad(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_kron_mode(SEXP, SEXP, SEXP);
SEXP xtx_min_nonzero(SEXP, SEXP);
SEXP xtx_all_finite(SEXP);
SEXP xtx_col_norms(SEXP);
SEXP xtx_in_range(SEXP, SEXP, SEXP);
SEXP xtx_codes_ok(SEXP, SEXP);
SEXP xtx_int_unique(SEXP);
SEXP xtx_str_unique(SEXP);
SEXP xtx_fe_gsum(SEXP, SEXP, SEXP, SEXP);
SEXP xtx_fe_crosstabs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_fe_struct_gsum(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP xtx_fe_code(SEXP, SEXP);
SEXP xtx_fe_recode(SEXP, SEXP);
SEXP xtx_fe_nested(SEXP, SEXP, SEXP);
SEXP xtx_fe_lookup(SEXP, SEXP);
/* fused.c */
SEXP xtx_ifelse(SEXP, SEXP, SEXP);
SEXP xtx_notnan(SEXP);
SEXP xtx_fe_moments(SEXP, SEXP);
SEXP xtx_unwrap_vars(SEXP);
SEXP xtx_levels_unused(SEXP, SEXP);
SEXP xtx_omp_max_threads(void);
/* error unless the integer vector codes is in 1..G (no NA), fe.c */
void xtx_check_codes(SEXP codes, int G);
#endif
