#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "xtx.h"
#ifdef _OPENMP
#include <omp.h>
#endif

SEXP xtx_omp_max_threads(void)
{
#ifdef _OPENMP
    return ScalarInteger(omp_get_max_threads());
#else
    return ScalarInteger(1);
#endif
}

static const R_CallMethodDef CallEntries[] = {
    {"xtx_setup", (DL_FUNC) &xtx_setup, 3},
    {"xtx_fill", (DL_FUNC) &xtx_fill, 3},
    {"xtx_jac_mult", (DL_FUNC) &xtx_jac_mult, 9},
    {"xtx_jac_crossprod", (DL_FUNC) &xtx_jac_crossprod, 11},
    {"xtx_jac_quad", (DL_FUNC) &xtx_jac_quad, 9},
    {"xtx_jac_sparse", (DL_FUNC) &xtx_jac_sparse, 6},
    {"xtx_jac_dense", (DL_FUNC) &xtx_jac_dense, 7},
    {"xtx_row_groups", (DL_FUNC) &xtx_row_groups, 8},
    {"xtx_jac_entries_mean", (DL_FUNC) &xtx_jac_entries_mean, 7},
    {"xtx_jac_finite", (DL_FUNC) &xtx_jac_finite, 7},
    {"xtx_jac_gram", (DL_FUNC) &xtx_jac_gram, 15},
    {"xtx_chol_pivot", (DL_FUNC) &xtx_chol_pivot, 8},
    {"xtx_struct_ok", (DL_FUNC) &xtx_struct_ok, 3},
    {"xtx_struct_info", (DL_FUNC) &xtx_struct_info, 2},
    {"xtx_struct_gram", (DL_FUNC) &xtx_struct_gram, 4},
    {"xtx_struct_xty", (DL_FUNC) &xtx_struct_xty, 4},
    {"xtx_struct_mult", (DL_FUNC) &xtx_struct_mult, 5},
    {"xtx_struct_quad", (DL_FUNC) &xtx_struct_quad, 8},
    {"xtx_kron_mode", (DL_FUNC) &xtx_kron_mode, 3},
    {"xtx_min_nonzero", (DL_FUNC) &xtx_min_nonzero, 2},
    {"xtx_all_finite", (DL_FUNC) &xtx_all_finite, 1},
    {"xtx_col_norms", (DL_FUNC) &xtx_col_norms, 1},
    {"xtx_in_range", (DL_FUNC) &xtx_in_range, 3},
    {"xtx_codes_ok", (DL_FUNC) &xtx_codes_ok, 2},
    {"xtx_int_unique", (DL_FUNC) &xtx_int_unique, 1},
    {"xtx_str_unique", (DL_FUNC) &xtx_str_unique, 1},
    {"xtx_fe_gsum", (DL_FUNC) &xtx_fe_gsum, 4},
    {"xtx_fe_crosstabs", (DL_FUNC) &xtx_fe_crosstabs, 6},
    {"xtx_fe_struct_gsum", (DL_FUNC) &xtx_fe_struct_gsum, 9},
    {"xtx_fe_code", (DL_FUNC) &xtx_fe_code, 2},
    {"xtx_fe_recode", (DL_FUNC) &xtx_fe_recode, 2},
    {"xtx_fe_nested", (DL_FUNC) &xtx_fe_nested, 3},
    {"xtx_fe_lookup", (DL_FUNC) &xtx_fe_lookup, 2},
    {"xtx_ifelse", (DL_FUNC) &xtx_ifelse, 3},
    {"xtx_notnan", (DL_FUNC) &xtx_notnan, 1},
    {"xtx_fe_moments", (DL_FUNC) &xtx_fe_moments, 2},
    {"xtx_unwrap_vars", (DL_FUNC) &xtx_unwrap_vars, 1},
    {"xtx_levels_unused", (DL_FUNC) &xtx_levels_unused, 2},
    {"xtx_omp_max_threads", (DL_FUNC) &xtx_omp_max_threads, 0},
    {NULL, NULL, 0}
};

void R_init_xtx(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
