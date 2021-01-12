#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ATAforecasting_ATAHoldoutForecast(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_AutoATA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_AutoATACore(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_AutoATACoreHoldhin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_AutoATACoreHoldout(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_AutoATADamped(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_AutoATADampedHoldhin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_AutoATADampedHoldout(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_AutoATAHoldhin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_AutoATAHoldout(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_inMASE(SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_inMASEholdin(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_meanIT(SEXP, SEXP);
extern SEXP _ATAforecasting_NaiveSD(SEXP, SEXP);
extern SEXP _ATAforecasting_NaiveSDholdin(SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_outMASE(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_ATAforecasting_ATAHoldoutForecast",   (DL_FUNC) &_ATAforecasting_ATAHoldoutForecast,   11},
    {"_ATAforecasting_AutoATA",              (DL_FUNC) &_ATAforecasting_AutoATA,              20},
    {"_ATAforecasting_AutoATACore",          (DL_FUNC) &_ATAforecasting_AutoATACore,          11},
    {"_ATAforecasting_AutoATACoreHoldhin",   (DL_FUNC) &_ATAforecasting_AutoATACoreHoldhin,   12},
    {"_ATAforecasting_AutoATACoreHoldout",   (DL_FUNC) &_ATAforecasting_AutoATACoreHoldout,   12},
    {"_ATAforecasting_AutoATADamped",        (DL_FUNC) &_ATAforecasting_AutoATADamped,        16},
    {"_ATAforecasting_AutoATADampedHoldhin", (DL_FUNC) &_ATAforecasting_AutoATADampedHoldhin, 17},
    {"_ATAforecasting_AutoATADampedHoldout", (DL_FUNC) &_ATAforecasting_AutoATADampedHoldout, 17},
    {"_ATAforecasting_AutoATAHoldhin",       (DL_FUNC) &_ATAforecasting_AutoATAHoldhin,       21},
    {"_ATAforecasting_AutoATAHoldout",       (DL_FUNC) &_ATAforecasting_AutoATAHoldout,       21},
    {"_ATAforecasting_inMASE",               (DL_FUNC) &_ATAforecasting_inMASE,                3},
    {"_ATAforecasting_inMASEholdin",         (DL_FUNC) &_ATAforecasting_inMASEholdin,          4},
    {"_ATAforecasting_meanIT",               (DL_FUNC) &_ATAforecasting_meanIT,                2},
    {"_ATAforecasting_NaiveSD",              (DL_FUNC) &_ATAforecasting_NaiveSD,               2},
    {"_ATAforecasting_NaiveSDholdin",        (DL_FUNC) &_ATAforecasting_NaiveSDholdin,         3},
    {"_ATAforecasting_outMASE",              (DL_FUNC) &_ATAforecasting_outMASE,               4},
    {NULL, NULL, 0}
};

void R_init_ATAforecasting(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
