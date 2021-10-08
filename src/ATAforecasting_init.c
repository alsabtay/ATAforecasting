#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ATAforecasting_ATAHoldoutForecast(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_SubATA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_SubATACore(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_SubATACoreHoldhin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_SubATACoreHoldout(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_SubATADamped(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_SubATADampedHoldhin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_SubATADampedHoldout(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_SubATAHoldhin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_SubATAHoldout(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_inMASE(SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_inMASEholdin(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_meanIT(SEXP, SEXP);
extern SEXP _ATAforecasting_NaiveSD(SEXP, SEXP);
extern SEXP _ATAforecasting_NaiveSDholdin(SEXP, SEXP, SEXP);
extern SEXP _ATAforecasting_outMASE(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_ATAforecasting_ATAHoldoutForecast",   (DL_FUNC) &_ATAforecasting_ATAHoldoutForecast,   11},
    {"_ATAforecasting_SubATA",              (DL_FUNC) &_ATAforecasting_SubATA,              21},
    {"_ATAforecasting_SubATACore",          (DL_FUNC) &_ATAforecasting_SubATACore,          12},
    {"_ATAforecasting_SubATACoreHoldhin",   (DL_FUNC) &_ATAforecasting_SubATACoreHoldhin,   13},
    {"_ATAforecasting_SubATACoreHoldout",   (DL_FUNC) &_ATAforecasting_SubATACoreHoldout,   12},
    {"_ATAforecasting_SubATADamped",        (DL_FUNC) &_ATAforecasting_SubATADamped,        17},
    {"_ATAforecasting_SubATADampedHoldhin", (DL_FUNC) &_ATAforecasting_SubATADampedHoldhin, 18},
    {"_ATAforecasting_SubATADampedHoldout", (DL_FUNC) &_ATAforecasting_SubATADampedHoldout, 17},
    {"_ATAforecasting_SubATAHoldhin",       (DL_FUNC) &_ATAforecasting_SubATAHoldhin,       22},
    {"_ATAforecasting_SubATAHoldout",       (DL_FUNC) &_ATAforecasting_SubATAHoldout,       21},
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
