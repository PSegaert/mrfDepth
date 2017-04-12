#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void adjprojout(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void HSDND(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void medcoupleC(void *, void *, void *, void *);
extern void projoutlyingness(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP dirOutl_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(bagplotf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(halfmed2d)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hsdep2)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hsdep3)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hsdepth_deepest)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(iso2hdw)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rdepth3)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rdepth4)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rdepthnd)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sweepmedres)(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"adjprojout",       (DL_FUNC) &adjprojout,       12},
  {"HSDND",            (DL_FUNC) &HSDND,            10},
  {"medcoupleC",       (DL_FUNC) &medcoupleC,        4},
  {"projoutlyingness", (DL_FUNC) &projoutlyingness, 15},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"dirOutl_cpp",     (DL_FUNC) &dirOutl_cpp,     7},
  {NULL, NULL, 0}
};



static const R_FortranMethodDef FortranEntries[] = {
  {"bagplotf",        (DL_FUNC) &F77_NAME(bagplotf),        18},
  {"halfmed2d",       (DL_FUNC) &F77_NAME(halfmed2d),        6},
  {"hsdep2",          (DL_FUNC) &F77_NAME(hsdep2),           8},
  {"hsdep3",          (DL_FUNC) &F77_NAME(hsdep3),          10},
  {"hsdepth_deepest", (DL_FUNC) &F77_NAME(hsdepth_deepest), 12},
  {"iso2hdw",         (DL_FUNC) &F77_NAME(iso2hdw),          9},
  {"rdepth3",         (DL_FUNC) &F77_NAME(rdepth3),         10},
  {"rdepth4",         (DL_FUNC) &F77_NAME(rdepth4),          9},
  {"rdepthnd",        (DL_FUNC) &F77_NAME(rdepthnd),         9},
  {"sweepmedres",     (DL_FUNC) &F77_NAME(sweepmedres),      7},
  {NULL, NULL, 0}
};

void R_init_mrfDepth(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}