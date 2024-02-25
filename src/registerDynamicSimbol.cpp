#include "RInterface.h"

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

static const R_CallMethodDef CallEntries[] = {
  {"RFnDist",       (DL_FUNC) &RFnDist,       6},
  {"Rfgrad",        (DL_FUNC) &Rfgrad,        6},
  {"RFnDist1",      (DL_FUNC) &RFnDist1,      7},
  {"Rfgrad1",       (DL_FUNC) &Rfgrad1,       7},
  {"Rfhess1",       (DL_FUNC) &Rfhess1,       7},

  {NULL, NULL, 0}
};

void R_init_HiDimDA(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
