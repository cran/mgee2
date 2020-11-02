#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void Cgetmgee2_i(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Cgetmgee2v_i(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Cgetordgee2_i(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"Cgetmgee2_i",   (DL_FUNC) &Cgetmgee2_i,   16},
    {"Cgetmgee2v_i",  (DL_FUNC) &Cgetmgee2v_i,  29},
    {"Cgetordgee2_i", (DL_FUNC) &Cgetordgee2_i, 13},
    {NULL, NULL, 0}
};

void R_init_mgee2(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
