/* The first half of this file is from */
/* ugee.h in the gee package */
#include <R.h>
// #include "S.h"
#include <Rmath.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R_ext/Applic.h> /* BLAS */
#include <R_ext/Linpack.h>
#include <inttypes.h>

#define PERMANENT 1
#define EPHEMERAL 0

typedef struct matrix
		{
		int nrows, ncols;
		double *data;
		int permanence;
		} MATRIX;
// class MATRIX
// {
//   int nrows, ncols;
//   double *data;
//   int permanence;
// } ;
#define ELREF( matp , s1, s2 ) ((matp)->data)+(s2)+((s1)*(matp->ncols))
#define MEL(X ,i, j) (*(ELREF( (X), (i), (j) ) ))
#define get_nelem( x ) (((x)->nrows) * ((x)->ncols))

#define malloc(n) S_alloc(n, 1)
#define calloc S_alloc
#define free(p) {p;}
#define cfree(p) {p;}
#define is_permanent( x ) (x)->permanence == PERMANENT
#define is_ephemeral( x ) (x)->permanence == EPHEMERAL
#define make_permanent( x ) (x)->permanence = PERMANENT;
#define make_ephemeral( x ) (x)->permanence = EPHEMERAL;
#define free_if_ephemeral( x ) if (is_ephemeral((x))) VC_GEE_destroy_matrix((x))

#define from_S( Sdblptr , Srowintptr , Scolintptr , Matptr ) \
Matptr = VC_GEE_create_matrix( (int)*Srowintptr, (int)*Scolintptr , EPHEMERAL ); \
{ \
int i, j, Scol, Srow; \
double *Sload; \
Scol = *Scolintptr; \
Srow = *Srowintptr; \
Sload = Sdblptr; \
for ( j = 0 ; j < Scol ; j++ ) \
	{ \
	for ( i = 0 ; i < Srow ; i++ ) \
		{ \
		MEL( Matptr , i , j ) = (double) * ( Sload ++ ); \
		} \
	} \
}

#define to_S( Matptr, Sdblptr ) \
{ \
int i, j; \
double *Sload; \
Sload = Sdblptr; \
for ( j = 0 ; j < Matptr->ncols ; j++ ) \
	{ \
	for ( i = 0 ; i < Matptr->nrows ; i++ ) \
		{ \
		* ( Sload ++ ) = MEL( Matptr , i , j ); \
		} \
	} \
}


/**************************************/
/*      Subfunctions Declarations     */
/**************************************/

static MATRIX *VC_GEE_create_matrix(int,int,int);
static MATRIX *VC_GEE_matcopy(MATRIX *);
static MATRIX *VC_GEE_extract_rows(MATRIX *, int, int);
static MATRIX *VC_GEE_matadd(MATRIX *, MATRIX *);
static MATRIX *VC_GEE_matsub(MATRIX *, MATRIX *);
static MATRIX *VC_GEE_matmult(MATRIX *, MATRIX *);
static MATRIX *VC_GEE_transp(MATRIX *);
static MATRIX *VC_GEE_col_1s(int);
static MATRIX *VC_GEE_matabs(MATRIX *);
static MATRIX *VC_GEE_matexp(MATRIX *);
static MATRIX *VC_GEE_px1_times_pxq(MATRIX *, MATRIX *);
static MATRIX *VC_GEE_pxq_divby_px1(MATRIX *, MATRIX *);
static MATRIX *VC_GEE_scalar_times_matrix(double, MATRIX *);
static MATRIX *VC_GEE_ident(int);
static MATRIX *VC_GEE_form_diag(MATRIX *);
static MATRIX *VC_GEE_extract_cols(MATRIX *,int,  int);
static MATRIX *VC_GEE_diag_as_vec(MATRIX *);
static MATRIX *VC_GEE_matsqrt(MATRIX *);
static MATRIX *VC_GEE_mat1over(MATRIX *);

static double VC_GEE_matmax(MATRIX *);
static double VC_GEE_elsum(MATRIX *);
static void VC_GEE_plug(MATRIX *,MATRIX *, int, int),
     VC_GEE_destroy_matrix(MATRIX *)
     ;
static int VC_GEE_split(MATRIX *,MATRIX *,MATRIX *[]),
     VC_GEE_nchanges(MATRIX *)
     ;


/* The following functions are written by Z. Chen */
/* ---------------------------------------------- */

static MATRIX *get_seq1(int, int),
     *get_rep_scalar(int, int),
     *get_rep(MATRIX *, int),
     *get_kronecker(MATRIX *, MATRIX *),
     *get_sum1row(MATRIX *),
     *get_sum2col(MATRIX *),
     *VC_GEE_matexpit(MATRIX *),
     *get_outer(MATRIX *, MATRIX *),
     *get_rbind(MATRIX *, MATRIX *),
     *get_cbind(MATRIX *, MATRIX *),
     *get_cholinv(MATRIX *),
     *matrix_subtract(MATRIX *,MATRIX *),
     *matrix_multiply(MATRIX *, MATRIX *),
     *get_matrix_row(MATRIX *, int)
     ;

static int get_rowindex(int, int, int, int, int, int);
static double get_max_reldif(MATRIX *, MATRIX *),
              get_1_colsum(MATRIX *, int),
              get_1_rowsum(MATRIX *, int);
static void get_mu_i(MATRIX *,MATRIX *,int,int),
     get_lambda_i(MATRIX *,MATRIX *,MATRIX *),
     get_dmu_i_dbetaT(MATRIX *,MATRIX *,int,int),
     get_dlambda_i_dbetaT(MATRIX *,MATRIX *,MATRIX *),
     get_bivar_cumuls_i(MATRIX *,
                        MATRIX *,
                        MATRIX *,
                        MATRIX *,
                        MATRIX *,
                        MATRIX *,
                        MATRIX *,
                        int ,
                        int ,
                        int ,
                        int ),
     get_bivar_marginals_i(MATRIX *,
                           MATRIX *,
                           MATRIX *,
                           MATRIX *,
                           MATRIX *,
                           MATRIX *,
                           MATRIX *,
                           MATRIX *,
                           MATRIX *,
                           MATRIX *,
                           MATRIX *,
                           MATRIX *,
                           MATRIX *,
                           int ,
                           int ,
                           int ,
                           int ),
     matrix_copyto(MATRIX *,MATRIX *),
     row_replace(MATRIX *, int, MATRIX *, int),
     col_replace(MATRIX *, int, MATRIX *, int),
     rows_plug(MATRIX *,int, int, MATRIX *, int),
     cols_plug(MATRIX *,int, int, MATRIX *, int),
     matrix_addto(MATRIX *, MATRIX *),
     get_matsub(MATRIX *,
                MATRIX *, MATRIX *),
     get_matadd(MATRIX *,
                MATRIX *, MATRIX *),
     get_mattransp(MATRIX *, MATRIX *),
     get_matmult(MATRIX *,
                 MATRIX *, MATRIX *),
     cholinv(MATRIX *, MATRIX *),
     add_outer_colvec_to(MATRIX *, MATRIX *),
     scalar_times_matrix(double,  MATRIX *),
     matrix_elem_mult(MATRIX *,MATRIX *),
     matrix_row_mult(MATRIX *,MATRIX *),
     matrix_col_mult(MATRIX *,MATRIX *),
     get_estfun(MATRIX *,
                MATRIX *,
                MATRIX *,
                MATRIX *,
                MATRIX *),
     get_dvd(MATRIX *,MATRIX *,MATRIX *,MATRIX *),
     fisherscoring(double,MATRIX *,MATRIX *,MATRIX *),
     set_zero(MATRIX *),
     get_dpXt_i_dvpT_l(MATRIX *,
                       int ,
                       MATRIX *,
                       MATRIX *,
                       MATRIX *,
                       MATRIX *,
                       MATRIX *,
                       MATRIX *,
                       MATRIX *,
                       int ,
                       int )
     ;

static MATRIX *get_dYtil_i_dgT(MATRIX *,
                               MATRIX *,
                               MATRIX *,
                               MATRIX *,
                               MATRIX *,
                               int ,
                               int ),
     *get_dZtil_i_dgT(MATRIX *,
                      MATRIX *,
                      int ,
                      int ),
     *form_matrix(double *, int, int, int);

void Cmgee2(double *,
            double *,
            double *,
            double *,
            double *,
            double *,
            double *,
            double *,
            int *,
            int *,
            int *,
            int *,
            int *,
            int *,
            double *,
            double *,
            double *,
            double *,
            int *,
            int *,
            double *),
     Cgetmgee2_i(double *,
                 double *,
                 double *,
                 double *,
                 double *,
                 double *,
                 int *,
                 int *,
                 int *,
                 int *,
                 int *,
                 double *,
                 double *,
                 double *,
                 double *,
                 double *),
     Cgetordgee2_i(double *,
                   double *,
                   double *,
                   double *,
                   int *,
                   int *,
                   int *,
                   int *,
                   double *,
                   double *,
                   double *,
                   double *,
                   double *)
     ;

void Cgetmgee2v_i(
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    int *,
    double *,
    double *,
    double *,
    double *,
    double *,
    double *
);
