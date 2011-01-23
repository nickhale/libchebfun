/*******************************************************************************
 * This file is part of libchebfun.
 * Coypright (c) 2010 The Chebfun Team.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 ******************************************************************************/

/* System-wide includes. */
#include <complex.h>

/* Local includes. */


/* Constants. */
#define fun_N0                          2
#define fun_defaultlen                  100


/* Flags. */
#define fun_flag_none                   0
#define fun_flag_init                   1
#define fun_flag_cplx                   2
#define fun_flag_haspts                 4
#define fun_flag_hasvals                8
#define fun_flag_hascoeffs              16


/* Error codes. */
#define fun_err_ok                      0
#define fun_err_null                    -1
#define fun_err_malloc                  -2
#define fun_err_fftw                    -3
#define fun_err_util                    -4
#define fun_err_fx                      -5
#define fun_err_uninit                  -6
#define fun_err_domain                  -7
#define fun_err_nyi                     -8
#define fun_err_lapack                  -9
#define fun_err_gnuplot                 -10
#define fun_err_converge                -11


/* Global, external variables. */
extern int fun_err;
extern const char *fun_err_msg[];


/**
 * This data structure represents a single polynomial approximant.
 */
struct fun {

    /** Bit-mask for flags. */
    unsigned int flags;
    
    /** Nr. of points in this fun. */
    unsigned int n;
    
    /** Size of points, vals and coeffs that were actually allocated. */
    unsigned int size;
    
    /** Scaling factor. */
    double scale;
    
    /** Ends of the fun. */
    double a, b;
    
    /** Pointer to an array of nodes. */
    double *points;
    
    /** Pointer to an array of funciton values at the nodes (real or complex). */
    union {
        double *real;
        double complex *cplx;
        } vals ;
        
    /** Pointer to an array of Chebyshev coefficients (real or complex). */
    union {
        double *real;
        double complex *cplx;
        } coeffs;

    };
    
    
/* Constant struct for virgin funs. */
extern const struct fun FUN_EMPTY;


/* Function declarations. */
int fun_add_const ( struct fun *A , double x , struct fun *B );
int fun_add ( struct fun *a , struct fun *b , struct fun *c );
int _fun_alloc ( struct fun *fun , unsigned int N );
int fun_clean ( struct fun *fun );
int fun_comp_fun ( struct fun *A , struct fun *op , struct fun *B );
int fun_comp_vec ( struct fun *A , int (*fx)( const double * , unsigned int , double * , void * ) , struct fun *B , void *data );
int fun_copy ( struct fun *fun , struct fun *fun2 );
int fun_create_coeffs( struct fun *fun , double *coeffs , double a , double b , unsigned int N );
int fun_create_nonadapt ( struct fun *fun , double (*fx)( double x , void * ) , double a , double b , unsigned int N , void *data );
int fun_create ( struct fun *fun , double (*fx)( double x , void * ) , double a , double b , void *data );
int fun_create_vals( struct fun *fun , double *vals , double a , double b , unsigned int N );
int fun_create_vec ( struct fun *fun , int (*fx)( const double * , unsigned int , double * , void * ) , double a , double b , void *data );
int fun_create_x ( struct fun *fun , double a , double b );
int fun_diff ( struct fun *f, struct fun *fp );
int fun_display ( struct fun *fun , FILE *out );
double fun_err_norm2 ( struct fun *A , struct fun *B );
double fun_err_norm_inf ( struct fun *A , struct fun *B );
double fun_eval ( struct fun *fun , double x );
double fun_eval_bary ( struct fun *fun , double x );
int fun_eval_bary_vec ( struct fun *fun , const double *x , unsigned int m , double *out );
double fun_eval_clenshaw ( struct fun *fun , double x );
int fun_eval_clenshaw_vec ( struct fun *fun , const double *x , unsigned int m , double *out );
int fun_eval_vec ( struct fun *fun , const double *x , unsigned int m , double *out );
int fun_gnuplot ( struct fun *f1 );
int fun_indef_integral ( struct fun *A , struct fun *B );
int fun_init( struct fun *fun , unsigned int N );
double fun_integrate ( struct fun *fun );
int fun_isequal ( struct fun *A , struct fun *B );
int fun_madd ( struct fun *A , double alpha , struct fun *B , double beta , struct fun *C );
int fun_max ( struct fun *fun , double *maxy , double *maxx );
int fun_minandmax ( struct fun *fun , double *miny , double *minx , double *maxy , double *maxx );
int fun_min ( struct fun *fun , double *miny , double *minx );
int fun_mul ( struct fun *A , struct fun *B , struct fun *C );
int fun_newdomain( struct fun *fun , double a , double b);
double fun_norm2 ( struct fun *fun );
double fun_norm_inf ( struct fun *fun );
int fun_points ( struct fun *fun );
int fun_poly ( struct fun *f1 , double *out );
int fun_prolong( struct fun *A , unsigned int N , struct fun *B );
int fun_rescale ( struct fun *fun );
void _fun_rescale ( struct fun *fun );
int fun_restrict ( struct fun *fun , double A , double B , struct fun *funout );
int fun_roots( struct fun *fun , double *roots  );
int fun_roots_unit ( struct fun *fun , double *roots , int sort );
int fun_roots_sort( struct fun *fun , double *roots );
int fun_scale ( struct fun *A , double w , struct fun *B );
int _fun_simplify ( struct fun *fun, double tol );
int fun_simplify ( struct fun *fun, double tol );
int fun_sub ( struct fun *A , struct fun *B , struct fun *C );
int fun_vals ( struct fun *fun );
