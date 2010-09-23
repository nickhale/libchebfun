/*******************************************************************************
 * This file is part of libchebfun.
 * Copyright (c) 2010 The Chebfun Team.
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
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <clapack.h>

/* Local includes. */
#include "errs.h"
#include "chebopts.h"
#include "util.h"
#include "fun.h"


/* Global variables. */
/** The value of the last error code. */
int fun_err = fun_err_ok;
/** Error messages for the (absolute) error codes. */
const char *fun_err_msg[] = {
    "All is well.",
    "An unexpected NULL-pointer was encountered.",
    "A call to malloc failed, probably due to insufficient memory.",
    "Something went wrong in a call to fftw.",
    "Something went wrong in a call to a util_* function.",
    "Something went wrong in a call to a user-supplied function.",
    "The fun has not been initialized.",
    "The funs do not span the same domain or operation requested outside of domain.",
    "Requested operation is not yet implemented." 
    "A call to LAPACK ended in a calamity." };
    
/* Define a macro to store the errors. */
#define error( id )     ( fun_err = errs_register( id , fun_err_msg[-id] , __LINE__ , __FILE__ ) )
    
    
/* Constant struct for virgin funs. */
const struct fun FUN_EMPTY = { 0 , 0 , 0.0 , 0.0 , 0.0 , NULL , { NULL } , { NULL } };


/**
 * @brief Compute the infinity norm of a fun on its domain.
 *
 * @param fun The #fun to take the norm of.
 *
 * @return The infinity norm of the fun or @c NAN if something went wrong
 *      (see #fun_err).
 */

double fun_norm_inf ( struct fun *fun ) {
	
	double miny, minx, maxy, maxx, norm;
    
    /* Bad apples? */
    if ( fun == NULL ) {
        error(fun_err_null);
        return NAN;
        }
    if ( !( fun->flags & fun_flag_init ) ) {
        error(fun_err_uninit);
        return NAN;
        }

    /* Call min/max on this fun. */
	if ( fun_minandmax( fun , &miny , &minx , &maxy , &maxx ) < 0 )
        return NAN;
        
    /* Get the larger of the two. */
	norm = -1.0*miny;
	if ( norm < maxy )
		norm = maxy;

    /* Well that was easy... */
    return norm;
    
	}


/**
 * @brief Compute the maximum of a fun on its domain.
 *
 * @param fun The #fun to be maximised.
 * @param maxy The maximum function value.
 * @param maxx The x co-ordinate for maxy.
 *
 * @return #fun_err_ok or < 0 if an error occured.
 */

int fun_max ( struct fun *fun , double *maxy , double *maxx ) {
	
	double miny, minx;

    
	if ( fun_minandmax ( fun , &miny , &minx , maxy , maxx ) < 0 )
        return error(fun_err);
	
    /* End on a good note. */
    return fun_err_ok;
	}

/**
 * @brief Compute the minimum of a fun on its domain.
 *
 * @param fun The #fun to be minimised.
 * @param miny The minimum function value.
 * @param minx The x co-ordinate for miny.
 *
 * @return #fun_err_ok or < 0 if an error occured.
 */

int fun_min ( struct fun *fun , double *miny , double *minx ) {
	
	double maxy, maxx;

	if ( fun_minandmax ( fun , miny , minx , &maxy , &maxx ) < 0 )
        return error(fun_err);
	
    /* End on a good note. */
    return fun_err_ok;
	}


/**
 * @brief Simultaneously compute the minimum and maxiumum of a fun on its domain.
 *
 * @param fun The #fun to be minimised.
 * @param miny The minimum function value.
 * @param minx The x co-ordinate for miny.
 * @param maxy The maximum function value.
 * @param maxx The x co-ordinate for maxy.
 *
 * @return #fun_err_ok or < 0 if an error occured.
 */

int fun_minandmax ( struct fun *fun , double *miny , double *minx , double *maxy , double *maxx ) {
	
    struct fun fp = FUN_EMPTY;
	double *roots, *vals, valR;
	int j;
	int nroots;
    
    /* Check for the usual bovine excrement. */
    if ( fun == NULL || miny == NULL || minx == NULL || maxx == NULL || maxy == NULL )
        return error(fun_err_null);
    if ( !( fun->flags & fun_flag_init ) )
        return error(fun_err_uninit);

	/* Compute the derivative */
	fun_diff( fun , &fp );

	/* Roots for the derivative */
    if ( ( roots = (double *)alloca( sizeof(double) * fun->n ) ) == NULL )
        return error(fun_err_malloc);
	if ( ( nroots = fun_roots( &fp , roots ) ) < 0 )
        return error(fun_err);

	/* Evaluate the function at these roots */
	if ( ( vals = (double *)alloca( sizeof(double) * nroots + 16 ) ) == NULL )
        return error(fun_err_malloc);
    vals = (double *)( (((size_t)vals) + 15 ) & ~15 );
	if ( fun_eval_clenshaw_vec ( fun , roots , nroots , vals ) < 0 )
        return error(fun_err);
	
	/* Find the maximum */
    /* TODO: call fun_eval_clenshaw_vec... */
	*miny = fun_eval_clenshaw( fun , fun->a ); 
	*minx = fun->a;
	for ( j = 0 ; j < nroots ; j++ ) {
		if ( vals[j] < *miny ) {
			*miny = vals[j];
			*minx = roots[j];
			}
		}
	valR = fun_eval_clenshaw( fun , fun->b );
	if ( valR < *miny ) {
		*miny = valR;
		*minx = fun->b;
		}

	/* Find the maximum */
	*maxy = fun_eval_clenshaw( fun , fun->a ); 
	*maxx = fun->a;
	for ( j = 0 ; j < nroots ; j++ ) {
		if ( vals[j] > *maxy ) {
			*maxy = vals[j];
			*maxx = roots[j];
			}
		}
	valR = fun_eval_clenshaw( fun , fun->b );
	if ( valR > *maxy ) {
		*maxy = valR;
		*maxx = fun->b;
		}

	/* Be free! */
	fun_clean( &fp );
	
    /* End on a good note. */
    return fun_err_ok;
    
	}


/**
 * @brief Copy data from one fun to another
 *
 * @param fun The #fun to be copied.
 * @param fun2 The #fun to copy into.
 *
 * @return #fun_err_ok or < 0 if an error occured.
 * 
 * If @c fun2 has not been initialized, or is too small, it will be
 * re-initialized (see #fun_init).
 */

int fun_copy ( struct fun *fun , struct fun *fun2 ) {

    /* Check inputs. */
    if ( fun == NULL || fun2 == NULL)
        return error(fun_err_null);
        
    /* If fun2 has not been initialized or is too small, re-init. */
    if ( !(fun2->flags & fun_flag_init) || fun2->n < fun->n )
        if ( fun_init( fun2 , fun->n ) < 0 )
            return error(fun_err);
            
    /* Copy the data */
    fun2->flags = fun->flags;
    fun2->n = fun->n;
    fun2->a = fun->a;
    fun2->b = fun->b;

    /* Copy the values. */
    memcpy( fun2->points , fun->points , sizeof(double) * fun->n );
    memcpy( fun2->vals.real , fun->vals.real , sizeof(double) * fun->n );
    memcpy( fun2->coeffs.real , fun->coeffs.real , sizeof(double) * fun->n );
        
    /* End on a good note. */
    return fun_err_ok;
    
    }

/**
 * @brief Restrict a fun to a subinterval of its domain
 *
 * @param fun The #fun which needs restricting.
 * @param A The new left endpoint.
 * @param B The new right endpoint.
 * @param funout The restricted fun.
 *
 * @return #fun_err_ok or < 0 if an error occured.
 *
 * Note that both @c A and @c B must be within the domain of @c fun.
 */

int fun_restrict ( struct fun *fun , double A , double B , struct fun *funout ) {
    
    double iba, *x, *out;
    int j;
    
    /* Check inputs. */
    if ( fun == NULL )
        return error(fun_err_null);
    if ( fun->a > A || fun->b < B )
        return error(fun_err_domain);
        
    /* Allocate some memory for the new evaluation points. */
    if ( ( x = (double *)alloca( sizeof(double) * (fun->n) ) ) == NULL )
        return error(fun_err_malloc);

    /* Note, writing in this form ensures the ends are mapped exactly, even with rounding errors. */
    iba = 1.0/(fun->b-fun->a);
    for (j = 0 ; j < fun->n ; j++ )
        //fun->points[j] = B * (fun->points[j] - fun->a) * iba + A * (fun->b - fun->points[j]) * iba;
        x[j] = B * (fun->points[j] + 1.0) * iba + A * (1.0 - fun->points[j]) * iba;

    if (fun == funout) {
        /* Allocate some storage for evaluations */
        if ( ( out = (double *)alloca( sizeof(double) * (fun->n) ) ) == NULL )
            return error(fun_err_malloc);
            
        /* Pedro: If you use clenshaw below, you don't neet to alloc out! */
    
        /* Get the restricted values. */
        //    fun_eval_clenshaw_vec ( fun , x , fun->n ,  fun->vals.real );
        fun_eval_vec ( fun , x , fun->n ,  out );
    
        /* Assign to funout */
        memcpy( funout->vals.real , out , sizeof(double) * fun->n );

        }
    else {
        /* Start by cleaning out funout */
        if ( !( funout->flags & fun_flag_init ) || ( funout->n < fun->n ) )
            if ( fun_init( funout , fun->n ) != fun_err_ok )
                return error(fun_err);
                
        /* Pass some variables to the new fun */
        funout->flags = fun->flags;

        /* Get the restricted values. */
        //    fun_eval_clenshaw_vec ( fun , x , fun->n ,  fun->vals.real );
        fun_eval_vec ( fun , x , fun->n ,  funout->vals.real );
        }
        
    /* Set the new limits. */
	funout->a = A;
	funout->b = B;

    /* Get the coefficients. */
    util_chebpoly ( funout->vals.real  , funout->n , funout->coeffs.real );

    /* Simplify and re-scale. */
// NEED TO PASS A SENSIBLE TOLERANCE
    if ( fun_simplify( funout , 0.0 ) < 0 )
        return error(fun_err);
    _fun_rescale( funout );
        
    /* End on a good note. */
    return fun_err_ok;
    
    }


/**
 * @brief Simplify a fun (remove trailing coefficients below tolerance).
 *      We simply call util_simplify. 
 *
 * @param fun The fun the be simplified.
 * @param tol Tolerance
 * @return The new length of the fun or < 0 on error.
 */

// AT THE MOMENT TOL IS IGNORED!
int fun_simplify ( struct fun *fun , double tol ) {
    
//	int j;
    unsigned int N;
//	double A05, B05;
    /* struct fun tmpfun = FUN_EMPTY; */

    /* Check inputs. */
    if ( fun == NULL )
        return error(fun_err_null);
    if ( !( fun->flags & fun_flag_init ) )
        return error(fun_err_uninit);

    /* Call util_simplify */
    N = util_simplify ( fun->points , fun->vals.real , fun->coeffs.real , fun->n , fun->b-fun->a , fun->scale , NULL );      

    /* Store in a new fun if simplify was successful. */
    if ( N < fun->n ) {
    
/*        _fun_alloc( &tmpfun , N );
//        util_chebpts( N, tmpfun.points );
        memcpy( tmpfun.points , fun->points , sizeof(double) * N );
        memcpy( tmpfun.vals.real , fun->vals.real , sizeof(double) * N );
        memcpy( tmpfun.coeffs.real , fun->coeffs.real , sizeof(double) * N );
        fun_clean( fun );
        _fun_alloc( fun , N );
        fun_copy ( &tmpfun , fun ); */
        
        /* Set the new size. */
        fun->n = N;
        
        /* Re-create the fun in the already-allocated memory. */
        util_chebpolyval( fun->coeffs.real , N , fun->vals.real );
        util_chebpts( N , fun->points );
//        util_chebptsAB( N , fun->points , fun->a , fun->b);
        /* Scale the points to the correct interval if needed.
		if ( fun->a != -1.0 || fun->b != 1.0 ) {
			A05 = 0.5*fun->a;
			B05 = 0.5*fun->b;
			for ( j = 0 ; j < N ; j++ )
				fun->points[j] = B05 * (1.0 + fun->points[j]) + A05 * (1.0 - fun->points[j]);
        	}
        */
        }

    /* Sweet. */
    return N;
    
    }

/**
 * @brief Compute the roots of a fun.
 *
 * @param fun The #fun to find roots of.
 * @param roots A pointer to an array of double of length at least
 *      equal to the length of @c fun - 1. This is where the roots will be stored.
 *
 * @return The number of roots found or < 0 if an error occured.
 */
 
int fun_roots( struct fun *fun , double *roots ) {
	
	int k, nroots;
	double scl;
    
    /* Check for the usual suspects. */
    if ( fun == NULL || roots == NULL )
        return error(fun_err_null);
    if ( !( fun->flags & fun_flag_init ) )
        return error(fun_err_uninit);
        
    /* Get scl once we're sure fun isn't bogus. */
    scl = 0.5 * (fun->b - fun->a);

    /* Call fun_roots_unit. */
	if ( ( nroots = fun_roots_unit( fun , roots ) ) < 0 )
        return error(fun_err);

    /* Scale the roots to the correct interval if needed. */
	if ( fun->a != -1.0 || fun->b != 1.0 )
	    for ( k = 0 ; k < nroots ; k++ )
		    roots[k] = ( 1.0 - roots[k] ) * scl + fun->a;

	return nroots;
    
	}


/**
 * @brief Compute the roots of a fun assuming it is on the unit interval.
 *
 * @param fun The #fun to find roots of.
 * @param roots A pointer to an array of double of length at least
 *      equal to the length of @c fun - 1. This is where the roots will be stored.
 *
 * @return The number of roots found or < 0 if an error occured.
 */
 
int fun_roots_unit ( struct fun *fun , double *roots ) {

    double *A, cN, c = -0.004849834917525, cL, cR, z, *rr, *ri, *work, tol, a, b;
    unsigned int N, nrootsL, nrootsR, nroots = 0;
    int j, k;
	char job = 'E', compz = 'N';
	int ilo, ihi, ldh, ldz, lwork, ok;
	const struct chebopts *opts;
    struct fun funL = FUN_EMPTY, funR = FUN_EMPTY;
    
    /* Check inputs. */
    if ( fun == NULL )
        return error(fun_err_null);

    if (fun->n < 101) {
    /* Solve the eigenvalue problem. */

        /* Remove small trailing coefficients */
        fun_rescale ( fun );
        N = fun->n-1;
        if ( fabs(fun->coeffs.real[N]) < 1e-14 * fun->scale)
        	while ( fabs(fun->coeffs.real[N]) < 1e-14 * fun->scale && N > 0)
           		N--;
        cN = -0.5 / fun->coeffs.real[N];
            
        /* Initialize the vector A.
           Allocation is done on the stack, so we don't need to worry about
           releasing these if the routine fails anywhere underway.
           Note that A is padded with 16 bytes so that it can be
           shifted and aligned to 16-byte boundaries. */
        if ( ( A = (double *)alloca( sizeof(double) * (N*N) + 16 ) ) == NULL )
            return error(fun_err_malloc);
        /* Shift A so that it is 16-byte aligned. */
        A = (double *)( (((size_t)A) + 15 ) & ~15 );
		/* Zero out A */
        bzero( A , sizeof(double) * ( N * N ) );

        /* Also initialize and shift the output and work vectors */
        if ( ( rr = (double *)alloca( sizeof(double) * (N) + 16 ) ) == NULL || 
			 ( ri = (double *)alloca( sizeof(double) * (N) + 16 ) ) == NULL ||
		     ( work = (double *)alloca( sizeof(double) * (N) + 16 ) ) == NULL )
            return error(fun_err_malloc);
        rr = (double *)( (((size_t)rr) + 15 ) & ~15 );
        ri = (double *)( (((size_t)ri) + 15 ) & ~15 );
        work = (double *)( (((size_t)work) + 15 ) & ~15 );
   
        /* Assign the coefficients to the first row */
        for (j = 0 ; j < N ; j++)
            A[N*(N-1-j)] = cN * fun->coeffs.real[j];
        A[N] += 0.5;
    
        /* Assign the super- and sub-diagonal */
        A[1] = 0.5;
        A[N+2] = 0.5;
        for (j = 2*N+1 ; j < N*(N-1) ; j += N+1) {
            A[j] = 0.5;
            A[j+2] = 0.5;
            }
		A[N*(N-1)-1] = 1.0;
        A[N*N-2] = 0.5;

        /* Call to LAPACK to solve eigenvalue problem. */
		ilo = 1; ihi = N; ldh = N; ldz = 1; lwork = N;
		dhseqr_( &job, &compz, &N , &ilo, &ihi, A, &ldh, rr, ri, &z, &ldz, work, &lwork, &ok);
		if ( ok < 0 )
			return error(fun_err_lapack);

		opts = &chebopts_default;
// THIS SHOULD INVOLVE A DECREASING HORZONTAL SCALE AS IN MATLAB/CHEBFUN
		tol = 100.0 * opts->eps;
        
        /* Count the number of valid roots and store them. */
		for (j = ok ; j < N ; j++)
			if (fabs(ri[j]) < tol && rr[j] >= -1.0-tol && rr[j] <= 1.0+tol)
                roots[nroots++] = rr[j];
                
		}
    else {
    /* Recurse. */

        /* Assume unit interval */
        if (fun->points[0] != -1.0 || fun->points[fun->n-1] != 1.0) {
			a = fun->a;
			b = fun->b;
            fun->a = -1.0;
            fun->b = 1.0;
            }

        /* Restrict to subintervals (with magic ChebConstant c) */
        fun_restrict( fun , -1.0, c , &funL );
        fun_restrict( fun , c , 1.0 , &funR );

		/* Reset the interval */
        fun->a = a;
        fun->b = b;

        /* Find roots recursively */
        nrootsR = fun_roots_unit( &funR ,  roots ); 
        nrootsL = fun_roots_unit( &funL ,  &( roots[ nrootsR ] ) );
        nroots = nrootsL + nrootsR;

		/* Adjust roots in output vector */
        k = 0;
        cR = 0.5 * (1.0 - c);
        for (j = 0 ; j < nrootsR ; j++, k++ )
            roots[k] = c + ( roots[j] + 1.0 ) * cR;
        cL = 0.5 * (1.0 + c); 
        for (j = 0 ; j < nrootsL ; j++, k++ )
            roots[k] = -1.0 + ( roots[k] + 1.0 ) * cL;

        fun_clean( &funL );
        fun_clean( &funR ); 
        
    }
        
    /* To the pub!. */
    return nroots;
    
    }


/**
 * @brief Compute the norm of a real valued function
 *
 * @param fun The #fun which needs re-scaling.
 *
 * @return #fun_err_ok or < 0 if an error occured.
 */
 
double fun_norm2_backup ( struct fun *fun ) {

    double norm;
    struct fun fun2 = FUN_EMPTY;

    /* Check for bad input. */
    if ( fun == NULL )
        return error(fun_err_null);
    if ( !(fun->flags & fun_flag_init ) )
        return error(fun_err_uninit);
        
    /* Compute the function ^2. */
    fun_mul( fun , fun , &fun2 );
    
    /* Compute the inner product. */
    norm = sqrt(fun_integrate( &fun2 ));
   
    /* Free up the created fun */
    fun_clean( &fun2);

    /* End on a good note. */
    return norm;
    
    }


/**
 * @brief Compute the norm of a real valued function
 *
 * @param fun The #fun which needs re-scaling.
 *
 * @return The 2-norm of @c fun or @c NAN if an error occured, in which
 *      case the error code is stored in #fun_err.
 */
 
double fun_norm2 ( struct fun *fun ) {

    double *c, w, val;
    int j, k;

    /* Check for bad input. */
    if ( fun == NULL ) {
        error(fun_err_null);
        return NAN;
        }
    if ( !(fun->flags & fun_flag_init ) ) {
        error(fun_err_uninit);
        return NAN;
        }
        
    /* Allocate c on the stack. */
    if ( ( c = (double *)alloca( sizeof(double) * ( 2 * fun->n - 1 ) ) ) == NULL ) {
        error(fun_err_malloc);
        return NAN;
        }
    bzero( c , sizeof(double) * ( 2 * fun->n - 1 ) );
        
    /* Fill c with the coefficients of fun.^2. */
    for ( j = 0 ; j < fun->n ; j++ ) {
        w = 0.5 * fun->coeffs.real[j] * fun->coeffs.real[j];
        c[2*j] += w;
        c[0] += w;
        for ( k = j+2 ; k < fun->n ; k += 2 ) {
            w = fun->coeffs.real[j] * fun->coeffs.real[k];
            c[j+k] += w;
            c[abs(j-k)] += w;
            }
        }
        
    /* Extract the integral from the coefficients. */
    val = 2.0 * c[0];
    for ( k = 2 ; k < 2*fun->n-1 ; k += 2 )
        val += 2.0 * c[k] / (1 - k*k);
    val = sqrt( val * 0.5 * (fun->b - fun->a) );
        
    /* End on a good note. */
    return val;
    
    }


/**
 * @brief Display a load of bumf about a fun. (Mostly for testing).
 *
 * @param fun The #fun to be displayed.
 *
 * @return #fun_err_ok or < 0 if an error occured.
 */
 
int fun_display ( struct fun *fun ) {

    int k;
    double xk, a05, b05;

    /* Check for bad input. */
    if ( fun == NULL )
        return error(fun_err_null);
    if ( !(fun->flags & fun_flag_init ) )
        return error(fun_err_uninit);
    
    a05 = 0.5*fun->a; 
    b05 = 0.5*fun->b; 
    /* Loop over the entries */
    for ( k = 0 ; k < fun->n ; k++ ) {
        xk = b05 * (fun->points[k] + 1.0) + a05 * (1.0 - fun->points[k]);
        printf("fun.points[%i]=%16.16e \tfun.vals[%i]=%16.16e \tfun.coeffs[%i]=%16.16e\n",
            k, xk, k, fun->vals.real[k], k, fun->coeffs.real[k]);
        }
    printf("\n");
    
    /* Huzzah! */
    return fun_err_ok;
    
    }


/**
 * @brief Re-computes and stores the scale of the given fun.
 *
 * @param fun The #fun which needs re-scaling.
 *
 * @return #fun_err_ok or < 0 if an error occured.
 */
 
int fun_rescale ( struct fun *fun ) {

    double scale = 0.0;
    int k;

    /* Check for bad input. */
    if ( fun == NULL )
        return error(fun_err_null);
    if ( !(fun->flags & fun_flag_init ) )
        return error(fun_err_uninit);
        
    /* Compute the scale. */
    for ( k = 0 ; k < fun->n ; k++ )
        if ( scale < fabs( fun->vals.real[k] ) )
            scale = fabs( fun->vals.real[k] );
    fun->scale = scale;
    
    /* End on a good note. */
    return fun_err_ok;
    
    }


/**
 * @brief Re-computes and stores the scale of the given fun.
 *
 * @param fun The #fun which needs re-scaling.
 *
 * Same as #fun_rescale, yet without error checking. Should be
 * used internally only!
 */
 
inline void _fun_rescale ( struct fun *fun ) {

    double scale = 0.0;
    int k;

    /* Compute the scale. */
    for ( k = 0 ; k < fun->n ; k++ )
        if ( scale < fabs( fun->vals.real[k] ) )
            scale = fabs( fun->vals.real[k] );
    fun->scale = scale;
    
    }


/**
 * @brief Allocate space for a fun to be filled with values.
 *
 * @param fun The #fun structure to be initialized.
 * @param N The size of the #fun.
 *
 * @return #fun_err_ok or < 0 on error.
 *
 */

int fun_init ( struct fun *fun , unsigned int N ) {

    /* Check for null and stuff. */
    if ( fun == NULL )
        return error(fun_err_null);
        
    /* Clean the fun, just to be safe. */
    if ( fun_clean( fun ) < 0 )
        return error(fun_err);
        
    /* Allocate the data inside the fun to the correct size. */
    _fun_alloc( fun , N);

    /* Set entries to zero. Is this needed? */
    bzero( fun->vals.real , sizeof(double) * N );
    bzero( fun->coeffs.real , sizeof(double) * N );
    
    /* Get the Chebyshev nodes. */
    if ( util_chebpts( N , fun->points ) < 0 )
        return error(fun_err_util);
        
    /* Jolly good! */
    return fun_err_ok;
    
    }

/**
 * @brief Allocate space for a fun to be filled with values.
 *     This differs slightly from fun_init in that it doesn't
 *     zero out entries or assign points.
 *
 * @param fun The #fun structure to be initialized.
 * @param N The size of the #fun.
 *
 * @return #fun_err_ok or < 0 on error.
 *
 */

int _fun_alloc ( struct fun *fun , unsigned int N ) {

    /* Allocate the data inside the fun to the correct size. */
    if ( posix_memalign( (void **)&(fun->vals.real) , 16 , sizeof(double) * N ) != 0 ||
         posix_memalign( (void **)&(fun->coeffs.real) , 16 , sizeof(double) * N ) != 0 ||
         ( fun->points = (double *)malloc( sizeof(double) * N ) ) == NULL )
        return error(fun_err_malloc);
    
    /* Write the data to the fun. */
    fun->n = N;
    fun->flags = fun_flag_init;
        
    /* Jolly good! */
    return fun_err_ok;
    
    }


/**
 * @brief Create a fun from a set of function values at Chebyshev points.
 *
 * @param fun The #fun structure to be created.
 * @param vals The values at the Chebyshev grid.
 * @param a The left hand end.
 * @param b The right hand end.
 * @param N The size of the #fun.
 *
 * @return #fun_err_ok or < 0 on error.
 *
 */

int fun_create_vals ( struct fun *fun , double *vals , double a , double b , unsigned int N ) {
    
    /* Check for null and stuff. */
    if ( fun == NULL || vals == NULL )
        return error(fun_err_null);
        
    /* Initialise a fun of the correct size. */
    fun_init( fun, N );

    /* Write the data to the fun. */
    memcpy( fun->vals.real , vals , sizeof(double) * N );
    fun->a = a;
    fun->b = b;
    
    /* Compute the coeffs from the values. */
    if ( util_chebpoly( vals , N , fun->coeffs.real ) < 0 )
        return error(fun_err_util);
        
    /* Update the scale */
    _fun_rescale( fun );
   
    /* Jolly good! */
    return fun_err_ok;
    
    }
    

/**
 * @brief Create a fun from a series of Chebyshev coefficients.
 *
 * @param fun The #fun structure to be created.
 * @param coeffs The Chebyshev coefficients.
 * @param a The left hand end.
 * @param b The right hand end.
 * @param N The size of the #fun.
 *
 * @return #fun_err_ok or < 0 on error.
 *
 */

int fun_create_coeffs ( struct fun *fun , double *coeffs , double a , double b , unsigned int N ) {
    
    /* Check for null and stuff. */
    if ( fun == NULL || coeffs == NULL )
        return error(fun_err_null);
        
    /* Initialise a fun of the correct size. */
    fun_init( fun, N );
   
    /* Write the data to the fun. */
    memcpy( fun->coeffs.real , coeffs , sizeof(double) * N );
    fun->a = a;
    fun->b = b;
    
    /* Compute the values from the coeffs. */
    if ( util_chebpolyval( coeffs , N , fun->vals.real ) < 0 )
        return error(fun_err_util);
        
    /* Update the scale */
    _fun_rescale( fun );
        
    /* Jolly good! */
    return fun_err_ok;
    
    } 

    
/**
 * @brief Differentiate a fun.
 *
 * @param f The #fun to be differentiates.
 * @param fp The differentiated #fun.
 * @return #fun_err_ok or < 0 if an error occurs.
 */

int fun_diff ( struct fun *f, struct fun *fp ) {
    int k, n;

    /* Check for null and stuff. */
    if ( f == NULL || fp == NULL ) {
        return error(fun_err_null);
        }
    if ( !( f->flags & fun_flag_init ) ) {
        return error(fun_err_uninit);
        }

    /* If fp is uninit or too small, init fp. */
    if ( !( fp->flags & fun_flag_init ) || fp->n < f->n - 1 )
        if ( fun_init( fp , f->n - 1 ) < 0 )
            return error(fun_err);

    /* Assign the end points */
    fp->a = f->a;
    fp->b = f->b;
  
    /* The case where f is a constant is trivial */
    if ( f->n == 1 ) {
        fp->vals.real[0] = 0.0;
        fp->coeffs.real[0] = 0.0;
        fp->scale = 0.0;
        fp->points[0] = 0.5*(f->b-f->a);
        return fun_err_ok;
        } 

    n = f->n;
    /* Apply the recurrence for the new coefficents */
    fp->coeffs.real[n-2] = (double)(2*(n-1)) * f->coeffs.real[n-1];
    fp->coeffs.real[n-3] = (double)(2*(n-2)) * f->coeffs.real[n-2];
    for ( k = n-4 ; k >= 0 ; k-- )
        fp->coeffs.real[k] = (double)(2*(k+1)) * f->coeffs.real[k+1] + fp->coeffs.real[k+2];
    fp->coeffs.real[0] = 0.5*fp->coeffs.real[0];

    /* Extract the vals from these coeffs. */
    if ( util_chebpolyval( fp->coeffs.real , fp->n , fp->vals.real ) < 0 )
         return fun_err_util;

    /* Scale values according to the domain */
    fun_scale( fp, 2.0/(fp->b-fp->a) , fp );

    fp->scale = 0.0;
    for ( k = 0 ; k < fp->n ; k++ )
        if ( fabs(fp->vals.real[k]) > fp->scale )
            fp->scale = fabs(fp->vals.real[k]);

    /* All is well. */
    return fun_err_ok;

    }    

     
     
/**
 * @brief Multiply two funs.
 *
 * @param A Pointer to a #fun.
 * @param B Pointer to a #fun.
 * @param C Pointer to a #fun in which to store the result of
 *      @c A * @c B. @c C can be equal to @c A or @c B. In any
 *      case, @c C should be either initialized or clean (see #fun_clean).
 * @return #fun_err_ok or < 0 if an error occurs.
 */
 
int fun_mul ( struct fun *A , struct fun *B , struct fun *C ) {

    int j, k, N;
    double w, *temp;

    /* Check for the usual nonsense... */
    if ( A == NULL || B == NULL || C == NULL )
        return error(fun_err_null);
    if ( !( A->flags & fun_flag_init ) || !( B->flags & fun_flag_init ) )
        return error(fun_err_uninit);
        
    /* Do both funs span the same domain? */
    if ( A->a != B->a || A->b != B->b )
        return error(fun_err_domain);
        
    /* Get the size of the new fun. */
    N = A->n + B->n - 1;
        
    /* Is C neither A or B? */
    if ( C != A && C != B ) {
    
        /* Clean-up C if needed. */
        if ( !( C->flags & fun_flag_init ) || C->n < N )
            if ( fun_init( C , N ) < 0 )
                return error(fun_err);
            
        /* Set some values. */
        C->a = A->a; C->b = A->b;
        
        /* Get the pointer to the coeffs. */
        temp = C->coeffs.real;
        
        }
        
    else {
            
        /* Allocate a buffer for the new coefficients. */
        if ( posix_memalign( (void **)&(temp) , 16 , sizeof(double) * N ) != 0 )
            return error(fun_err_malloc);
        bzero( temp , sizeof(double) * (A->n + B->n - 1) );
        
        /* Allocate the new values. */
        if ( C->vals.real != NULL )
            free( C->vals.real );
        if ( posix_memalign( (void **)&(C->vals.real) , 16 , sizeof(double) * N ) != 0 )
            return error(fun_err_malloc);
            
        /* Create the new points. */
        if ( C->points != NULL )
            free( C->points );
        if ( ( C->points = util_chebpts_alloc( N ) ) == NULL )
            return error(fun_err_malloc);
                
        }

    /* Collect the coefficients from A and B. */
    for ( j = 0 ; j < A->n ; j++ )
        for ( k = 0 ; k < B->n ; k++ ) {
            w = 0.5 * A->coeffs.real[j] * B->coeffs.real[k];
            temp[j+k] += w;
            temp[abs(j-k)] += w;
            }

    /* Set the length of C. */
    C->n = N;
    
    /* Set the coefficients. */
    if ( C->coeffs.real != temp )
        free( C->coeffs.real );
    C->coeffs.real = temp;
    
    /* Compute the new values. */
    if ( util_chebpolyval( C->coeffs.real , C->n , C->vals.real ) < 0 )
        return error(fun_err_util);
                                
    /* Get the scale of the new fun. */
    _fun_rescale( C );
            
    /* If we haven't bombed-out yet, then all is well... */
    return fun_err_ok;
            
    }


/**
 * @brief Integrate a fun over its domain.
 *
 * @param fun The #fun to be integrated.
 * @return the value of the intagral or @c NAN if an error was encountered.
 */

double fun_integrate ( struct fun *fun ) {

    double val;
    int k;

    /* Check for null and stuff. */
    if ( fun == NULL ) {
        error(fun_err_null);
        return NAN;
        }
    if ( !( fun->flags & fun_flag_init ) ) {
        error(fun_err_uninit);
        return NAN;
        }

    /* Sum the coefficients with the correct weights. */
    val = 2.0 * fun->coeffs.real[0];
    for ( k = 2 ; k < fun->n ; k += 2 )
        val += 2.0 * fun->coeffs.real[k] / (1 - k*k);
    val *= 0.5 * (fun->b - fun->a);

    /* All is well. */
    return val;

    }
     
     
/**
 * @brief Scale the given fun by a scalar value.
 *
 * @param A Pointer to a #fun.
 * @param w Weight with which to scale the #fun @c A.
 * @param B Pointer to a #fun in which to store the result. This can
 *      be the same value as @c A.
 * @return #fun_err_ok or < 0 if an error occurs.
 */

int fun_scale ( struct fun *A , double w , struct fun *B ) {

    int k;

    /* Check for the usual problems... */
    if ( A == NULL || B == NULL )
        return error(fun_err_null);
    if ( !( A->flags & fun_flag_init ) )
        return error(fun_err_uninit);
        
    /* Are A and B one and the same? */
    if ( A == B ) {
    
        /* Run through the coeffs and vals and scale them all. */
        for ( k = 0 ; k < A->n ; k++ ) {
            A->vals.real[k] *= w;
            A->coeffs.real[k] *= w;
            }
            
        /* Adjust the scale of A. */
        A->scale *= fabs(w);
        
        }
        
    /* A and B are not the same fun. */
    else {
    
        /* Start by cleaning out B */
        if ( fun_clean( B ) != fun_err_ok )
            return error(fun_err);
            
        /* Copy the values from A to B */
        B->flags = A->flags;
        B->a = A->a; B->b = A->b;
        B->n = A->n;
        B->scale = fabs(w) * A->scale;
        
        /* Allocate memory in B */
        if ( posix_memalign( (void **)&(B->vals.real) , 16 , sizeof(double) * B->n ) != 0 ||
             posix_memalign( (void **)&(B->coeffs.real) , 16 , sizeof(double) * B->n ) != 0 ||
             ( B->points = (double *)malloc( sizeof(double) * B->n ) ) == NULL )
            return error(fun_err_malloc);
            
        /* Copy the scaled values from A. */
        for ( k = 0 ; k < B->n ; k++ ) {
            B->points[k] = A->points[k];
            B->vals.real[k] = w * A->vals.real[k];
            B->coeffs.real[k] = w * A->coeffs.real[k];
            }
    
        }
        
    /* We made it! */
    return fun_err_ok;
        
    }


/**
 * @brief Add two funs.
 * 
 * @param A Pointer to a #fun.
 * @param B Pointer to a #fun.
 * @param C Pointer to a #fun in which to store the result of
 *      @c A + @c B. @c C can be equal to @c A or @c B. In any
 *      case, @c C should be either initialized or clean (see #fun_clean).
 * @return #fun_err_ok or < 0 if an error occurs.
 *
 * This function is just a wrapper for #fun_madd.
 */
 
int fun_add ( struct fun *A , struct fun *B , struct fun *C ) {

    return fun_madd( A , 1.0 , B , 1.0 , C );

    }


/**
 * @brief Add two weighted funs.
 * 
 * @param A Pointer to a #fun.
 * @param alpha Scalar with which to scale @c A.
 * @param B Pointer to a #fun.
 * @param beta Scalar with which to scale @c B.
 * @param C Pointer to a #fun in which to store the result of
 *      @c alpha*A + @c beta*B. @c C can be equal to @c A or @c B. In any
 *      case, @c C should be either initialized or clean (see #fun_clean).
 * @return #fun_err_ok or < 0 if an error occurs.
 *
 * The inputs @c A and @c B must be defined on the same domain.
 * If @c A and @c B have the same length/degree, their coefficients and
 * values are added into @c C.
 * If they are of different lengths, the coefficients of the shorter
 * are added to those of the longer and the values are re-computed using
 * an inverse DCT (see #util_chebpolyval).
 */
 
int fun_madd ( struct fun *A , double alpha , struct fun *B , double beta , struct fun *C ) {

    int k;
    struct fun *a, *b;
    double wa, wb, *temp;

    /* Check for the usual problems... */
    if ( A == NULL || B == NULL || C == NULL )
        return error(fun_err_null);
    if ( !( A->flags & fun_flag_init ) || !( B->flags & fun_flag_init ) )
        return error(fun_err_uninit);
        
    /* Do both funs span the same domain? */
    if ( A->a != B->a || A->b != B->b )
        return error(fun_err_domain);
        
    /* Re-name a and b such that |a| >= |b|. */
    if ( A->n >= B->n ) {
        a = A; b = B; wa = alpha; wb = beta;
        }
    else {
        a = B; b = A; wa = beta; wb = alpha;
        }
        
    /* Is C neither a nor b? */
    if ( C != a && C != b ) {
    
        /* Init C if needed */
        if ( !( C->flags & fun_flag_init ) || C->n < a->n )
            if ( fun_init( C , a->n ) < 0 )
                return error(fun_err);
            
        /* Set the domain. */
        C->a = a->a; C->b = a->b;
            
        /* Copy the points from a. */
        memcpy( C->points , a->points , sizeof(double) * C->n );

        /* Merge the coefficients from a and b. */
        for ( k = 0 ; k < b->n ; k++ )
            C->coeffs.real[k] = wa * a->coeffs.real[k] + wb * b->coeffs.real[k];
        for ( k=k ; k < a->n ; k++ )
            C->coeffs.real[k] = wa * a->coeffs.real[k];
                
        /* Are a and b the same length? */
        if ( a->n == b->n )
        
            /* Merge the vals from a and b. */
            for ( k = 0 ; k < C->n ; k++ )
                C->vals.real[k] = wa * a->vals.real[k] + wb * b->vals.real[k];
                
        /* Different lengths in a and b. */
        else
        
            /* Extract the vals from these coeffs. */
            if ( util_chebpolyval( C->coeffs.real , C->n , C->vals.real ) < 0 )
                return fun_err_util;
                
        /* As of here, this fun is initialized. */
        C->flags |= fun_flag_init;
            
        }
        
    /* Is C == a? */
    else if ( C == a ) {
    
        /* Are a and b the same length? */
        if ( a->n == b->n ) {
        
            /* Merge the vals and coeffs from b into a. */
            for ( k = 0 ; k < C->n ; k++ ) {
                C->vals.real[k] = wa * C->vals.real[k] + wb * b->vals.real[k];
                C->coeffs.real[k] = wa * C->coeffs.real[k] + wb * b->coeffs.real[k];
                }
                
            }
            
        /* Different lengths in a and b. */
        else {
    
            /* Merge the coefficients from b into a. */
            for ( k = 0 ; k < b->n ; k++ )
                C->coeffs.real[k] = wa * C->coeffs.real[k] + wb * b->coeffs.real[k];
            for ( k = k ; k < a->n ; k++ )
                C->coeffs.real[k] = wa * C->coeffs.real[k];

            /* Extract the vals from these coeffs. */
            if ( util_chebpolyval( C->coeffs.real , C->n , C->vals.real ) < 0 )
                return fun_err_util;
                
            }
            
        }
        
    /* Otherwise, C == b. */
    else {
    
        /* Are a and b the same length? */
        if ( a->n == b->n ) {
        
            /* Merge the vals and coeffs from a into b. */
            for ( k = 0 ; k < C->n ; k++ ) {
                C->vals.real[k] = wa * a->vals.real[k] + wb * C->vals.real[k];
                C->coeffs.real[k] = wa * a->coeffs.real[k] + wb* C->coeffs.real[k];
                }
                
            }
            
        /* Different lengths in a and b. */
        else {
        
            /* Free the old points, re-allocate them and copy them from a. */
            free( C->points );
            if ( ( C->points = (double *)malloc( sizeof(double) * a->n ) ) == NULL )
                return error(fun_err_malloc);
            for ( k = 0 ; k < a->n ; k++ )
                C->points[k] = a->points[k];
                
            /* Allocate memory for the new coeffs and merge from a and b. */
            if ( posix_memalign( (void **)&(temp) , 16 , sizeof(double) * a->n ) != 0 )
                return error(fun_err_malloc);
            for ( k = 0 ; k < b->n ; k++ )
                temp[k] = wa * a->coeffs.real[k] + wb * b->coeffs.real[k];
            for ( k = k ; k < a->n ; k++ )
                temp[k] = wa * a->coeffs.real[k];
                
            /* Free and replace the coeffs from b. */
            free( C->coeffs.real );
            C->coeffs.real = temp;
            
            /* Set the new length. */
            C->n = a->n;
    
            /* Extract the vals from the merged coeffs. */
            free( C->vals.real );
            if ( ( C->vals.real = util_chebpolyval_alloc( C->coeffs.real , C->n ) ) == NULL )
                return error(fun_err_util);
                
            }
    
        }

    /* Get the scale of the new fun. */
    _fun_rescale( C );
            
    /* If we haven't bombed-out yet, then all is well... */
    return fun_err_ok;
            
    }
    
    
/**
 * @brief Evaluates the given fun at the given node @c x using the
 *      Clenshaw algorithm.
 *
 * @param fun The #fun to be evaluated.
 * @param x The position at which to evaluate @c fun.
 * @return The #fun evaluated at @c x or @c NaN if an error occured
 *      (see #fun_err).
 *
 * @sa fun_eval, fun_eval_vec, fun_eval_clenshaw_vec
 */
 
double fun_eval_clenshaw ( struct fun *fun , double x ) {

    int k;
    double yn = 0.0, ynp1 = 0.0, ynp2, m, ih;

    /* Check for nonsense. */
    if ( fun == NULL ) {
        error(fun_err_null);
        return NAN;
        }
    if ( !( fun->flags & fun_flag_init ) ) {
        error(fun_err_uninit);
        return NAN;
        }
        
    /* Check for the simplest cases. */
    if ( fun->n == 0 )
        return 0.0;
    else if ( fun->n == 1 )
        return fun->coeffs.real[0];
        
    /* Map x to the interval [-1,1]. */
    m = ( fun->a + fun->b ) * 0.5;
    ih = 2.0 / ( fun->b - fun->a );
    x = ( x - m ) * ih;
        
    /* Evaluate the recurrence. */
    for ( k = fun->n-1 ; k >= 0 ; k-- ) {
        ynp2 = ynp1; ynp1 = yn;
        yn = fun->coeffs.real[k] + 2 * x * ynp1 - ynp2;
        }
        
    /* Return the result. */
    return yn - x * ynp1;
        
    }
    
    
/**
 * @brief Evaluates the given fun at the given vector of nodes @c x using
 *      the Clenshaw algorithm.
 *
 * @param fun The #fun to be evaluated.
 * @param x A pointer to an array of double values containing the positions
 *      at which to evaluate @c fun.
 * @param m The number of positions in @c x.
 * @param out A pointer to an array of doubles of length @c m in which the
 *      values of @c fun at @c x will be stored.
 * @return #fun_err_ok or < 0 if an error occured.
 *
 * @sa fun_eval_clenshaw, fun_eval, fun_eval_vec
 */
 
 
int fun_eval_clenshaw_vec ( struct fun *fun , double *x , unsigned int m , double *out ) {

    int j, k;
    double yn, ynp1 = 0.0, ynp2, xj, mi, ih;

    /* Check for nonsense. */
    if ( fun == NULL || x == NULL || out == NULL )
        return error(fun_err_null);
    if ( !( fun->flags & fun_flag_init ) )
        return error(fun_err_uninit);
        
    /* Check for the simplest cases. */
    if ( fun->n < 2 ) {
        if ( fun->n == 0 )
            for ( k = 0 ; k < m ; k++ )
                out[k] = 0.0;
        else
            for ( k = 0 ; k < m ; k++ )
                out[k] = fun->coeffs.real[0];
        return fun_err_ok;
        }
        
    /* Get the centre and width of the interval. */
    mi = ( fun->a + fun->b ) * 0.5;
    ih = 2.0 / ( fun->b - fun->a );
        
    /* Loop over the input values. */
    for ( j = 0 ; j < m ; j++ ) {
        
        /* Map x[j] back to the [-1,1] interval. */
        xj = ( x[j] - mi ) * ih;
    
        /* Init the recurrence. */
        ynp1 = 0.0; yn = 0.0;

        /* Evaluate the recurrence. */
        for ( k = fun->n-1 ; k >= 0 ; k-- ) {
            ynp2 = ynp1; ynp1 = yn;
            yn = fun->coeffs.real[k] + 2 * xj * ynp1 - ynp2;
            }

        /* store the result. */
        out[j] = yn - xj * ynp1;
        
        }
        
    /* If all went well... */
    return fun_err_ok;
        
    }
    
    
/**
 * @brief Evaluates the given fun at the given node @c x using barycentric
 *      interpolation.
 *
 * @param fun The #fun to be evaluated.
 * @param x The position at which to evaluate @c fun.
 * @return The #fun evaluated at @c x or @c NaN if an error occured
 *      (see #fun_err).
 *
 * @sa fun_eval_clenshaw, fun_eval, fun_eval_clenshaw_vec
 */
 
double fun_eval ( struct fun *fun , double x ) {

    int k;
    double w, u = 0.0, v = 0.0, m, ih;
    
    /* Check for nonsense. */
    if ( fun == NULL ) {
        error(fun_err_null);
        return NAN;
        }
    if ( !( fun->flags & fun_flag_init ) ) {
        error(fun_err_uninit);
        return NAN;
        }
        
    /* Check for the simplest cases. */
    if ( fun->n == 0 )
        return 0.0;
    else if ( fun->n == 1 )
        return fun->vals.real[0];
        
    /* Map x to the interval [-1,1]. */
    m = ( fun->a + fun->b ) * 0.5;
    ih = 2.0 / ( fun->b - fun->a );
    x = ( x - m ) * ih;
        
    /* Do the barycentric form in place. */
    /* Do the first node separately due to the half-weight. */
    if ( x == fun->points[0] )
        return fun->vals.real[0];
    w = 0.5 / ( x - fun->points[0] );
    u = fun->vals.real[0] * w;
    v = w;

    /* Do the interior nodes. */
    for ( k = 1 ; k < fun->n-1 ; k++ ) {
        if ( x == fun->points[k] )
            return fun->vals.real[k];
        w = (double)( 1 - 2 * ( k & 1 ) ) / ( x - fun->points[k] );
        u += fun->vals.real[k] * w;
        v += w;
        }
        
    /* Do the last node separately due to the half-weight. */
    if ( x == fun->points[k] )
        return fun->vals.real[k];
    w = ( 0.5 - ( k & 1 ) ) / ( x - fun->points[k] );
    u += fun->vals.real[k] * w;
    v += w;

    /* Return the fraction. */
    return u / v;

    }


/**
 * @brief Evaluates the given fun at the given vector of nodes @c x using
 *      barycentric interpolation.
 *
 * @param fun The #fun to be evaluated.
 * @param x A pointer to an array of double values containing the positions
 *      at which to evaluate @c fun.
 * @param m The number of positions in @c x.
 * @param out A pointer to an array of doubles of length @c m in which the
 *      values of @c fun at @c x will be stored.
 * @return #fun_err_ok or < 0 if an error occured.
 *
 * @sa fun_eval_clenshaw, fun_eval_vec, fun_eval_clenshaw_vec
 */
 
int fun_eval_vec ( struct fun *fun , double *x , unsigned int m , double *out ) {

    int j, k;
    double w, u, v, xj, mi, ih;
    
    /* Check for nonsense. */
    if ( fun == NULL || x == NULL || out == NULL )
        return error(fun_err_null);
    if ( !( fun->flags & fun_flag_init ) )
        return error(fun_err_uninit);
        
    /* Check for the simplest cases. */
    if ( fun->n < 2 ) {
        if ( fun->n == 0 )
            for ( k = 0 ; k < m ; k++ )
                out[k] = 0.0;
        else
            for ( k = 0 ; k < m ; k++ )
                out[k] = fun->vals.real[0];
        return fun_err_ok;
        }
        
    /* Get the centre and width of the interval. */
    mi = ( fun->a + fun->b ) * 0.5;
    ih = 2.0 / ( fun->b - fun->a );
        
    /* For each element of x... */
    for ( j = 0 ; j < m ; j++ ) {
    
        /* Init u and w. */
        u = 0.0; w = 0.0;
        
        /* Map x[j] back to the [-1,1] interval. */
        xj = ( x[j] - mi ) * ih;
    
        /* Do the barycentric form in place. */
        /* Do the first node separately due to the half-weight. */
        if ( xj == fun->points[0] ) {
            out[j] = fun->vals.real[0];
            continue;
            }
        w = 0.5 / ( xj - fun->points[0] );
        u = fun->vals.real[0] * w;
        v = w;
        
        /* Do the interior nodes. */
        for ( k = 1 ; k < fun->n-1 ; k++ ) {
            if ( x[j] == fun->points[k] )
                break;
            w = (double)( 1 - 2 * ( k & 1 ) ) / ( xj - fun->points[k] );
            u += fun->vals.real[k] * w;
            v += w;
            }
            
        /* Did the loop terminate ok? */
        if ( k == fun->n-1 ) {
        
            /* Get the last node/weight. */
            if ( xj == fun->points[k] ) {
                out[j] = fun->vals.real[k];
                continue;
                }
            w = ( 0.5 - ( k & 1 ) ) / ( xj - fun->points[k] );
            u += fun->vals.real[k] * w;
            v += w;
            
            /* Store the result. */
            out[j] = u / v;
            
            }
            
        /* Otherwise, store the result where the loop failed. */
        else
            out[j] = fun->vals.real[k];
        
        }
        
    /* All is well... */
    return fun_err_ok;

    }


/**
 * @brief Release all memory associated with a #fun.
 *
 * @param fun The #fun to be cleaned-up.
 * @return #fun_err_ok or < 0 on error.
 */
 
int fun_clean ( struct fun *fun ) {

    /* Check for null. */
    if ( fun == NULL )
        return error(fun_err_null);
        
    /* Check if this fun was initialized at all. */
    if ( fun->flags & fun_flag_init ) {
        free( fun->points );
        free( fun->vals.real );
        free( fun->coeffs.real );
        }
        
    /* Re-set some values. */
    fun->flags = fun_flag_none;
    fun->n = 0;
    fun->scale = 1.0;
    fun->points = NULL;
    fun->vals.real = NULL;
    fun->coeffs.real = NULL;
        
    /* All is well. */
    return fun_err_ok;

    }

/**
 * @brief Constructs a @a fun from a real scalar-valued function.
 *
 * @param fun The #fun structure to be initialized.
 * @param fx A pointer to the target function. @a fx takes two parameters,
 *      a double @a x and an optional void pointer and returns the value
 *      of the function at @a x.
 * @param a
 * @param b The left and right boundaries of the interval over which
 *      @a f is to be approximated.
 * @param N The desired length of the fun.
 * @param data An pointer to additional data that will be passed
 *      to @a fx at every call.
 * @return #fun_err_ok or < 0 on error. 
 */

int fun_create_nonadapt ( struct fun *fun , double (*fx)( double x , void * ) , double a , double b , unsigned int N , void *data ) {
    
	int j;
	double a05, b05;

    /* Check inputs. */
    if ( fun == NULL || fx == NULL )
        return error(fun_err_null);

	/* initialise the fun */
	fun_init( fun , N );
	fun->a = a;
	fun->b = b;

	a05 = 0.5*fun->a; 
	b05 = 0.5*fun->b;

	/* Evaluate the function */
	for ( j = 0 ; j < N ; j++ )
    	fun->vals.real[j] = (fx)( b05 * (fun->points[j] + 1.0) + a05 * (1.0 - fun->points[j]) , data );

    /* Compute the coeffs from the values. */
    if ( util_chebpoly( fun->vals.real , N , fun->coeffs.real ) < 0 )
        return error(fun_err_util);
       
    /* Update the scale */
    _fun_rescale( fun );
        
    /* All is well. */
    return fun_err_ok;
    
    }


/**
 * @brief Constructs a @a fun from a real scalar-valued function.
 *
 * @param fun The #fun structure to be initialized.
 * @param fx A pointer to the target function. @a fx takes two parameters,
 *      a double @a x and an optional void pointer and returns the value
 *      of the function at @a x, e.g. for the function @f$f(x)=\sin(4\pi x)@f$:
 *      @code
double myfun ( double x , void *data ) {

    return sin( 4.0 * M_PI * x );

    }
        @endcode
 *      The void pointer is set to the parameter
 *      @a data and can be used to pass additional data to the function,
 *      e.g. for the function @f$f(x) = \sin( \omega \pi x )@f$:
 *      @code
double myfun ( double x , void *data ) {

    double omega = *((double *)data);

    return sin( omega * M_PI * x );

    }
        @endcode
 *      where @c fun_create is then called as
 *      @code
 struct fun f;
 double omega = 4.0;
 fun_create( &f , &myfun , -1.0 , 1.0 , &chebopts_default , &omega );
        @endcode
 * @param a
 * @param b The left and right boundaries of the interval over which
 *      @a f is to be approximated.
 * @param opts A pointer to a #chebopts structure containing the
 *      parameters that will be used to construct the fun. If
 *      this parameter is set to @a NULL, the default parameters
 *      #chebopts_default will be used.
 * @param data An pointer to additional data that will be passed
 *      to @a fx at every call.
 * @return #fun_err_ok or < 0 on error.
 *
 * This function only wraps the function @a fx into a vector-valued
 * function an calls #fun_create_vec.
 *
 * @sa fun_create_vec
 */
 
int fun_create ( struct fun *fun , double (*fx)( double x , void * ) , double a , double b , const struct chebopts *opts , void *data ) {

    /* Check inputs. */
    if ( fun == NULL || fx == NULL )
        return error(fun_err_null);
        
    /* declare a wrapper function for fun_create_vec. */
    int fx_wrapper ( const double *x , unsigned int N , double *out , void *data ) {
        int i;
        if ( x == NULL || fx == NULL )
            return error(fun_err_null);
        for ( i = 0 ; i < N ; i++ )
            out[i] = (*fx)( x[i] , data );
        return fun_err_ok;
        }
        
    /* Call the vectorized version with the wrapper function. */
    return fun_create_vec( fun , &fx_wrapper , a , b , opts , data );

    }
    
    
/**
 * @brief Constructs a @a fun from a real vector-valued function.
 *
 * @param fun The #fun structure to be initialized.
 * @param fx A pointer to the target function. @a fx takes four parameters,
 *      a pointer to an array of doubles @a x , an integer @a N, a pointer
 *      to an array of doubles @a out and an optional void pointer.
 *      The function @a fx should evaluate the target function at the
 *      @a N values in @a x and store the results in the pre-allocated
 *      array @a out, e.g. for the function @f$f(x)=\sin(4\pi x)@f$:
 *      @code
int myfun ( const double *x , unsigned int N , double *out , void *data ) {

    int k;

    for ( k = 0 ; k < N ; k++ )
        out[k] = sin( 4.0 * M_PI * x[k] );

    return 0;

    }
        @endcode
 *      It is assumed that the function will return a negative value
 *      if an error occurs.
 *      Note that this may be more efficient than using the scalar
 *      version as the compiler may be able to vectorize
 *      operations inside the loop.
 *      The void pointer is set to the parameter
 *      @a data and can be used to pass additional data to the function
 *      (see #fun_create).
 * @param a
 * @param b The left and right boundaries of the interval over which
 *      @a f is to be approximated.
 * @param opts A pointer to a #chebopts structure containing the
 *      parameters that will be used to construct the fun. If
 *      this parameter is set to @a NULL, the default parameters
 *      #chebopts_default will be used.
 * @param data An pointer to additional data that will be passed
 *      to @a fx at every call.
 * @return #fun_err_ok or < 0 on error.
 *
 * @sa fun_create
 */
 
int fun_create_vec ( struct fun *fun , int (*fx)( const double * , unsigned int , double * , void * ) , double a , double b , const struct chebopts *opts , void *data ) {

    double *x, *v, *coeffs, scale = 0.0;
    double m = (a + b) * 0.5, h= (b - a) * 0.5;
    unsigned int N;
    static double *xi = NULL;
    static int nr_xi = 0;
    int k, stride, N_new;
    
    /* Check inputs. */
    if ( fun == NULL || fx == NULL )
        return error(fun_err_null);
        
    /* If no options were specified (NULL), use the default options. */
    if ( opts == NULL )
        opts = &chebopts_default;
        
    /* Initialize the vectors x, v and coeffs to the maximum length.
       Allocation is done on the stack, so we don't need to worry about
       releasing these if the routine fails anywhere underway.
       Note that v and coeffs are padded with 16 bytes so that they can be
       shifted and aligned to 16-byte boundaries. */
    if ( ( x = (double *)alloca( sizeof(double) * (opts->maxdegree + 1) ) ) == NULL ||
         ( v = (double *)alloca( sizeof(double) * (opts->maxdegree + 1) + 16 ) ) == NULL ||
         ( coeffs = (double *)alloca( sizeof(double) * (opts->maxdegree + 1) + 16 ) ) == NULL )
        return error(fun_err_malloc);
        
    /* Shift v and coeffs such that they are 16-byte aligned. */
    v = (double *)( (((size_t)v) + 15 ) & ~15 );
    coeffs = (double *)( (((size_t)coeffs) + 15 ) & ~15 );
    
    /* Init the length N. */
    N = opts->minsamples;
        
    /* Make sure the nodes have been pre-allocated. */
    if ( !(opts->flags & chebopts_flag_resampling) &&
         ( nr_xi < opts->maxdegree + 1 || (nr_xi - 1) % (N - 1) != 0 ) ) {
    
        /* Clean up old nodes. */
        if ( xi != NULL )
            free(xi);
            
        /* Set the nr of nodes. */
        for ( nr_xi = opts->minsamples ; nr_xi < opts->maxdegree + 1 ; nr_xi = 2*nr_xi - 1 );
            
        /* Allocate the nodes. */
        if ( ( xi = (double *)malloc( sizeof(double) * nr_xi ) ) == NULL )
            return error(fun_err_malloc);
            
        /* Fill the nodes. */
        if ( util_chebpts( nr_xi , xi ) < 0 )
            return error(fun_err_util);
            
        }
        
        
    /* Main loop. */
    while ( 1 ) {
    
        /* If we are re-sampling, get the nodes for arbitrary N. */
        if ( opts->flags & chebopts_flag_resampling ) {
        
            /* Get the N chebpts. */
            if ( util_chebpts( N , x ) < 0 )
                return error(fun_err_util);
                
            /* Scale them to the correct interval. */
            for ( k = 0 ; k < N ; k++ )
                x[k] = m + h * x[k];
                
            /* Evaluate fx at the N nodes. */
            if ( (*fx)( x , N , v , data ) < 0 )
                return error(fun_err_fx);
                
            /* Update the scale. */
            for ( k = 0 ; k < N ; k++ )
                if ( fabs(v[k]) > scale )
                    scale = fabs(v[k]);
                
            }
            
        /* Otherwise, get the missing nodes from xi. */
        else {
        
            /* Set the stride for this N. */
            stride = (nr_xi - 1) / (N - 1);
        
            /* If this is the first go, just copy all the nodes. */
            if ( N == opts->minsamples ) {
            
                /* Pick the N nodes out of xi. */
                for ( k = 0 ; k < N ; k++ )
                    x[k] = m + h * xi[ k * stride ];
                    
                /* Evaluate fx at the N nodes. */
                if ( (*fx)( x , N , v , data ) < 0 )
                    return error(fun_err_fx);
                    
                /* Update the scale. */
                for ( k = 0 ; k < N ; k++ )
                    if ( fabs(v[k]) > scale )
                        scale = fabs(v[k]);
                
                }
                
            /* Not the first run, skip nodes from previous go. */
            else {
            
                /* Pick-out the interlacing nodes from xi. */
                for ( k = 0 ; k < N / 2 ; k++ )
                    x[k] = m + h * xi[ (2 * k + 1 ) * stride ];
                    
                /* Evaluate fx at the (N-1)/2 new nodes and store
                    the result in the second half of x. */
                if ( (*fx)( x , N / 2 , &(x[N/2]) , data ) < 0 )
                    return error(fun_err_fx);
                    
                /* Update the scale. */
                for ( k = 0 ; k < N/2 ; k++ )
                    if ( fabs(x[N/2+k]) > scale )
                        scale = fabs(x[N/2+k]);
                
                /* Fill the old and new values of fx into v. */
                for ( k = N/2 ; k >= 0 ; k-- ) {
                    v[2*k] = v[k];
                    v[2*k-1] = x[N/2+k-1];
                    }
                    
                /* Correct the values of x for util_simplify. */
                for ( k = 0 ; k < N ; k++ )
                    x[k] = m + h * xi[ k * stride ];
                    
                }
        
            }
            
        
        /* Compute the coeffs from the values. */
        if ( util_chebpoly( v , N , coeffs ) < 0 )
            return error(fun_err_util);
        
        
        /* Check convergence of the coefficients. */
        if ( ( N_new = util_simplify( x , v , coeffs , N , 2*h , scale , opts ) ) < 0 )
            return fun_err_util;
            
        /* TODO: Sampletest? */
            
        /* Did this converge? */
        if ( N_new < N ) {
            N = N_new;
            break;
            }
            
        /* No convergence, adjust N. */
        if ( ( opts->flags & chebopts_flag_resampling ) && 
             ( N < 64 ) )
            N = (int)( ( N - 1 ) * M_SQRT2 + 1 ) | 1;
        else
            N = 2 * N - 1;
        
        } /* Main loop. */
        
        
    /* Allocate the data inside the fun to the correct size. */
    if ( posix_memalign( (void **)&(fun->vals.real) , 16 , sizeof(double) * N ) != 0 ||
         posix_memalign( (void **)&(fun->coeffs.real) , 16 , sizeof(double) * N ) != 0 ||
         ( fun->points = (double *)malloc( sizeof(double) * N ) ) == NULL )
        return error(fun_err_malloc);
    
    /* Write the data to the fun. */
    memcpy( fun->vals.real , v , sizeof(double) * N );
    memcpy( fun->coeffs.real , coeffs , sizeof(double) * N );
    memcpy( fun->points , x , sizeof(double) * N );
    fun->scale = scale;
    fun->a = a;
    fun->b = b;
    fun->n = N;
    fun->flags = fun_flag_init;
    
    /* If nothing bad happened until here, we're done! */
    return fun_err_ok;

    }
