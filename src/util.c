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
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <fftw3.h>

/* Local includes. */
#include "errs.h"
#include "util.h"


/* Global constants. */
/** The value of the last error code. */
int util_err = util_err_ok;
/** Error messages for the (absolute) error codes. */
const char *util_err_msg[] = {
    "All is well.",
    "An unexpected NULL-pointer was encountered.",
    "A call to malloc failed, probably due to insufficient memory." };

/* Define a macro to store the errors. */
#define error( id )     ( util_err = errs_register( id , util_err_msg[-id] , __LINE__ , __FILE__ ) )
    

/**
 * @brief Checks if a given interpolation is converged or not and
 *  returns the expected number of points necessary for its
 *  representation.
 *
 * @param x Pointer to an array of double values containing the
 *      interpolation nodes.
 * @param v Pointer to an array of double values containing the
 *      function values at the nodes.
 * @param coeffs Pointer to an array of double values containing
 *      the Chebyshev coefficients of the interpolation.
 * @param N Number of nodes.
 * @param hscale Horizontal scale of the interval.
 * @param vscale Vertical scale of the interval, i.e. @f$\max_k |v_k|@f$.
 * @param opts Pointer to a #chebopts structure or @c NULL to use the
 *      default parameters #chebopts_default.
 * @return The expected number of points necessary for the accurate
 *      representation of the function or < 0 on error. If the
 *      interpolation is not converged, @a N is returned.
 */
 
/* TODO: add sampletest here? */
 
int util_simplify ( double *x , double *v , double *coeffs , unsigned int N , double hscale , double vscale , double eps ) {

    int k, tail;
    double temp, tail_max, diff_max;
    
    /* The usual checks and balances. */
    if ( x == NULL || v == NULL || coeffs == NULL )
        return error(util_err_null);
        
    /* Get the expected tail length. */
    if ( N < 4 )
        tail = N;
    else if ( N < 26 )
        tail = 3;
    else
        tail = ( N - 1 ) / 8;

    /* Get the maximum tail magnitude. */
    diff_max = 0.0;
    for ( k = 1 ; k < N ; k++ ) {
        temp = fabs( v[k] - v[k-1] );
        if ( fabs( x[k] - x[k-1] ) > hscale * DBL_EPSILON )
            temp /= fabs( x[k] - x[k-1] );
        else
            temp /= hscale * DBL_EPSILON;
        if ( temp > diff_max )
            diff_max = temp;
        }
    diff_max *= ( hscale / vscale );
    tail_max = pow( tail , 2.0/3.0 );
    if ( tail_max < diff_max )
        tail_max = diff_max;
    if ( tail_max > 1.0e12 )
        tail_max = 1.0e12;
    tail_max *= DBL_EPSILON;
    if ( tail_max < eps )
        tail_max = eps;
        
    /* Get the index of the last coefficient |coeffs[k]| < tail_max. */
    for ( k = N-1 ; k >= 0 && fabs(coeffs[k]) < tail_max ; k-- );
    
    /* Is the tail long enough? */
    if ( N - k > tail ) {
    
        /* For now, just use k+1... */
        /* TODO: actually implement the strategy from simplify.m! */
        N = k+1;
        
        /* Get the new values. */
        if ( util_chebpolyval( coeffs , N , v ) < 0 )
            return util_err;
            
        /* Get the new nodes. */
        if ( util_chebpts( N , x ) < 0 )
            return util_err;
    
        }
        
    /* Otherwise, don't change anything. */
    return N;
        
    }
    
    
/**
 * @brief Evaluate a polynomial given by its Chebyshev coefficients at the
 *      Chebyshev nodes using a Discrete Cosine Transform (DCT).
 *
 * @param coeffs Pointer to an array of real Chebyshev coefficients.
 * @param N Number of Chebyshev coefficients.
 * @param v Pointer to an array of doubles of length at least @a N
 *      where the results will be stored.
 * @return #util_err_ok or < 0 on error.
 * @sa util_chebpolyval_alloc
 */
 
int util_chebpolyval ( double *coeffs , unsigned int N , double *v ) {

    int k;
    fftw_plan plan;
    
    /* Check sanity of inputs. */
    if ( coeffs == NULL || v == NULL )
        return error(util_err_null);
        
    /* Mind the edges... */
    coeffs[0] *= 2.0;
    coeffs[N-1] *= 2.0;

    /* Note, as of here, that fftw's routines do not produce error codes! */
    /* Create a plan for the DCT. */
    plan = fftw_plan_r2r_1d( N , coeffs , v , FFTW_REDFT00 , FFTW_ESTIMATE );
    
    /* Execute our devious plan. */
    fftw_execute(plan);
    
    /* Destroy the old plan (Mwahaha!). */
    fftw_destroy_plan(plan);
    
    /* Mind the edges... */
    coeffs[0] *= 0.5;
    coeffs[N-1] *= 0.5;

    /* Scale the function values back. */
    for ( k = 0 ; k < N ; k++ )
        v[k] *= 0.5;
    
    /* If all went well... */
    return util_err_ok;

    }
    
    
/**
 * @brief Allocates and fills an array with the values of a polynomial
 *      given by its Chebyshev coefficients at the Chebyshev nodes using
 *      a Discrete Cosine Transform.
 *
 * @param coeffs Pointer to an array of real Chebyshev coefficients.
 * @param N Number of Chebyshev coefficients.
 * @return A pointer to the newly-allocated array or @c NULL on error
 *      (error code in #util_err).
 * @sa util_chebpolyval
 */
 
double *util_chebpolyval_alloc ( double *coeffs , unsigned int N ) {

    double *v;

    /* Check the inputs. */
    if ( coeffs == NULL ) {
        error(util_err_null);
        return NULL;
        }
        
    /* Allocate the return value array. */
    if ( posix_memalign( (void **)&(v) , 16 , sizeof(double) * N ) != 0 ) {
        error(util_err_malloc);
        return NULL;
        }
        
    /* Just call util_chebpoly and check the result. */
    if ( util_chebpolyval( coeffs , N , v ) < 0 ) {
        free(v);
        return NULL;
        }
    else
        return v;

    }
    

/**
 * @brief Fills an array with @a N Chebyshev points on [A,B].
 *
 * @param N Number of Chebyshev points to create.
 * @param x Pointer to an array of double values of length at least @a N.
 * @param A Left endpoint.
 * @param B Right endpoint.
 * @return #util_err_ok or < 0 on error.
 *
 * This function is just a wrapper for util_chebpts.
 *
 * @sa util_chebpts_alloc, util_chebpts
 * TODO: don't call chebpts first, do it all in one loop, will be faster
 *      since the whole thing can be pipelined.
 */
 
int util_chebptsAB ( unsigned int N , double *x , double A, double B ) {

    int j;
	double A05 = 0.5*A, B05 = 0.5*B;

    /* Check inputs. */
    if ( x == NULL )
        return error(util_err_null);
        
    /* Call util_chebpts on [-1 1] */
	if ( util_chebpts( N , x ) < 0 )
		return error(util_err_null);

	/* Scale to [A B] */
	for ( j = 0 ; j < N ; j++ )
		x[j] = A05 * (1.0 + x[j]) + B05 * (1.0 - x[j]);
        
    /* If all went well... */
    return util_err_ok;

    }

/**
 * @brief Fills an array with @a N Chebyshev points on [-1,1].
 *
 * @param N Number of Chebyshev points to create.
 * @param x Pointer to an array of double values of length at least @a N.
 * @return #util_err_ok or < 0 on error.
 *
 * This function computes the values \f$x_i = \cos\left( \pi  i/(N-1) \right)\f$
 * for \f$i=0\dots N-1\f$.
 * Note that to symmetry is preserverd artificially and that the computed
 * values may not match those produced by @a cos(..) exactly.
 *
 * @sa util_chebpts_alloc
 */
 
int util_chebpts ( unsigned int N , double *x ) {

    int i;
    double w = M_PI / (N - 1);

    /* Check inputs. */
    if ( x == NULL )
        return error(util_err_null);
        
    /* Fill the values (fist half) */
    for ( i = 0 ; i < (N+1)/2 ; i++ )
        x[i] = cos( i * w );
        
    /* Fill the remaining values, exploiting symmetry. */
    for ( i = (N+1)/2 ; i < N ; i++ )
        x[i] = -x[ N - i - 1 ];
        
    /* If all went well... */
    return util_err_ok;

    }
    

/**
 * @brief Allocates and fills an array of @a N Chebyshev points on [-1,1].
 *
 * @param N Number of Chebyshev points to create.
 * @return A pointer to the newly-allocated array or @c NULL on error
 *      (error code in #util_err).
 *
 * This function computes the values \f$x_i = \cos\left( \pi  i/(N-1) \right)\f$
 * for \f$i=0\dots N-1\f$.
 * 
 * Note that to symmetry is preserverd artificially and that the returned
 * values may not match those produced by @a cos(..) exactly.
 *
 * @sa util_chebpts
 */
 
double * util_chebpts_alloc ( unsigned int N ) {

    double *x;
    
    /* Allocate the array. */
    if ( ( x = (double *)malloc( sizeof(double) * N ) ) == NULL ) {
        error(util_err_malloc);
        return NULL;
        }
        
    /* Call util_chebpts to fill x and check if anything went wrong. */
    if ( util_chebpts( N , x ) < 0 ) {
        free(x);
        return NULL;
        }
    else
        return x;

    }


/**
 * @brief Computes the Chebyshev coefficients of a real-valued polynomial
 *      sampled at the Chebyshev nodes.
 *
 * @param v Pointer to an array of doubles containing the polynomial
 *      sampled at the @a N Chebyshev nodes.
 * @param N Number of nodes.
 * @param c Pointer to an array of doubles of length at least @a N
 *      where the Chebyshev coefficients will be stored.
 * @return #util_err_ok or < 0 on error.
 *
 * This function computes the coefficients
 *
 *      \f[ c_k = \frac{1}{2N-2}\left[v_0 + (-1)^kv_{N-1} + 2 \sum_{j=1}^{N-2} v_j \cos(\pi j k / (N-1))\right] \f]
 *
 * using the Discrete Cosine Transform (DCT) provided by FFTW.
 *
 * Note that for the computation of the DCT to exploit SIMD instructions
 * (e.g. SSE2 or AltiVec), the pointers in the input and output arrays should
 * be aligned to 16 bytes (see posix_memalign).
 *
 * @sa util_chebpoly_alloc
 */
 
int util_chebpoly ( double *v , unsigned int N , double *c ) {

    int k;
    double w = 1.0 / ( N - 1 );
    fftw_plan plan;
    
    /* Check the inputs. */
    if ( v == NULL || c == NULL )
        return error(util_err_null);
        
    /* Note, as of here, that fftw's routines do not produce error codes! */
    /* Create a plan for the DCT. */
    plan = fftw_plan_r2r_1d( N , v , c , FFTW_REDFT00 , FFTW_ESTIMATE );
    
    /* Execute our devious plan. */
    fftw_execute(plan);
    
    /* Destroy the old plan (Mwahaha!). */
    fftw_destroy_plan(plan);
    
    /* Mind the edges. */
    c[0] *= 0.5;
    c[N-1] *= 0.5;
    
    /* Scale the coefficients. */
    for ( k = 0 ; k < N ; k++ )
        c[k] *= w;
    
    /* If all went well... */
    return util_err_ok;

    }


/**
 * @brief Allocates and fills an array of the Chebyshev coefficients
 *      of a real-valued polynomial sampled at the Chebyshev nodes.
 *
 * @param v Pointer to an array of doubles containing the polynomial
 *      sampled at the @a N Chebyshev nodes.
 * @param N Number of nodes.
 * @return A pointer to the newly-allocated array or @c NULL on error
 *      (error code in #util_err).
 *
 * This function computes the coefficients
 *
 *      \f[ c_k = \frac{1}{2N-2}\left[v_0 + (-1)^kv_{N-1} + 2 \sum_{j=1}^{N-2} v_j \cos(\pi j k / (N-1))\right] \f]
 *
 * using the Discrete Cosine Transform (DCT) provided by FFTW.
 *
 * Note that for the computation of the DCT to exploit SIMD instructions
 * (e.g. SSE2 or AltiVec), the pointers to the input array should
 * be aligned to 16 bytes (see posix_memalign).
 *
 * @sa util_chebpoly
 */
 
double * util_chebpoly_alloc ( double *v , unsigned int N ) {

    double *c;
    
    /* Check the inputs. */
    if ( v == NULL ) {
        error(util_err_null);
        return NULL;
        }
        
    /* Allocate the return value array. */
    if ( posix_memalign( (void **)&(c) , 16 , sizeof(double) * N ) != 0 ) {
        error(util_err_malloc);
        return NULL;
        }
        
    /* Just call util_chebpoly and check the result. */
    if ( util_chebpoly(v,N,c) < 0 ) {
        free(c);
        return NULL;
        }
    else
        return c;

    }

