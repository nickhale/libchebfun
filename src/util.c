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
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <fftw3.h>
#ifdef __SSE2__
    #include <emmintrin.h>
#endif

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
#define error( id )     ( util_err = errs_register( id , util_err_msg[-id] , __LINE__ , __FUNCTION__ , __FILE__ ) )


/**
 * @brief Evaluate a polynomial using Barycentric interpolation
 *
 * @param vals A pointer to an array of values.
 * @param points A pointer to an array of nodes (may be @c NULL).
 * @param N The number of coefficients.
 * @param x The position at which to evaluate the polynomial.
 * @return The polynomial evaluated at @c x or @c NaN if an error occured
 *      (see #util_err).
 *
 * @note Even when @c points are supplied, they are assumed to be the
 * Chebyshev nodes on @f$[-1,1]@f$!
 *
 * @sa util_clenshaw_vec
 */
 
double util_bary_real ( const double *vals , const double *points , unsigned int N , double x ) {

    int k;
    double w, u = 0.0, v = 0.0;
    
    /* Check for nonsense. */
    if ( N != 0 && vals == NULL ) {
        error(util_err_null);
        return NAN;
        }

    /* Check for the simplest cases. */
    if ( N == 0 )
        return 0.0;
    else if ( N == 1 )
        return vals[0];
        
    /* Get the Chebyshev points, if needed. */
    if ( points == NULL ) {
        if ( ( points = (double *)alloca( sizeof(double) * N ) ) == NULL ) {
            error(util_err_malloc);
            return NAN;
            }
        util_chebpts( N , (double *)points );
        }
        
    /* Do the barycentric form in place. */
    /* Do the first node separately due to the half-weight. */
    if ( x == points[0] )
        return vals[0];
    w = 0.5 / ( x - points[0] );
    u = vals[0] * w;
    v = w;

    /* Do the interior nodes. */
    for ( k = 1 ; k < N-1 ; k++ ) {
        if ( x == points[k] )
            return vals[k];
        w = (double)( 1 - 2 * ( k & 1 ) ) / ( x - points[k] );
        u += vals[k] * w;
        v += w;
        }
        
    /* Do the last node separately due to the half-weight. */
    if ( x == points[k] )
        return vals[k];
    w = ( 0.5 - ( k & 1 ) ) / ( x - points[k] );
    u += vals[k] * w;
    v += w;

    /* Return the fraction. */
    return u / v;

    }
    
    
/**
 * @brief Evaluate a polynomial using Barycentric interpolation
 *
 * @param vals A pointer to an array of values.
 * @param points A pointer to an array of nodes (may be @c NULL).
 * @param N The number of coefficients.
 * @param x A pointer to the positions at which to evaluate the polynomial.
 * @param M The number of points in @c x.
 * @return The polynomial evaluated at @c x or @c NaN if an error occured
 *      (see #util_err).
 *
 * @note Even when @c points are supplied, they are assumed to be the
 * Chebyshev nodes on @f$[-1,1]@f$!
 *
 * @sa util_bary_real, util_clenshaw_vec
 */
 
double util_bary_vec_real ( const double *vals , const double *points , unsigned int N , const double *x , double *out , unsigned int M ) {

    int j, k, start = 0;
    double w, u = 0.0, v = 0.0, xj;
    #ifdef __SSE2__
    union {
        double f[2];
        unsigned long long i[2];
        __v2df v;
        } w2, u2, v2, xj2, pts2;
    #endif
    
    /* Check for nonsense. */
    if ( N != 0 && vals == NULL )
        return error(util_err_null);

    /* Check for the simplest cases. */
    if ( N == 0 ) {
        for ( k = 0 ; k < M ; k++ )
            out[k] = 0.0;
        return util_err_ok;
        }
    else if ( N == 1 ) {
        for ( k = 0 ; k < M ; k++ )
            out[k] = vals[0];
        return util_err_ok;
        }
        
    /* Get the Chebyshev points, if needed. */
    if ( points == NULL ) {
        if ( ( points = (double *)alloca( sizeof(double) * N ) ) == NULL )
            return error(util_err_malloc);
        util_chebpts( N , (double *)points );
        }
        
    #ifdef __SSE2__
    
        /* Check if x and out are properly aligned. */
        if ( ((size_t)x & 15) == 0 && ((size_t)out & 15) == 0 ) {
        
            /* Loop through the x in chunks of two. */
            for ( start = 0 ; start < (int)M-3 ; start += 2 ) {
            
                /* Get the pair for xj. */
                xj2.v = _mm_load_pd( &x[start] );
            
                /* Do the barycentric form in place. */
                /* Do the first node separately due to the half-weight. */
                pts2.v = _mm_set1_pd( points[0] );
                w2.v = _mm_div_pd( _mm_set1_pd( 0.5 ) , _mm_sub_pd( xj2.v , pts2.v ) );
                u2.v = _mm_mul_pd( _mm_set1_pd( vals[0] ) , w2.v );
                v2.v = w2.v;

                /* Do the interior nodes. */
                for ( k = 1 ; k < N-1 ; k++ ) {
                    pts2.v = _mm_set1_pd( points[k] );
                    w2.v = _mm_div_pd( _mm_set1_pd( (double)( 1 - 2 * ( k & 1 ) ) ) , _mm_sub_pd( xj2.v , pts2.v ) );
                    u2.v = _mm_add_pd( u2.v , _mm_mul_pd( _mm_set1_pd( vals[k] ) , w2.v ) );
                    v2.v = _mm_add_pd( v2.v , w2.v );
                    }

                /* Do the last node separately due to the half-weight. */
                pts2.v = _mm_set1_pd( points[k] );
                w2.v = _mm_div_pd( _mm_set1_pd( (double)( 0.5 - ( k & 1 ) ) ) , _mm_sub_pd( xj2.v , pts2.v ) );
                u2.v = _mm_add_pd( u2.v , _mm_mul_pd( _mm_set1_pd( vals[k] ) , w2.v ) );
                v2.v = _mm_add_pd( v2.v , w2.v );
                
                /* Store the results. */
                _mm_store_pd( &out[start] , _mm_div_pd( u2.v , v2.v ) );

                }

            }
    
    #endif
    
    /* Otherwise, no fancy SSE-stuff or get leftovers. */
    for ( j = start ; j < M ; j++ ) {
    
        /* Do the barycentric form in place. */
        /* Do the first node separately due to the half-weight. */
        xj = x[j];
        w = 0.5 / ( xj - points[0] );
        u = vals[0] * w;
        v = w;

        /* Do the interior nodes. */
        for ( k = 1 ; k < N-1 ; k++ ) {
            w = (double)( 1 - 2 * ( k & 1 ) ) / ( xj - points[k] );
            u += vals[k] * w;
            v += w;
            }

        /* Do the last node separately due to the half-weight. */
        w = ( 0.5 - ( k & 1 ) ) / ( xj - points[k] );
        u += vals[k] * w;
        v += w;
        
        /* Store the result. */
        out[j] = u / v;
        
        }
        
    /* Run through out and look for non-numerical values. */
    for ( j = 0 ; j < M ; j++ )
        if ( !isfinite( out[j] ) ) {
        
            /* Map x[j] back to the corresponding point. */
            out[j] = vals[ (int)round( acos( x[j] ) / M_PI * (N-1) ) ];
                
            }

    /* Happy end. */
    return util_err_ok;

    }
    
    
/**
 * @brief Evaluate a polynomial using the Clenshaw algorithm.
 *
 * @param coeffs A pointer to an array of coefficients.
 * @param N The number of coefficients.
 * @param x The position at which to evaluate the polynomial.
 * @return The polynomial evaluated at @c x or @c NaN if an error occured
 *      (see #util_err).
 *
 * @sa util_clenshaw_vec
 */
 
double util_clenshaw ( double *coeffs , unsigned int N , double x ) {

    int k;
    double yn = 0.0, ynp1 = 0.0, ynp2, x2 = 2*x;

    /* Check for nonsense. */
    if ( coeffs == NULL ) {
        error(util_err_null);
        return NAN;
        }
        
    /* Check for the simplest cases. */
    if ( N == 0 )
        return 0.0;
    else if (N == 1 )
        return coeffs[0];
        
    /* Evaluate the recurrence. */
    for ( k = N-1 ; k >= 0 ; k-- ) {
        ynp2 = ynp1; ynp1 = yn;
        yn = coeffs[k] + x2 * ynp1 - ynp2;
        }
        
    /* Return the result. */
    return yn - x * ynp1;
        
    }
    
    
/**
 * @brief Evaluate a polynomial at a vector of points using the
 *      Clenshaw algorithm.
 *
 * @param coeffs A pointer to an array of coefficients.
 * @param N The number of coefficients.
 * @param x A pointer to an array of @c M double values containing the positions
 *      at which to evaluate the polynomial.
 * @param M The number of positions in @c x.
 * @param out A pointer to an array of doubles of length @c M in which the
 *      values of the polynomial at @c x will be stored.
 * @return #util_err_ok or < 0 if an error occured (see #util_err).
 *
 * @sa util_clenshaw
 */

int util_clenshaw_vec ( double *coeffs , unsigned int N , double *x , unsigned int M , double *out ) {

    int j, k, start = 0;
    #ifdef __SSE2__
        __v2df vxj_1, vyn_1, vynp1_1, vynp2_1;
        __v2df vxj_2, vyn_2, vynp1_2, vynp2_2;
        __v2df c, half = _mm_set1_pd( 0.5 ), two = _mm_set1_pd( 2.0 );
    #endif
    double yn, ynp1, ynp2, x2;

    /* Check for nonsense. */
    if ( coeffs == NULL || x == NULL || out == NULL )
        return error(util_err_null);
        
    /* Check for the simplest cases. */
    if ( N < 2 ) {
        if ( N == 0 )
            for ( k = 0 ; k < M ; k++ )
                out[k] = 0.0;
        else
            for ( k = 0 ; k < M ; k++ )
                out[k] = coeffs[0];
        return util_err_ok;
        }
        
    /* Loop over the input values. */
    #ifdef __SSE2__
    /* Are all the values aligned correctly? */
    if ( ((size_t)x) % 16 == 0 && ((size_t)out) % 16 == 0 ) {
    
        /* Take two chunks at a time. */
        for ( start = 0 ; start < (int)M-3 ; start += 4 ) {

            /* Init the recurrence. */
            vxj_1 = _mm_load_pd( &(x[start]) );
            vxj_2 = _mm_load_pd( &(x[start+2]) );
            vxj_1 = _mm_mul_pd( vxj_1 , two );
            vxj_2 = _mm_mul_pd( vxj_2 , two );
            vynp1_1 = _mm_setzero_pd();
            vynp1_2 = _mm_setzero_pd();
            vyn_1 = _mm_setzero_pd();
            vyn_2 = _mm_setzero_pd();

            /* Evaluate the recurrence. */
            for ( k = N-1 ; k >= 0 ; k-- ) {
                vynp2_1 = vynp1_1; vynp1_1 = vyn_1;
                vynp2_2 = vynp1_2; vynp1_2 = vyn_2;
                c = _mm_set1_pd( coeffs[k] );
                vyn_1 = _mm_add_pd( c , _mm_sub_pd( _mm_mul_pd( vxj_1 , vynp1_1 ) , vynp2_1 ) );
                vyn_2 = _mm_add_pd( c , _mm_sub_pd( _mm_mul_pd( vxj_2 , vynp1_2 ) , vynp2_2 ) );
                }

            /* store the result. */
            _mm_store_pd( &(out[start]) , _mm_sub_pd( vyn_1 , _mm_mul_pd( _mm_mul_pd( half , vxj_1 ) , vynp1_1 ) ) );
            _mm_store_pd( &(out[start+2]) , _mm_sub_pd( vyn_2 , _mm_mul_pd( _mm_mul_pd( half , vxj_2 ) , vynp1_2 ) ) );

            }
            
        /* Two or more left? */
        if ( start < (int)M-1 ) {

            /* Init the recurrence. */
            vxj_1 = _mm_load_pd( &(x[start]) );
            vxj_1 = _mm_mul_pd( vxj_1 , two );
            vynp1_1 = _mm_setzero_pd();
            vyn_1 = _mm_setzero_pd();

            /* Evaluate the recurrence. */
            for ( k = N-1 ; k >= 0 ; k-- ) {
                vynp2_1 = vynp1_1; vynp1_1 = vyn_1;
                vyn_1 = _mm_add_pd( _mm_set1_pd( coeffs[k] ) , _mm_sub_pd( _mm_mul_pd( vxj_1 , vynp1_1 ) , vynp2_1 ) );
                }

            /* store the result. */
            _mm_store_pd( &(out[start]) , _mm_sub_pd( vyn_1 , _mm_mul_pd( _mm_mul_pd( half , vxj_1 ) , vynp1_1 ) ) );

            }
            
        }
    #endif
      
            
    /* Otherwise, no SSE-stuff. */
    for ( j = start ; j < M ; j++ ) {

        /* Init the recurrence. */
        x2 = 2 * x[j];
        ynp1 = 0.0; yn = 0.0;

        /* Evaluate the recurrence. */
        for ( k = N-1 ; k >= 0 ; k-- ) {
            ynp2 = ynp1; ynp1 = yn;
            yn = coeffs[k] + x2 * ynp1 - ynp2;
            }

        /* store the result. */
        out[j] = yn - 0.5 * x2 * ynp1;

        }
        
            
    /* If all went well... */
    return util_err_ok;
        
    }
    
    
/**
 * @brief construct the Chebyshev differentiation matrix.
 *
 * @param N Size of the matrix.
 *
 * @return Pointer to the matrix (stored as a vector)..
 */

double *util_diffmat ( unsigned int N ) {

    double *x, *w, *D, s, val, xj, wj;
    int j, k, sgn = 1.0; 
    

    /* Trivial case */
    if ( N == 1 ) {
        if ( ( D = (double *)malloc( sizeof(double) * (N * N) + 16 ) ) == NULL ) {
            error(util_err_malloc);
            return NULL;
            }
        D[0] = 0.0;
        return D;
        }

    /* Allocate the memory */
    if ( ( x = (double *)alloca( sizeof(double) * N ) ) == NULL ||
         ( w = (double *)alloca( sizeof(double) * N ) ) == NULL ||
         ( D = (double *)malloc( sizeof(double) * (N * N) + 16 ) ) == NULL ) {
        error(util_err_malloc);
        return NULL;
        }
    if ( posix_memalign( (void **)(&D) , 16 , sizeof(double) * N * N ) != 0 ) {
        error(util_err_malloc);
        return NULL;
        }

    /* Get the Chebyshev points. */
    if ( util_chebpts ( N , x ) != 0 ) {
        error(util_err_malloc);
        return NULL;
        }

    /* Get the weights. */
    for ( k = 0 ; k < N ; k++ ) {
        w[k] = sgn;
        sgn = -sgn;
        }
    w[0] = 0.5*w[0];
    w[N-1] = 0.5*w[N-1];

    /* Make D. */
    for ( j = 0 ; j < N ; j++ ) {
        s = 0.0; wj = w[j]; xj = x[j];
        for ( k = 0; k < N; k++ ) {
            if ( k == j )
                continue;
            val =  w[k]/( ( x[k]-xj ) * wj );
            D[k+j*N] = val;
            s += val;
            }
        /* Negative sum trick */
        D[j*(N+1)] = -s;
        }

    /* print
    for ( j = 0 ; j < N ; j++ ) {
        for ( k = 0; k < N ; k++ ) {
            printf("%8.8e ",D[k+j*N]);
            }
        printf("\n");
        }
    */

    /* Weee! */
    return D;

    }

/**
 * @brief Checks if a given interpolation is converged or not and
 *  returns the expected number of points necessary for its
 *  representation.
 *
 * @param v Pointer to an array of double values containing the
 *      function values at the Chebyshev nodes.
 * @param coeffs Pointer to an array of double values containing
 *      the Chebyshev coefficients of the interpolation.
 * @param N Number of nodes.
 * @param hscale Horizontal scale of the interval.
 * @param vscale Vertical scale of the interval, i.e. @f$\max_k |v_k|@f$.
 * @param eps The relative tolerace for which to test.
 * @return The expected number of points necessary for the accurate
 *      representation of the function or < 0 on error. If the
 *      interpolation is not converged, @a N is returned.
 */
 
/* TODO: add sampletest here? */
 
int util_simplify ( double *v , double *coeffs , unsigned int N , double hscale , double vscale , double eps ) {

    int k, tail;
    double temp, tail_max, diff_max, dx, cos_new, pin, cos_last;
    
    /* The usual checks and balances. */
    if ( v == NULL || coeffs == NULL )
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
    pin = M_PI / (N - 1); cos_last = 1.0;
    for ( k = 1 ; k < N ; k++ ) {
        temp = fabs( v[k] - v[k-1] );
        cos_new = cos( k * pin );
        dx = cos_last - cos_new;
        cos_last = cos_new;
        if ( dx > 2 * DBL_EPSILON )
            temp /= dx;
        else
            temp /= 2 * DBL_EPSILON;
        if ( temp > diff_max )
            diff_max = temp;
        }
    diff_max /= vscale;
    tail_max = pow( tail , 2.0/3.0 );
    if ( tail_max < diff_max )
        tail_max = diff_max;
    if ( tail_max > 1.0e12 )
        tail_max = 1.0e12;
    tail_max *= DBL_EPSILON;
    if ( tail_max < eps )
        tail_max = eps;
    tail_max *= vscale;
        
    /* Get the index of the last coefficient |coeffs[k]| < tail_max. */
    for ( k = N-1 ; k >= 0 && fabs(coeffs[k]) < tail_max ; k-- );
    
    /* Is the function actually just zero? */
    if ( k < 0 ) {
    
        /* Set N to 1. */
        N = 1;
        
        /* Set the coefficients and values. */
        coeffs[0] = 0.0;
        v[0] = 0.0;
    
        }
    
    /* Is the tail long enough? */
    else if ( N - k > tail ) {
    
        /* For now, just use k+1... */
        /* TODO: actually implement the strategy from simplify.m! */
        N = k+1;
        
        /* Get the new values. */
        if ( util_chebpolyval( coeffs , N , v ) < 0 )
            return error(util_err);
            
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
        
    /* Catch the trivial case of N=1. */
    if ( N == 1 ) {
    
        /* Set v. */
        v[0] = coeffs[0];
        
        }
        
    /* Case for N > 1. */
    else {
        
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
        
        }
    
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
 
int util_chebptsAB ( unsigned int N , double A, double B , double *x ) {

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
		x[j] = A05 * (1.0 - x[j]) + B05 * (1.0 + x[j]);
        
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
        
    /* Catch the trivial case of N=1. */
    if ( N == 1 ) {
    
        c[0] = v[0];
        
        }
        
    /* Otherwise, call fftw. */
    else {
        
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
            
        }
    
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

