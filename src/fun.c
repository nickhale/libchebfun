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
#include <float.h>
#include <string.h>
#include <clapack.h>
#include <fftw3.h>

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
    "Requested operation is not yet implemented.",
    "A call to LAPACK ended in a calamity.",
    "A call to GNUPLOT wasn't very cool." };
    
/* Define a macro to store the errors. */
#define error( id )     ( fun_err = errs_register( id , fun_err_msg[-id] , __LINE__ , __FUNCTION__ , __FILE__ ) )
    
    
/* Constant struct for virgin funs. */
const struct fun FUN_EMPTY = { 0 , 0 , 0 , 0.0 , 0.0 , 0.0 , NULL , { NULL } , { NULL } };


/**
 * @brief Constructs a #fun of x on the interval a b.
 *
 * @param fun The #fun structure to be initialized.
 * @param a
 * @param b The left and right boundaries of the interval.
 * @return #fun_err_ok or < 0 on error. 
 */

int fun_create_x ( struct fun *fun , double a , double b ) {

    /* Check inputs. */
    if ( fun == NULL )
        return error(fun_err_null);

	/* Initialise the fun */
	fun_init( fun , 2 );
	fun->a = a;
	fun->b = b;

	/* Assign the vals. */
    fun->vals.real[0] = b;
    fun->vals.real[1] = a;

    /* Assign the coeffs. */
    fun->coeffs.real[0] = 0.5*(b+a);
    fun->coeffs.real[1] = 0.5*(b-a);

    /* Assign the points. */
    fun->points[0] = 1.0;
    fun->points[1] = -1.0;
       
    /* Update the scale */
    a = fabs(a); b = fabs(b);
    if ( a > b )
        fun->scale = a;
    else
        fun->scale = b;
        
    /* All is well. */
    return fun_err_ok;
    
    }
    
    
/**
 * @brief Compose a #fun with another #fun
 *
 * @param A The input #fun.
 * @param op A second #fun to compose @c A with.
 * @param B The output #fun, may be the same as @c A.
 *
 * @return #fun_err_ok or < 0 on error (see #fun_err).
 *
 * @sa fun_comp_vec
 */
 
int fun_comp_fun ( struct fun *A , struct fun *op , struct fun *B ) {

    /* Wrapper function. */
    int thefun ( const double *x , unsigned int N , double *v , void *data ) {
        return fun_eval_vec( op , x , N , v );
        }

    /* Check inputs. */
    if ( A == NULL || op == NULL || B == NULL )
        return error(fun_err_null);
    if ( !( A->flags & fun_flag_init ) || !( op->flags & fun_flag_init ) )
        return error(fun_err_uninit);
        
    /* Refer to fun_comp_vec. */
    if ( fun_comp_vec( A , &thefun , B , NULL ) < 0 )
        return error(fun_err);
    else
        return fun_err_ok;

    }


/**
 * @brief Compose a #fun wth a function.
 *
 * @param A The input #fun.
 * @param op The function to compose with.
 * @param B The output fun (which maybe the same as A).
 * @param data Pointer to arbitrary data to be passed to @c op.
 *
 * @return #fun_err_ok or < 0 on error (see #fun_err).
 *
 * Note, fun_comp does not deal with the case of an op
 * of two funs, i.e., op(A,B).
 *
 * This function takes the role of MATLAB/@fun/growfun
 * for when there are 4 input arguments. The op(A,B)
 * case mentioned about would be with 5 inputs.
 *
 * Furthernote, fun_comp does not yet allow "resampling
 * off" (although neither does the MATLAB implementation).
 */

int fun_comp_vec ( struct fun *A , int (*op)( const double * , unsigned int , double * , void * ) , struct fun *B , void *data ) {
    
    double *x, *v, *coeffs, scale = 0.0;
    double a, b;
    double m, h;
    unsigned int N;
    int k, N_new;
    
    /* Routine checks of sanity. */
    if ( ( A == NULL ) || ( op == NULL ) )
        return error(fun_err_null);
    if ( !( A->flags & fun_flag_init ) )
        return error(fun_err_uninit);
        
    /* Get the domain and scalings */
    a = A->a;
    b = A->b;    
    m = (a + b) * 0.5;
    h = (b - a) * 0.5;
        
    /* Initialize the vectors x, v and coeffs to the maximum length.
       Allocation is done on the stack, so we don't need to worry about
       releasing these if the routine fails anywhere underway.
       Note that v and coeffs are padded with 16 bytes so that they can be
       shifted and aligned to 16-byte boundaries. */
    if ( ( x = (double *)alloca( sizeof(double) * (chebopts_opts->maxdegree + 1) ) ) == NULL ||
         ( v = (double *)alloca( sizeof(double) * (chebopts_opts->maxdegree + 1) + 16 ) ) == NULL ||
         ( coeffs = (double *)alloca( sizeof(double) * (chebopts_opts->maxdegree + 1) + 16 ) ) == NULL )
        return error(fun_err_malloc);

    /* Shift v and coeffs such that they are 16-byte aligned. */
    v = (double *)( (((size_t)v) + 15 ) & ~15 );
    coeffs = (double *)( (((size_t)coeffs) + 15 ) & ~15 );
    
    /* Init the length N. */
    N = chebopts_opts->minsamples;
    
    /* Main loop. */
    while ( 1 ) {  

        /* Evaluate A at N chebyshev points and put the values in v. */
        if ( N > A->n ) {
            memcpy( coeffs, A->coeffs.real , sizeof(double) * A->n );
            bzero( &(coeffs[A->n]) , sizeof(double) * (N - A->n) );
            if ( util_chebpolyval( coeffs , N , x ) < 0 )
                return error(fun_err_util);
            }
        else if ( A->n > N ) {
            if ( util_chebpts( N , x ) < 0 )
                return error(fun_err_util);
            if ( fun_eval_clenshaw_vec( A , x , N , x ) < 0 )
                return error(fun_err);
            }
        else
            memcpy( x , A->vals.real , sizeof(double) * N );
                
        /* If resampling is on or this is the first run, just evaluate the op. */
        if ( N == chebopts_opts->minsamples || chebopts_opts->flags & chebopts_flag_resampling ) {
        
            /* Evaluate the op. */
            if ( (*op)( x , N , v , data ) < 0 )
                return error(fun_err_fx);
                
            }
            
        /* Otherwise, we only need to re-evaluate at the interlacing points. */
        else {
        
            /* Wrap up the points. */
            for ( k = 0 ; k < N/2 ; k++ )
                x[k] = x[2*k+1];
                
            /* Evaluate op at the points and store the result in the second half of v. */
            if ( (*op)( x , N/2 , &(x[N/2]) , data ) < 0 )
                return error(fun_err_fx);
                
            /* Unwrap the results into v. */
            for ( k = N/2 ; k > 0 ; k-- ) {
                v[2*k] = v[k];
                v[2*k-1] = x[N/2+k-1];
                }
                
            }

        /* Update the scale. */
        for ( k = 0 ; k < N ; k++ )
            if ( fabs(v[k]) > scale )
                scale = fabs(v[k]);

        /* Compute the coeffs from the values. */
        if ( util_chebpoly( v , N , coeffs ) < 0 )
            return error(fun_err_util);
            
        /* Get a fresh set of points. */
        if ( util_chebpts( N , x ) < 0 )
            return error(fun_err_util);

        /* Check convergence of the coefficients. */
        if ( ( N_new = util_simplify( x , v , coeffs , N , 2*h , scale , chebopts_opts->eps ) ) < 0 )
            return fun_err_util;

        /* TODO: Sampletest? */

        /* Did this converge? */
        if ( N_new < N ) {
            N = N_new;
            break;
            }

        /* No convergence, adjust N. */
        if ( ( chebopts_opts->flags & chebopts_flag_resampling ) && 
             ( N < 64 ) )
            N = (int)( ( N - 1 ) * M_SQRT2 + 1 ) | 1;
        else
            N = 2 * N - 1;

        } /* Main loop. */
        
    /* Allocate the data inside the fun to the correct size. */
    if ( !( B->flags & fun_flag_init ) || B->size < N ) {
        if ( _fun_alloc( B , N ) < 0 )
            return error(fun_err);
        }
    else
        B->n = N;
    
    /* Write the data to the fun. */
    memcpy( B->vals.real , v , sizeof(double) * N );
    memcpy( B->coeffs.real , coeffs , sizeof(double) * N );
    /* Get the N chebpts. */
    if ( util_chebpts( N , B->points ) < 0 )
        return error(fun_err_util);
    B->scale = scale;
    B->a = a;
    B->b = b;

    /* We're outta here! */
    return fun_err_ok;
    
    }


/**
 * @brief Compute the polynomial coefficients of a #fun.
 *
 * @param f1 The input #fun.
 * @param out A vector to contain the coefficients.
 *
 * @return #fun_err_ok or < 0 on error (see #fun_err).
 */

int fun_poly ( struct fun *f1 , double *out ) {

    int j, k, n, lwr;
    int *tn, *tn1, *tn2, *tmp;
    int *fac;
    double alpha, beta, binom, tmpout, *powa, *powb;

    /* Routine checks of sanity. */
    if ( f1 == NULL )
        return error(fun_err_null);
    if ( !( f1->flags & fun_flag_init ) )
        return error(fun_err_uninit);

    n = f1->n;

    if ( n == 1 )  /* Trivial constant case */
        out[0] = f1->coeffs.real[0];

    else if ( n == 2 ) {  /* Trivial linear case */
        out[0] = f1->coeffs.real[0];
        out[1] = f1->coeffs.real[1];
        if ( f1->a != -1.0 || f1->b != 1.0 ) {
            /* Constants for rescaling */
            alpha = 2.0 / (f1->b - f1->a);
            beta = - (f1->b + f1->a) / (f1->b - f1->a);
            out[0] = out[0]+beta*out[1];
            out[1] *= alpha;
            }
        } /* End trivial linear case */

    else {  /* General case */
        
        /* Allocate some work memory */
        if ( ( ( tn = (int *)alloca( sizeof(int) * n ) ) == NULL ) ||
             ( ( tn1 = (int *)alloca( sizeof(int) * n ) ) == NULL ) ||
             ( ( tn2 = (int *)alloca( sizeof(int) * n ) ) == NULL ) )
            return error(fun_err_malloc);

        /* Zero it out */
        bzero( tn , sizeof(int) * n );
        bzero( tn1 , sizeof(int) * n );
        bzero( tn2 , sizeof(int) * n );

        /* Initialise */
        tn1[0] = 0;   tn1[1] = 1;  tn2[0] = 1;
        out[1] = f1->coeffs.real[0];
        out[0] = f1->coeffs.real[1];

        /* The unit interval case */
        lwr = 2; 
        for ( j = 2 ; j < n ; j++ ) {
            if ( lwr <= j+2 )
                lwr = j-2;

            /* Update the matrix row */
            tn[0] = -tn2[0];
            for ( k = 1 ; k < j+1 ; k++ )
                tn[k] = 2*tn1[k-1]-tn2[k];

            /* Update the coefficients */
            for ( k = j ; k > 0 ; k-- )
                out[k] = f1->coeffs.real[j]*(double)tn[j-k] + out[k-1];
            out[0] = f1->coeffs.real[j]*tn[j];

            /* Juggle memory for storage */
            tmp = tn2;  tn2 = tn1;  tn1 = tn;  tn = tmp;

            }

        /* Rescale for arbitrary intervals */
        if ( f1->a != -1.0 || f1->b != 1.0 ) {
        
            /* Pre-compute the factorials. */
            if ( ( fac = (int *)alloca( sizeof(int) * n ) ) == NULL )
                return error(fun_err_malloc);
            for ( fac[0] = 1 , k = 1 ; k < n ; k++ )
                fac[k] = fac[k-1] * k;
            
            /* Constants for rescaling */
            alpha = 2.0 / (f1->b - f1->a);
            beta = - (f1->b + f1->a) / (f1->b - f1->a);
            if ( ( ( powa = (double *)alloca( sizeof(double) * n ) ) == NULL ) ||
                 ( ( powb = (double *)alloca( sizeof(double) * n ) ) == NULL ) )
                return error(fun_err_malloc);
            powa[0] = 1.0; powb[0] = 1.0;
            for ( k = 1 ; k < n ; k++ ) {
                powa[k] = powa[k-1] * alpha;
                powb[k] = powb[k-1] * beta;
                }
            
            /* Compute binomial coefficients */
            for ( j = 0 ; j < n ; j++ ) {
                for ( k = j ; k < n ; k++ ) {
                
                    /* Get the binomial number */
                    binom = ((double)fac[k]) / ( fac[k-j] * fac[j] );
    
                    /* Apply the update */
                    if ( k == j )
                        out[n-1-j] *= binom * powa[j];
                    else
                        out[n-1-j] += out[n-1-k] * binom * powb[k-j] * powa[j];
                    }
                }
            } /* End rescale for arbitrary intervals */

        /* Reverse the ordering */
        for ( k = 0 ; k < n / 2  ; k++ ) {
            tmpout = out[k];
            out[k] = out[n-1-k];
            out[n-1-k] = tmpout;
            }
        
        }/* End general case */

    /* Sweet. */
    return fun_err_ok; 
    
    }


/**
 * @brief Plot a #fun with gnuplot.
 *
 * @param f1 The input #fun.
 *
 * @return #fun_err_ok or < 0 on error (see #fun_err).
 */
 
int fun_gnuplot ( struct fun *f1 ) {

    int k, npts = 1000;
    double tk, xk, vk, scl1, scl2;
    FILE *lines, *marks, *pipe;

    /* Routine checks of sanity. */
    if ( f1 == NULL )
        return error(fun_err_null);
    if ( !( f1->flags & fun_flag_init ) )
        return error(fun_err_uninit);

    /* Scaling for the interval [a b] */
    scl1 = (f1->a+f1->b) * 0.5;
    scl2 = (f1->b-f1->a) * 0.5;

    /*** LINES ***/
    /* Open the output file. */
    if ( ( lines = fopen( "fun_gnuplot.dump" , "w" ) ) == NULL ) {
        printf("fun_test: unable to create fun_gnuplot.dump.\n");
        return error(fun_err_gnuplot);
        }
    /* Get the lines */
    for ( k = 0 ; k < npts ; k++ ) {
        tk = (2.0 * k) / (double)(npts - 1) - 1.0;
        xk = scl1 + scl2 * tk;
        fprintf(lines," %.20e %.20e\n", xk , fun_eval( f1 , xk ));
        }     
    /* Close the output file. */
    fclose(lines);

    /*** MARKS ***/
    /* Open the output file. */
    if ( ( marks = fopen( "fun_gnuplot2.dump" , "w" ) ) == NULL ) {
        printf("fun_gnuplot: unable to create fun_gnuplot2.dump.\n");
        return error(fun_err_gnuplot);
        }
    /* Get the marks */
    for ( k = 0 ; k < f1->n ; k++ ) {
        vk = scl1 + scl2 * f1->points[k];
        fprintf(marks," %.20e %.20e\n", vk , f1->vals.real[k] );
        }
    /* Close the output file. */
    fclose(marks);

    /* Fire-up gnuplot. */
    if ( ( pipe = popen( "gnuplot -persist" , "w" ) ) == NULL ) {
        printf("fun_gnuplot: unable to create a pipe to gnuplot.\n");
        return error(fun_err_gnuplot);
        }
    
    /* Pipe the data */
    fprintf( pipe , "set term wxt 1\np 'fun_gnuplot.dump' u 1:2 w lines title \"f\", 'fun_gnuplot2.dump' u 1:2 w points ps 2 title \"\"\n" );

    /* Close the pipe */
    pclose(pipe);

    return fun_err_ok; 
    
    }


/**
 * @brief Add a constant value to a #fun.
 *
 * @param A The input #fun.
 * @param x The constant to add to @c A.
 * @param B A #fun in which to store the result, may also be @c A.
 *
 * @return #fun_err_ok or < 0 on error (see #fun_err).
 */
 
/* TODO: Should crash if A->n == 0! */
 
int fun_add_const ( struct fun *A , double x , struct fun *B ) {

    int k;

    /* Routine checks of sanity. */
    if ( A == NULL || B == NULL)
        return error(fun_err_null);
    if ( !( A->flags & fun_flag_init ) )
        return error(fun_err_uninit);

    /* We have to copy A into B if A != B. */
    if ( A != B ) {
    
        /* Make a copy of A in B. */
        if ( fun_copy( A , B ) < 0 )
            return error(fun_err);
            
        }

    /* Shift the 0th coefficient. */
    B->coeffs.real[0] += x;
    
    /* Shift the vals and re-compute the scale while we're at it. */
    B->scale = 0.0;
    for ( k = 0 ; k < B->n ; k++ ) {
        B->vals.real[k] += x;
        if ( fabs( B->vals.real[k] ) > B->scale )
            B->scale = fabs( B->vals.real[k] );
        }
    
    /* And so to bed... */
    return fun_err_ok;

    }


/**
 * @brief Check if two funs are identical.
 *
 * @param A The first input #fun.
 * @param B the second input #fun.
 *
 * @return 1 if A==B, 0 if A!=B, or < 0 on error (see #fun_err).
 */
 
int fun_isequal ( struct fun *A , struct fun *B ) {

    int k;

    /* Routine checks of sanity. */
    if ( A == NULL || B == NULL)
        return error(fun_err_null);
    if ( !( A->flags & B->flags & fun_flag_init ) )
        return error(fun_err_uninit);

    /* Trivial case. */
    if ( A == B )
        return 1;

    /* Perform some basic checks. */
    if ( (A->n != B->n) || (A->a != B->a) || (A->b != B->b) )
        return 0;

    /* Check values. */
    for ( k = 0 ; k < A->n ; k++ )
        if ( A->vals.real[k] != B->vals.real[k] )
            return 0;        
    
    /* We made it! */
    return 1;

    }


/**
 * @brief Prolong a fun to a different number of Chebyshev points.
 *
 * @param A The #fun to prolong.
 * @param N The desired length.
 * @param B The prolonged #fun. This can
 *      be the same value as @c A.
 *
 * @return #fun_err_ok or < 0 if an error occured.
 */

int fun_prolong ( struct fun *A , unsigned int N , struct fun *B ) {
	
	int j, k, m;
    double *c;

    /* Bad apples? */
    if ( A == NULL || B == NULL)
        return error(fun_err_null);
    if ( !( A->flags & fun_flag_init ) )
        return error(fun_err_uninit);

    /* The difference in lengths. */
	m = N - A->n;
    
	/* Trivial case m == 0 */
	if (m == 0) {
    
        /* Copy A to B if A is not already B. */
		if ( B != A && fun_copy( A , B ) < 0 )
            return error(fun_err);
            
        /* Leave the room, now. */
		return fun_err_ok;
        
		}

    /* Treatment if B != A. */
    if (B != A) {
    
        /* Start by cleaning out B */
        if ( fun_init( B , N ) != fun_err_ok )
            return error(fun_err);
            
        /* Copy the values from A to B */
        B->flags = A->flags;
        B->a = A->a; B->b = A->b;
        B->scale = A->scale;
        
        /* Trivial constant case */
        if ( A->n == 1 ) {
        
            /* Just copy the vals and we're done. */
            B->coeffs.real[0] = A->coeffs.real[0];
            for ( k = 1 ; k < N ; k++ )
                B->vals.real[k] = A->vals.real[0];
                
            }
            
        /* Otherwise, re-create the vals. */
        else {
        
            /* Copy the coefficients from A to B (pad zeros). */
            if ( N >= A->n ){
                memcpy( B->coeffs.real , A->coeffs.real , sizeof(double) * A->n );
                bzero( &(B->coeffs.real[A->n]) , sizeof(double) * (N - A->n) );
                }
            
            /* Restricting - need to consider aliasing */
            else {
                bzero( B->coeffs.real , sizeof(double) * N );
                for ( k = 0 ; k < A->n ; k += 2*(N-2) )
                    B->coeffs.real[0] += A->coeffs.real[k];
                for ( j = 1 ; j < N-1 ; j++ ) {
                    for ( k = j ; k < A->n ; k += 2*(N-2) )
                        B->coeffs.real[j] += A->coeffs.real[k];
                    for ( k = 2*(N-2)-j+2 ; k < A->n ; k += 2*(N-2) )
                        B->coeffs.real[j] += A->coeffs.real[k];
                    }
                for ( k = N-1 ; k < A->n ; k += 2*(N-2) )
                    B->coeffs.real[N-1] += A->coeffs.real[k];
                }
                
            /* Get the vals at the new points. */
            if ( util_chebpolyval( B->coeffs.real , N , B->vals.real ) < 0 )
                return error(fun_err_util);
            
            }
            
        }
        
    /* Treatment if B == A and the function needs to grow. */
    else if ( N > A->size ) {
    
        /* Allocate a new coeffs array and fill it with the old coeffs. */
        if ( posix_memalign( (void **)&c , 16 , sizeof(double) * N ) != 0 )
            return error(fun_err_malloc);
        memcpy( c , A->coeffs.real , sizeof(double) * A->n );
        bzero( &(c[A->n]) , sizeof(double) * (N - A->n) );
        free( B->coeffs.real );
        B->coeffs.real = c;
        
        /* Generate the new points. */
        free( B->points );
        if ( ( B->points = util_chebpts_alloc( N ) ) == NULL )
            return error(fun_err_util);
    
        /* Re-allocate the vals and generate the new values with the FFT. */
        free( B->vals.real );
        if ( ( B->vals.real = util_chebpolyval_alloc( B->coeffs.real , N ) ) == NULL )
            return error(fun_err_util);
            
        /* Set the new length. */
        B->n = N;
        B->size = N;

        }
        
    /* Otherwise, B == A and the function is shrinking. */
    else {
    
        /* Adjust the coefficients. */
        if ( ( c = (double *)alloca( sizeof(double) * N ) ) == NULL )
            return error(fun_err_malloc);
        bzero( c , sizeof(double) * N );
        for ( k = 0 ; k < A->n ; k += 2*(N-2) )
            c[0] += A->coeffs.real[k];
        for ( j = 1 ; j < N-1 ; j++ ) {
            for ( k = j ; k < A->n ; k += 2*(N-2) )
                c[j] += A->coeffs.real[k];
            for ( k = 2*(N-2)-j+2 ; k < A->n ; k += 2*(N-2) )
                c[j] += A->coeffs.real[k];
            }
        for ( k = N-1 ; k < A->n ; k += 2*(N-2) )
            c[N-1] += A->coeffs.real[k];
        memcpy( A->coeffs.real , c , sizeof(double) * N );
                    
        /* Re-create the function values. */
        if ( util_chebpolyval( A->coeffs.real , N , B->vals.real ) < 0 )
            return error(fun_err_util);
            
        /* Re-create the points. */
        if ( util_chebpts( N , B->points ) < 0 )
            return error(fun_err_util);
    
        /* Set the new length. */
        B->n = N;

        }
        
    return fun_err_ok;

    }           

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
 * @brief Redefine the domain of a fun.
 *
 * @param fun The #fun to move.
 * @param a The new left point.
 * @param b The new right point.
 *
 * @return #fun_err_ok or < 0 if an error occured
 */

int fun_newdomain ( struct fun *fun , double a , double b ) {
	
    /* Bad apples? */
    if ( fun == NULL ) {
        return error(fun_err_null);
        }
    if ( !( fun->flags & fun_flag_init ) ) {
        return error(fun_err_uninit);
        }

    /* A quick check. */
    if ( b < a )
        return error(fun_err_domain);

    /* Assign the new domain. */
    fun->a = a;
    fun->b = b;

    /* Well that was easy... */
    return fun_err_ok;
    
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
    if ( !(fun2->flags & fun_flag_init) || fun2->size < fun->n )
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
 * Note that both @c A and @c B must be within the domain of @c fun, 
 * and that the output fun will NOT be 'simplified'. 
 */

int fun_restrict ( struct fun *fun , double A , double B , struct fun *funout ) {
    
    double *x;
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
    for (j = 0 ; j < fun->n ; j++ )
        //fun->points[j] = B * (fun->points[j] - fun->a) * iba + A * (fun->b - fun->points[j]) * iba;
        x[j] = B * 0.5 * (fun->points[j] + 1.0) + A * 0.5 * (1.0 - fun->points[j]);
        
    if (fun == funout) {

        /* Get the restricted values. */
        fun_eval_clenshaw_vec ( fun , x , fun->n ,  funout->vals.real );

        }
    else {
    
        /* Start by cleaning out funout */
        if ( !( funout->flags & fun_flag_init ) || ( funout->size < fun->n ) ) {
            if ( fun_init( funout , fun->n ) != fun_err_ok )
                return error(fun_err);
            }
        else {
            funout->n = fun->n;
            memcpy( funout->points , fun->points , sizeof(double) * fun->n );
            }
                
        /* Pass some variables to the new fun */
        funout->flags = fun->flags;

        /* Get the restricted values. */
        fun_eval_clenshaw_vec ( fun , x , fun->n , funout->vals.real );
        
        }
        
    /* Set the new limits. */
	funout->a = A;
	funout->b = B;

    /* Get the coefficients. */
    if ( util_chebpoly( funout->vals.real  , funout->n , funout->coeffs.real ) < 0 )
        return error(fun_err_util);

    /* Simplify */
// NEED TO PASS A SENSIBLE TOLERANCE
//    if ( fun_simplify( funout , 1e-15 ) < 0 )
//        return error(fun_err);

    /* Re-scale. */
    _fun_rescale( funout );
        
    /* End on a good note. */
    return fun_err_ok;
    
    }


/**
 * @brief Simplify a fun (remove trailing coefficients below tolerance).
 *
 * @param fun The fun the be simplified.
 * @param tol Tolerance
 * @return The new length of the fun or < 0 on error.
 *
 * This wraps _fun_simplify (internal) and adds points and values!
 */

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
    N = _fun_simplify ( fun , tol );      

    /* Store newdata simplify was successful. */
    if ( N < fun->n && N > 0) {
        
        /* Re-create the fun in the already-allocated memory. */
        util_chebpolyval( fun->coeffs.real , N , fun->vals.real );
        util_chebpts( N , fun->points );

        }

    /* Sweet. */
    return N;

    }


/**
 * @brief The 2-norm difference of two functions
 *
 * @param A The first #fun.
 * @param B The second #fun.
 * @return The 2-norm of their difference.
 */

double fun_err_norm2 ( struct fun *A , struct fun *B ) {

    double norm;
    struct fun C = FUN_EMPTY;

    /* Check inputs. */
    if ( A == NULL || B == NULL )
        return error(fun_err_null);
    if ( !( A->flags & B->flags & fun_flag_init ) )
        return error(fun_err_uninit);

    /* If A == B this must be zero */
    if ( A == B )
        return 0.0; 

    /* Subtract A and B */
    fun_madd ( A , 1.0 , B , -1.0 , &C );  

    /* Take the 2-norm of the difference. */
    norm = fun_norm2 ( &C ); 

    /* Clean C */
    fun_clean( &C );

    /* Sweet. */
    return norm;

    }


/**
 * @brief The inf-norm difference of two functions
 *
 * @param A The first #fun.
 * @param B The second #fun.
 * @return The inf-norm of their difference.
 */

double fun_err_norm_inf ( struct fun *A , struct fun *B ) {

    double norm;
    struct fun C = FUN_EMPTY;

    /* Check inputs. */
    if ( A == NULL || B == NULL )
        return error(fun_err_null);
    if ( !( A->flags & B->flags & fun_flag_init ) )
        return error(fun_err_uninit);

    /* If A == B this must be zero */
    if ( A == B )
        return 0.0;   

    /* Subtract A and B */
    fun_madd ( A , 1.0 , B , -1.0 , &C );  

    /* Take the inf-norm of the difference. */
    norm = fun_norm_inf ( &C ); 

    /* Clean C */
    fun_clean( &C );

    /* Sweet. */
    return norm;

    }


    
/**
 * @brief Simplify a fun (remove trailing coefficients below tolerance).
 *      Note that values or points are not assigned, only coefficients.
 *
 * @param fun The fun the be simplified.
 * @param tol Tolerance
 * @return The new length of the fun or < 0 on error.
 *
 * We simply call util_simplify.  Should be used internally only!
 */

int _fun_simplify ( struct fun *fun , double tol ) {
    
    unsigned int N;

    /* Call util_simplify */
    if ( ( N = util_simplify( fun->points , fun->vals.real , fun->coeffs.real , fun->n , fun->b-fun->a , fun->scale , tol ) ) < 0 )
        return error(fun_err_util);

    /* Set the new size. */
    if ( N < fun->n )
        fun->n = N;  
    
    /* Sweet. */
    return N;
    
    }


/**
 * @brief Compute the sorted roots of a fun.
 *
 * @param fun The #fun to find roots of.
 * @param roots A pointer to an array of double of length at least
 *      equal to the length of @c fun - 1. This is where the roots will be stored.
 *
 * @return The number of roots found or < 0 if an error occured.
 * 
 * @sa fun_roots 
 */
 
int fun_roots_sort( struct fun *fun , double *roots ) {
	
	int k, nroots;
	double m, h;
    
    /* Check for the usual suspects. */
    if ( fun == NULL || roots == NULL )
        return error(fun_err_null);
    if ( !( fun->flags & fun_flag_init ) )
        return error(fun_err_uninit);

    /* Call fun_roots_unit (with the sort flag). */
	if ( ( nroots = fun_roots_unit( fun , roots , 1 ) ) < 0 )
        return error(fun_err);
        
    /* Get m and h for the interval. */
    m = 0.5 * (fun->a + fun->b);
    h = 0.5 * (fun->b - fun->a);

    /* Scale the roots to the correct interval if needed. */
	if ( fun->a != -1.0 || fun->b != 1.0 )
	    for ( k = 0 ; k < nroots ; k++ )
		    roots[k] = m + roots[k] * h;

	return nroots;
    
	}

/**
 * @brief Compute the roots of a fun.
 *
 * @param fun The #fun to find roots of.
 * @param roots A pointer to an array of double of length at least
 *      equal to the length of @c fun - 1. This is where the roots will be stored.
 *
 * @return The number of roots found or < 0 if an error occured.
 * 
 * @sa fun_roots 
 *
 * Note, if the roots need to be sorted, it is more efficient to 
 * call fun_roots_sort, as they are they sorted during calculation.
 */
 
int fun_roots( struct fun *fun , double *roots ) {
	
	int k, nroots;
	double m, h;
    
    /* Check for the usual suspects. */
    if ( fun == NULL || roots == NULL )
        return error(fun_err_null);
    if ( !( fun->flags & fun_flag_init ) )
        return error(fun_err_uninit);

    /* Call fun_roots_unit (without the sort flag). */
	if ( ( nroots = fun_roots_unit( fun , roots , 0 ) ) < 0 )
        return error(fun_err);
        
    /* Get m and h for the interval. */
    m = 0.5 * (fun->a + fun->b);
    h = 0.5 * (fun->b - fun->a);

    /* Scale the roots to the correct interval if needed. */
	if ( fun->a != -1.0 || fun->b != 1.0 )
	    for ( k = 0 ; k < nroots ; k++ )
		    roots[k] = m + roots[k] * h;

	return nroots;
    
	}


/**
 * @brief Compute the roots of a fun assuming it is on the unit interval.
 *
 * @param fun The #fun to find roots of.
 * @param roots A pointer to an array of double of length at least
 *      equal to the length of @c fun - 1. This is where the roots will be stored.
 * @param sort Determines whether the roots should be sorted or not. 
 *      sort == 1 will sort, sort == 0 will not.
 *
 * @return The number of roots found or < 0 if an error occurred.
 */

/* TODO : Implement a quicksort for sorting the roots.
 * TODO: Changed the minimum nr. of coefficients for the eigenvalue solver
 *          to N <= 20 until we get a decent lapack to link against.
 */
 
int fun_roots_unit ( struct fun *fun , double *roots , int sort ) {

    const int split = 20;
    double c = -0.004849834917525, tail_max, temp, *v, hscl;
    int nroots;
    int j, k;
    static double *Tleft = NULL, *Tright = NULL;
    static int sizeT = 0, skipT;
    fftw_plan plan;
    fftw_r2r_kind kind = FFTW_REDFT00;
    
    
    /* Recursive routine that works only on the coefficients. */
    int fun_roots_unit_rec ( double *coeffs , unsigned int N , double h , double *roots , int sort ) {
    
        int j, k;
	    int ilo, ihi, ldh, ldz, lwork, ok;
        int Nl, Nr, nrootsL = 0, nrootsR = 0, nroots = 0;
	    char job = 'E', compz = 'N';
        double tol, cN, w, temp;
        double *cleft, *cright, *rr, *ri, *work, *A, z, *v;
        
        /* Is the degree small enough for a colleague matrix approach? */
        if ( N <= split ) {
        
            /* Set the horizontal tolerance for roots inside the interval. */
            tol = 100.0 * h * DBL_EPSILON;
            
            /* Adjust N and compute scaling cN. */
        	for ( N -= 1 ; fabs(fun->coeffs.real[N]) < 1e-14 * fun->scale && N > 0 ; N-- );
            cN = -0.5 / coeffs[N];

            /* Trivial case of N==0. */
            if ( N == 0 ) {
                if ( coeffs[0] == 0.0 ) {
                    roots[0] = 0.0;
                    return 1;
                    }
                else
                    return 0;
                }
                
            /* Trivial linear case (N==1). */
            if ( N == 1 ) {
                roots[0] = -coeffs[0] / coeffs[1];
                if ( roots[0] > -(1+tol) && roots[0] < (1+tol) )
                    return 1;
                else
                    return 0;
                }
        
            /* Initialize the matrix A. */
            if ( ( A = (double *)alloca( sizeof(double) * (N*N) + 16 ) ) == NULL )
                return error(fun_err_malloc);
            A = (double *)( (((size_t)A) + 15 ) & ~15 );
            bzero( A , sizeof(double) * (N*N) );

            /* Also initialize and shift the output and work vectors */
            if ( ( rr = (double *)alloca( sizeof(double) * N + 16 ) ) == NULL || 
			     ( ri = (double *)alloca( sizeof(double) * N + 16 ) ) == NULL ||
		         ( work = (double *)alloca( sizeof(double) * N + 16 ) ) == NULL )
                return error(fun_err_malloc);
            rr = (double *)( (((size_t)rr) + 15 ) & ~15 );
            ri = (double *)( (((size_t)ri) + 15 ) & ~15 );
            work = (double *)( (((size_t)work) + 15 ) & ~15 );
            
            /* Construct the colleague matrix. */
            /* Do N = 2 case by hand. */
            if ( N == 2 ) { 
            
                A[0] = cN * coeffs[1];
                A[1] = 1.0;
                A[2] = cN * coeffs[0] + 0.5;

                }
                
            /* General case. */ 
            else {
            
                /* Assign the coefficients to the first row */
                for (j = 0 ; j < N ; j++)
                    A[N*(N-1-j)] = cN * coeffs[j];
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
                }

            /* Call to LAPACK to solve eigenvalue problem. */
		    ilo = 1; ihi = N; ldh = N; ldz = 1; lwork = N;
		    dhseqr_( &job, &compz, &N , &ilo, &ihi, A, &ldh, rr, ri, &z, &ldz, work, &lwork, &ok);
		    if ( ok < 0 )
			    return error(fun_err_lapack);

            /* Count the number of valid roots and store them. */
            /* Note, imag tol in Matlab is 0.5*tol ... */
		    for (j = ok ; j < N ; j++)
			    if ( fabs(ri[j]) < tol && rr[j] >= -(1.0+tol) && rr[j] <= (1.0+tol) )
                    roots[nroots++] = rr[j];  

            /* Sort if required. */
            if ( sort == 1 )
                for ( j = 0 ; j < nroots ; j++ )
                    for ( k = j + 1 ; k < nroots ; k++ )
                        if ( roots[j] > roots[k] ) {
                            temp = roots[j];
                            roots[j] = roots[k];
                            roots[k] = temp;
                            }
                              
            }       
            
        /* Otherwise, split and recurse. */
        else {
        
            /* Allocate arrays for cleft and cright. */
            if ( ( cleft = (double *)alloca( sizeof(double) * N + 16 ) ) == NULL ||
                 ( cright = (double *)alloca( sizeof(double) * N + 16 ) ) == NULL )
                return error(fun_err_malloc);
            cleft = (double *)( (((size_t)cleft) + 15 ) & ~15 );
            cright = (double *)( (((size_t)cright) + 15 ) & ~15 );
            bzero( cleft , sizeof(double) * N );
            bzero( cright , sizeof(double) * N );
            
            /* If N is small enough, compute both halves with Tleft and Tright. */
            if ( N <= sizeT ) {
        
                /* Multiply coeffs with Tleft and Tright to get cleft and cright. */
                for ( k = 0 ; k < N ; k++ )
                    for ( j = 0 ; j <= k ; j++ ) {
                        cleft[j] += coeffs[k] * Tleft[k*skipT+j];
                        cright[j] += coeffs[k] * Tright[k*skipT+j];
                        }
                        
                }

            /* Otherwise, restrict both halves "by hand". */
            else {

                /* Allocate the temporary vector for the function values. */
                if ( ( v = (double *)alloca( sizeof(double) * N + 16 ) ) == NULL )
                    return error(fun_err_malloc);
                v = (double *)( (((size_t)v) + 15 ) & ~15 );

                /* Get the coefficients on the left. */
                if ( util_chebptsAB( N , -1.0 , c , cleft ) < 0 ||
                     util_clenshaw_vec( coeffs , N , cleft , N , v ) < 0 ||
                     util_chebpoly( v , N , cleft ) < 0 )
                    return error(fun_err_util);

                /* Get the coefficients on the right. */
                if ( util_chebptsAB( N , c , 1.0 , cright ) < 0 ||
                     util_clenshaw_vec( coeffs , N , cright , N , v ) < 0 ||
                     util_chebpoly( v , N , cright ) < 0 )
                    return error(fun_err_util);

                }
            
            /* Scale the lengths of cleft and cright. */
            for ( Nl = N ; Nl > 0 && fabs(cleft[Nl-1]) < tail_max ; Nl-- );
            for ( Nr = N ; Nr > 0 && fabs(cright[Nr-1]) < tail_max ; Nr-- );

            /* Recurse on the left if Nl > 0. */
            if ( Nl > 0 ) {
                if ( ( nrootsL = fun_roots_unit_rec( cleft , Nl , 2.0*hscl/(c+1.0) , roots , sort ) ) < 0 )
                    return error(fun_err);
                w = 0.5 * (1.0 + c);
                for ( nroots = 0 ; nroots < nrootsL ; nroots++ )
                    roots[nroots] = -1.0 + w * (roots[nroots] + 1.0);
                }
            /* Recurse on the right if Nr > 0. */
            if ( Nr > 0 ) {
                if ( ( nrootsR = fun_roots_unit_rec( cright , Nr , 2.0*hscl/(1.0-c) , &(roots[nrootsL] ) , sort ) ) < 0 )
                    return error(fun_err);
                w = 0.5 * (1.0 - c);
                for ( nroots = nrootsL ; nroots < nrootsL + nrootsR ; nroots++ )
                    roots[nroots] = c + w * (roots[nroots] + 1.0);
                }
            
            }
    
        /* Bail in a good way. */
        return nroots;
                
        }
    
    
    /* Check inputs. */
    if ( fun == NULL )
        return error(fun_err_null);
        
    /* Deal with this stupid case. */
    if ( fun->n == 0 ) 
        return 0;
        
    /* Check if the coefficient matrices Tleft and Tright have been
        computed and are the right size. */
    if ( fun->n > split && Tleft == NULL ) {
    
        /* printf("fun_roots_unit: building transformation matrices... "); fflush(stdout); */
    
        /* Free any old instance of Tleft and Tright. */
        if ( Tleft != NULL )
            free(Tleft);
        if ( Tright != NULL )
            free(Tright);
            
        /* Set the size. */
        /* for ( sizeT = 2 ; sizeT < fun->n ; sizeT = 2*sizeT - 1 ); */
        sizeT = 513;
        
        /* Set the stride for the Tleft and Tright matrices such as to
            preserve alignment. */
        skipT = (sizeT + 3) & ~3;
            
        /* Allocate the new matrices. */
        if ( posix_memalign( (void**)(&Tleft) , 16 , sizeof(double) * sizeT * skipT ) != 0 ||
             posix_memalign( (void**)(&Tright) , 16 , sizeof(double) * sizeT * skipT ) != 0 )
            return error(fun_err_malloc);
        if ( posix_memalign( (void**)(&v) , 16 , sizeof(double) * sizeT * skipT ) != 0 )
            return error(fun_err_malloc);
            
        /* Compute the left matrix. */
        for ( k = 0 ; k < sizeT ; k++ )
            v[k] = 1.0;
        if ( util_chebptsAB( sizeT , -1.0 , c , &(v[skipT]) ) < 0 )
            return error(fun_err_util);
        for ( k = 2 ; k < sizeT ; k++ )
            for ( j = 0 ; j < sizeT ; j++ )
                v[k*skipT+j] = 2 * v[skipT+j] * v[(k-1)*skipT+j] - v[(k-2)*skipT+j];
                
        /* Transform v to Tleft. */
        plan = fftw_plan_many_r2r( 1 , &sizeT , sizeT , v , &sizeT , 1 , skipT , Tleft , &sizeT , 1 , skipT , &kind , FFTW_ESTIMATE );
        fftw_execute(plan);
        
        /* Clean up the edge coefficients. */
        for ( k = 0 ; k < sizeT ; k++ )
            Tleft[k*skipT] *= 0.5;
        Tleft[(skipT+1)*(sizeT-1)] *= 0.5;
        
        /* Compute the right matrix. */
        if ( util_chebptsAB( sizeT , c , 1.0 , &(v[skipT]) ) < 0 )
            return error(fun_err_util);
        for ( k = 2 ; k < sizeT ; k++ )
            for ( j = 0 ; j < sizeT ; j++ )
                v[k*skipT+j] = 2 * v[skipT+j] * v[(k-1)*skipT+j] - v[(k-2)*skipT+j];
                
        /* Transform v to Tleft. */
        /* plan = fftw_plan_many_r2r( 1 , &sizeT , sizeT , v , &sizeT , 1 , skipT , Tright , &sizeT , 1 , skipT , &kind , FFTW_ESTIMATE ); */
        fftw_execute_r2r( plan , v , Tright );
        fftw_destroy_plan(plan);
        
        /* Clean up the edge coefficients. */
        for ( k = 0 ; k < sizeT ; k++ )
            Tright[k*skipT] *= 0.5;
        Tright[(skipT+1)*(sizeT-1)] *= 0.5;
        
        /* Scale both matrices. */
        temp = 1.0 / (sizeT - 1);
        for ( k = 0 ; k < skipT * sizeT ; k++ ) {
            Tleft[k] *= temp;
            Tright[k] *= temp;
            }
        
        /* Clean up v. */
        free(v);
        
        /* printf("done.\n"); fflush(stdout); */
        
        }
        
    /* Find out what the coefficient cutoff tail_max is for the global fun. */
    tail_max = 0.0;
    for ( k = 1 ; k < fun->n ; k++ ) {
        temp = fabs( fun->vals.real[k] - fun->vals.real[k-1] );
        if ( fabs( fun->points[k] - fun->points[k-1] ) > 2 * DBL_EPSILON )
            temp /= fabs( fun->points[k] - fun->points[k-1] );
        else
            temp /= 2.0 * DBL_EPSILON;
        if ( temp > tail_max )
            tail_max = temp;
        }
    if ( tail_max > 1.0e12 )
        tail_max = 1.0e12;
    tail_max *= DBL_EPSILON;

    /* Get the horizontal scale (following Matlab defn). */
    hscl = fabs(fun->a);
    if ( fabs(fun->b) > hscl )
        hscl = fabs(fun->b);

    /* Call the recursion. */
    if ( ( nroots = fun_roots_unit_rec( fun->coeffs.real , fun->n , 2.0*hscl/(fun->b - fun->a) , roots , sort ) ) < 0 )
        return error(fun_err);
        
    /* To the pub!. */
    return nroots;
    
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
 * @param out A pointer to a @c FILE to which output will be directed.
 *      If in doubt, just use @c stdout.
 *
 * @return #fun_err_ok or < 0 if an error occured.
 */
 
int fun_display ( struct fun *fun , FILE *out ) {

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
        fprintf(out,"fun.points[%i]=%16.16e \tfun.vals[%i]=%16.16e \tfun.coeffs[%i]=%16.16e\n",
            k, xk, k, fun->vals.real[k], k, fun->coeffs.real[k]);
        }
    fprintf(out,"\n");
    
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
        
    /* Allocate the data inside the fun to the correct size, if needed. */
    if ( !( fun->flags & fun_flag_init ) || fun->size < N )
        _fun_alloc( fun , N );
    else
        fun->n = N;

    /* Set entries to zero. Is this needed?
    bzero( fun->vals.real , sizeof(double) * N );
    bzero( fun->coeffs.real , sizeof(double) * N ); */
    
    /* Get the Chebyshev nodes. */
    if ( util_chebpts( N , fun->points ) < 0 )
        return error(fun_err_util);
    
    /* Set the init flag. */    
    fun->flags = fun_flag_init;
    
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

    /* Clean the fun, just to be safe. */
    if ( fun_clean( fun ) < 0 )
        return error(fun_err);

    /* Allocate the data inside the fun to the correct size. */
    if ( posix_memalign( (void **)&(fun->vals.real) , 16 , sizeof(double) * N ) != 0 ||
         posix_memalign( (void **)&(fun->coeffs.real) , 16 , sizeof(double) * N ) != 0 ||
         ( fun->points = (double *)malloc( sizeof(double) * N ) ) == NULL )
        return error(fun_err_malloc);

    /* Write the data to the fun. */
    fun->n = N;
    fun->size = N;
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
    if ( N == 1 )
        fun->coeffs.real[0] = vals[0];
    else {
        if ( util_chebpoly( vals , N , fun->coeffs.real ) < 0 )
            return error(fun_err_util);
        }
        
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
    if ( N == 1 )
        fun->vals.real[0] = fun->coeffs.real[0];
    else {
        if ( util_chebpolyval( coeffs , N , fun->vals.real ) < 0 )
        return error(fun_err_util);
        }
        
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
    if ( !( fp->flags & fun_flag_init ) || fp->size < f->n - 1 )
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
    
    /* Apply the recurrence for the new coefficents */
    n = f->n;
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
        if ( !( C->flags & fun_flag_init ) || C->size < N )
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
        
        /* Is the output fun large enough? */
        if ( C->size < N ) {
        
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
                
            /* Set the new size. */
            C->size = N;
                
            }
            
        /* Otherwise, just need to re-create the points. */
        else {
            if ( util_chebpts( N , C->points ) < 0 )
                return error(fun_err_util);
            }
                
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
 * @return the value of the integral or @c NAN if an error was encountered.
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
 * @brief Indefinite integral of a fun.
 *
 * @param A The #fun to be integrated.
 * @param B The indefinite integral.
 * @return #fun_err_ok or < 0 if an error occurs.
 */

int fun_indef_integral ( struct fun *A , struct fun *B ) {

    int k, n, sgn;
    double *c, scl;

    /* Check for the usual problems... */
    if ( A == NULL || B == NULL )
        return error(fun_err_null);
    if ( !( A->flags & fun_flag_init ) )
        return error(fun_err_uninit);
        
    /* Get the size */
    n = A->n;

    /* Are A and B one and the same? */
    if ( A == B ) {
    
        /* Allocate a temporary array for the new coefficients. */
        if ( ( c = (double *)alloca( sizeof(double) * A->n ) ) == NULL )
            return error(fun_err_malloc);
        memcpy( c , A->coeffs.real , sizeof(double) * A->n );
            
        /* Re-allocate the data in A. */
        if ( fun_init( B , A->n + 1 ) < 0 )
            return error(fun_err);
        
        }
        
    /* A and B are not the same fun. */
    else {
    
        /* Start by cleaning out B if needed. */
        if ( !( B->flags & fun_flag_init ) || B->size < A->n + 1 ) {
            if ( fun_init( B , A->n + 1 ) < 0 )
                return error(fun_err);
            }
        else {
            B->n = A->n + 1;
            if ( util_chebpts( B->n , B->points ) < 0 )
                return error(fun_err_util);
            }
            
        /* Set the bounds of B. */
        B->a = A->a;
        B->b = A->b;

        /* Set c. */
        c = A->coeffs.real;
        
        }

    /* Compute the coefficients of the integral. */
    if ( n == 1 ) {
        B->coeffs.real[1] = c[0];
        }
    else if (n == 2) {
        B->coeffs.real[1] = c[0];
        B->coeffs.real[2] = 0.25 * c[1];
        }
    else {
        /* Recurrence */
        B->coeffs.real[1] = c[0] - 0.5*c[2];
        for ( k = 0 ; k < n-3 ; k++ )
            B->coeffs.real[k+2] = 0.5 * ( c[k+1] - c[k+3] ) / (k + 2);
        B->coeffs.real[n-1] = 0.5 * c[n-2] / (n - 1);
        B->coeffs.real[n] = 0.5 * c[n-1] / n;
        }

    /* Get the constant term */
    sgn = 1;
    B->coeffs.real[0] = 0.0;
    for ( k = 1 ; k < n+1 ; k++ ) {
        B->coeffs.real[0] += sgn * B->coeffs.real[k];
        sgn = -sgn;
        }

    /* Scaling for non-standard intervals */
    scl = 0.5 * (A->b - A->a);
    if ( scl != 1.0 )
        for ( k = 0 ; k < n+1 ; k++ )
            B->coeffs.real[k] *= scl;

    /* Compute the new values. */
    if ( util_chebpolyval( B->coeffs.real , B->n , B->vals.real ) < 0 )
        return error(fun_err_util);

    /* Get the scale of the new fun. */
    _fun_rescale( B );

    /* We made it! */
    return fun_err_ok;
        
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
    
        /* Make sure B is the right size. */
        if ( !( B->flags & fun_flag_init ) || ( B->size < A->n ) )
            if ( fun_init( B , A->n ) != fun_err_ok )
                return error(fun_err);
            
        /* Copy the values from A to B */
        B->flags = A->flags;
        B->a = A->a; B->b = A->b;
        B->n = A->n;
        B->scale = fabs(w) * A->scale;
        
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
 * @brief Subtract two funs.
 * 
 * @param A Pointer to a #fun.
 * @param B Pointer to a #fun.
 * @param C Pointer to a #fun in which to store the result of
 *      @c A - @c B. @c C can be equal to @c A or @c B. In any
 *      case, @c C should be either initialized or clean (see #fun_clean).
 * @return #fun_err_ok or < 0 if an error occurs.
 *
 * This function is just a wrapper for #fun_madd.
 */
 
int fun_sub ( struct fun *A , struct fun *B , struct fun *C ) {

    return fun_madd( A , 1.0 , B , -1.0 , C );

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
 
/* TODO: Case where C == b and b->size >= a->n! */
 
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
        if ( !( C->flags & fun_flag_init ) || C->size < a->n )
            if ( fun_init( C , a->n ) < 0 )
                return error(fun_err);
            
        /* Set the domain. */
        C->n = a->n;
        C->a = a->a; C->b = a->b;
            
        /* Copy the points from a. */
        memcpy( C->points , a->points , sizeof(double) * C->n );

        /* Merge the coefficients from a and b. */
        for ( k = 0 ; k < b->n ; k++ )
            C->coeffs.real[k] = wa * a->coeffs.real[k] + wb * b->coeffs.real[k];
        for ( k = b->n ; k < a->n ; k++ )
            C->coeffs.real[k] = wa * a->coeffs.real[k];
                
        /* Are a and b the same length? */
        if ( a->n == b->n ) {
        
            /* Merge the vals from a and b. */
            for ( k = 0 ; k < C->n ; k++ )
                C->vals.real[k] = wa * a->vals.real[k] + wb * b->vals.real[k];
                
            }
                
        /* Different lengths in a and b. */
        else {
        
            /* Extract the vals from these coeffs. */
            if ( util_chebpolyval( C->coeffs.real , C->n , C->vals.real ) < 0 )
                return fun_err_util;
                
            }
                
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
            C->size = a->n;
    
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
 * @brief Evaluates the given fun at the given node @c x.
 *
 * @param fun The #fun to be evaluated.
 * @param x The position at which to evaluate @c fun.
 * @return The #fun evaluated at @c x or @c NaN if an error occured
 *      (see #fun_err).
 *
 * This function delegates the input to either #fun_eval_bary or
 * #fun_eval_clenshaw depending on the flags in #chebopts_opts or
 * what values are present in the #fun.
 *
 * @sa fun_eval_bary, fun_eval_vec, fun_eval_clenshaw_vec
 */
 
double fun_eval ( struct fun *fun , double x ) {

    /* Check for nonsense. */
    if ( fun == NULL ) {
        error(fun_err_null);
        return NAN;
        }
    if ( !( fun->flags & fun_flag_init ) ) {
        error(fun_err_uninit);
        return NAN;
        }
        
    /* Re-direct according to what's specified in chebopts. */
    if ( ( chebopts_opts->flags & chebopts_flag_evalbary ) ||
         ( fun->coeffs.real == NULL ) )
        return fun_eval_bary( fun , x );
    else
        return fun_eval_clenshaw( fun , x );

    }
    

/**
 * @brief Evaluates the given fun at the given vector of nodes @c x.
 *
 * @param fun The #fun to be evaluated.
 * @param x A pointer to an array of double values containing the positions
 *      at which to evaluate @c fun.
 * @param m The number of positions in @c x.
 * @param out A pointer to an array of doubles of length @c m in which the
 *      values of @c fun at @c x will be stored.
 * @return #fun_err_ok or < 0 if an error occured.
 *
 * This function delegates the input to either #fun_eval_bary_vec or
 * #fun_eval_clenshaw_vec depending on the flags in #chebopts_opts or
 * what values are present in the #fun.
 *
 * @sa fun_eval_clenshaw, fun_eval, fun_eval_bary, fun_eval_bary_vec, fun_eval_clenshaw_vec
 */
 
int fun_eval_vec ( struct fun *fun , const double *x , unsigned int m , double *out ) {
    
    /* Check for nonsense. */
    if ( fun == NULL || x == NULL || out == NULL )
        return error(fun_err_null);
    if ( !( fun->flags & fun_flag_init ) )
        return error(fun_err_uninit);
        
    /* Re-direct according to what's specified in chebopts. */
    if ( ( chebopts_opts->flags & chebopts_flag_evalbary ) ||
         ( fun->coeffs.real == NULL ) )
        return fun_eval_bary_vec( fun , x , m , out );
    else
        return fun_eval_clenshaw_vec( fun , x , m , out );

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

int fun_eval_clenshaw_vec ( struct fun *fun , const double *x , unsigned int m , double *out ) {

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
 * @sa fun_eval_clenshaw, fun_eval, fun_eval_clenshaw_vec, fun_eval_bary_vec
 */
 
double fun_eval_bary ( struct fun *fun , double x ) {

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
 * @sa fun_eval_clenshaw, fun_eval_bary, fun_eval_vec, fun_eval_clenshaw_vec
 */
 
int fun_eval_bary_vec ( struct fun *fun , const double *x , unsigned int m , double *out ) {

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
	    fun->vals.real[j] = (*fx)( b05 * (fun->points[j] + 1.0) + a05 * (1.0 - fun->points[j]) , data );

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
 fun_create( &f , &myfun , -1.0 , 1.0 , omega );
        @endcode
 * @param a
 * @param b The left and right boundaries of the interval over which
 *      @a f is to be approximated.
 * @param data An pointer to additional data that will be passed
 *      to @a fx at every call.
 * @return #fun_err_ok or < 0 on error.
 *
 * This function only wraps the function @a fx into a vector-valued
 * function an calls #fun_create_vec.
 *
 * @sa fun_create_vec
 */
 
int fun_create ( struct fun *fun , double (*fx)( double x , void * ) , double a , double b , void *data ) {

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
    return fun_create_vec( fun , &fx_wrapper , a , b , data );

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
 * @param data An pointer to additional data that will be passed
 *      to @a fx at every call.
 * @return #fun_err_ok or < 0 on error.
 *
 * @sa fun_create
 */
 
int fun_create_vec ( struct fun *fun , int (*fx)( const double * , unsigned int , double * , void * ) , double a , double b , void *data ) {

    double *x, *v, *coeffs, scale = 0.0;
    double m = (a + b) * 0.5, h= (b - a) * 0.5;
    unsigned int N;
    static double *xi = NULL;
    static int nr_xi = 0;
    int k, stride, N_new;
    
    /* Check inputs. */
    if ( fun == NULL || fx == NULL )
        return error(fun_err_null);
        
    /* Initialize the vectors x, v and coeffs to the maximum length.
       Allocation is done on the stack, so we don't need to worry about
       releasing these if the routine fails anywhere underway.
       Note that v and coeffs are padded with 16 bytes so that they can be
       shifted and aligned to 16-byte boundaries. */
    if ( ( x = (double *)alloca( sizeof(double) * (chebopts_opts->maxdegree + 1) ) ) == NULL ||
         ( v = (double *)alloca( sizeof(double) * (chebopts_opts->maxdegree + 1) + 16 ) ) == NULL ||
         ( coeffs = (double *)alloca( sizeof(double) * (chebopts_opts->maxdegree + 1) + 16 ) ) == NULL )
        return error(fun_err_malloc);
        
    /* Shift v and coeffs such that they are 16-byte aligned. */
    v = (double *)( (((size_t)v) + 15 ) & ~15 );
    coeffs = (double *)( (((size_t)coeffs) + 15 ) & ~15 );
    
    /* Init the length N. */
    N = chebopts_opts->minsamples;
        
    /* Make sure the nodes have been pre-allocated. */
    if ( !(chebopts_opts->flags & chebopts_flag_resampling) &&
         ( nr_xi < chebopts_opts->maxdegree + 1 || (nr_xi - 1) % (N - 1) != 0 ) ) {
    
        /* Clean up old nodes. */
        if ( xi != NULL )
            free(xi);
            
        /* Set the nr of nodes. */
        for ( nr_xi = chebopts_opts->minsamples ; nr_xi < chebopts_opts->maxdegree + 1 ; nr_xi = 2*nr_xi - 1 );
            
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
        if ( chebopts_opts->flags & chebopts_flag_resampling ) {
        
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
            if ( N == chebopts_opts->minsamples ) {
            
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
        if ( ( N_new = util_simplify( x , v , coeffs , N , 2*h , scale , chebopts_opts->eps ) ) < 0 )
            return fun_err_util;
            
        /* TODO: Sampletest? */
            
        /* Did this converge? */
        if ( N_new < N ) {
            N = N_new;
            break;
            }
            
        /* No convergence, adjust N. */
        if ( ( chebopts_opts->flags & chebopts_flag_resampling ) && 
             ( N < 64 ) )
            N = (int)( ( N - 1 ) * M_SQRT2 + 1 ) | 1;
        else
            N = 2 * N - 1;
        
        } /* Main loop. */
        
        
    /* Allocate the data inside the fun to the correct size. */
    if ( !( fun->flags & fun_flag_init ) || fun->n < N ) {
        if ( _fun_alloc( fun , N ) < 0 )
            return error(fun_err);
        }
    else
        fun->n = N;
    
    /* Write the data to the fun. */
    memcpy( fun->vals.real , v , sizeof(double) * N );
    memcpy( fun->coeffs.real , coeffs , sizeof(double) * N );
    memcpy( fun->points , x , sizeof(double) * N );
    fun->scale = scale;
    fun->a = a;
    fun->b = b;
    
    /* If nothing bad happened until here, we're done! */
    return fun_err_ok;

    }
