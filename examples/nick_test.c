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
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>

/* Local includes. */
#include "chebopts.h"
#include "util.h"
#include "fun.h"


/**
 * @brief A simple scalar function.
 */
 
double myfun ( double x , void *data ) {

    return sin( 4.0 * M_PI * x );

    }


/**
 * @brief A simple scalar function with a parameter.
 */
 
double myfun_param ( double x , void *data ) {

    double omega = *((double *)data);

    return sin( omega * M_PI * x );

    }


/**
 * @brief A simple vectorized function.
 */
 
int myfun_vec ( const double *x , unsigned int N , double *out , void *data ) {
    
    int k;

    for ( k = 0 ; k < N ; k++ )
        out[k] = sin( 4.0 * x[k] + 1.0);

    return 0;

    }


/**
 * @brief Derivative of the above.
 */

int myfun_vec_prime ( const double *x , unsigned int N , double *out , void *data ) {
    
    int k;

    for ( k = 0 ; k < N ; k++ )
        out[k] = 4.0 * cos( 4.0 * x[k] + 1.0);

    return 0;

    }


/**
 * @brief A simple vectorized function with a parameter.
 */
 
int myfun_vec_param ( const double *x , unsigned int N , double *out , void *data ) {
    
    int k;
    double omega = *((double *)data);

    for ( k = 0 ; k < N ; k++ )
        out[k] = sin( omega * M_PI * x[k] );

    return 0;

    }


/**
 * @brief Create a simple #fun.
 */
 
int main ( int argc , char *argv[] ) {

    struct fun f1 = FUN_EMPTY, f4 = FUN_EMPTY, fp = FUN_EMPTY, fp2 = FUN_EMPTY;
    struct chebopts opts;
    int k, res;
    
    /* Get a copy of the default options. */
    memcpy( &opts , &chebopts_default , sizeof(struct chebopts) );
    /* opts.flags |= chebopts_flag_resampling; */
    
    
    /* Initialize the fun f1 (vector real). */
    if ( ( res = fun_create_vec( &f1 , &myfun_vec , -1.0 , 1.0 , &opts , NULL ) ) < 0 )
        printf("fun_test: fun_create_vec failed with fun_err=%i (%s).\n",
            fun_err, fun_err_msg[-fun_err]);

    /* Test 2 norm */
    printf("||f1||_2 = %e\n", fun_norm2( &f1 ));

    /* Test differentiation */
    fun_init( &fp , f1.n-1 );
    printf("fun_test: funp.n=%u\n",fp.n);
    fun_diff( &f1, &fp );
    for ( k = 0 ; k < fp.n ; k++ )
        printf("funp_test: fun.points[%i]=%e \tfun.vals[%i]=%e \tfun.coeffs[%i]=%e\n",
            k, fp.points[k], k, fp.vals.real[k], k, fp.coeffs.real[k]);
    fun_create_vec( &fp2 , &myfun_vec_prime , -1.0 , 1.0 , &opts , NULL );
    fun_scale ( &fp2 , -1.0 , &fp2);
    fun_add ( &fp , &fp2 , &fp2 );
    printf("||fp-fp2||_2 = %e\n", fun_norm2( &fp2 ));
    fun_clean( &fp );
    fun_clean( &fp2 );

    /* Test construction from vals */
    fun_create_vals( &f4 , f1.vals.real , -1.0 , 1.0 , f1.n );
//    fun_display( &f4 );
    fun_scale ( &f4 , -1.0 , &f4);
    fun_add ( &f1 , &f4 , &f4 );
    printf("||f-f4||_2 = %e\n", fun_norm2( &f4 ));
    fun_clean( &f4 );

    /* Test construction from coeffs */
    fun_create_coeffs( &f4 , f1.coeffs.real , -1.0 , 1.0 , f1.n );
//    fun_display( &f4 );
    fun_scale ( &f4 , -1.0 , &f4);
    fun_add ( &f1 , &f4 , &f4 );
    printf("||f-f4||_2 = %e\n", fun_norm2( &f4 ));
    fun_clean( &f4 );
  
    /* Clean-up the fun. */
    if ( fun_clean( &f1 ) < 0 )
        printf("fun_test: FUN_EMPTY failed with fun_err=%i (%s).\n",
            fun_err, fun_err_msg[-fun_err]);   
            
    /* All is well. */
    return 0;

    }
