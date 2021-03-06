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
        out[k] = sin( 4.0 * M_PI * x[k] + 1 );

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

    struct fun f1 = FUN_EMPTY, f2 = FUN_EMPTY, f3 = FUN_EMPTY;
    double omega = 5.0;
    int k, res;
    FILE *out, *pipe;
    
    
    /* Initialize the fun f1 (vector real). */
    if ( ( res = fun_create_vec( &f1 , &myfun_vec , -1.0 , 1.0 , NULL ) ) < 0 )
        printf("fun_test: fun_create_vec failed with fun_err=%i (%s).\n",
            fun_err, fun_err_msg[-fun_err]);
            
    /* Initialize the fun f2 (vector reall, param). */
    if ( ( res = fun_create_vec( &f2 , &myfun_vec_param , -1.0 , 1.0 , &omega ) ) < 0 )
        printf("fun_test: fun_create_vec failed with fun_err=%i (%s).\n",
            fun_err, fun_err_msg[-fun_err]);
            
    /* Create f3 = 0.5*f1 + f2. */
    if ( ( res = fun_madd( &f1 , 0.5 , &f2 , 1.0 , &f2 ) ) < 0 )
        printf("fun_test: fun_madd failed with fun_err=%i (%s).\n",
            fun_err, fun_err_msg[-fun_err]);
    /* if ( ( res = fun_mul( &f1 , &f2 , &f3 ) ) < 0 )
        printf("fun_test: fun_mul failed with fun_err=%i (%s).\n",
            fun_err, fun_err_msg[-fun_err]); */
            
            
    /* Open the output file. */
    if ( ( out = fopen( "fun_test.dump" , "w" ) ) == NULL ) {
        printf("fun_test: unable to create fun_test.dump.\n");
        return -1;
        }
    
    /* Show what we got. */
    printf("fun_test: f1.n=%u, f1.coeffs.real[0]=%e\n",f1.n,f1.coeffs.real[0]);
    printf("fun_test: f2.n=%u\n",f2.n);
    printf("fun_test: f3.n=%u\n",f3.n);
    printf("fun_test: int of f3 is %e\n",fun_integrate(&f3));
    for ( k = 0 ; k < 500 ; k++ )
        fprintf(out," %.20e %.20e %.20e %.20e \n", (2.0 * k) / 499 - 1 ,
        fun_eval_clenshaw( &f1 , (2.0 * k)/499 - 1 ) , 
        fun_eval_clenshaw( &f2 , (2.0 * k)/499 - 1 ) , 
        fun_eval_clenshaw( &f3 , (2.0 * k)/499 - 1 ) );
            
    /* Clean-up the fun. */
    if ( fun_clean( &f1 ) < 0 )
        printf("fun_test: FUN_EMPTY failed with fun_err=%i (%s).\n",
            fun_err, fun_err_msg[-fun_err]);
            
    /* Close the output file. */
    fclose(out);
    
    
    /* Fire-up gnuplot. */
    if ( ( pipe = popen( "gnuplot -persist" , "w" ) ) == NULL ) {
        printf("fun_test: unable to create a pipe to gnuplot.\n");
        return -1;
        }
    fprintf( pipe , "set term wxt 0\np 'fun_test.dump' u 1:2 w lp title \"f1\"\n" );
    fprintf( pipe , "set term wxt 1\np 'fun_test.dump' u 1:3 w lp title \"f2\"\n" );
    fprintf( pipe , "set term wxt 2\np 'fun_test.dump' u 1:4 w lp title \"f3\"\n" );
    pclose(pipe);
    
            
    /* All is well. */
    return 0;

    }
