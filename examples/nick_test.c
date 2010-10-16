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
 * @brief A simple (short) vectorized function.
 */
 
int myfun_vec ( const double *x , unsigned int N , double *out , void *data ) {
    
    int k;

    for ( k = 0 ; k < N ; k++ ) {
        out[k] = x[k];
        }

    return 0;

    }


/**
 * @brief A simple (short) non-vectorized function.
 */
 
double myfun ( const double x , void *data ) {
    
    return cos( x ) + sin( x );
//        out[k] = 4.0 + x[k]*(3.0 + x[k]*(2.0 + x[k]));

    }


/**
 * @brief A simple (short) vectorized function.
 */
 
int myfun2_vec ( const double *x , unsigned int N , double *out) {
    
    int k;

    for ( k = 0 ; k < N ; k++ ) {
        out[k] = cos( M_PI * x[k] ) + sin( M_PI * x[k] );
//        out[k] = 4.0 + x[k]*(3.0 + x[k]*(2.0 + x[k]));
        }

    return 0;

    }

/**
 * @brief Do some tests.
 */
 
int main ( int argc , char *argv[] ) {

    struct fun f1 = FUN_EMPTY, f2 = FUN_EMPTY;
    int k, res, nroots;
    double *p;
    
    /* Initialize the fun f1 (vector real).  */
    fun_create_vec( &f1 , &myfun_vec , -1.0 , 1.0 , NULL );

//    fun_create_vec( &f2 , &myfun , 0 , 1.0 ,  NULL );

	/* Test cumsum.
    fun_display( &f1 , stdout );
    fun_indef_integral( &f1, &f2 );
    fun_display( &f2 , stdout ); */

    /* Test isequal 
    printf("error = %e\n", fun_err_norm_inf( &f1, &f2 ));
    printf("equal = %i\n", fun_isequal( &f1, &f2 )); */

    /* Test plotting
    fun_plot( &f1 ); */ 

//    fun_display( &f1 , stdout );
//    fun_restrict( &f1, 0.9 , 1.0 , &f1 );
    fun_display( &f1 , stdout );

    /* test roots
    p = (double *)alloca( sizeof(double) * f1.n );
    nroots = fun_roots( &f1, p );
    printf("nroots = %i\n",nroots);
    for ( k = 0 ; k < nroots ; k++ )  
        printf("%16.16e\n",p[k]); */
        
    /* test comp */
    fun_comp_vec( &f1 , &myfun2_vec , &f2 );
    fun_display( &f2 , stdout );
    fun_plot( &f2 );

    

    /* Test poly
    if ( ( p = (double *)alloca( sizeof(double) * f1.n ) ) == NULL )
        return error(fun_err_malloc);
    fun_poly( &f1, p); 
    for ( k = 0 ; k < f1.n ; k++ )  
        printf("%16.16e\n",p[k]);
    printf("\n"); */
    
	/* Clean */
	fun_clean( &f1 );
    fun_clean( &f2 );
	
    /* All is well. */
    return 0;

    }
