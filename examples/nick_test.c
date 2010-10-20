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
#include "errs.h"


/**
 * @brief Do some tests.
 */
 
int main ( ) {

    struct fun x = FUN_EMPTY, f = FUN_EMPTY;
    double *roots;
    int nroots; 

    int thefun ( const double *x , unsigned int N , double *v , void *data ) {
        int k;
        for ( k = 0 ; k < N ; k++ )
            v[k] = x[k]*x[k]-0.5;
        return 0;
        }

    // Make x.
    fun_x( &x, -1.0 , 1.0 );
//    fun_display( &x , stdout );
//    fun_plot( &x );

    // Compose.
    fun_comp_vec( &x , &thefun , &f );
//    fun_plot( &f );
    fun_display( &f , stdout );

    // Compute roots
    roots = (double *)alloca( sizeof(double) * f.n );
    nroots = fun_roots( &f , roots );
  
	// Clean.
	fun_clean( &x );
	fun_clean( &f );
	
    // All is well.
    return 0;

    }

/**
 * @brief Do some tests.
 */
 
/* ExtremeExtrema
int main ( ) {

    struct fun x = FUN_EMPTY, f = FUN_EMPTY;

    // The function cos(x)*sin(exp(x)).
    int thefun ( const double *x , unsigned int N , double *v , void *data ) {
        int k;
        for ( k = 0 ; k < N ; k++ )
            v[k] = cos(x[k])*sin(exp(x[k]));
        return 0;
        }

    // Make x.
    fun_x( &x, 0.0 , 6.0 );
//    fun_display( &x , stdout );
//    fun_plot( &x );

    // Compose.
    fun_comp_vec( &x , &thefun , &f );
//    fun_plot( &f );

	// Clean.
	fun_clean( &x );
	fun_clean( &f );
	
    // All is well.
    return 0;

    }
*/

/* BesselRoots 
int main ( ) {

    struct fun x = FUN_EMPTY, f = FUN_EMPTY;
    double *roots;
    int nroots;

    // The Bessel function.
    int thefun ( const double *x , unsigned int N , double *v ) {
        int k;
        for ( k = 0 ; k < N ; k++ )
            v[k] = j0(x[k]);
        return 0;
        }

    // Make x.
    fun_x( &x, 0.0 , 100.0 );
//    fun_display( &x , stdout );
    if ( fun_plot( &x ) < 0 )
        errs_dump( stdout );
    

    // Compose.
    fun_comp_vec( &x , &thefun , &f );
    if ( fun_plot( &f ) < 0 )
        errs_dump( stdout );

    // Find the roots.
    nroots = fun_roots( &f , roots );
    printf("The are %i roots of the Bessel function in [%9.9e %9.9e]\n",nroots,f.a,f.b);
//    for ( k = 0; k < nroots ; k++ )
//        printf( "%16.16e\n",roots[k]);

	// Clean.
	fun_clean( &x );
	fun_clean( &f );

    // Make x again.
    fun_x( &x, 1000000.0 , 1001000.0 );
//    fun_display( &x , stdout );
//    fun_plot( &x );

    // Compose.
    if ( fun_comp_vec( &x , &thefun , &f ) < 0 )
        errs_dump( stdout );
//    fun_display( &f , stdout );
//    if ( fun_plot( &f ) < 0 )
//        errs_dump( stdout );
//    printf("The length of f is %i\n",f.n);
    
    // Find the roots.
    nroots = fun_roots( &f , roots );
    printf("The are %i roots of the Bessel function in [%9.9e %9.9e]\n",nroots,f.a,f.b);

	// Clean.
	fun_clean( &x );
	fun_clean( &f );
	
    // All is well.
    return 0;

    }
*/

