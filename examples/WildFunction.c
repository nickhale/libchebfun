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

/* WildFunction */
int main ( int argc , char *argv[] ) {

    struct fun f = FUN_EMPTY, s = FUN_EMPTY;
    double sum, maxx, maxy;
    int j, k, runs = 10;

    // sin(pi*x).
    int sinpix ( const double *x , unsigned int N , double *v , void *data ) {
        int k;
        for ( k = 0 ; k < N ; k++ )
            v[k] = sin( M_PI * x[k] );
        return 0;
        }
    // x^4.
    int pow4 ( const double *x , unsigned int N , double *v , void *data ) {
        int k;
        double tmp;
        for ( k = 0 ; k < N ; k++ ) {
            tmp = x[k]*x[k];
            v[k] = tmp*tmp;
            }
        return 0;
        }
    int thefun ( const double *x , unsigned int N , double *v , void *data ) {
        int k;
        double tmp;
        for ( k = 0 ; k < N ; k++ ) {
            tmp = x[k]*x[k];
            v[k] = (0.75)*(1.0-2.0*tmp*tmp);
            }
        return 0;
        }
        
    /* Set the default to use bary. */
    /* chebopts_opts->flags |= chebopts_flag_evalbary; */
        
    /* Did the user specify a number of runs? */
    if ( argc > 1 )
        runs = atoi(argv[1]);
        
    /* Main loop. */
    for ( k = 0 ; k < runs ; k++ ) {

        // Make starting funs.
        fun_create_vec( &f , &sinpix , -1.0 , 1.0 , NULL );
        fun_copy( &f , &s);

        // Iterative definition.
        for ( j = 1 ; j <= 15 ; j++ ) { 
    //        f = (3/4)*(1 - 2*f.^4);
    /*
            fun_comp_vec( &f , &pow4 , &f , NULL );   // f -> f^4
            fun_scale( &f , -2.0 , &f );              // f -> -2*f
            fun_add_const( &f , 1.0 , &f );           // f -> f+1
            fun_scale( &f , 3.0/4.0 , &f );           // f -> 3/4*f
    */
            fun_comp_vec( &f , &thefun , &f , NULL ); // f -> (3/4)*(1-2*f^4)
    //        s = s + f;
            fun_add( &s , &f , &s );
            }

        // Plot.
        /* if ( fun_gnuplot( &s ) < 0 )
            errs_dump( stdout ); */

        // Integrate s
        sum = fun_integrate( &s );
        printf("length(s) = %i\n", s.n);
        printf("sum(s) = %16.16e\n", sum);
        printf("err from MATLAB/Chebfun = %16.16e\n", fabs(sum-15.265483825826742));

        // Maximum of s
        fun_max( &s , &maxy , &maxx );
        printf("max(s) = %16.16e\n", maxy);
        printf("err from MATLAB/Chebfun = %16.16e\n", fabs(maxy-9.295487828334499));

	    // Clean.
	    fun_clean( &f );
	    fun_clean( &s );
        
        }
	
    // All is well.
    return 0;

    }

