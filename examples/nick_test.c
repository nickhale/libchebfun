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
 
int myfun ( const double *x , unsigned int N , double *out , void *data ) {
    
    int k;

    for ( k = 0 ; k < N ; k++ ) {
//        out[k] = cos( M_PI * x[k] ) + sin( M_PI * x[k] );
        out[k] = 1.0+0*x[k];
        }

    return 0;

    }

/**
 * @brief Do some tests.
 */
 
int main ( int argc , char *argv[] ) {

    struct fun f1 = FUN_EMPTY, f2 = FUN_EMPTY;
    int k, res, nroots;
    
    /* Initialize the fun f1 (vector real). */
    fun_create_vec( &f1 , &myfun , -1.0 , 1.0 , NULL );

	/* Test cumsum. */
    fun_display( &f1 , stdout );
    fun_indef_integral( &f1, &f2 );
    fun_display( &f2 , stdout );

	/* Clean f1 */
	fun_clean( &f1 );
    fun_clean( &f2 );
	
    /* All is well. */
    return 0;

    }
