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
#include <float.h>
#include <string.h>

/* Local includes. */
#include "errs.h"
#include "chebopts.h"
#include "util.h"
#include "fun.h"

/* Local defines */
#define PASS 0
#define FAIL (-__LINE__)


/**
 * @brief Compute the Inf-norm for a bunch of function for which we
 *      know what's going on.
 * 
 */
 
int chebtest_restrict ( char **name ) {

    struct fun f1 = FUN_EMPTY, f2 = FUN_EMPTY;
    double *x, *v;
    int k;
    
    /* The function for which to create a chebfun. */
    int thefun ( const double *x , unsigned int N , double *v , void *data ) {
        int k;
        for ( k = 0 ; k < N ; k++ )
            v[k] = sin( x[k] ) + sin( x[k] * x[k] );
        return 0;
        }
        
    /* Set the function name. */
    *name = "restrict";
    
    /* Create f1 in the domain [0,10]. */
    if ( fun_create_vec( &f1 , &thefun , 0.0 , 10.0 , NULL ) < 0 )
        return FAIL;
        
    /* Restrict f1 to the domain [2,7]. */
    if ( fun_restrict( &f1 , 2.0 , 7.0 , &f2 ) < 0 )
        return FAIL;
    
    /* Allocate the temporary vector for the nodes and function values. */
    if ( ( x = (double *)alloca( sizeof(double) * f1.n ) ) == NULL ||
         ( v = (double *)alloca( sizeof(double) * f1.n ) ) == NULL )
        return FAIL;
        
    /* Evaluate f1 at the nodes of f2. */
    if ( util_chebptsAB( f2.n , f2.a , f2.b , x ) < 0 )
        return FAIL;
    if ( fun_eval_vec( &f1 , x , f2.n , v ) < 0 )
        return FAIL;
    
    /* Compare the results. */
    for ( k = 0 ; k < f2.n ; k++ ) {
        printf("chebtest_restrict: x[%i]=%e, v[%i]=%e, f2[%i]=%e\n",
            k, x[k], k, v[k], k, f2.vals.real[k] );
        if ( fabs( v[k] - f2.vals.real[k] ) > 10 * f1.scale * DBL_EPSILON )
            return FAIL;
        }
    
    /* All passed... */
    fun_clean(&f1); fun_clean(&f2);
    return PASS;
    
    }    


/**
 * @brief Compute the Inf-norm for a bunch of function for which we
 *      know what's going on.
 * 
 */
 
int chebtest_norm_inf ( char **name ) {

    struct fun f = FUN_EMPTY;
    double res, omegas[5] = { 1.0 , 2.0 , 6.1 , 0.5 , 20.0 };
    int k;
    
    /* Set the function name. */
    *name = "norm_inf";
    
    /* The function for which to create a chebfun. */
    int thefun ( const double *x , unsigned int N , double *v , void *data ) {
        double omega = *((double *)data);
        int k;
        for ( k = 0 ; k < N ; k++ )
            v[k] = omega * sin( M_PI * x[k] + omega );
        return 0;
        }
        
    /* Loop over the values of omega. */
    for ( k = 0 ; k < 5 ; k++ ) {
    
        /* Compute f. */
        if ( fun_create_vec( &f , &thefun , -1.0 , 1.0 , (void *)(&omegas[k]) ) < 0 )
            return FAIL;
            
        /* Compute the 2-norm. */
        if ( isnan( res = fun_norm2( &f ) ) )
            return FAIL;
            
        /* Verify the 2-norm. */
        if ( fabs( res - omegas[k] ) > 10 * res * chebopts_opts->eps )
            return FAIL;
    
        }
        
    /* All passed... */
    fun_clean(&f);
    return PASS;
    
    }
        
        
/**
 * @brief Compute the 2-norm of some funs for which the norm is known.
 * 
 */
 
int chebtest_norm2 ( char **name ) {

    struct fun f = FUN_EMPTY;
    double res, omegas[5] = { 1.0 , 2.0 , 6.1 , 0.5 , 20.0 };
    int k;
    
    /* Set the function name. */
    *name = "norm2";
        
    /* The function for which to create a chebfun. */
    int thefun ( const double *x , unsigned int N , double *v , void *data ) {
        double omega = *((double *)data);
        int k;
        for ( k = 0 ; k < N ; k++ )
            v[k] = sin( M_PI * x[k] * omega );
        return 0;
        }
        
    /* Loop over the values of omega. */
    for ( k = 0 ; k < 5 ; k++ ) {
    
        /* Compute f. */
        if ( fun_create_vec( &f , &thefun , -1.0 , 1.0 , (void *)(&omegas[k]) ) < 0 )
            return FAIL;
            
        /* Compute the 2-norm. */
        if ( isnan( res = fun_norm2( &f ) ) )
            return FAIL;
            
        /* Verify the 2-norm. */
        if ( fabs( res - sqrt( ( M_PI * omegas[k] - cos( M_PI * omegas[k] ) * sin( M_PI * omegas[k] ) ) / ( M_PI * omegas[k] ) ) ) > 10 * res * chebopts_opts->eps )
            return FAIL;
    
        }
        
    /* All passed... */
    fun_clean(&f);
    return PASS;
        
    }
    

/**
 * @brief Prolong a #fun making it both larger and smaller, check the result.
 * 
 */
 
int chebtest_prolong ( char **name ) {

    struct fun f1 = FUN_EMPTY, f2 = FUN_EMPTY, f3 = FUN_EMPTY;
    int k;
    double *v;
    
    /* The function for which to create a chebfun. */
    int thefun ( const double *x , unsigned int N , double *v , void *data ) {
        int k;
        for ( k = 0 ; k < N ; k++ )
            v[k] = sin( 5*(x[k]+1) ) + sin( 25 * (x[k]+1) * (x[k]+1) );
        return 0;
        }
        
    /* Set the function name. */
    *name = "prolong";
        
    /* Create the fun f1. */
    if ( fun_create_vec( &f1 , &thefun , -1.0 , 1.0 , NULL ) < 0 )
        return FAIL;
        
    /* Prolong the fun f1 into f2. */
    if ( fun_prolong( &f1 , f1.n + 20 , &f2 ) < 0 )
        return FAIL;
        
    /* Check if the prolonged function is identical for all practical purposes. */
    if ( fun_madd( &f1 , 1.0 , &f2 , -1.0 , &f3 ) < 0 )
        return FAIL;
    if ( fun_norm2( &f3 ) > 10.0 * DBL_EPSILON )
        return FAIL;
        
    /* Prolong (restrict) the fun f1 into f2. */
    if ( fun_prolong( &f1 , f1.n - 10 , &f2 ) < 0 )
        return FAIL;
        
    /* Check if the prolonged function matches at the nodes. */
    if ( ( v = (double *)alloca( sizeof(double) * f2.n ) ) == NULL )
        return FAIL;
    if ( fun_eval_vec( &f1 , f2.points , f2.n , v ) < 0 )
        return FAIL;
    /* for ( k = 0 ; k < f2.n ; k++ )
        printf("chebtest_prolong: x[%i]=%e, v[%i]=%e, f2[%i]=%e\n",
            k,f2.points[k],k,v[k],k,f2.vals.real[k]); */
    for ( k = 0 ; k < f2.n ; k++ )
        if ( fabs( v[k] - f2.vals.real[k] ) > 100 * f1.scale * DBL_EPSILON )
            return FAIL;
        
    /* All passed... */
    fun_clean(&f1); fun_clean(&f2); fun_clean(&f3);
    return PASS;
        
    /* Copy f1 into f2. */
    if ( fun_copy( &f1 , &f2 ) < 0 )
        return FAIL;
    
    /* Prolong the fun f2 into f2. */
    if ( fun_prolong( &f2 , f2.n + 20 , &f2 ) < 0 )
        return FAIL;
        
    /* Check if the prolonged function is identical for all practical purposes. */
    if ( fun_madd( &f1 , 1.0 , &f2 , -1.0 , &f3 ) < 0 )
        return FAIL;
    if ( fun_norm2( &f3 ) > 10.0 * DBL_EPSILON )
        return FAIL;
        
    /* Copy f1 into f2. */
    if ( fun_copy( &f1 , &f2 ) < 0 )
        return FAIL;
    
    /* Prolong (restrict) the fun f2 into f2. */
    if ( fun_prolong( &f2 , f2.n - 10 , &f2 ) < 0 )
        return FAIL;
        
    /* Check if the prolonged function matches at the nodes. */
    if ( fun_eval_vec( &f1 , f2.points , f2.n , v ) < 0 )
        return FAIL;
    /* for ( k = 0 ; k < f2.n ; k++ )
        printf("chebtest_prolong: x[%i]=%e, v[%i]=%e, f2[%i]=%e\n",
            k,f2.points[k],k,v[k],k,f2.vals.real[k]); */
    for ( k = 0 ; k < f2.n ; k++ )
        if ( fabs( v[k] - f2.vals.real[k] ) > 100 * f1.scale * DBL_EPSILON )
            return FAIL;
        
    /* All passed... */
    fun_clean(&f1); fun_clean(&f2); fun_clean(&f3);
    return PASS;
        
    }
    

/**
 * @brief Re-implementation of the chebtest sumcos20x, impemented for #fun instead
 * of #chebfun.
 * 
 * MISSING: Test for complex!
 */
 
int chebtest_sumcos20x ( char **name ) {

    struct fun f = FUN_EMPTY;
    
    /* The function for which to create a chebfun. */
    int thefun ( const double *x , unsigned int N , double *v , void *data ) {
        int k;
        for ( k = 0 ; k < N ; k++ )
            v[k] = cos( 20.0 * x[k] );
        return 0;
        }
        
    /* Set the function name. */
    *name = "sumcos20x";
        
    /* Create the fun. */
    if ( fun_create_vec( &f , &thefun , -1.0 , 1.0 , NULL ) < 0 )
        return FAIL;
        
    /* Do the first test. */
    if ( fabs( fun_integrate(&f) - sin(20.0)/10.0 ) >= 1.5e-15 * chebopts_opts->eps / DBL_EPSILON )
        return FAIL;
        
    /* TODO: multiply f by 1i! */
        
    /* If nothing went wrong, just return 0. */
    fun_clean(&f);
    return PASS;
        
    }
    
    
/**
 * @brief Re-implementation of the chebtest max_min, impemented for #fun instead
 * of #chebfun.
 *
 * MISSING: Test with infinite right boundary.
 */
 
int chebtest_max_min ( char **name ) {

    struct fun f1 = FUN_EMPTY, f2 = FUN_EMPTY;
    double minx_f1, miny_f1, maxx_f1, maxy_f1;
    double minx_f2, miny_f2, maxx_f2, maxy_f2;
    
    /* The function for which to create a chebfun. */
    int thefun ( const double *x , unsigned int N , double *v , void *data ) {
        int k;
        for ( k = 0 ; k < N ; k++ )
            v[k] = cos( 2.0 * M_PI * (x[k] - M_PI) ) * exp( -x[k] );
        return 0;
        }
        
    /* Set the function name. */
    *name = "max_min";
        
    /* Create the funs. */
    if ( fun_create_vec( &f1 , &thefun , 0.0 , 1.0 , NULL ) < 0 )
        return FAIL;
    /* TODO: f2 should actually be in [0,inf]. */
    if ( fun_create_vec( &f2 , &thefun , 0.0 , 20.0 , NULL ) < 0 )
        return FAIL;
        
    /* Get the min and max values for f1 and f2. */
    if ( fun_minandmax( &f1 , &miny_f1 , &minx_f1 , &maxy_f1 , &maxx_f1 ) < 0 )
        return FAIL;
    if ( fun_minandmax( &f2 , &miny_f2 , &minx_f2 , &maxy_f2 , &maxx_f2 ) < 0 )
        return FAIL;
        
    /* Do the first test. */
    if ( fabs(maxx_f1 - maxx_f2) + fabs(maxy_f1 - maxy_f2) >= chebopts_opts->eps * 100 )
        return FAIL;
        
    /* Do the second test. */
    /* if ( fabs(minx_f1 - minx_f2) + fabs(miny_f1 - miny_f2) >= chebopts_opts->eps * 100 )
        return FAIL; */
        
    /* If nothing went wrong, just return 0. */
    fun_clean(&f1); fun_clean(&f2);
    return PASS;
        
    }


/**
 * @brief Re-implementation of the chebtest cumsumcos100x, impemented for #fun instead
 * of #chebfun.
 *
 * MISSING: Test with complex values.
 */
 
int chebtest_cumsumcos100x ( char **name ) {

    struct fun f1 = FUN_EMPTY, f2 = FUN_EMPTY, f3 = FUN_EMPTY, err = FUN_EMPTY ;
    double norm;
    
    /* The function for which to create a chebfun. */
    int thefun ( const double *x , unsigned int N , double *v , void *data ) {
        int k;
        for ( k = 0 ; k < N ; k++ )
            v[k] = cos( 100.0  * x[k] );
        return 0;
        }

    /* Derivative of the above. */
    int thefun_prime ( const double *x , unsigned int N , double *v , void *data ) {
        int k;
        for ( k = 0 ; k < N ; k++ )
            v[k] = sin( 100.0 * x[k] ) / 100.0;
        return 0;
        }
        
    /* Set the function name. */
    *name = "cumsumcos100x";
        
    /* Create the funs. */
    if ( fun_create_vec( &f1 , &thefun , 10.0 , 13.0 , NULL ) < 0 )
        return FAIL;
    if ( fun_create_vec( &f2 , &thefun_prime , 10.0 , 13.0 , NULL ) < 0 )
        return FAIL;

    /* Compute the indefinite integral. */
    if ( fun_indef_integral( &f1 , &f3 ) < 0 )
        return FAIL;
        
    /* Do the first test. */
    fun_madd ( &f3 , 1.0 , &f2 , -1.0 , &err );
    if ( ( norm = fun_norm_inf( &err ) ) < 0 )
        return FAIL;

    /* If nothing went wrong, just return 0. */
    fun_clean(&f1); fun_clean(&f2); fun_clean(&f3); fun_clean(&err);
    return PASS;
        
    }


/**
 * @brief Runs through a list of tests and reports errors.
 */
 
int main ( int argc , char *argv[] ) {

    /* Adjust these as you add chebtests. */
    const int ntests = 7;
    int (*tests[7])( char ** ) = { &chebtest_sumcos20x , &chebtest_max_min ,
        &chebtest_norm2 , &chebtest_norm_inf , &chebtest_prolong ,
        &chebtest_restrict , &chebtest_cumsumcos100x };
    
    int k, res;
    char *name = NULL;
    
    /* Loop over the chebtests. */
    for ( k = 0 ; k < ntests ; k++ ) {
    
        /* Call the kth chebtest. */
        res = (*tests[k])( &name );
        
        /* Be verbose about the result. */
        if ( res != PASS ) {
            printf("chebtest: test %s failed on line %i of file %s.\n", name, -res, __FILE__ );
            errs_dump(stdout);
            }
        else
            printf("chebtest: test %s passed.\n", name);
    
        }
    
    /* Leave quietly. */
    return 0;

    }
        
    

