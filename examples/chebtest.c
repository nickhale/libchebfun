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


/*
 * These are the chebtests that we would like to implement...
 * adaptapply.m
 * adapt_minsamples.m
 * ad.m
 * ad_vs_diff_1.m
 * ad_vs_diff_2.m
 * airy_extrema.m
 * aliasing.m
 * applycumsumop.m
 * applydiffop.m
 * bcchange.m
 * besseljextrema.m
 * besseljroots.m
 * blowupscl.m
 * blowuptest.m
 * breakpoints.m
 * build.m
 * bvp1.m
 * bvp2.m
 * bvpsimplifytest.m
 * bvptest.m
 * callpref.m
 * ceiltest.m
 * cftest.m
 * chebdomain.m
 * cheboptests.m
 * chebpadetest.m
 * chebpolytest.m
 * chebtest_report.txt
 * circplate.m
 * complexrotation.m
 * composetest.m
 * convspline.m
 * cumsumcos100x.m
 * cumsum_fracexps.m
 * cumsumunbnd.m
 * diffexp.m
 * diffexpsmaps.m
 * diffinf.m
 * diffsingmaps.m
 * divide.m
 * eigsho.m
 * elementary.m
 * ellipj_ode.m
 * ellipjtest.m
 * evalcomplex.m
 * exact_endpoints.m
 * exps_ctor.m
 * exps_cumsum.m
 * exps_diff.m
 * exps_sum.m
 * extracting_roots.m
 * falknerskan.m
 * findexps_test.m
 * fraccalctest.m
 * functionals.m
 * funscale.m
 * gmrestest.m
 * grammatrix.m
 * guide2tests.m
 * infimps.m
 * interp1test.m
 * intops.m
 * inverseq.m
 * invtest.m
 * ivp1.m
 * ivptestcomplex.m
 * ivptest.m
 * ivp_ty_test.m
 * jacptstest.m
 * lebesguetest.m
 * legptstest.m
 * lntAD.m
 * mapends_nan_inf.m
 * mathieu.m
 * matrixnorms.m
 * max_and_imps.m
 * maxdegree.m
 * maxdoubleroots.m
 * maxtest.m
 * mfun_integrate.m
 * misclnttests1.m
 * nanavoidance.m
 * normtests.m
 * operarith.m
 * orrsommerfeld.m
 * orthosincos.m
 * outerprod.m
 * plotrow.m
 * plusquasimatrix.m
 * pointevals.m
 * polyfittest.m
 * polytest.m
 * qrtest.m
 * ratinterptest.m
 * realize_eye.m
 * remeztest.m
 * resampletest.m
 * residuetest.m
 * rest.m
 * restrictimps.m
 * restrict_roots.m
 * restrictscl.m
 * rootspol.m
 * scaleinvariance2.m
 * scaleinvariance.m
 * scaletest.m
 * scribbles.m
 * sinx.m
 * smallintervals.m
 * splittingtest.m
 * sqrt_test.m
 * std_test.m
 * stringinput.m
 * subspacetest.m
 * suminftest.m
 * sumtest.m
 * systemapply.m
 * systemeig.m
 * systemexpm.m
 * systemsolve1.m
 * systemsolve2.m
 * unbndpolys.m
 * vectornorms.m
 * webexamples.m
 */
 
 
/**
 * @brief Prolong an fun making it both larger and smaller, check the result.
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
    if ( fun_norm2( &f3 ) > 2 * DBL_EPSILON )
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
        
    /* Copy f1 into f2. */
    if ( fun_copy( &f1 , &f2 ) < 0 )
        return FAIL;
    
    /* Prolong the fun f2 into f2. */
    if ( fun_prolong( &f2 , f2.n + 20 , &f2 ) < 0 )
        return FAIL;
        
    /* Check if the prolonged function is identical for all practical purposes. */
    if ( fun_madd( &f1 , 1.0 , &f2 , -1.0 , &f3 ) < 0 )
        return FAIL;
    if ( fun_norm2( &f3 ) > 2 * DBL_EPSILON )
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
    if ( fabs(minx_f1 - minx_f2) + fabs(miny_f1 - miny_f2) >= chebopts_opts->eps * 100 )
        return FAIL;
        
    /* If nothing went wrong, just return 0. */
    fun_clean(&f1); fun_clean(&f2);
    return PASS;
        
    }
    
    
/**
 * @brief Runs through a list of tests and reports errors.
 */
 
int main ( int argc , char *argv[] ) {

    /* Adjust these as you add chebtests. */
    const int ntests = 3;
    int (*tests[3])( char ** ) = { &chebtest_sumcos20x , &chebtest_max_min ,
        &chebtest_prolong };
    
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
        
    

