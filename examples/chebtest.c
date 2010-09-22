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
#include "chebopts.h"
#include "util.h"
#include "fun.h"


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
 * max_min.m
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
 * @brief Re-implementation of the chebtest sumcos20x, impemented for #fun instead
 * of #chebfun.
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
    if ( fun_create_vec( &f , &thefun , -1.0 , 1.0 , NULL , NULL ) < 0 )
        return -__LINE__;
        
    /* Do the first test. */
    if ( fabs( fun_integrate(&f) - sin(20.0)/10.0 ) >= 1.5e-15 * chebopts_default.eps / DBL_EPSILON )
        return -__LINE__;
        
    /* If nothing went wrong, just return 0. */
    return 0;
        
    }
    
    
/**
 * @brief Runs through a list of tests and reports errors.
 */
 
int main ( int argc , char *argv[] ) {

    /* Adjust these as you add chebtests. */
    const int ntests = 1;
    const int (*tests[1])( char ** ) = { &chebtest_sumcos20x };
    
    int k, res;
    char *name = NULL;
    
    /* Loop over the chebtests. */
    for ( k = 0 ; k < ntests ; k++ ) {
    
        /* Call the kth chebtest. */
        res = (*tests[k])( &name );
        
        /* Be verbose about the result. */
        if ( res < 0 )
            printf("chebtest: test %s failed on line %i of file %s.\n", name, -res, __FILE__ );
        else
            printf("chebtest: test %s passed.\n", name);
    
        }
    
    /* Leave quietly. */
    return 0;

    }
        
    

