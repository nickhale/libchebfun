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


/* Local includes. */


/* Constants. */


/* Error codes. */
#define util_err_ok                      0
#define util_err_null                    -1
#define util_err_malloc                  -2
#define util_err_fftw                    -3


/* Global, external variables. */
extern int util_err;
extern const char *util_err_msg[];


/* Functions. */
int util_chebptsAB ( unsigned int N , double A , double B , double *x );
int util_chebpts ( unsigned int N , double *x );
double * util_chebpts_alloc ( unsigned int N );
int util_chebpoly ( double *v , unsigned int N , double *c );
double * util_chebpoly_alloc ( double *v , unsigned int N );
int util_chebpolyval ( double *coeffs , unsigned int N , double *v );
double *util_chebpolyval_alloc ( double *coeffs , unsigned int N );
int util_simplify ( double *v , double *coeffs , unsigned int N , double hscale , double vscale , double eps );
double * util_diffmat ( unsigned int N );
double util_clenshaw ( double *coeffs , unsigned int N , double x );
int util_clenshaw_vec ( double *coeffs , unsigned int N , double *x , unsigned int M , double *out );
double util_bary_real ( const double *vals , const double *points , unsigned int N , double x );
double util_bary_vec_real ( const double *vals , const double *points , unsigned int N , const double *x , double *out , unsigned int M );
