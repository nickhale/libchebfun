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


/* Flags. */
#define chebopts_flag_splitting         1
#define chebopts_flag_resampling        2
#define chebopts_flag_sampletest        4
#define chebopts_flag_blowup            8
#define chebopts_flag_extrapolate       16
#define chebopts_flag_polishroots       32

#define chebopts_flag_compcoeffs        64
#define chebopts_flag_compvals          128
#define chebopts_flag_comppts           256
#define chebopts_flag_evalbary          512


/* Error codes. */
#define chebopts_err_ok                 0


/* Global, external variables. */
extern struct chebopts *chebopts_opts;
extern struct chebopts chebopts_default;


/**
 * The @a chebopts struct contains information on the options
 * and parameters used when constructing and manipulating
 * @a chebfuns.
 */
struct chebopts {

    /** Flags. */
    unsigned int flags;
    
    /** Minimum nr. of samples. */
    unsigned int minsamples;
    
    /** Maximum polynomial degree. */
    unsigned int maxdegree;
    
    /** Maximum nr. of points. */
    unsigned int maxlength;
    
    /** Degree as of which a chebfun should be split. */
    unsigned int splitdegree;
    
    /** Default domain. */
    double a, b;
    
    /** Relative tolerance. */
    double eps;

    };

