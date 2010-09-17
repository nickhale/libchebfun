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
#include <float.h>

/* Local includes. */
#include "chebopts.h"


/**
 * @brief Default #chebopts settings.
 * 
 * This #chebopts struct is used as a default when no options are provided
 * by the user.
 */
const struct chebopts chebopts_default = {
    0 ,             /* flags */
    9 ,             /* minsamples */
    65536 ,         /* maxdegree */
    6000 ,          /* maxlength */
    128 ,           /* splitdegree */
    -1.0 ,          /* a */
    1.0 ,           /* b */
    DBL_EPSILON     /* eps */
    };


