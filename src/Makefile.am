
# This file is part of libchebfun.
# Coypright (c) 2010 The Chebfun Team.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Automake stuff
AUTOMAKE_OPTIONS=gnu

# Add the debug flag to the whole thing
AM_CFLAGS = -g -O3 $(SIMD_FLAGS) $(CFLAGS)
# AM_CFLAGS = -g -O2 -Wall -msse2

# Build the libchebfun library
lib_LTLIBRARIES = libchebfun.la
libchebfun_la_SOURCES = chebopts.c fun.c util.c errs.c

# List required headers
include_HEADERS = chebopts.h fun.h util.h

# Assign a "safe" version number
libchebfun_la_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) -version-info 0:0:0

