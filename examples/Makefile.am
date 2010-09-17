
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

# Add the source directory to CPPFLAGS
AM_CFLAGS = -g -I../src

# Set-up the library
bin_PROGRAMS = fun_test test_pedro test_nick

# Sources for fun_test
fun_test_SOURCES = fun_test.c
fun_test_LDADD = ../src/.libs/libchebfun.a

# Sources for test_pedro
test_pedro_SOURCES = test_pedro.c
test_pedro_LDADD = ../src/.libs/libchebfun.a

# Sources for test_nick
test_nick_SOURCES = test_nick.c
test_nick_LDADD = ../src/.libs/libchebfun.a


