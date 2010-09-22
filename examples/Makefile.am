
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

# Add the source directory and debug to CFLAGS
AM_CFLAGS = -g -I../src

# Set-up the library
bin_PROGRAMS = fun_test pedro_test nick_test chebtest

# Sources for fun_test
fun_test_SOURCES = fun_test.c
fun_test_LDADD = ../src/.libs/libchebfun.a

# Sources for pedro_test
pedro_test_SOURCES = pedro_test.c
pedro_test_LDADD = ../src/.libs/libchebfun.a

# Sources for nick_test
nick_test_SOURCES = nick_test.c
nick_test_LDADD = ../src/.libs/libchebfun.a

# Sources for chebtest
chebtest_SOURCES = chebtest.c
chebtest_LDADD = ../src/.libs/libchebfun.a


