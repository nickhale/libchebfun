libchebfun
==========

Unmaintained C implementation of some Chebfun basics.

This repo contains some C files which Pedro Gonnet and Nick Hale began working on in 2010. 

I can no longer get the code to compile, but it does still contain a number of C implementations of methods similar to those used in Chebfun V4 and Chebfuin V5. 

It might be useful to anyone hoping to implement a C or C++ version of core Chebfun functionality.


Installation and requirements
=============================

Below are the last known set of instructions which produced a successful compile:


```
        cd libchebfun
        ./autogen.sh
        ./configure
        make
```
If the configure script assumes you have LAPACK, BLAS and FFTW
installed. If it can't find the former two on its own, you can use the
options "--with-lapack=/path/to/your/lapack" and
"--with-blas=/path/to/your/blas".
There are some tests and examples in the "examples folder, e.g.
```
        cd examples
        ./chebtest
```
Note that things may go horribly wrong, since I just added AVX and FMA
instructions that I could not really test yet. If these fail, try
re-configuring with "--with-gcc-arch=pentiumpro"
