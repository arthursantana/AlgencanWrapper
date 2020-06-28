# Algencan has to be compiled with -fPIC
# You only need to add it to FFLAGS in the Makefile
ar -x libalgencan.a; gcc -shared -o ../libalgencan.so *.o -lgfortran
clang -fpic -shared -Wl, -all_load libalgencan.a -o ../libalgencan.dylib
rm *.o
