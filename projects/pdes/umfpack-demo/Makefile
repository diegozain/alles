# without the Makefile, compilation is:
# 
# gcc -c -Wall -pedantic -std=c99 -O2 xmalloc.c
# gcc -c -Wall -pedantic -std=c99 -O2 sparse.c
# gcc -c -Wall -pedantic -std=c99 -O2 umfpack-simple.c -lumfpack
# 
# ./umfpack-simple

# declare TARGET variable to handle output
# 
TARGET_SIMPLE = umfpack-simple
TARGET_COOL = umfpack-cool

# declare OFILES variable to handle names
OFILES_SIMPLE = umfpack-simple.o ../../../src/c/xmalloc.o ../../../src/c/sparse.o
OFILES_COOL = umfpack-cool.o ../../../src/c/xmalloc.o

# declare LIBS variable to handle libs
LIBS = -L../../../../SuiteSparse/lib -lumfpack

# choose compiler
# 
CC = gcc

# set CFLAGS
# 
CFLAGS = -Wall -pedantic -std=c99 -O2

# compile SIMPLE 
# $(TARGET_SIMPLE) : $(OFILES_SIMPLE)
# 	$(CC) $(OFILES_SIMPLE) -o $@ $(LIBS)
# compile COOL
$(TARGET_COOL) : $(OFILES_COOL)
	$(CC) $(OFILES_COOL) -o $@ $(LIBS)
# clean
# 
clean:
	rm -f $(OFILES_SIMPLE)
	rm -f $(OFILES_COOL)

# gcc -MM output. Generated after running:
# 
# CFIL = umfpack-simple.c xmalloc.c sparse.c
# 
# $ gcc -MM $(CFIL) >>Makefile
# 
xmalloc.o: ../../../src/c/xmalloc.c ../../../src/c/xmalloc.h
sparse.o: ../../../src/c/sparse.c ../../../src/c/sparse.h ../../../src/c/array.h ../../../src/c/xmalloc.h
umfpack-simple.o: umfpack-simple.c ../../../src/c/sparse.h ../../../src/c/array.h ../../../src/c/xmalloc.h ../../../../SuiteSparse/UMFPACK/Include/umfpack.h
umfpack-cool.o: umfpack-cool.c ../../../src/c/array.h ../../../src/c/xmalloc.h ../../../../SuiteSparse/UMFPACK/Include/umfpack.h