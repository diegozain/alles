# without the Makefile, compilation is:
# 
# gcc -c -Wall -pedantic -std=c99 -O2 array.c
# gcc -c -Wall -pedantic -std=c99 -O2 xmalloc.c
# gcc -c -Wall -pedantic -std=c99 -O2 haar.c
# gcc -c -Wall -pedantic -std=c99 -O2 haar-test.c
# 
# ./haar-test

# declare TARGET variable to handle output
# 
TARGET = haar-test

# declare OFILES variable to handle names
# 
OFILES = haar-test.o ../../../src/c/haar.o ../../../src/c/xmalloc.o

# declare LIBS variable to handle libs
# 
# LIBS = -L../path/to/libs

# choose compiler
# 
CC = gcc

# set CFLAGS
# 
CFLAGS = -Wall -pedantic -std=c99 -O2

# compile
# 
$(TARGET) : $(OFILES)
	$(CC) $(OFILES) -o $@

# clean
# 
clean:
	rm -f $(OFILES)

# gcc -MM output. Generated after running:
# 
# CFIL = haar-test.c array.c haar.c
# 
# $ gcc -MM $(CFIL) >>Makefile
# 
haar-test.o: haar-test.c ../../../src/c/haar.h ../../../src/c/array.h ../../../src/c/xmalloc.h
xmalloc.o: ../../../src/c/xmalloc.c ../../../src/c/xmalloc.h
array.o: ../../../src/c/array.h
haar.o: ../../../src/c/haar.c ../../../src/c/haar.h
