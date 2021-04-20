# without the Makefile, compilation is:
# 
# gcc -c -Wall -pedantic -std=c99 -O2 xmalloc.c
# gcc -c -Wall -pedantic -std=c99 -O2 sparse.c
# gcc -c -Wall -pedantic -std=c99 -O2 umfpack-simple.c -lumfpack
# 
# ./umfpack-simple

# declare TARGET variable to handle output
# 
TARGET = umfpack-simple

# declare OFILES variable to handle names
# 
OFILES = umfpack-simple.o ../../../src/c/xmalloc.o ../../../src/c/sparse.o

# declare LIBS variable to handle libs
# 
LIBS = -L../../../../SuiteSparse/lib -lumfpack

# choose compiler
# 
CC = gcc

# set CFLAGS
# 
CFLAGS = -Wall -pedantic -std=c99 -O2

# compile
# 
$(TARGET) : $(OFILES)
	$(CC) $(OFILES) -o $@ $(LIBS)

# clean
# 
clean:
	rm -f $(OFILES)

# gcc -MM output. Generated after running:
# 
# CFIL = umfpack-simple.c xmalloc.c sparse.c
# 
# $ gcc -MM $(CFIL) >>Makefile
# 
# NOTE: I dont know what all these local/include files are doing here lol
umfpack-simple.o: umfpack-simple.c ../../../src/c/sparse.h ../../../src/c/array.h ../../../src/c/xmalloc.h ../../../../SuiteSparse/UMFPACK/Include/umfpack.h \
  /usr/local/include/SuiteSparse_config.h \
  /usr/local/include/umfpack_symbolic.h \
  /usr/local/include/umfpack_numeric.h \
  /usr/local/include/umfpack_solve.h \
  /usr/local/include/umfpack_free_symbolic.h \
  /usr/local/include/umfpack_free_numeric.h \
  /usr/local/include/umfpack_defaults.h \
  /usr/local/include/umfpack_qsymbolic.h \
  /usr/local/include/umfpack_wsolve.h \
  /usr/local/include/umfpack_triplet_to_col.h \
  /usr/local/include/umfpack_col_to_triplet.h \
  /usr/local/include/umfpack_transpose.h \
  /usr/local/include/umfpack_scale.h \
  /usr/local/include/umfpack_get_lunz.h \
  /usr/local/include/umfpack_get_numeric.h \
  /usr/local/include/umfpack_get_symbolic.h \
  /usr/local/include/umfpack_save_numeric.h \
  /usr/local/include/umfpack_load_numeric.h \
  /usr/local/include/umfpack_save_symbolic.h \
  /usr/local/include/umfpack_load_symbolic.h \
  /usr/local/include/umfpack_get_determinant.h \
  /usr/local/include/umfpack_report_status.h \
  /usr/local/include/umfpack_report_info.h \
  /usr/local/include/umfpack_report_control.h \
  /usr/local/include/umfpack_report_matrix.h \
  /usr/local/include/umfpack_report_triplet.h \
  /usr/local/include/umfpack_report_vector.h \
  /usr/local/include/umfpack_report_symbolic.h \
  /usr/local/include/umfpack_report_numeric.h \
  /usr/local/include/umfpack_report_perm.h \
  /usr/local/include/umfpack_timer.h /usr/local/include/umfpack_tictoc.h \
  /usr/local/include/amd.h /usr/local/include/umfpack_global.h
xmalloc.o: ../../../src/c/xmalloc.c ../../../src/c/xmalloc.h
sparse.o: ../../../src/c/sparse.c ../../../src/c/sparse.h ../../../src/c/array.h ../../../src/c/xmalloc.h
