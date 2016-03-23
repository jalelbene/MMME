CC=g++

ITINCLUDE3=-I/home/s/it++3.7.0/linux/usr.local.include/it++3.7.0/

ITINCLUDE = ${ITINCLUDE3}

ITLIB=-L/usr/lib/ -L/usr/local/lib/


INCLUDE_PATH=-I./
LDFLAGS=${ITLIB}
LIBRARY=-lstdc++  -lm -lblas -llapack  -lg2c -lgslcblas
WARN = -Wno-deprecated 
CPPFLAGS= -ffast-math -funroll-loops ${WARN} ${INCLUDE_PATH}
DEFINES = -DHAVE_LAPACK -DHAVE_CBLAS -Wall

.cpp.o: ; ${CC} ${DEFINES} ${INCLUDE_PATH} ${WARN} -g -c $*.cpp

OBJ= matrice.o pendule.o
all:	${OBJ}
	${CC} ${CPPFLAGS}  ${DEFINES}  ${OBJ} -o pendule++ ${LDFLAGS}
clean:
	/bin/rm -rf *.o *.exe *~ pendule++
