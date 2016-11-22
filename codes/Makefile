
#CC=icc
CC=g++
CFLAGS=-I. -Wno-write-strings
OBJECTS = allocate.o data_term.o em.o mpm.o ndmulti_e.o qc_badEO.o qc_badIO.o qc_doubleSeg.o qc_freeEndsConnection.o qc_freeSeg.o qc_singleSeg.o qc_transitionP.o qualityCandy.o randlib.o tiff.o util.o

default: ChannelMPP 


ChannelMPP: main.o $(OBJECTS)
	$(CC) -o ChannelMPP main.o $(OBJECTS)


main.o: main.cpp allocate.h random.h tiff.h neck.h ndmulti_e.h QualityCandy.h
	$(CC) $(CFLAGS) -c main.cpp 

allocate.o: allocate.cpp allocate.h
	$(CC) $(CFLAGS) -c allocate.cpp 

data_term.o: data_term.cpp tiff.h allocate.h util.h
	$(CC) $(CFLAGS) -c data_term.cpp 

em.o: em.cpp em.h mpm.h
	$(CC) $(CFLAGS) -c em.cpp 

mpm.o: mpm.cpp em.h mpm.h
	$(CC) $(CFLAGS) -c mpm.cpp 

ndmulti_e.o: neck.h ndmulti_e.h random.h QualityCandy.h ndmulti_e.cpp
	$(CC) $(CFLAGS) -c ndmulti_e.cpp

qc_badEO.o: QualityCandy.h qc_badEO.cpp
	$(CC) $(CFLAGS) -c qc_badEO.cpp

qc_badIO.o: QualityCandy.h qc_badIO.cpp
	$(CC) $(CFLAGS) -c qc_badIO.cpp

qc_doubleSeg.o: QualityCandy.h qc_doubleSeg.cpp
	$(CC) $(CFLAGS) -c qc_doubleSeg.cpp

qc_freeEndsConnection.o: QualityCandy.h tiff.h qc_freeEndsConnection.cpp
	$(CC) $(CFLAGS) -c qc_freeEndsConnection.cpp

qc_freeSeg.o: QualityCandy.h qc_freeSeg.cpp
	$(CC) $(CFLAGS) -c qc_freeSeg.cpp

qc_singleSeg.o: QualityCandy.h qc_singleSeg.cpp
	$(CC) $(CFLAGS) -c qc_singleSeg.cpp

qc_transitionP.o: QualityCandy.h tiff.h
	$(CC) $(CFLAGS) -c qc_transitionP.cpp

qualityCandy.o: tiff.h allocate.h neck.h qualityCandy.cpp
	$(CC) $(CFLAGS) -c qualityCandy.cpp

randlib.o: randlib.h randlib.cpp
	$(CC) $(CFLAGS) -c randlib.cpp

tiff.o: tiff.cpp tiff.h
	$(CC) $(CFLAGS) -c tiff.cpp

util.o: util.h util.cpp
	$(CC) $(CFLAGS) -c util.cpp

clean: 
	$(RM) count *.o *~