# Use this Makefile with make

# Executable name
CMD = zstruct.exe

# -------- description of DFLAGS ---------------


# -------- Define environmental variable C_COMPILER -----------
# Make sure it is defined
#          ifeq ($(strip$(FORTRAN_COMPILER)),)
# Otherwise you can define it here also by uncommenting next line
 FC = icpc -openmp -I$(MKLROOT)/include
# FC = g++ -fopenmp -I$(MKLROOT)/include
DFLAGS =  #-Define the cpp flags to be used
#DFLAGS =  #-Define the cpp flags to be used
OFLAGS =  # optimization

#Intel parallel openmp (only w/icpc compiler)
#LINKERFLAGS =  -L$(MKLROOT)/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
LINKERFLAGS =  -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread 
# MAC OS linkers
#LINKERFLAGS = -lm -framework Accelerate



#
# Implicit rules for handling src files
#  ( uses static pattern rules, see info for make )
.c.o:
	$(FC) -c -g $(DFLAGS) -Wimplicit $<
.cpp.o:
	$(FC) -c -g $(DFLAGS) $<

OBJECTS = zstruct.o main.o icoord.o pTable.o stringtools.o utils.o mem.o print.o mopac.o dft.o opt.o mm_grad.o lsq.o rxndb.o rxnftr.o atom.o align.o iccomp.o write.o knnr.o rtype.o tm.o nbo.o

$(CMD) : $(OBJECTS)
	$(FC) $(DEBUG_FLAGS) $(OFLAGS) $(LINKERFLAGS) $(OBJECTS)  -o ./$(CMD)

clean:
	/bin/rm -f *.o *.i *.mod *.exe a.out make.log

cleano:
	rm -f *.o *.i

depend :
	g++ -MM *.cpp *.c >> Makefile 

# DO NOT DELETE created with g++ -MM *.cpp *.c
main.o: main.cpp zstruct.h icoord.h lsq.h rxndb.h rxnftr.h align.h
lsq.o: lsq.cpp lsq.h utils.h
rxndb.o: rxndb.cpp rxndb.h utils.h rxnftr.h knnr.h rtype.h
rxnftr.o: rxnftr.cpp rxnftr.h utils.h
align.o: align.cpp align.h icoord.h utils.h
paths.o: paths.cpp zstruct.h icoord.h
pgsm.o: pgsm.cpp pgsm.h constants.h
dft.o: dft.cpp dft.h constants.h
iso.o: iso.cpp zstruct.h icoord.h geombasis.h constants.h utils.h mopac.h
atom.o: atom.cpp atom.h icoord.h
zstruct.o: zstruct.cpp zstruct.h utils.h constants.h print.h lsq.h rxndb.h rxnftr.h align.h dft.h rtype.h nbo.h
icoord.o: icoord.cpp icoord.h zstruct.h
mm_grad.o: mm_grad.cpp icoord.h
mopac.o: mopac.cpp mopac.h constants.h
geombasis.o: geombasis.cpp geombasis.h icoord.h
mem.o: mem.cpp icoord.h
opt.o: opt.cpp icoord.h
pTable.o: pTable.cpp pTable.h
print.o: print.cpp icoord.h
stringtools.o: stringtools.cpp stringtools.h
utils.o: utils.cpp utils.h 
iccomp.o: iccomp.cpp zstruct.h utils.h
write.o: write.cpp zstruct.h utils.h constants.h
knnr.o: knnr.cpp knnr.h utils.h
rtype.o: rtype.cpp rtype.h utils.h stringtools.h
tm.o: tm.cpp zstruct.h constants.h
gaussian.o: gaussian.h gaussian.cpp constants.h utils.h stringtools.h pTable.h
nbo.o: nbo.h nbo.cpp zstruct.h
