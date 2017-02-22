CC = g++
FC = gfortran

all: FortranFromC.o splineLib.a

FortranFromC.o: FortranFromC.cpp
	$(CC) -fpic -c FortranFromC.cpp $(FLAGS)

splineLib.a: css.o dchold.o dlv.o dsetup.o ifind.o sortm.o
	ar -r splineLib.a css.o dchold.o dlv.o dsetup.o ifind.o sortm.o

css.o: css.f
	$(FC) -c css.f $(FLAGS)

dchold.o: dchold.f
	$(FC) -c dchold.f $(FLAGS)

dlv.o: dlv.f
	$(FC) -c dlv.f $(FLAGS)

dsetup.o: dsetup.f
	$(FC) -c dsetup.f $(FLAGS)

ifind.o: ifind.f
	$(FC) -c ifind.f $(FLAGS)

sortm.o: sortm.f
	$(FC) -c sortm.f $(FLAGS)

clean:
	rm -rf *o *a