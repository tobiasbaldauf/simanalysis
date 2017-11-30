
source=power_halo
source=bispect

AUX1 = eval
AUX2= read_snap

F90 = gfortran
F90FLAGS = -O3 -mcmodel=medium -fopenmp -ffixed-line-length-none -ffree-line-length-none


mod1 = -I${INCLUDE}
mods = ${mod1}

lib1 = -lfftw3_threads -lfftw3 -lm -lpthread
libs = ${lib1}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

makeall: aux prime

prime: ${source}.f90
	${F90} ${mods} ${F90FLAGS} ${source}.f90 -o ${source}.exe ${AUX1}.o ${AUX2}.o ${libs}

aux: 
	${F90} ${F90FLAGS} ${mods} -c ${AUX1}.f90 ${libs}
	${F90} ${F90FLAGS} ${mods} -c ${AUX2}.f90

clean:
	rm -f *.o
	rm -f *.exe
	rm -f *.mod
