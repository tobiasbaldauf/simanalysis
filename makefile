#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/cosmos/users/dc-bald1/lib/HDF5-GNU/lib
source=power_halo
source=bispect_2lpt
source=bispect_mmh_init

AUX1 = eval
AUX2= read_snap

F90 = gfortran
F90FLAGS = -O3 -mcmodel=medium -fopenmp -ffixed-line-length-none -ffree-line-length-none


HDF5 = /nfs/software/apps/hdf5-1.8.13-intel-15.0.4-serial
HDF5 = /home/cosmos/users/dc-bald1/lib/HDF5-GNU

mod1 = -I${INCLUDE}
mod2 = -I${HDF5}/include
mods = ${mod1} ${mod2}

lib1 = -lfftw3_threads -lfftw3 -lm -lpthread

LIBZ    = ${HDF5}/lib/libz.a
LIBSZ   = ${HDF5}/lib/libsz.a
lib2 = -L${HDF5}/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lhdf5_hl -lz

libs = ${lib1} ${lib2}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

makeall: aux prime

prime: ${source}.f90
	${F90} ${mods} ${F90FLAGS} ${source}.f90 -o ${source}.exe ${AUX1}.o ${AUX2}.o ${libs}

aux: 
	${F90} ${F90FLAGS} ${mods} -c ${AUX1}.f90 ${libs}
	${F90} ${F90FLAGS} ${mods} -c ${AUX2}.f90

remap:
	${F90} ${F90FLAGS} remap.f90 -o remap.exe
clean:
	rm -f *.o
	rm -f *.exe
	rm -f *.mod
