CFLAGS = -I/opt/system/app/intel/2022.3.1/mpi/2021.7.1/bin/mpiifort

#FC              = mpifrtpx
FC              = mpiifort
#FC              = ftn

#FFLAGS          = -Qt -x3000 -Kfast,visimpact,nomfunc,noalias=s,fsimple,prefetch_indirect,prefetch_strong,parallel,optmsg=2,array_private
#FFLAGS          =  -O3 -no-prec-div -static -fp-model fast=2 -align array64byte -qopt-prefetch=5 -qopt-report=5
FFLAGS          =

LINKER          = $(FC)

#SUB             = subpack.o

OBJ0            = parareal_params.o
OBJ1            = global_variables.o
OBJ2            = constants.o
OBJ3            = flags.o
OBJ6            = OML_currents.o
OBJ7            = parareal_mpi.o
OBJ8            = run_parareal_mpi.o
OBJ9            = run_seq1.o
OBJ10           = run_seq2.o

PROGRAM0        = parareal_mpi_exe
PROGRAM1        = seq1_exe
PROGRAM2        = seq2_exe

all:            $(PROGRAM0) $(PROGRAM1) $(PROGRAM2)

# parareal_params
$(OBJ0):        %.o : %.f90
		$(FC) $(CFLAGS) $(FFLAGS) -c $<

# global_variables
$(OBJ1):        %.o : %.f90 
		$(FC) $(CFLAGS) $(FFLAGS) -c $<

# constants
$(OBJ2):        %.o : %.f90 
		$(FC) $(CFLAGS) $(FFLAGS) -c $<

# flags   
$(OBJ3):        %.o : %.f90 
		$(FC) $(CFLAGS) $(FFLAGS) -c $<

# OML_currents
$(OBJ6):        %.o : %.f90 $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3)
		$(FC) $(CFLAGS) $(FFLAGS) -c $<

# parareal_mpi
$(OBJ7):        %.o : %.f90 $(OBJ1) $(OBJ3) $(OBJ6)
		$(FC) $(CFLAGS) $(FFLAGS) -c $<

# run_parareal_mpi
$(OBJ8):        %.o : %.f90 $(OBJ0) $(OBJ1) $(OBJ3) $(OBJ6) $(OBJ7)
		$(FC) $(CFLAGS) $(FFLAGS) -c $<

# run_seq1
$(OBJ9):        %.o : %.f90 $(OBJ1) $(OBJ6)
		$(FC) $(CFLAGS) $(FFLAGS) -c $<
   
# run_seq2
$(OBJ10):        %.o : %.f90 $(OBJ1) $(OBJ6)
		$(FC) $(CFLAGS) $(FFLAGS) -c $<

# make execution program
$(PROGRAM0):    $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ6) $(OBJ7) $(OBJ8)
		$(LINKER) $(CFLAGS) $(FFLAGS) -o $@ $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ6) $(OBJ7) $(OBJ8)

$(PROGRAM1):    $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ6) $(OBJ7) $(OBJ9)
		$(LINKER) $(CFLAGS) $(FFLAGS) -o $@ $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ6) $(OBJ7) $(OBJ9)

$(PROGRAM2):    $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ6) $(OBJ7) $(OBJ10)
		$(LINKER) $(CFLAGS) $(FFLAGS) -o $@ $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ6) $(OBJ7) $(OBJ10) 

clean:
		rm -f *.o *.mod *.lst *.optrpt *.dat *.out *.err *.png $(PROGRAM0) $(PROGRAM1) $(PROGRAM2)
