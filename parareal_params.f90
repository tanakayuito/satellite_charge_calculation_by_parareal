!>
!! Module that reads in and provides parameters for the simulation.
!!
MODULE parareal_params

IMPLICIT NONE

!> If TRUE, the solution will be written out. Otherwise only timings are generated.
LOGICAL :: do_io_first_coarse
LOGICAL :: do_io_fine
LOGICAL :: do_io_cor_coarse

!> If TRUE, Parareal will generate a number of messages during runtime.
LOGICAL ::be_verbose

!conductors


!> Final simulation time
DOUBLE PRECISION :: Tstart,Tend

!> Time steps *per timeslice* for the fine integrator
INTEGER :: N_fine
INTEGER :: N_coarse

INTEGER :: Niter


CONTAINS

  !>
  !! Read parameter namelist from input file provided as first command line argument

END MODULE parareal_params
