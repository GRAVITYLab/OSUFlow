!---------------------------------------------------------------------------
!
! diy Fortran example
!
! Tom Peterka
! Argonne National Laboratory
! 9700 S. Cass Ave.
! Argonne, IL 60439
! tpeterka@mcs.anl.gov
!
! (C) 2011 by Argonne National Laboratory.
! See COPYRIGHT in top-level directory.
!
!--------------------------------------------------------------------------
!
program example

  ! set some default parameters
  dim = 3;
  tot_blocks = 64;

  ! start MPI
  call MPI_INIT(ierror);

  ! start DIY
   call DIY_begin(dim, tot_blocks);

  ! do some work

  ! end DIY
  call DIY_end();

  ! end MPI
   call MPI_Finalize(ierror);

end program example
