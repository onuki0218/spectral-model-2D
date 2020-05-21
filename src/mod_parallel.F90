module mod_parallel
  use, intrinsic :: iso_c_binding
  use mod_common
  use misk
  implicit none

  public parallel_init
  public parallel_finalize

contains

! *****************************************************************************

  subroutine parallel_init
    implicit none
    integer :: provided_mpi_thread

    call MPI_Init_thread(MPI_THREAD_FUNNELED, provided_mpi_thread, ierr)

    call assert(provided_mpi_thread >= MPI_THREAD_FUNNELED,  &
      &  "ERROR: Provided MPI threading less than required.")

    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

    call MPI_Comm_size(MPI_COMM_WORLD, N_process, ierr)

    return 
  end subroutine parallel_init

! *****************************************************************************

  subroutine parallel_finalize
    implicit none

    call MPI_Finalize(ierr)

    return 
  end subroutine parallel_finalize

end module mod_parallel
