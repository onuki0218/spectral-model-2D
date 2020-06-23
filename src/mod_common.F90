module mod_common
  ! use, intrinsic :: iso_c_binding
  use OMP_LIB
  implicit none
  public
  include 'mpif.h'
  ! include 'fftw3-mpi.f03'

! system variables
  integer :: time_begin_system_count
  integer :: time_system_count
  integer :: system_count_per_second
  real(8) :: time_begin_system
  real(8) :: time_system
  integer :: record_factor = 4
    ! depends on compiler
    ! = 4 for gfortran or ifort with option "-assume byterecl"
    ! = 1 for ifort without option
 
 ! IO unit
  integer, parameter :: unit_stderr = 0
  integer, parameter :: unit_stdin = 5
  integer, parameter :: unit_stdout = 6

  integer, parameter :: unit_setting = 10
  integer, parameter :: unit_log = 11

  integer, parameter :: unit_real = 22
  integer, parameter :: unit_budget = 23
  integer, parameter :: unit_budget_csv = 24
 
! flags
  logical, save :: flag_restart = .False.
 
! constants
  real(8), parameter :: PI = 3.1415926535897932384d0
  real(8), parameter :: PI_2 = PI * 2
  ! real(8), parameter :: PI_half = PI / 2
  ! complex(8), parameter :: UI = (0.0d0, 1.0d0)
 
! dimensions - basic parameters
  integer :: NK
  integer :: NL
  integer :: NK_truncate
  integer :: NL_truncate
 
! dimensions - coordinates
  real(8), save :: length_domain_X = PI
  real(8), save :: length_domain_Y = PI
  real(8) :: Del_X
  real(8) :: Del_Y
  real(8) :: wavenumber_truncate
  integer, allocatable :: K(:)
  integer, allocatable :: L(:)
  integer, allocatable :: mask(:,:)  ! for dealiasing
 
! MPI parameters
  integer :: N_process
  integer :: my_rank
  integer :: ierr
 
! Variable output
  integer :: N_output_file
  integer, save :: interval_write_variable = 16

! Time setting
  real(8), save :: time_current = 0.0d0
  real(8) :: time_end
  real(8) :: Del_time  ! change each step
  real(8), save :: Del_time_max = 0.001d0
  ! real(8), save :: CFL_factor = 1.83d0
  real(8), save :: CFL_factor = 0.1d0

! Time output
  real(8) :: time_write = 0.01d0
  real(8), save :: count_time_write = 0.0d0
  integer, save :: step_current = 0
  integer, save :: file_number_current = 0
  integer, save :: interval_write_budget = 1
  logical, save :: flag_write = .False.
  logical, save :: flag_file_number_header = .True.

! - basic
  ! real(8) :: viscosity
  real(8), save :: Rossby_radius_re = 0.0d0

! parameter: initial condition
  real(8) :: total_energy = 1.0d0
  real(8) :: alpha = - 1.0d0
  real(8) :: circulation = 0.0d0

! parameter: distortion
  logical, save :: flag_distortion = .False.
  real(8), save :: aspect_initial = 2.0d0
  real(8) :: tau = 1.0d0

! variables - basic
  real(8), allocatable :: Q(:,:)
  real(8), allocatable :: H(:,:)
  ! complex(8), allocatable :: mu(:,:)

! variables - energy
  ! real(8), allocatable :: E(:,:)
  real(8), allocatable :: production(:,:)
 
contains

! *****************************************************************************
 
   subroutine common_dimensions
     implicit none
     integer :: i
 
! allocate variables
     allocate(Q(1:NL-1, 1:NK-1))
     allocate(H(1:NL-1, 1:NK-1))
     ! allocate(mu(1:NL-1, 1:NK-1))
     ! allocate(E(1:NL-1, 1:NK-1))
     allocate(Production(1:NL-1, 1:NK-1))
 
     allocate(K(1:NK-1))
     allocate(L(1:NL-1))
     allocate(mask(1:NL-1, 1:NK-1))

! set coordinates
     Del_X = length_domain_X / NK
     Del_Y = length_domain_Y / NL

     do i = 1, NK-1
       K(i) = i
     end do

     do i = 1, NL-1
       L(i) = i
     end do

! Topography
     H(:,:) = 0.0d0

     return
   end subroutine common_dimensions

 end module mod_common
