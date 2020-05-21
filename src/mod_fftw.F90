module mod_fftw
  use, intrinsic :: iso_c_binding
  use OMP_LIB
  use mod_common
  use misk
  implicit none
  private
  include 'fftw3.f03'

  integer :: ierr_fftw

  public fftw_init
  public fftw_forward_SS
  public fftw_backward_SS
  public fftw_backward_CS
  public fftw_backward_SC

contains

! *****************************************************************************

  subroutine fftw_init
    implicit none

    ierr_fftw = fftw_init_threads()
    call assert(ierr_fftw /= 0, "Error: FFTW initialization failed.")
    call FFTW_PLAN_WITH_NTHREADS(omp_get_num_threads())

    return
  end subroutine fftw_init

! *****************************************************************************

  subroutine fftw_forward_SS(Q_in, Q_out)
    implicit none
  
    real(8), intent(in) :: Q_in(1:NL-1, 1:NK-1)
    real(8), intent(out) :: Q_out(1:NL-1, 1:NK-1)

    integer(c_int) :: NK_c
    integer(c_int) :: NL_c
    type(c_ptr) :: plan
    real(c_double), allocatable :: Q1(:,:)
    real(c_double), allocatable :: Q2(:,:)

    NK_c = NK
    NL_c = NL
    allocate(Q1(1:NL_c-1, 1:NK_c-1))
    allocate(Q2(1:NL_c-1, 1:NK_c-1))

    plan = fftw_plan_r2r_2d(NK_c-1, NL_c-1, Q1, Q2,  &
      &  FFTW_RODFT00, FFTW_RODFT00, FFTW_MEASURE)

    Q1(:,:) = Q_in(:,:)
    call fftw_execute_r2r(plan, Q1, Q2)
    call fftw_destroy_plan(plan)
    Q_out(:,:) = Q2(:,:) / (2 * NK * NL)

    return
  end subroutine fftw_forward_SS

! *****************************************************************************

  subroutine fftw_backward_SS(Q_in, Q_out)
    implicit none
    real(8), intent(in) :: Q_in(1:NL-1, 1:NK-1)
    real(8), intent(out) :: Q_out(1:NL-1, 1:NK-1)

    type(c_ptr) :: plan
    integer(c_int) :: NK_c
    integer(c_int) :: NL_c
    real(c_double), allocatable :: Q1(:,:)
    real(c_double), allocatable :: Q2(:,:)

    NK_c = NK
    NL_c = NL
    allocate(Q1(1:NL_c-1, 1:NK_c-1))
    allocate(Q2(1:NL_c-1, 1:NK_c-1))

    plan = fftw_plan_r2r_2d(NK_c-1, NL_c-1, Q1, Q2,  &
      &  FFTW_RODFT00, FFTW_RODFT00, FFTW_MEASURE)

    Q1(:,:) = Q_in(:,:)
    call fftw_execute_r2r(plan, Q1, Q2)
    call fftw_destroy_plan(plan)
    Q_out(:,:) = Q2(:,:) / 2

    return
  end subroutine fftw_backward_SS

! *****************************************************************************

  subroutine fftw_backward_CS(Q_in, Q_out)
    implicit none
    real(8), intent(in) :: Q_in(1:NL-1, 1:NK-1)
    real(8), intent(out) :: Q_out(1:NL-1, 1:NK-1)

    type(c_ptr) :: plan
    integer(c_int) :: NK_c
    integer(c_int) :: NL_c
    real(c_double), allocatable :: Q1(:,:)
    real(c_double), allocatable :: Q2(:,:)

    NK_c = NK
    NL_c = NL
    allocate(Q1(0:NL_c, 1:NK_c-1))
    allocate(Q2(0:NL_c, 1:NK_c-1))

    plan = fftw_plan_r2r_2d(NK_c-1, NL_c+1, Q1, Q2,  &
      &  FFTW_RODFT00, FFTW_REDFT00, FFTW_MEASURE)

    Q1(:,:) = 0.0d0
    Q1(1:NL_c-1, 1:NK_c-1) = Q_in(:,:)
    call fftw_execute_r2r(plan, Q1, Q2)
    call fftw_destroy_plan(plan)

    Q_out(:,:) = Q2(1:NL_c-1, 1:NK_c-1) / 2

    return
  end subroutine fftw_backward_CS

! *****************************************************************************

  subroutine fftw_backward_SC(Q_in, Q_out)
    implicit none
    real(8), intent(in) :: Q_in(1:NL-1, 1:NK-1)
    real(8), intent(out) :: Q_out(1:NL-1, 1:NK-1)
    type(c_ptr) :: plan
    integer(c_int) :: NK_c
    integer(c_int) :: NL_c
    real(c_double), allocatable :: Q1(:,:)
    real(c_double), allocatable :: Q2(:,:)

    NK_c = NK
    NL_c = NL
    allocate(Q1(1:NL_c-1, 0:NK_c))
    allocate(Q2(1:NL_c-1, 0:NK_c))

    plan = fftw_plan_r2r_2d(NK_c+1, NL_c-1, Q1, Q2,  &
      &  FFTW_REDFT00, FFTW_RODFT00, FFTW_MEASURE)

    Q1(:,:) = 0.0d0
    Q1(1:NL_c-1, 1:NK_c-1) = Q_in(:,:)
    call fftw_execute_r2r(plan, Q1, Q2)
    call fftw_destroy_plan(plan)

    Q_out(:,:) = Q2(1:NL_c-1, 1:NK_c-1) / 2

    return
  end subroutine fftw_backward_SC

end module mod_fftw
