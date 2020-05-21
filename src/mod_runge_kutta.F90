module mod_runge_kutta
  use OMP_LIB
  use mod_common
  use mod_governing_equation
  use mod_fftw
  use mod_dealias
  use misk
  implicit none

  public runge_kutta_3

contains

! *****************************************************************************

  subroutine runge_kutta_3
    implicit none
    real(8), parameter :: factor_A = 8.0d0 / 15
    real(8), parameter :: factor_B = 5.0d0 / 12
    real(8), parameter :: factor_C = 17.0d0 / 60
    real(8), parameter :: factor_D = 3.0d0 / 4
    real(8), parameter :: factor_E = 5.0d0 / 12

    real(8), parameter :: factor_A_rcp = 1.0d0 / factor_A
    real(8), parameter :: factor_B_rcp = 1.0d0 / factor_B
    real(8), parameter :: factor_ABC = factor_A + factor_B - factor_C
                        ! = 2.0d0 / 3.0

    real(8) :: time_A
    real(8) :: time_B
    real(8) :: Del_time_short

    real(8) :: term_A(1:NL-1, 1:NK-1)
    real(8) :: term_B(1:NL-1, 1:NK-1)

! Runge-Kutta step 1
    time_A = time_current
    Del_time_short = Del_time * factor_A
    time_B = time_current + Del_time_short

    call eqn_nonlinear(term_A, time_A)
    Q(:,:) = Q(:,:) + term_A(:,:) * factor_A
    call dealias_truncate

! Runge-Kutta step 2
    time_A = time_B
    Del_time_short = Del_time * factor_ABC
    time_B = time_current + Del_time_short

    call eqn_nonlinear(term_B, time_A)
    Q(:,:) = Q(:,:) + term_B(:,:) * factor_B - term_A(:,:) * factor_C
    call dealias_truncate

! Runge-Kutta step 3
    time_A = time_B
    Del_time_short = Del_time
    time_B = time_current + Del_time_short

    call eqn_nonlinear(term_A, time_A)
    Q(:,:) = Q(:,:) + term_A(:,:) * factor_D - term_B(:,:) * factor_E
    call dealias_truncate

    return
  end subroutine runge_kutta_3

end module mod_runge_kutta
