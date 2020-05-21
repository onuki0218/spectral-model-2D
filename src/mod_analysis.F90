module mod_analysis
  use mod_common
  use mod_fftw
  use mod_governing_equation
  use misk
  implicit none

  public analysis_cal_energy

contains

! *****************************************************************************

  subroutine analysis_cal_energy(time)
    implicit none
    real(8), intent(in) :: time
    real(8) :: aspect
    real(8) :: mu(1:NL-1, 1:NK-1)

    call eqn_aspect_ratio(aspect=aspect, time=0.0d0)
    call eqn_mu(mu, aspect)

    E(:,:) = Q(:,:)**2 / mu(:,:) * 0.5d0

    return
  end subroutine analysis_cal_energy

end module mod_analysis
