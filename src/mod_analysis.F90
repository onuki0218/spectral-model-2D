module mod_analysis
  use mod_common
  use mod_fftw
  use mod_governing_equation
  use misk
  implicit none

  public analysis_cal_energy

contains

! *****************************************************************************

  subroutine analysis_cal_energy(E, time)
    implicit none
    real(8), intent(out) :: E(1:NL-1, 1:NK-1)
    real(8), intent(in) :: time
    real(8) :: aspect
    real(8) :: mu(1:NL-1, 1:NK-1)

    call eqn_aspect_ratio(aspect=aspect, time=time)
    call eqn_mu(mu, aspect)
    E(:,:) = Q(:,:)**2 / mu(:,:) * 0.5d0

    return
  end subroutine analysis_cal_energy

  subroutine analysis_average(var_ave, var, time)
    implicit none
    real(8), intent(out) :: var_ave(1:NL-1, 1:NK-1)
    real(8), intent(in) :: var(1:NL-1, 1:NK-1)
    real(8), intent(in) :: time

    call MPI_Allreduce(var(1,1), var_ave(:,:), (NK-1)*(NL-1),  &
      &  MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    var_ave(:,:) = var_ave(:,:) / N_process

    return
  end subroutine analysis_average

end module mod_analysis
