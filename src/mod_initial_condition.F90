module mod_initial_condition
  use mod_common
  use mod_fftw
  use mod_dealias
  use mod_governing_equation
  use mod_analysis
  use misk
  implicit none

  public initial_condition_default

contains

! *****************************************************************************

  subroutine initial_condition_default
    implicit none
    real(8) :: random_homogeneous(1:2, 1:NL-1, 1:NK-1)
    real(8) :: random_gauss(1:NL-1, 1:NK-1)
    integer, allocatable :: seed(:)
    integer :: is, seedsize

    real(8) :: E(1:NL-1, 1:NK-1)
    real(8) :: ave_E(1:NL-1, 1:NK-1)
    real(8) :: aspect
    real(8) :: mu(1:NL-1, 1:NK-1)

    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    do is = 1, seedsize
      call system_clock(count=seed(is))
      seed(is) = seed(is) * (my_rank + 1)
    end do
    call random_seed(put=seed(:))

    call random_number(random_homogeneous)
    random_gauss(:,:) = sqrt(-2 * log(random_homogeneous(1,:,:)))  &
      &  * cos(2 * PI * random_homogeneous(2,:,:))

    call eqn_aspect_ratio(aspect=aspect, time=0.0d0)
    call eqn_mu(mu, aspect)

    Q(:,:) = random_gauss(:,:) * mu(:,:) ** alpha
    call dealias_truncate

    call analysis_cal_energy(E, time=0.0d0)
    call analysis_average(ave_E, E, time=0.0d0)
    Q(:,:) = Q(:,:) / sqrt(sum(ave_E))

    production(:,:) = 0.0d0

    return
  end subroutine initial_condition_default  

end module mod_initial_condition
