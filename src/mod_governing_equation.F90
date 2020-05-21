module mod_governing_equation
  use mod_common
  ! use mod_wave
  use mod_fftw
  use misk
  implicit none

  public eqn_nonlinear
  public eqn_Del_time
  public eqn_aspect_ratio
  public eqn_mu

contains

! *****************************************************************************

  subroutine eqn_nonlinear(term, time)
    implicit none
    real(8), intent(inout) :: term(1:NL-1, 1:NK-1)
    real(8), intent(in) :: time

    real(8) :: aspect
    real(8) :: a_sqrt
    real(8) :: a_sqrt_re

    real(8) :: P(1:NL-1, 1:NK-1)
    real(8) :: U_real(1:NL-1, 1:NK-1)
    real(8) :: V_real(1:NL-1, 1:NK-1)
    real(8) :: dQdX_real(1:NL-1, 1:NK-1)
    real(8) :: dQdY_real(1:NL-1, 1:NK-1)

    real(8) :: mu(1:NL-1, 1:NK-1)
    real(8) :: work(1:NL-1, 1:NK-1)

    integer :: ik, il

    call eqn_aspect_ratio(aspect=aspect, aspect_sqrt=a_sqrt,  &
      &  aspect_sqrt_re=a_sqrt_re, time=time)
    call eqn_mu(mu, aspect)

    P(:,:) = (h(:,:) - Q(:,:)) / mu(:,:)

! U
    do ik = 1, NK-1
      do il = 1, NL-1
        work(il,ik) = - PI_2 * a_sqrt * L(il) * P(il, ik)
      end do
    end do
    call fftw_backward_cs(work, U_real)

! V
    do ik = 1, NK-1
      do il = 1, NL-1
        work(il,ik) = PI_2 * a_sqrt_re * K(ik) * P(il, ik)
      end do
    end do
    call fftw_backward_sc(work, V_real)

! dQ/dx
    do ik = 1, NK-1
      do il = 1, NL-1
        work(il,ik) = PI_2 * a_sqrt_re * K(ik) * Q(il, ik)
      end do
    end do
    call fftw_backward_sc(work, dQdX_real)

! dQ/dy
    do ik = 1, NK-1
      do il = 1, NL-1
        work(il,ik) = PI_2 * a_sqrt * L(il) * Q(il, ik)
      end do
    end do
    call fftw_backward_cs(work, dQdY_real)

    work(:,:) = -(U_real(:,:) * dQdX_real(:,:) + V_real(:,:) * dQdY_real(:,:))
    call fftw_forward_SS(work, term)
    term(:,:) = term(:,:) * Del_time

    return
  end subroutine eqn_nonlinear
  
! *****************************************************************************

  subroutine eqn_Del_time(time)
    implicit none
    real(8), intent(in) :: time
    real(8) :: P(1:NL-1, 1:NK-1)
    real(8) :: U_real(1:NL-1, 1:NK-1)
    real(8) :: V_real(1:NL-1, 1:NK-1)
    real(8) :: U_dot_K(1:NL-1, 1:NK-1)
    real(8) :: U_dot_K_max
    real(8) :: mu(1:NL-1, 1:NK-1)
    real(8) :: work(1:NL-1, 1:NK-1)

    real(8) :: aspect
    real(8) :: a_sqrt
    real(8) :: a_sqrt_re
    real(8) :: K_max
    real(8) :: L_max
    real(8) :: Del_time_CFL

    integer :: ik, il

    call eqn_aspect_ratio(aspect=aspect, aspect_sqrt=a_sqrt,  &
      &  aspect_sqrt_re=a_sqrt_re, time=time)
    call eqn_mu(mu, aspect)
    K_max = NK_truncate * PI * a_sqrt_re
    L_max = NL_truncate * PI * a_sqrt

    P(:,:) = (H(:,:) - Q(:,:)) / mu(:,:)

! U
    do ik = 1, NK-1
      do il = 1, NL-1
        work(il,ik) = - PI_2 * a_sqrt * L(il) * P(il, ik)
      end do
    end do
    call fftw_backward_cs(work, U_real)

! V
    do ik = 1, NK-1
      do il = 1, NL-1
        work(il,ik) = PI_2 * a_sqrt_re * K(ik) * P(il, ik)
      end do
    end do
    call fftw_backward_sc(work, V_real)

    U_dot_K(:,:) = abs(U_real(:,:)) * K_max + abs(V_real(:,:)) * L_max
    U_dot_K_max = maxval(U_dot_K)

    Del_time_CFL = CFL_factor / U_dot_K_max
    Del_time = min(Del_time_max, Del_time_CFL)

    if (Del_time >= time_write - count_time_write) then
      Del_time = time_write - count_time_write
      flag_write = .True.
    end if

    ! call print_main("Del_time", Del_time)

    return
  end subroutine eqn_Del_time

! *****************************************************************************

  subroutine eqn_aspect_ratio(aspect, aspect_sqrt, aspect_sqrt_re, time)
    implicit none
    real(8), intent(out) :: aspect
    real(8), intent(out), optional :: aspect_sqrt
    real(8), intent(out), optional :: aspect_sqrt_re
    real(8), intent(in) :: time

    if (flag_distortion) then
      aspect = (aspect_initial + 1.0d0 / aspect_initial) / 2  &
        &  + (aspect_initial - 1.0d0 / aspect_initial) / 2  &
        &  * cos(2 * PI * time / tau)
    else
      aspect = 1.0d0
    end if

    if (present(aspect_sqrt)) then
      aspect_sqrt = sqrt(aspect)
      if (present(aspect_sqrt_re)) then
        aspect_sqrt_re = 1.0d0 / aspect_sqrt
      end if
    end if

    return
  end subroutine eqn_aspect_ratio

! *****************************************************************************

  subroutine eqn_mu(mu, aspect)
    implicit none
    real(8), intent(out) :: mu(1:NL-1, 1:NK-1)
    real(8), intent(in) :: aspect

    integer :: ik, il

    do ik = 1, NK-1
      do il = 1, NL-1
        mu(il,ik) = PI**2 * K(ik)**2 / aspect + PI**2 * L(il)**2 * aspect  &
          &  + Rossby_radius_re ** 2
      end do
    end do

    return
  end subroutine eqn_mu

end module mod_governing_equation
