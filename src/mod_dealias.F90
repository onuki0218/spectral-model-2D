module mod_dealias
  use, intrinsic :: iso_c_binding
  use mod_common
  implicit none

  public dealias_prepare_two_three
  public dealias_truncate

contains

! *****************************************************************************

  subroutine dealias_prepare_two_three
    implicit none
    integer :: ik, il

    mask(1:NL-1, 1:NK-1) = 1
    mask(NL_truncate:NL-1, :) = 0
    mask(:, NK_truncate:NK-1) = 0

    return 
  end subroutine dealias_prepare_two_three

! *****************************************************************************

  subroutine dealias_truncate
    implicit none

    Q(:,:) = Q(:,:) * mask(:,:)
  
    return
  end subroutine dealias_truncate

end module mod_dealias
