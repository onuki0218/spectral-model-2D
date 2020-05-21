module misk
    use mod_common
    implicit none

    public assert
    public print_main

  contains

! *****************************************************************************

  subroutine assert(flag, message)
      implicit none
      logical, intent(in) :: flag
      character(*), intent(in) :: message

      if (.not. flag) then
        write(unit_stderr, *) message
        stop
      end if

      return
    end subroutine assert

! *****************************************************************************

    subroutine print_main(message, value)
      implicit none
      character(*), intent(in) :: message
      real(8), intent(in), optional :: value

      if (my_rank == 0) then
        write(unit_stdout, *) message
        if (present(value)) write(unit_stdout, *) value
      end if

      return
    end subroutine print_main

  end module misk
