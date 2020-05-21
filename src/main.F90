program main
    use mod_common
    use mod_core
    implicit none
  
    call init
    do while (time_current < time_end)
      call step
      call time_increment
    end do
    call finalize
  
    stop
  end program main
