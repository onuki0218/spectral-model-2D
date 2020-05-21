module mod_core
  use mod_common
  ! use mod_wave
  use mod_io
  use mod_runge_kutta
  use mod_governing_equation
  use mod_initial_condition
  use mod_dealias
  use mod_parallel
  use mod_fftw
  use mod_log
  use misk
  implicit none

  public init
  public step
  public time_increment
  public finalize

contains

! *****************************************************************************

  subroutine init
    implicit none
    character(100) :: message

    call log_time_begin

    call io_get_environment("setting", file_setting)
    call io_get_environment("dir_name_result", dir_name_result)

    call io_read_setting
    call common_dimensions
    call dealias_prepare_two_three
  
! Parallel 
    call parallel_init
    call print_main("Starting with setting file: "//file_setting)
    call fftw_init

    if (my_rank == 0) call log_setting

    if (flag_restart) then 
      time_current = file_number_current * time_write

      call io_get_environment("dir_name_input", dir_name_input)
      call io_read_all
      call io_write_all

      if (my_rank == 0) call log_file_number

      write(message, '("Integration is restarting at file number ", i4.4)')  &
        &  file_number_current
      call print_main(message)

    else

      file_number_current = 0
      time_current = 0.0d0
      call initial_condition_default
      
      call io_write_all

      if (my_rank == 0) call log_file_number

      write(message, '("Integration is starting from a defalut condition")')
      call print_main(message)

    end if

    call log_time_energy

    return 
  end subroutine init 

! *****************************************************************************

  subroutine step
    implicit none

    call eqn_Del_time(time_current)

    call runge_kutta_3
    call analysis_cal_energy(time=0.0d0)

    return
  end subroutine step

! *****************************************************************************

  subroutine time_increment
    implicit none
    character(100) :: message

    step_current = step_current + 1
    time_current = time_current + Del_time
    count_time_write = count_time_write + Del_time

    if (flag_write) then 

      file_number_current = file_number_current + 1

      if (mod(file_number_current, interval_write_budget) == 0) then
        call io_write_real_sum
        write(message, '("*** Output budgets ***")')
        call print_main(trim(message))
      end if

      if (mod(file_number_current, interval_write_variable) == 0) then
        call io_write_all
        write(message, '("*** Output variables ***")')
        call print_main(trim(message))
      end if

      flag_write = .False.
      count_time_write = 0.0d0
      if (my_rank == 0) call log_file_number

    end if

    call log_time_energy

    return
  end subroutine time_increment
  
! *****************************************************************************

  subroutine finalize
    implicit none
    character(100) :: message

    ! call fftw_finalize_parallel
    call parallel_finalize
    write(message,  &
      &  '("Integration is completed with MPI error code:", i4.4)') ierr
    call print_main(message)

    return
  end subroutine finalize  

end module mod_core
