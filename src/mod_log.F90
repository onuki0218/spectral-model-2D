module mod_log
  use, intrinsic :: iso_c_binding
  use mod_common
  use mod_io
  use mod_analysis
  use misk
  implicit none

  private

  public log_time_begin
  public log_time
  public log_energy
  public log_time_energy
  public log_file_number
  public log_setting

  contains

! *****************************************************************************

  subroutine log_time_begin
    implicit none

    if (my_rank == 0) then
      call cpu_time(time_begin_system)
      call system_clock(count=time_begin_system_count,  &
        &  count_rate=system_count_per_second)
    endif

    return 
  end subroutine log_time_begin

! *****************************************************************************

  subroutine log_time
    implicit none
  
    if (my_rank == 0) then
  
      call cpu_time(time_system)
      call system_clock(count=time_system_count)
  
      write(unit_stdout,*) 'STEP:', step_current
      write(unit_stdout,*) '  Calculation time:', time_current
      write(unit_stdout,*) '    Until', time_end
      write(unit_stdout,*) '  Real time:',  &
        &  (time_system_count - time_begin_system_count)  &
        &  / system_count_per_second
      write(unit_stdout,*) '  CPU time:', time_system - time_begin_system

    endif
   
    return 
  end subroutine log_time
  
! *****************************************************************************

  subroutine log_energy
    implicit none
    real(8) :: E(1:NL-1, 1:NK-1)
    real(8) :: sum_E

    call analysis_cal_energy(E, time=time_current)
    sum_E = sum(E(:,:))
    if (my_rank == 0) then
      write(unit_stdout,*) '  Total energy:', sum_E
      write(unit_stdout,*) '  Total enstrophy:', sum(Q(:,:)**2) / 2
    end if

    return 
  end subroutine log_energy

! *****************************************************************************

  subroutine log_time_energy
    implicit none

    call log_time
    call log_energy

    return
  end subroutine log_time_energy

! *****************************************************************************

  subroutine log_file_number
    implicit none
    character(90) :: file_tmp
    character(90) :: message
    integer :: access

    file_tmp = trim(dir_name_result)//trim(file_number_out)

    if (file_number_current == 0 .or. access(file_tmp, " ") /= 0) then

      open(unit=unit_log, file=file_tmp,  &
        &  status='replace', form='formatted')

    else

      open(unit=unit_log, file=file_tmp, status='old',  &
        &  form='formatted', position='append')

    end if

    if (flag_file_number_header) then
      write(unit_log, *) file_setting
      flag_file_number_header = .False.
    end if

    write(unit_log,'(i7)') file_number_current

    close(unit_log)

    return
  end subroutine log_file_number

! *****************************************************************************

  subroutine log_setting
    implicit none
    character(90) :: file_tmp

    file_tmp = trim(dir_name_result)//trim(file_setting)//'.out'

    open(unit=unit_log, file=file_tmp, status='replace', form='formatted')
 
    write(unit_log, '(A)') "&"//trim(file_setting)

    write(unit_log, '("NK=                      ",I5)') NK
    write(unit_log, '("NL=                      ",I5)') NL
    write(unit_log, '("NK_truncate=              ",I5)') NK_truncate
    write(unit_log, '("N_process=               ",I5)') N_process
    write(unit_log, '("N_output_file=           ",I5)') N_output_file
    write(unit_log, '("length_domain_X=         ",e16.8)') length_domain_X
    write(unit_log, '("length_domain_Y=         ",e16.8)') length_domain_Y
    write(unit_log, '("time_end=                ",e16.8)') time_end
!    write(unit_log, '("CFL_factor=              ",e16.8)') CFL_factor
    write(unit_log, '("time_write=              ",e16.8)') time_write
    write(unit_log, '("interval_write_budget=   ",I5)') interval_write_budget
    write(unit_log, '("interval_write_variable= ",I5)') interval_write_variable
!    write(unit_log, '("viscosity=               ",e16.8)') viscosity

    write(unit_log, '(A)') "/"
    write(unit_log, *)

    close(unit_log)
  
  end subroutine log_setting

end module mod_log
