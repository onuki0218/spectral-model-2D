module mod_io
  use mod_common
  use mod_analysis
  use misk
  implicit none
  
  public
  
! directory
  character(60), save :: dir_name_result = '../data/data_default'
  character(60), save :: dir_name_input = '../data/data_default'
  
! control file
  character(30), save :: dir_setting = '../data/setting/'
  character(30), save :: file_setting = 'setting'
  
! Energy budget file
!  character(30), save :: file_budget_out = 'budget.out'
!  character(30), save :: file_budget_csv = 'budget.csv'
  
! log file
  character(30), save :: file_number_out = 'file_number.out'
  
  public io_get_environment
  public io_read_setting
  public io_write_real_sum
  public io_read_all
  public io_write_real
  public io_read_real

contains
  
! *****************************************************************************
  
  subroutine io_get_environment(name, value)
    implicit none
    character(*), intent(in) :: name
    character(*), intent(out) :: value
    integer len, status
    intrinsic get_environment_variable
  
    call get_environment_variable(name, status=status, length=len,  &
      &  value=value)
    call assert(status == 0, "Error for environment variable "//name)
    value = trim(value)
  
    return
  end subroutine io_get_environment

  subroutine io_read_setting
    implicit none
    character(90) :: file_tmp
  
    namelist /nml_dimension/ NK, NL
    namelist /nml_time/ time_end, CFL_factor, Del_time_max
    namelist /nml_initial/ alpha
    namelist /nml_output/ time_write, interval_write_variable
    namelist /nml_budget/ interval_write_budget
    namelist /nml_distortion/ flag_distortion, aspect_initial, tau
    namelist /nml_restart/ flag_restart, file_number_current

    file_tmp = trim(dir_setting)//trim(file_setting)//'.in'

    open(unit=unit_setting, file=file_tmp, status='old', action='read')

    read(unit=unit_setting, nml=nml_dimension)
    read(unit=unit_setting, nml=nml_time)
    read(unit=unit_setting, nml=nml_output)
    read(unit=unit_setting, nml=nml_budget)
    read(unit=unit_setting, nml=nml_distortion)
    read(unit=unit_setting, nml=nml_restart)

    close(unit_setting)  
  
! setting dimensions
    NK_truncate = int(NK * 2.0d0 / 3)
    NL_truncate = int(NK * 2.0d0 / 3)

		return
  end subroutine io_read_setting

! *****************************************************************************

  subroutine io_write_real_sum
    implicit none

    real(8) :: E(1:NL-1, 1:NK-1)
    real(8) :: sum_E
    ! real(8) :: sum_nl_transfer_E
    real(8) :: sum_production
    ! real(8) :: sum_dissipation_E  
    integer :: N_column
    character(90) :: file_tmp
    integer :: access

    call analysis_cal_energy(E, time=time_current)
    sum_E = sum(E(:,:))
		sum_production = sum(production)

    write(file_tmp,'(A, A, i4.4, ".csv")')  &
      &  trim(dir_name_result), 'budget', my_rank

    N_column = 3 !

    if (file_number_current == 0 .or. access(file_tmp, " ") /= 0) then

      open(unit=unit_budget_csv, file=file_tmp,  &
        &  status='replace', form='formatted')

      write(unit_budget_csv, '(A)') 'T,E,PE'

    else

      open(unit=unit_budget_csv, file=file_tmp,  &
        &  status='old', form='formatted', position='append')

    end if

    write(unit_budget_csv, '(e16.8, ",", e16.8, ",", e16.8)')  &
      &  time_current,          &
      &  sum_E,                 &
!      &  sum_nl_transfer_E,    &
      &  sum_production
!      &  sum_dissipation_E

    close(unit_budget_csv)

    write(file_tmp,'(A, A, i4.4, ".out")')  &
      &  trim(dir_name_result), 'budget', my_rank

    open(unit=unit_budget, file=file_tmp,  &
      &  status='unknown', form='unformatted', access='direct',  &
      &  recl=record_factor * N_column)

    write(unit_budget, rec=file_number_current+1)  &
      &  real(time_current),              &
      &  real(sum_E),                     &
!      &  real(sum_nl_transfer_E),         &
      &  real(sum_production)
!      &  real(sum_dissipation_E)

    close(unit_budget)

    return
  end subroutine io_write_real_sum

! *****************************************************************************

  subroutine io_write_all
    implicit none
  
    call io_write_real(dir_name_result, 'Q')
    call io_write_real(dir_name_result, 'PE')
    call io_write_energy_ave
    
    return
  end subroutine io_write_all

! *****************************************************************************

  subroutine io_read_all
    implicit none

    call io_read_real(dir_name_input, 'Q')
    call io_read_real(dir_name_input, 'PE')

    return
  end subroutine io_read_all

! *****************************************************************************

  subroutine io_write_real(dir_name, var_name)
    implicit none
    character(*), intent(in) :: dir_name
    character(*), intent(in) :: var_name
    real(4), allocatable :: R_data(:,:)  ! Single precision
    integer :: N_record_length
    integer,dimension(MPI_STATUS_SIZE) :: ista
    character(50) :: file_tmp
    integer :: i

		N_record_length = NK_truncate * NL_truncate
		allocate(R_data(1:NL_truncate, 1:NK_truncate))
!
    if (var_name == 'Q') R_data(:,:) = real(Q(1:NL_truncate, 1:NK_truncate))
		if (var_name == 'PE')  &
		  &  R_data(:,:) = real(Production(1:NL_truncate, 1:NK_truncate))

    write(file_tmp,'(A, A, i4.4, "_", i6.6, ".out")')  &
      &  trim(dir_name), trim(var_name), my_rank, file_number_current

    open(unit=unit_real, file=file_tmp, access='direct',  &
      &  status='unknown', form='unformatted',  &
      &  recl=record_factor * N_record_length)

    write(unit_real, rec=1) R_data

		close(unit_real)

    deallocate(R_data)

    return
  end subroutine io_write_real

! *****************************************************************************

  subroutine io_read_real(dir_name, var_name)
    implicit none
    character(*), intent(in) :: dir_name
    character(*), intent(in) :: var_name
    real(8), allocatable :: R_all(:,:)
    real(4), allocatable :: R_data(:,:)

    integer :: N_record_length
    integer,dimension(MPI_STATUS_SIZE) :: ista
    character(50) :: file_tmp
    integer :: i

    N_record_length = NK_truncate * NL_truncate
		allocate(R_all(1:NL-1, 1:NK-1))
		allocate(R_data(1:NL_truncate, 1:NK_truncate))

    write(file_tmp,'(A, A, i4.4, "_", i6.6, ".out")')  &
      &  trim(dir_name), trim(var_name), my_rank, file_number_current

    open(unit=unit_real, file=file_tmp, access='direct',  &
      &  status='unknown', form='unformatted',  &
      &  recl=record_factor * N_record_length)

    read(unit_real, rec=1) R_data

		close(unit_real)

		R_all(:,:) = 0.0d0
		R_all(1:NL_truncate, 1:NK_truncate) = R_data(:,:)
    if (var_name == 'Q') Q(:,:) = R_all(:,:)
		if (var_name == 'PE') Production(:,:) = R_all(:,:)
!
    deallocate(R_data)
!
    return
  end subroutine io_read_real

! *****************************************************************************

  subroutine io_write_energy_ave()
    implicit none
    real(8) :: E(1:NL-1,1:NK-1)
    real(8) :: E_ave(1:NL-1,1:NK-1)
    real(4), allocatable :: R_data(:,:)  ! Single precision
    integer :: N_record_length
    integer,dimension(MPI_STATUS_SIZE) :: ista
    character(50) :: file_tmp
    integer :: i

!
    call analysis_cal_energy(E, time=time_current)
    call analysis_average(E_ave, E, time=time_current)

    if (my_rank == 0) then

      N_record_length = NK_truncate * NL_truncate
      allocate(R_data(1:NL_truncate, 1:NK_truncate))
      R_data(:,:) = real(E_ave(1:NL_truncate, 1:NK_truncate))

      write(file_tmp,'(A, A, i6.6, ".out")')  &
        &  trim(dir_name_result), 'E', file_number_current
      open(unit=unit_real, file=file_tmp, access='direct',  &
        &  status='unknown', form='unformatted',  &
        &  recl=record_factor * N_record_length)
      write(unit_real, rec=1) R_data
      close(unit_real)
      
      deallocate(R_data)

    end if

    return
  end subroutine io_write_energy_ave

end module mod_io
