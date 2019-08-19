module temperature
  use types
  use cmd_io
  implicit none

contains
  ! --------------------Public Contents-----------------------------------------
  ! SUBROUTINE: adjust_T(Total_MC_Sweeps)
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  !
  ! ----------------------------------------------------------------------------
  subroutine adjust_T(n_sweeps)
    ! SUBROUTINE adjust_T(Total_MC_Sweeps)
    integer :: n_sweeps                                     !Total number of MC sweeps
    select case(trim(cool_opt))
    case('constant')
       T = T
    case('linear')
       T  = T - T_param
    case('exp')
       T = T*T_param
    case('tlinear')
       T = T0 - T_param*n_sweeps
    case('tlog')
       T = T_param/log(real(n_sweeps + 1))
    case('texp')
       T = T0*T_param**n_sweeps
    case('default')
       write(*,'(A80)') 'Failure in calculating new Temperature. Invalid cooling schedule option.'
       stop
    end select
    if (T<min_T) then
       T = min_T
    end if
  end subroutine adjust_T
  ! -----------------------------------------------------------------------------

  ! function calc_T_intervals() result(r)
  !   ! FUNCTION: calc_T_intervals()
  !   ! Calculated the number of temperature intervals required to move from T0 to Tmin.
  !   ! Cooling schedules are presently built w.r.t. number of sweeps not number of Temperature steps
  !   ! Excpet in the case of contsant cooling
  !   integer :: r                                            !Number of Temperature steps to minimum

  !   select case(trim(cool_opt))
  !   case('constant')
  !      r = int(-1.*(min_T - T0)/T_param)
  !   case('linear')
  !      r = int((-1.*(min_T - T0)/T_param)/max_sweeps)
  !   case('log')
  !      r = int((exp(T_param/min_T) - 1)/max_sweeps)
  !   case('exp')
  !      r = int((log(T0/min_T)/log(T_param))/max_sweeps)
  !   case('default')
  !      write(*,'(A80)') 'Failure in calculating number of temperature intervals. Invalid cooling schedule option.'
  !      stop
  !   end select
  ! end function calc_T_intervals
  ! -----------------------------------------------------------------------------
end module temperature
