! LAST EDIT: Phil Maffettone 2017-03-18
! New options will need to be added in time
! Some options are outdated and unused. Should be removed for final
! TODO: Program in threshold arg
! TODO: Option for using more than Ca in distance list. Already partly in dist_list.f90
! TODO: Fix cmd_out to refelct MC sweeps and Metroplois moves differently
module cmd_io
  ! Module for retaining desired inputs and outputs of user
  ! All variables are saved to be accessed by any using programs/modules
  ! Depends on filename.pin user input
  use types

  implicit none

  character(len=80), save :: run_title                      !User supplied run title
  character(len=120), save :: structure_fname               !Structural data filename
  character(len=120), save :: exp_fname                     !Experimental data filename
  character(len=80), save :: init_config                    !Initial adjustment to configuration
  character(len=10), save :: exp_opt                        !Experimental Option (XRD or NEUTRON)
  character(len=10), save :: cs_opt                         !Chemical shift option (DANGLE,TALOS,NONE)
  character(len=10), save :: cool_opt                       !Cooling schedule speicification(linear,log,exp)
  character(len=80), save :: cs_dir                         !Chemical shift files directory
  character(len=80), allocatable, save :: output_fname(:)   !Set of output filenames
  character(len=80), save :: dist_flag                      !Flag for inclusion in dist list
  character(len=80), save :: pdf_flag                       !Flag for PDF flavor
  character(len=80), save :: pdf_scat                       !Flag for scattering length
  logical, save :: incl_omega                               !Whether to bend omega
  logical, save :: bbonly                                   !Backbone data to be used exclusively
  logical, save :: annealing                                !Whether to anneal
  logical, save :: cs_only                                  !Chemical shifts only, don't include ramachandran
  logical, save :: dist_gen                                 !Quick generation of distance list then quit
  logical, save :: dl                                       !True for distance list data, false for PDF data
  integer, save :: xye_file                                 !Index in output_fname, 0 for none
  integer, save :: pdb_file                                 !Index in output_fname, 0 for none
  integer, save :: pot_file                                 !Index in output_fname, 0 for none
  integer, save :: ram_file                                 !Index in output_fname, 0 for none
  integer, save :: time_file                                !Index in output_fname, 0 for none
  integer, save :: dist_file                                !Index in output_fname, 0 for none
  integer, save :: therm_file                               !Index in output_fname, 0 for none
  integer, save :: max_sweeps                               !Maximum number of MC sweeps before finish
  integer, save :: n_snap                                   !No of printed snapshots per run in maxsweeps limit (distlist/pdb)
  integer, save :: snap_met_moves                           !Number of metropolis moves between snapshots
  integer, save :: tot_accept                               !Total accepted Metropolis moves
  integer, save :: snap_accept                              !Metropolis Moves accepted in snapshot interval
  integer, save :: mobile_range(2)                          !Range of residues considered mobile during MC
  integer, save :: tau_calc_sweeps                          !Number of sweeps used to calculate autocorrelation
  integer, save :: autocorr_save                            !Number of covariences saved when calculating tau
  real(dp), save :: threshold                               !Chi threshold for satisfaction
  real(dp), save :: T                                       !Monte Carlo Temperature
  real(dp), save :: T0                                      !Initial simulation temperature
  real(dp), save :: max_T                                   !Maximum annealing temperature
  real(dp), save :: min_T                                   !Minimum annealing temperature
  real(dp), save :: T_param                                 !Temperature parameter for cooling schedule
  real(dp), save :: NMR_weight                              !Weighting on NMR-based sampling dist (<0 for all NMR)
  real(dp), save :: pdf_rmin                                !PDF calculation minimum radius
  real(dp), save :: pdf_rmax                                !PDF calculation maximum radius
  real(dp), save :: pdf_dr                                  !PDF bin size
contains

  ! --------------------Public Contents-----------------------------------------
  ! SUBROUTINE: get_input(filename)
  ! SUBROUTINE: write_pot(filename,counter,[continuation,chisquared])
  ! SUBROUTINE: write_time(filename,total_count,[chisquared])
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  !
  ! ----------------------------------------------------------------------------

  subroutine get_input(filename)
    ! SUBROUTINE get_input(filename)
    ! Retreives input from .pin file type to run RMCParty
    ! In general, the keywords are kept in switch statements, and any settings are expected
    ! Directly below keywords.
    character(len=80), intent(in) :: filename               !Filename
    integer :: unit                                         !File unit
    integer :: info                                         !Error information. 0 for success
    character(len=80) :: errmsg                             !Error message
    character(len=120) :: line                              !Read line
    integer :: i_file                                       !Output file counter
    character(len=80) :: output(5)                          !Local array of outputs, max size 5
    character(len=1) :: buffer                              !Character by character buffer
    integer :: i_line                                       !Line counter
    integer :: i                                            !Dummy integer
    real :: power                                           !Power measure for log based inputs

    info = 0
    open(newunit=unit,file=trim(filename), status='old',IOSTAT = info,IOMSG=errmsg)

    ! Quitting immediately if input is unsuccessful
    if(info /=0) then
       write(*,*) errmsg
       stop
    end if

    ! Initializing default values
    bbonly = .true.
    annealing = .false.
    dist_gen = .false.
    cs_only = .false.
    dl = .true.
    NMR_weight = 1.0_dp
    max_T = 1.0_dp
    min_T = 0.0_dp
    T=0.0_dp
    max_sweeps = 1000000
    n_snap = 0
    tau_calc_sweeps = 1000
    autocorr_save = 10
    threshold = 0._dp
    xye_file = 0
    pdb_file = 0
    pot_file = 0
    ram_file = 0
    time_file = 0
    dist_file = 0
    mobile_range = 0
    init_config = 'n/a'
    structure_fname = 'N/A'
    exp_fname = 'N/A'
    exp_opt = 'N/A'
    run_title = 'Run_Title'
    pdf_flag = 'none'
    pdf_rmin=0.
    pdf_rmax=100.
    pdf_dr=0.5
    pdf_scat='XRD'
    if(allocated(output_fname))deallocate(output_fname)

    i_line = 0
    do                                                      !Reading until EOF
       info = 0
       read(unit,'(A120)',iostat = info) line
       i_line = i_line+1
       if (info /= 0) exit
       line = adjustl(line)                                 !Ignoring leading blanks
       if(line(1:1) == '!') cycle
       call lower_case(line)

       ! Select case for keywords, continuous if/then for options
       select case(trim(line))
       case('title')                                        !Line after 'title' is title string
          read(unit,'(A80)',iostat = info) run_title
          i_line=i_line+1
          if(info/=0) run_title = 'Syntax Error in Title'

       case('structure input')
          read(unit,'(A120)',iostat=info)line
          i_line=i_line+1
          if(info==0) then
             structure_fname = trim(adjustl(line))
          else
             write(*,"('Bad input at line ',i3,' for structure input.')") i_line
          end if

       case('initial configuration')                        !Line after is adjustment to init config
          read(unit,'(A120)',iostat = info) line
          i_line=i_line+1
          if(info/=0)then
             init_config = 'n/a'
             write(*,"('Bad input at line ',i3,' for initial configuration.')") i_line
          end if
          call lower_case(line)
          init_config = adjustl(line)
          if (trim(init_config) == 'distance generation') then
             dist_gen = .true.
             init_config = 'n/a'
          end if

       case('begin experiment')                             !Section for experimental data
          do
             info = 0
             read(unit,'(A120)',iostat=info) line
             i_line=i_line+1
             if(info/=0) then
                write(*,"('Bad input at line ',i3,' for experimental data specifications.')")i_line
                info=0
                exit
             end if
             line = adjustl(line)
             call lower_case(line)
             if(line(1:1) == '!') cycle
             if(trim(line) == 'pdf') then
                dl = .false.
             else if (trim(line) == 'dl') then
                dl = .true.
             else if(trim(line) == 'input') then                 !Line after is source file
                read(unit,'(A120)',iostat=info) exp_fname
                i_line=i_line+1
                if(info/=0) then
                   write(*,"('Bad input at line ',i3, ' for experimental data input.')")i_line
                   info=0
                   exit
                end if
             else if(trim(line) == 'sidechain')then
                bbonly = .false.
             else if(trim(line) == 'xrd' .or. trim(line) == 'neutron') then
                exp_opt = trim(line)
             else if(trim(line) == 'flag') then
                read(unit,'(A80)', iostat=info) dist_flag
                 if(info/=0) then
                   write(*,*)'Bad string for distance list flag'; info=0; dist_flag='CA'
                end if
                call lower_case(dist_flag)
             else if(trim(line) == 'end experiment') then
                exit
             else
                write(*,"('Unrecognized input at line ',i3, ' for experimental data specifications.')")i_line
                exit
             end if
          end do

       case('begin pdf')
          ! Section for PDF output details
          do
             info = 0
             read(unit,'(A120)',iostat=info) line
             i_line=i_line+1
             if(info/=0) then
                write(*,"('Bad input at line ',i3,' for experimental data specifications.')")i_line
                info=0
                exit
             end if
             line = adjustl(line)
             call lower_case(line)
             if(line(1:1) == '!') then
                cycle
             else if(line(1:4) == 'flag') then
                pdf_flag = line(5:len(trim(line)))
             else if(line(1:4) == 'scat') then
                pdf_scat = line(5:len(trim(line)))
             else if(line(1:5) == 'r_max') then
                read(line(6:len_trim(line)),*,iostat=info) power
                if(info/=0)then
                   write(*,*) 'Bad value for r_max'; pdf_rmax = 50.; info=0
                else
                   pdf_rmax=power
                end if
             else if(line(1:5) == 'r_min') then
                read(line(6:len_trim(line)),*,iostat=info) power
                if(info/=0)then
                   write(*,*) 'Bad value for r_min'; pdf_rmin = 0.; info=0
                else
                   pdf_rmin=power
                end if
             else if(line(1:2) == 'dr') then
                read(line(3:len_trim(line)),*,iostat=info) power
                if(info/=0)then
                   write(*,*) 'Bad value for dr'; pdf_rmin = 0.; info=0
                else
                   pdf_dr=power
                end if
             else if(trim(line) == 'end pdf') then
                exit
             else
                write(*,"('Unrecognized input at line ',i3, ' for experimental data specifications.')")i_line
                exit
             end if
          end do

       case('begin cs')
          do
             info = 0
             read(unit,'(A120)',iostat=info) line
             i_line=i_line+1
             if(info/=0) then
                write(*,"('Bad input at line ',i3,' for chemical shift specifications.')")i_line
                info=0
                exit
             end if
             line = adjustl(line)
             if(line(1:1) == '!') cycle
             call lower_case(line)
             if(line(1:3) == 'dir') then
                cs_dir = trim(adjustl(line(4:len(line))))
                cycle
             end if
             if(trim(line) == 'dangle') then
                cs_opt = 'dangle'
             else if(trim(line) == 'talos') then
                cs_opt = 'talos'
             else if(trim(line) == 'none') then
                cs_opt = 'none'
             else if(trim(line) == 'only') then
                cs_only = .true.
             else if(line(1:6) == 'weight') then
                read(line(7:len_trim(line)),*,iostat=info) NMR_weight
                if(info/=0) then
                   write(*,*) 'Bad value for NMR based probability weighting'; NMR_weight = 1.0_dp; info=0
                end if
             else if(trim(line) == 'end cs') then
                exit
             else
                write(*,"('Unrecognized input at line ',i3, ' for chemical shift specifications.')")i_line
                exit
             end if
          end do

       case('begin mc')                                     !Section for Monte Carlo Parameters
          do
             info = 0
             read(unit,'(A120)',iostat=info) line
             i_line=i_line+1
             if(info/=0) then
                write(*,*) 'Bad input at line ', i_line, ' for MC specifications.'
                exit
             end if
             line = trim(adjustl(line))
             if(line(1:1) == '!') then
                cycle
             else if(line(1:9) == 'maxsweeps') then
                read(line(11:len_trim(line)),*,iostat=info) power
                if(info/=0)then
                   write(*,*) 'Bad value for maxsweeps'; max_sweeps = 0; info=0
                else
                   max_sweeps = int(10**power)
                end if
             else if(line(1:8) == 'autocorr') then
                read(line(9:len_trim(line)),*, iostat=info) power
                if(info/=0)then
                   write(*,*) 'Bad value for number of autocorrelation sweeps'; tau_calc_sweeps=1000; info=0
                else
                   tau_calc_sweeps = int(power)
                end if
             else if(line(1:8) == 'autosave') then
                read(line(9:len_trim(line)),*, iostat=info) power
                if(info/=0)then
                   write(*,*) 'Bad value for number of autocorrelation sweeps'; autocorr_save=10; info=0
                else
                   autocorr_save = int(power)
                end if
             else if(line(1:8) == 'snapshot') then
                read(line(9:len_trim(line)),*, iostat=info) power
                if(info/=0)then
                   write(*,*) 'Bad value for number of program snapshots'; n_snap=0; info=0
                else
                   n_snap = int(power)
                end if
             else if(line(1:6) == 'mobile') then
                i=8
                do
                   read(line(i:i),"(A1)") buffer
                   if(buffer ==',') then
                      read(line(7:i-1),*) mobile_range(1)
                      read(line(i+1:len_trim(line)),*) mobile_range(2)
                      exit
                   else
                      i = i + 1
                   end if
                end do
             else if(line(1:9) == 'threshold') then
                read(line(10:len_trim(line)),*,iostat=info) threshold
                if(info/=0) then
                   write(*,*)'Bad value for threshold'; threshold = 0._dp; info=0
                end if
             else if(line(1:3) == 'end') then
                exit
             else
                write(*,"('Unrecognized input at line ',i3, ' for Monte Carlo Parameters.')")i_line
                exit
             end if
          end do

       case('begin temperature')                            !Section for Temperature Control Parameters
          do
             info = 0
             read(unit,'(A120)',iostat=info) line
             i_line=i_line+1
             if(info/=0) then
                write(*,*) 'Bad input at line ', i_line, ' for temperature specifications.'
                exit
             end if
             line = trim(adjustl(line))
             call lower_case(line)
             if(line(1:1) == '!') then
                cycle
             else if(line(1:9) == 'annealing') then
                annealing = .true.
             else if(line(1:8) == 'constant') then
                cool_opt = 'constant'
             else if(line(1:6) == 'linear') then
                cool_opt = 'linear'
                read(line(7:len_trim(line)),*,iostat=info) T_param
                if(info/=0)then
                   write(*,*) 'Bad value for linear cooling schedule parameter'; T_param=1._dp; info=0
                end if
             else if(line(1:3) == 'exp') then
                cool_opt = 'exp'
                read(line(4:len_trim(line)),*,iostat=info) T_param
                if(info/=0)then
                   write(*,*) 'Bad value for exponential cooling schedule parameter'; T_param=1._dp; info=0
                end if
                if(T_param <= 0. .or. T_param >= 1.)then
                   write(*,*) 'Poor value for exponential cooling parameter; set to 0.1'; T_param =0.1_dp
                end if
             else if(line(1:7) == 'tlinear') then
                cool_opt = 'tlinear'
                read(line(8:len_trim(line)),*,iostat=info) T_param
                if(info/=0)then
                   write(*,*) 'Bad value for constant cooling schedule parameter'; T_param=1._dp; info=0
                end if
             else if(line(1:4) == 'tlog') then
                cool_opt = 'tlog'
                read(line(5:len_trim(line)),*,iostat=info) T_param
                if(info/=0)then
                   write(*,*) 'Bad value for logarithmic cooling schedule parameter'; T_param=1._dp; info=0
                end if
             else if(line(1:4) == 'texp') then
                cool_opt = 'texp'
                read(line(5:len_trim(line)),*,iostat=info) T_param
                if(info/=0)then
                   write(*,*) 'Bad value for exponential cooling schedule parameter'; T_param=1._dp; info=0
                end if
                if(T_param > 0.)then
                   write(*,*) 'Positive value for exponential cooling parameter, set negative'; T_param =T_param*(-1.)
                end if
             else if(line(1:2) == 't0') then
                read(line(3:len_trim(line)),*,iostat=info) T
                if(info/=0)then
                   write(*,*) 'Bad value for initial temperature, T0'; T=0._dp; info=0
                end if
                max_T = T
             else if(line(1:4) == 'maxt')then
                read(line(5:len_trim(line)),*,iostat=info)max_T
                if(info/=0)then
                   write(*,*) 'Bad value for maxT.'; max_T=0.0_dp; info=0
                end if
             else if(line(1:4) == 'mint')then
                read(line(5:len_trim(line)),*,iostat=info)min_T
                if(info/=0)then
                   write(*,*) 'Bad value for minT.'; min_T=0.0_dp; info=0
                end if
             else if(line(1:3) == 'end') then
                exit
             else
                write(*,"('Unrecognized input at line ',i3, ' for Monte Carlo Temperature Parameters.')")i_line
                exit
             end if
          end do

       case('begin output')                                 !Section for output demands
          i_file = 1
          do
             info = 0
             read(unit,'(A120)',iostat=info) line
             i_line=i_line+1
             if(info/=0) then
                write(*,"('Bad input at line ', i3, ' for Output specifications.')")i_line
                info=0
                exit
             end if
             line = trim(adjustl(line))
             if(line(1:1) == '!')then
                cycle
             else if(line(1:3) == 'pdf') then
                xye_file = i_file
                read(line(4:80),'(A77)',iostat = info) output(i_file)
                if(info/=0) then
                   write(*,"('Bad filename for xye at line ', i3)")i_line; info=0
                end if
                i_file = i_file+1
             else if(line(1:3) == 'pdb') then
                pdb_file = i_file
                read(line(4:80),'(A77)',iostat = info) output(i_file)
                if(info/=0) then
                   write(*,"('Bad filename for pdb at line ', i3)")i_line; info=0
                end if
                i_file = i_file+1
             else if(line(1:3) == 'pot') then
                pot_file = i_file
                read(line(4:80),'(A77)',iostat = info) output(i_file)
                if(info/=0)then
                   write(*,"('Bad filename for pot at line ', i3)")i_line; info=0
                end if
                i_file = i_file+1
             else if(line(1:3) == 'ram') then
                ram_file = i_file
                read(line(4:80),'(A77)',iostat=info)output(i_file)
                if(info/=0)then
                   write(*,"('Bad filename for ram at line ', i3)")i_line; info=0
                end if
                i_file = i_file+1
             else if(line(1:4) == 'time') then
                time_file = i_file
                read(line(5:80),'(A76)',iostat=info) output(i_file)
                if(info/=0)then
                   write(*,"('Bad filename for time at line ', i3)")i_line; info=0
                end if
                i_file = i_file+1
             else if(line(1:4) == 'dist') then
                dist_file = i_file
                read(line(5:80),'(A76)',iostat=info) output(i_file)
                if(info/=0)then
                   write(*,"('Bad filename for dist at line ', i3)")i_line; info=0
                end if
                i_file = i_file+1
             else if(line(1:6) == 'thermo') then
                therm_file = i_file
                read(line(7:80),'(A76)',iostat=info) output(i_file)
                if(info/=0)then
                   write(*,"('Bad filename for thermo at line ', i3)")i_line; info=0
                end if
                i_file = i_file+1
             else if(line(1:3) == 'end') then
                exit
             else
                write(*,"('Unrecognized input at line ',i3,' for Outputs.')")i_line
             end if
          end do
          allocate(output_fname(i_file-1))
          do i=1,size(output_fname)
             output_fname(i) = adjustl(output(i))
          end do

       case('')
         cycle

       case default
          write(*, "('Unclassifiable statement at line ',i3,'.')")i_line
       end select
    end do


    ! System of Checks
    if(max_sweeps==0) then
       write(*,*) 'Bad MC parameters, program terminated.'
       write(*,"('Maximum MC Sweeps = ',i2.2)")max_sweeps
       stop
    end if

    if (.not. allocated(output_fname)) allocate(output_fname(0))

    if(max_T <=0.0 .or. min_T <0.0 .or. T<0.0) then
       write(*,*) 'Bad MC temperature parameters, program terminated.'
       write(*,"('Maximum temperature = ',f8.3)") max_T
       write(*,"('Minimum temperature = ',f8.3)") min_T
       write(*,"('Starting temperature = ',f8.3)") T
       stop
    end if

    if(annealing)then
       if(T>max_T .or. T<min_T .or. T <=0.)then
          T = max_T
          write(*,"('Poor value for initial annealing Temperature; reset to maximum. T = ',f8.3)") T
       end if
    end if

    ! Finalizing calculated variables
    if(min_T==0.) min_T=0.00005_dp
    if(T==0.)T=0.00005_dp
    T0 = T

    close(unit)
  end subroutine get_input
  ! ----------------------------------------------------------------------------

  subroutine write_time(filename,total_count,chi,p)
    ! SUBROUTINE: write_time(filename,total_count,[chisquared])
    character(len=*), intent(in) :: filename                !Filename
    integer, intent(in) :: total_count                      !MC Counter for tracking
    real(dp), intent(in), optional :: chi                   !Chi_squared value from most recent run
    real(dp), intent(in),optional :: p                      !Prior probability from most recent run
    logical, save :: FirstCall = .true.                     !Indicator for first write call
    integer :: info                                         !Error flag (0 for success)
    character(len=80) :: error                              !Error message
    integer :: unit                                         !File handle

    if(FirstCall) then
       open(newunit=unit,file=trim(filename),status='replace',IOSTAT=info,IOMSG=error)
       if(info/=0) then
          write(*,*)'Error in writing ',trim(filename)
          write(*,*) error
       end if
       FirstCall = .false.
       write(unit,"('       Sweeps,   Temperature,   Chi^2, log10(liklihood)')")
       write(unit,"('-----------------------------------------------------------------------------')")
       if (present(chi) .and. present(p)) then
          write(unit,"(I10,',',ES15.3,',',ES12.3,',',F8.3)") total_count, T, chi,p
       else if (present(chi)) then
          write(unit,"(I10,',',ES15.3,',',ES12.3)") total_count, T, chi
       else if(present(p)) then
          write(unit,"(I10,',',ES15.3,',',F8.3)") total_count, T, p
       else
          write(unit,"(I10,',',ES15.3,',',F8.3)") total_count, T, 0.0
       end if
    else
       open(newunit=unit,file=trim(filename),status ='old',position='append',IOSTAT=info,IOMSG=error)
       if(info/=0) then
          write(*,*)'Error in writing ',trim(filename)
          write(*,*) error
       end if
       if (present(chi) .and. present(p)) then
          write(unit,"(I10,',',ES15.3,',',ES12.3,',',F8.3)") total_count, T, chi,p
       else if (present(chi)) then
          write(unit,"(I10,',',ES15.3,',',ES12.3)") total_count, T, chi
       else if(present(p)) then
          write(unit,"(I10,',',ES15.3,',',F8.3)") total_count, T, p
       else
          write(unit,"(I10,',',ES15.3,',',F8.3)") total_count, T, 0.0
       end if
    end if
    close(unit)
  end subroutine write_time
  !-----------------------------------------------------------------------------

end module cmd_io
