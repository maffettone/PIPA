! LAST EDIT: Phil Maffettone 2017-04-06
! TODO: CHECK for superflluous variables from version 0.0.
! TODO: Fix start and finish pdf print in main and cmd_io
! CONTAINS TEST CODE; REMOVE FOR PROPER USE.
! Presently collects data for calculating acceptance percentages between snpashots and across total
! without actually calculating and outputing these values.
program main

  use types
  use maths
  use class_proteins
  use pdb_io
  use ramachandran
  use cmd_io
  use exp_io
  use angle_manip
  use dist_list
  use random
  use likelihood
  use CSprby
  use temperature
  use lists
  use pair_dist_funcs
  use scd

  implicit none
  type(protein) :: prot                                     !Protein
  character(len=80) :: filename                             !Generic filename
  character(len=80) :: input_filename                       !Input filename
  character(len=80) :: out_dir                              !Ouput director
  character(len=80) :: error                                !Error message
  real(dp), allocatable :: trip_prby(:,:,:)                 !Ramachandran triplet proby for all angle combos
  real(dp), allocatable :: cum_prby(:,:,:)                  !Ramachandran cumulative proby from trip
  real(dp), allocatable :: sim_mat(:,:)                     !Simulated distance matrix
  real(dp), allocatable :: sim_mat_new(:,:)                 !Simulated dist_mat after change
  real(dp), allocatable :: dat(:)                           !Experimental distance list
  real(dp), allocatable :: sim(:)                           !Simulated distance list
  real(dp), allocatable :: sim_new(:)                       !Simulated dist_list after change
  real(dp) :: chi                                           !Chi squared value
  real(dp) :: chi_new                                       !Chi squared value after change
  real(dp) :: ran                                           !Random number [0,1]
  real(dp) :: phi_change                                    !Phi change retained for rejection
  real(dp) :: psi_change                                    !Psi change retained for rejection
  real(dp) :: w                                             !Transition probability
  real(dp) :: new_angles(2)                                 !New angles from the sampling (degrees)
  real(dp) :: timing(2)                                     !Start and finish time
  real(dp) :: prior_prby(2)                                 !Storing liklihoods from baysian info
  integer :: total_sweep_count                              !Total Monte Carlo Counter
  integer :: autocorr_sweeps                                !Max num of sweeps to calculate autocorrelation
  integer :: snap_counter                                   !Snapshot counter for file name
  integer :: mc_count                                       !Counter for mc sweeps after equilibration
  integer :: MM_per_sweep                                   !Metropolis moves per sweep
  integer :: tau                                            !Autocorrelation time
  integer :: pro_length                                     !Protein length not including heteroatoms
  integer :: i_move                                         !Residue index to shift
  integer :: info                                           !Error flag (0 for success)
  integer :: i,j,k                                          !Dummy integers
  logical :: accepted                                       !Whether the move was accepted
  ! Variables for Autocorrelation
  type(dp_list) :: eSave                                    !Linked list saving 'energy' values
  integer :: nSave                                          !Number of saved values for autocorrelation
  integer :: nCorr                                          !Number of values accumulated in autocorr sums
  real(dp),allocatable :: cee(:)                            !Energy correlation sums (size of nSave+1)
  real(dp) :: eItem                                         !Energy item off list
  real(dp) :: eAv                                           !Expectation value of the energy

  ! ---------------------Global Data defined in input file------------------------
  ! ------------------------------------------------------------------------------
  ! character(len=80), save :: run_title                      !User supplied run title
  ! character(len=120), save :: structure_fname               !Structural data filename
  ! character(len=120), save :: exp_fname                     !Experimental data filename
  ! character(len=80), save :: init_config                    !Initial adjustment to configuration
  ! character(len=10), save :: exp_opt                        !Experimental Option (XRD or NEUTRON)
  ! character(len=10), save :: cs_opt                         !Chemical shift option (DANGLE,TALOS,NONE)
  ! character(len=10), save :: cool_opt                       !Cooling schedule speicification(linear,log,exp)
  ! character(len=80), save :: cs_dir                         !Chemical shift files directory
  ! character(len=80), allocatable, save :: output_fname(:)   !Set of output filenames
  ! character(len=80), save :: dist_flag                      !Flag for inclusion in dist list
  ! character(len=80), save :: pdf_flag                       !Flag for PDF flavor
  ! character(len=80), save :: pdf_scat                       !Flag for scattering length
  ! logical, save :: incl_omega                               !Whether to bend omega
  ! logical, save :: bbonly                                   !Backbone data to be used exclusively
  ! logical, save :: annealing                                !Whether to anneal
  ! logical, save :: cs_only                                  !Chemical shifts only, don't include ramachandran
  ! logical, save :: dist_gen                                 !Quick generation of distance list then quit
  ! logical, save :: dl                                       !True for distance list data, false for PDF data
  ! integer, save :: xye_file                                 !Index in output_fname, 0 for none
  ! integer, save :: pdb_file                                 !Index in output_fname, 0 for none
  ! integer, save :: pot_file                                 !Index in output_fname, 0 for none
  ! integer, save :: ram_file                                 !Index in output_fname, 0 for none
  ! integer, save :: time_file                                !Index in output_fname, 0 for none
  ! integer, save :: dist_file                                !Index in output_fname, 0 for none
  ! integer, save :: therm_file                               !Index in output_fname, 0 for none
  ! integer, save :: max_sweeps                               !Maximum number of MC sweeps before finish
  ! integer, save :: n_snap                                   !No of printed snapshots per run in maxsweeps limit (distlist/pdb)
  ! integer, save :: snap_met_moves                           !Number of metropolis moves between snapshots
  ! integer, save :: tot_accept                               !Total accepted Metropolis moves
  ! integer, save :: snap_accept                              !Metropolis Moves accepted in snapshot interval
  ! integer, save :: mobile_range(2)                          !Range of residues considered mobile during MC
  ! integer, save :: tau_calc_sweeps                          !Number of sweeps used to calculate autocorrelation
  ! integer, save :: autocorr_save                            !Number of covariences saved when calculating tau
  ! real(dp), save :: threshold                               !Chi threshold for satisfaction
  ! real(dp), save :: T                                       !Monte Carlo Temperature
  ! real(dp), save :: T0                                      !Initial simulation temperature
  ! real(dp), save :: max_T                                   !Maximum annealing temperature
  ! real(dp), save :: min_T                                   !Minimum annealing temperature
  ! real(dp), save :: T_param                                 !Temperature parameter for cooling schedule
  ! real(dp), save :: NMR_weight                              !Weighting on NMR-based sampling dist (<0 for all NMR)
  ! real(dp), save :: pdf_rmin                                !PDF calculation minimum radius
  ! real(dp), save :: pdf_rmax                                !PDF calculation maximum radius
  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------
  ! TEMP VAR
  ! -----------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------


  ! Initialization
  call init_random_seed()
  call cpu_time(timing(1))
  info = 0
  chi = 0._dp
  chi_new = 0._dp
  i = command_argument_count()
  if (i<1) then
     write(*,*)'Please start the program by writing ./RMCParty [input file name]'
     write(*,*)''
     stop
  else
     call get_command_argument(1,input_filename)
  end if
  call get_input(input_filename)                            !Retreieve input fname.pin

  call read_pdb(structure_fname,prot,bbonly,info,error)          !Reads initial structure
  if(info/=0) then
     write(*,*) error
     write(*,*) structure_fname
     stop
  end if

  ! Option for distance list and quit
  if (dist_gen) then
     ! Make directories and output distance list and PDF
     out_dir = './Output/'//trim(adjustl(run_title))
     call system('mkdir -p '//trim(out_dir))
     out_dir = './Output/'//trim(adjustl(run_title))//'/Data'
     call system('mkdir -p '//trim(out_dir))
     sim_mat = dist_mat(prot,dist_flag)
     allocate(sim(size(sim_mat,1)*(size(sim_mat,1)-1)/2))
     call sort_dl(sim_mat,sim)
     out_dir = './Output/'//trim(adjustl(run_title))
     write(filename,"('distancegen.dat')")
     filename = trim(out_dir)//'/Data/'//trim(filename)
     call write_dist_list(filename,sim)
     write(filename,"('_genpdf.csv')")
     filename = trim(out_dir)//'/Data/'//trim(output_fname(xye_file))//trim(filename)
     if (pdf_flag=='none') then
        call write_pdf(prot,filename)
     else
        call write_pdf(prot,filename,pdf_flag,pdf_rmax,pdf_dr,pdf_rmin,pdf_scat)
     end if
     write(*,'(A)') "Distance list and PDF generated and program terminated."
     write(*,'(A)') "Located at directory",trim(out_dir)//'/Data/'
     call exit()
  end if

  ! Sets mobile range to whole length if not or poorly set prior
  pro_length = prot%n_res
  do
     if(prot%pep(pro_length)%heteroatom) then
        pro_length = pro_length - 1
     else
        exit
     end if
  end do
  if (mobile_range(1)<= 0) mobile_range(1) = 1
  if (mobile_range(2)<=0 .or. mobile_range(2)>pro_length) mobile_range(2) = pro_length
  call Triplet_Prby(prot%pep,trip_prby)               !Finds ramachandran triplets
  if (cs_only) then
     call Uniform_Triplet_Prby(prot%pep, trip_prby)
  end if
  call incorpCS(trip_prby,prot%n_res)                 !Adds in chemical shift to proby
  prior_prby(1) = rama_likely(prot%pep,trip_prby)
  call initial_config(prot,trip_prby)                 !Potentially randomizes distribution
  cum_prby = cumulative_prby(trip_prby)               !Creates cumulative distribution
  if (dl) then
     call read_dist_list(exp_fname,dat,info,error)    !Reads in experimental data
  else
     call read_exp(prot,dist_flag,exp_fname,dat,info,error)
  end if
  if(info/=0) then
     write(*,*) info,error
     stop
  end if


  !Makes directory for outputs
  if(size(output_fname)>0) then
     out_dir = './Output/'//trim(adjustl(run_title))
     call system('mkdir -p '//trim(out_dir))
     if (xye_file /= 0) then
        out_dir = './Output/'//trim(adjustl(run_title))//'/Data'
        call system('mkdir -p '//trim(out_dir))
     end if
     if (dist_file /= 0) then
        out_dir = './Output/'//trim(adjustl(run_title))//'/Distances'
        call system('mkdir -p '//trim(out_dir))
     end if
     if (pdb_file/=0)then
        out_dir = './Output/'//trim(adjustl(run_title))//'/PDB'
        call system('mkdir -p '//trim(out_dir))
     end if
     if (ram_file/=0)then
        out_dir = './Output/'//trim(adjustl(run_title))//'/Ram_Plots'
        call system('mkdir -p '//trim(out_dir))
        filename = trim(out_dir)//'/'//trim(output_fname(ram_file))
        call write_triplet(trip_prby,filename)
     end if
     out_dir = './Output/'//trim(adjustl(run_title))
  end if

  ! Building simulated correlation
  sim_mat = dist_mat(prot,dist_flag)
  allocate(sim(size(sim_mat,1)*(size(sim_mat,1)-1)/2))
  allocate(sim_new(size(sim_mat,1)*(size(sim_mat,1)-1)/2))
  allocate(sim_mat_new(size(sim_mat,1),size(sim_mat,2)))
  call sort_dl(sim_mat,sim)
  chi = chisq(dat,sim)

  ! Initializing parameters
  tot_accept = 0
  snap_accept = 0
  snap_counter = 0
  accepted = .false.
  total_sweep_count = 0
  snap_met_moves = 0
  MM_per_sweep = mobile_range(2) - mobile_range(1)+1
  nSave = autocorr_save
  allocate(cee(nSave+1))
  if (tau_calc_sweeps>200000) then
     autocorr_sweeps = tau_calc_sweeps
  else
     autocorr_sweeps = 200000
  end if

  ! Initializing output files
  if (.not. dl) then
     filename = trim(out_dir)//'/PDF2Distances.dat'
     call write_dist_list(filename,dat)
  end if
  if(dist_file/=0) then
     write(filename,"('_',i3.3,'.dat')")0
     filename = trim(out_dir)//'/Distances/'//trim(output_fname(dist_file))//trim(filename)
     call write_dist_list(filename,sim)
  end if
  if(pdb_file/=0)then
     write(filename,"('_',i3.3,'.pdb')")0
     filename = trim(out_dir)//'/PDB/'//trim(output_fname(pdb_file))//trim(filename)
     call write_pdb(filename,prot,bbonly)
  end if
  if(time_file/=0) then
     write(filename,"('.csv')")
     filename = trim(out_dir)//'/'//trim(output_fname(time_file))//trim(filename)
     prior_prby(2) = rama_likely(prot%pep,trip_prby)
     call write_time(filename,total_sweep_count,chi=chi,p=prior_prby(2))
  end if
  if(xye_file/=0) then
     write(filename,"('_startpdf.csv')")
     filename = trim(out_dir)//'/Data/'//trim(output_fname(xye_file))//trim(filename)
     if (pdf_flag=='none') then
        call write_pdf(prot,filename)
     else
        call write_pdf(prot,filename,pdf_flag,pdf_rmax,pdf_dr,pdf_rmin,pdf_scat)
     end if
  end if



  T_loop: do
     eAv = 0.
     cee = 0.
     nCorr = 0
     Autocorr: do k=1,autocorr_sweeps/tau_calc_sweeps
        tau_calc: do i=1,tau_calc_sweeps
           total_sweep_count = total_sweep_count + 1
           call mcsweep(MM_per_sweep)
           if(eSave%len() == nSave) then
              eAv = eAv + chi
              nCorr = nCorr + 1
              cee(1) = cee(1) + chi**2
              call eSave%iter_restart()
              do j=2,nSave+1
                 call eSave%iter_next(eItem)
                 cee(j) = cee(j) + eItem*chi
              end do
              call eSave%pop_back()
           end if
           call eSave%push_front(chi)

           ! Storing data w.r.t time every 100 sweeps
           if (total_sweep_count == max_sweeps) then
              exit T_loop
           else if (mod(total_sweep_count,100) == 0 .and. time_file/=0) then
              write(filename,"('.csv')")
              filename = trim(out_dir)//'/'//trim(output_fname(time_file))//trim(filename)
              prior_prby(2) = rama_likely(prot%pep,trip_prby)
              call write_time(filename,total_sweep_count,chi=chi,p=prior_prby(2))
           end if

           Snapshots: if(mod(total_sweep_count,max_sweeps/n_snap) == 0) then
              snap_counter = snap_counter + 1
              if(dist_file/=0) then
                 write(filename,"('_',i3.3,'.dat')")snap_counter
                 filename = trim(out_dir)//'/Distances/'//trim(output_fname(dist_file))//trim(filename)
                 call write_dist_list(filename,sim,info,error)
                 if (info/=0) then
                    write(*,*) error
                    stop
                 end if
              end if
              if(pdb_file/=0)then
                 write(filename,"('_',i3.3,'.pdb')")snap_counter
                 filename = trim(out_dir)//'/PDB/'//trim(output_fname(pdb_file))//trim(filename)
                 call write_pdb(filename,prot,bbonly)
              end if
              if(xye_file/=0) then
                 write(filename,"('_',i3.3,'.csv')")snap_counter
                 filename = trim(out_dir)//'/Data/'//trim(output_fname(xye_file))//trim(filename)
                 if (pdf_flag=='none') then
                    call write_pdf(prot,filename)
                 else
                    call write_pdf(prot,filename,pdf_flag,pdf_rmax,pdf_dr,pdf_rmin,pdf_scat)
                 end if
              end if
              snap_accept = 0
              snap_met_moves = 0
           end if Snapshots
        end do tau_calc
        ! Checking on Tau
        tau = 0
        call ComputeAutocorrelationTime()
        if (tau<1) cycle Autocorr
        if (20*tau < nCorr) then
           ! Calc thermo data, and skip over extra equilibration
           if (therm_file/=0) then
              call calc_thermo(nCorr*4)
              if (total_sweep_count >= max_sweeps) exit T_loop
           end if
           go to 100
        end if
     end do Autocorr

     ! MC_equilibrate if autocorrelaiton didn't meet the 20*tau threshold
     if (tau<1) then
        tau =1
        write(*,'(A17,ES15.3)') 'Temperature = ',T
        write(*,'(A100)') 'Warning: Autocorrelation time less than 1 sweep. Consider increasing time used for calculation'
     end if
     MC_equilibrate: do mc_count = 1,20*tau
        total_sweep_count = total_sweep_count + 1
        call mcsweep(MM_per_sweep)
        ! Storing data w.r.t time every 100 sweeps
        if (total_sweep_count == max_sweeps) then
           exit T_loop
        else if (mod(total_sweep_count,100) == 0 .and. time_file/=0) then
           write(filename,"('.csv')")
           filename = trim(out_dir)//'/'//trim(output_fname(time_file))//trim(filename)
           prior_prby(2) = rama_likely(prot%pep,trip_prby)
           call write_time(filename,total_sweep_count,chi=chi,p=prior_prby(2))
        end if
        Snapshot: if(mod(total_sweep_count,max_sweeps/n_snap) == 0) then
           snap_counter = snap_counter + 1
           if(dist_file/=0) then
              write(filename,"('_',i3.3,'.dat')")snap_counter
              filename = trim(out_dir)//'/Distances/'//trim(output_fname(dist_file))//trim(filename)
              call write_dist_list(filename,sim,info,error)
              if (info/=0) then
                 write(*,*) error
                 stop
              end if
           end if
           if(pdb_file/=0)then
              write(filename,"('_',i3.3,'.pdb')")snap_counter
              filename = trim(out_dir)//'/PDB/'//trim(output_fname(pdb_file))//trim(filename)
              call write_pdb(filename,prot,bbonly)
           end if
           if(xye_file/=0) then
              write(filename,"('_',i3.3,'.csv')")snap_counter
              filename = trim(out_dir)//'/Data/'//trim(output_fname(xye_file))//trim(filename)
              if (pdf_flag=='none') then
                 call write_pdf(prot,filename)
              else
                 call write_pdf(prot,filename,pdf_flag,pdf_rmax,pdf_dr,pdf_rmin,pdf_scat)
              end if
           end if
           snap_accept = 0
           snap_met_moves = 0
        end if Snapshot
     end do MC_equilibrate
     ! Storing Thermo Data
     if (therm_file/=0) then
        call calc_thermo((nCorr+mc_count-2)*4)
        if (total_sweep_count >= max_sweeps) exit T_loop
     end if

100  call eSave%finalize()
     call adjust_T(total_sweep_count)
     if (T<=min_T) exit T_loop
  end do T_loop


  ! Outputing final information
  if(dist_file/=0) then
     write(filename,"('_final')")
     filename = trim(out_dir)//'Distances/'//trim(output_fname(dist_file))//trim(filename)
     call write_dist_list(filename,sim)
  end if
  if(pdb_file/=0)then
     write(filename,"('_final.pdb')")
     filename = trim(out_dir)//'/PDB/'//trim(output_fname(pdb_file))//trim(filename)
     call write_pdb(filename,prot,bbonly)
  end if
  if(time_file/=0) then
     write(filename,"('.csv')")
     filename = trim(out_dir)//'/'//trim(output_fname(time_file))//trim(filename)
     prior_prby(2) = rama_likely(prot%pep,trip_prby)
     call write_time(filename,total_sweep_count,chi=chi,p=prior_prby(2))
  end if
  if (xye_file /= 0) then
     write(filename,"('_finalpdf.csv')")
     filename = trim(out_dir)//'/Data/'//trim(output_fname(xye_file))//trim(filename)
     if (pdf_flag=='none') then
        call write_pdf(prot,filename)
     else
        call write_pdf(prot,filename,pdf_flag,pdf_rmax,pdf_dr,pdf_rmin,pdf_scat)
     end if
  end if
  call decorate(prot);
  write(filename,"('_decorated.pdb')")
  filename = trim(out_dir)//'/PDB/'//trim(output_fname(pdb_file))//trim(filename)
  call write_pdb(filename,prot)

  call cpu_time(timing(2))
  write(*,*)'Simulation time = ', timing(2)-timing(1)
  if (allocated(trip_prby))deallocate(trip_prby)
  if (allocated(cum_prby))deallocate(cum_prby)
  if (allocated(dat))deallocate(dat)
  if(allocated(sim))deallocate(sim)
  if(allocated(sim_new))deallocate(sim_new)
  if(allocated(sim_mat))deallocate(sim_mat)
  if(allocated(sim_mat_new))deallocate(sim_mat_new)
  call exit()


contains

  subroutine calc_thermo(sweeps)
    ! Calculates thermodynamic quantities if desired
    ! Does not recaclculate integrated autocorrelation time
    integer, intent(in) :: sweeps                           !Number of sweeps for equilibration
    integer :: i                                            !Dummy int
    real(dp) :: E, E2                                       !'Thermodynamic' Energy and Energy^2
    E=0.
    E2=0.
    do i=1, sweeps
       call mcsweep(MM_per_sweep)
       E = E + chi
       E2 = E2 + chi**2
       total_sweep_count = total_sweep_count+1
    end do
    E = E/sweeps
    E2 = E2/sweeps
    filename = trim(out_dir)//'/'//trim(output_fname(therm_file))//'.xyz'
    call write_thermo(filename,T,E,E2,tau,info,error)
  end subroutine calc_thermo
  ! -----------------------------------------------------------------------------
  
  subroutine mcsweep(MM_per_sweep)
    integer, intent(in) :: MM_per_sweep                     !Metropolis moves per MC sweep
    integer :: i                                            !Dummy counter
    do i=1,MM_per_sweep
       snap_met_moves = snap_met_moves+1
       ! Copy the old for future comparison
       sim_mat_new = sim_mat
       sim_new = sim

       ! Choosing which residue to change
       call random_number(ran)
       i_move = ceiling(ran * (mobile_range(2)-mobile_range(1)+1) + mobile_range(1)-1)

       ! Sampling new angles (changing to radians in range)
       call random_number(ran)
       new_angles= sample_cum(cum_prby,ran,i_move)
       phi_change = (new_angles(1) - get_dihedral(prot%pep(i_move),'phi'))*pi/180._dp
       psi_change = (new_angles(2) - get_dihedral(prot%pep(i_move),'psi'))*pi/180._dp
       if (phi_change > pi) then
          phi_change = phi_change - 2._dp*pi
       else if (phi_change < -1*pi) then
          phi_change = phi_change + 2._dp*pi
       end if
       if(psi_change > pi) then
          psi_change = psi_change - 2._dp*pi
       else if(psi_change < -1.*pi) then
          psi_change = psi_change + 2._dp*pi
       end if

       ! Building and bending new protein
       if(i_move /=1)then
          call protein_bend(prot,i_move,phi_change,'phi',info,error)
          if(info/=0)then
             write(*,*)error
             info=0
          end if
       end if
       if(i_move /= pro_length)then
          call protein_bend(prot,i_move,psi_change,'psi',info,error)
          if(info/=0)then
             write(*,*)error
             info=0
          end if
       end if

       ! Calculating new distance list and then X^2, then accepting move based on criterion
       accepted = .false.
       call change_dist_mat(sim_mat_new,prot,i_move,dist_flag)
       call sort_dl(sim_mat_new,sim_new)
       chi_new = chisq(dat,sim_new)
       w = exp(-1._dp*(chi_new-chi)/(T))
       call random_number(ran)
       if (ran < w) accepted = .true.

       if(accepted) then
          chi = chi_new
          sim_mat = sim_mat_new
          sim = sim_new
          tot_accept = tot_accept+1
          snap_accept = snap_accept+1
       else
          if(i_move /=1)then
             phi_change = -1._dp * phi_change
             call protein_bend(prot,i_move,phi_change,'phi',info,error)
             if(info/=0)then
                write(*,*)error
                info=0
             end if
          end if
          if(i_move /= pro_length)then
             psi_change = -1._dp * psi_change
             call protein_bend(prot,i_move,psi_change,'psi',info,error)
             if(info/=0)then
                write(*,*)error
                info=0
             end if
          end if
       end if
    end do
  end subroutine mcsweep

  subroutine ComputeAutocorrelationTime()
    ! Local subroutine which overwrites variables in the main
    real(dp) :: av
    real(dp) :: c0
    real(dp) :: tau_loc
    integer :: n
    integer :: i
    n = nCorr
    av = eAv / n
    c0 = cee(1)/n - av*av
    tau_loc = 0._dp
    do i=2,nSave+1
       tau_loc =  tau_loc + (cee(i)/n - av*av)/c0
    end do
    tau = ceiling(tau_loc)
  end subroutine ComputeAutocorrelationTime
end program main
