! LAST EDIT: Phil Maffettone 2017-03-20
! Used for making one off calculations of pair distribution functions
! TODO: Create finalizing functions to remove objects from mem
! TODO: Create options for partial PDFS in public routine
module pair_dist_funcs

  use types
  use maths
  use class_proteins

  implicit none

  private
  public :: write_pdf

  type partial_rcf
     ! PARTIAL RADIAL COUNTING FUNCTIONS, N_ab(r)
     ! List of histograms according to spec_pair, where N_ab = N_ba
     real(dp), allocatable :: dat(:,:)
     real(dp) :: r_min
     real(dp) :: r_max
     real(dp) :: dr
     integer :: n_bin
     character(len=2) ::  spec_pair(2,15)
  end type partial_rcf

  type partial_rdf
     ! PARTIAL RADIAL DISTRIBUTION FUNCITONS, R_ab(r)
     ! Square matrix of histograms, where R_ab /= R_ba
     real(dp),allocatable :: dat(:,:,:)
     real(dp) :: r_min,r_max, dr
     integer :: n_bin
     character(len=2) :: spec_pair(2,5,5)
  end type partial_rdf

  type rho_r
     ! Billinge style number density
     real(dp),allocatable :: dat(:,:)
     real(dp),allocatable :: b_ij(:,:)
     real(dp),allocatable :: partials(:,:,:)
     real(dp) :: b_avg
     real(dp) :: rho0
  end type rho_r

contains
  
  ! --------------------Public Contents-----------------------------------------
  ! SUBROUTINE: write_pdf(protein,filename,[flag,r_max,dr,r_min,bbonly,info])
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  ! FUNCTION: partial_rcf_(protein,r_max, dr,r_min)
  ! FUNCITON: partial_rdf_(partial_rcf,protein)
  ! FUNCTION rho_(partial_rdf,protein, [scattering factor option])
  ! FUNCTION: search_spec_pair(atom_symbol(2))
  ! SUBROUTINE: mask_spec_pair(mask1,mask2,atom_symbol(2),protein)
  ! ----------------------------------------------------------------------------

  subroutine write_pdf(prot,fname,flag,r_max,dr,r_min,scat,info)
    ! SUBROUTINE: write_pdf(protein,filename,[flag,r_max,dr,r_min,bbonly,info])
    ! Outer routine, to write pdf flavor according to flag and filename
    ! No explicit request for bbonly needed as contained in the protein object

    type(partial_rcf) :: rcf                                !Radial counting function
    type(partial_rdf) :: rdf                                !Radial distribution functino
    type(rho_r) :: rho                                      !Density function
    type(protein), intent(in) :: prot                       !Protein
    character(len=*),intent(in) :: fname                    !Filename or root
    character(len=*),intent(in),optional :: flag            !Flavor set
    character(len=*),intent(inout),optional :: scat         !Scattering length type
    real(dp),intent(in),optional :: r_max                   !Maximum r for PDF
    real(dp),intent(in),optional :: dr                      !r spacing for PDF
    real(dp),intent(in),optional :: r_min                   !Minumum r for PDF
    integer,intent(inout),optional :: info                  !Info for errors
    character(len=80) :: flag_loc                           !Local flag
    real(dp) :: r_max_loc                                   !Local r_max
    real(dp) :: dr_loc                                      !Local dr
    real(dp) :: r_min_loc                                   !Local r_min
    real(dp), allocatable :: func(:,:)                      !Desired functinon data for printing
    integer :: info_loc                                     !Local info
    integer :: unit                                         !File unit
    integer :: i,j

    ! Setting initials and default values
    info_loc=0
    if(present(flag)) then
       flag_loc = trim(adjustl(flag))
    else
       flag_loc = 'total pdf'
    end if
    if(present(r_max)) then
       r_max_loc = r_max
    else
       r_max_loc = 50.
    end if
    if(present(dr)) then
       dr_loc = dr
    else
       dr_loc = 0.1
    end if
    if(present(r_min)) then
       r_min_loc = r_min
    else
       r_min_loc = 0.0
    end if

    ! Checks for failure
    if (r_min_loc < 0.) then
       info_loc = 1
    end if
    if (r_min_loc >= r_max_loc) then
       info_loc = 2
    end if
    if (dr_loc>r_max_loc-r_min_loc) then
       info_loc = 3
    end if
    if (present(info) .and. info_loc /= 0) then
       info = info_loc
       return
    end if


    ! Constructing function from potential objects
    select case(trim(flag_loc))
    case('rcf')
       ! Sum of partial radial counting functions, similar to full distance list
       rcf = partial_rcf_(prot,r_max_loc,dr_loc,r_min_loc)
       allocate(func(size(rcf%dat,1),2))
       do i=1,size(rcf%dat,2)
          func(:,2) = func(:,2) + rcf%dat(:,i)
       end do
       do i=1,size(func,1)
          func(i,1) = rcf%r_min + (i-0.5)*rcf%dr
       end do
    case('density')
       rcf = partial_rcf_(prot,r_max_loc,dr_loc,r_min_loc)
       rdf = partial_rdf_(rcf,prot)
       if(present(scat)) then
          rho = rho_(rdf,prot,scat)
       else
          rho = rho_(rdf,prot)
       end if
       allocate(func(size(rho%dat,1),size(rho%dat,2)))
       func = rho%dat
    case('gpdf')
       rcf = partial_rcf_(prot,r_max_loc,dr_loc,r_min_loc)
       rdf = partial_rdf_(rcf,prot)
       if(present(scat)) then
          rho = rho_(rdf,prot,scat)
       else
          rho = rho_(rdf,prot)
       end if
       allocate(func(size(rho%dat,1),size(rho%dat,2)))
       func = rho%dat
       do i=1,size(func,1)
          func(i,2) = (func(i,2) - rho%rho0)*4*pi*func(i,1)
       end do
    case('d(r)')
       rcf = partial_rcf_(prot,r_max_loc,dr_loc,r_min_loc)
       rdf = partial_rdf_(rcf,prot)
       if(present(scat)) then
          rho = rho_(rdf,prot,scat)
       else
          rho = rho_(rdf,prot)
       end if
       allocate(func(size(rho%dat,1),size(rho%dat,2)))
       func = rho%dat
       do i=1,size(func,1)
          func(i,2) = (func(i,2) - rho%rho0)*4*pi*func(i,1)
       end do
       func(:,2) = func(:,2)*rho%b_avg**2
    case('g(r)')
       rcf = partial_rcf_(prot,r_max_loc,dr_loc,r_min_loc)
       rdf = partial_rdf_(rcf,prot)
       if(present(scat)) then
          rho = rho_(rdf,prot,scat)
       else
          rho = rho_(rdf,prot)
       end if
       allocate(func(size(rho%dat,1),size(rho%dat,2)))
       func = rho%dat
       do i=1,size(func,1)
          func(i,2) = (func(i,2) - rho%rho0)*4*pi*func(i,1)
       end do
       func(:,2) = func(:,2)*rho%b_avg**2
       do i=1,size(func,1)
          func(i,2) = func(i,2)/(4*pi*rho%rho0*func(i,1))
       end do
    end select

    ! Writing output of desired function
    open(newunit=unit,file=trim(fname),status='replace',IOSTAT = info_loc)
    if (present(info) .and. info_loc /= 0) then
       info = info_loc
       return
    end if
    if (present(flag)) write(unit,*) flag
    do i=1,size(func,1)
       write(unit,*) func(i,1),',',func(i,2)
    end do
    close(unit)
  end subroutine write_pdf
  ! -----------------------------------------------------------------------------
  
  function partial_rcf_(prot,r_max,dr,r_min) result(r)
    ! FUNCTION: partial_rcf_(protein,r_max, dr,[r_min])
    ! Constructor for partial radial counting funtion
    type(protein), intent(in) :: prot                       !Protein
    real(dp), intent(in) :: r_max                           !Maximum r for PDF
    real(dp),intent(in) :: r_min                            !Minimum r for PDF
    real(dp), intent(in) :: dr                              !r spacing for PDF
    type(partial_rcf) :: r                                  !Returned RCF
    real(dp) :: dist                                        !Distance between atoms
    real(dp) :: dat_range                                   !Range of experimental data
    character(len=2) :: sym(2)                              !Atomic pair
    logical,allocatable  :: mask1(:),mask2(:)               !Local masks for looping atoms
    integer :: n_bin                                        !Number of bins in histogram
    integer :: n_pair                                       !Number of pairs
    integer :: bin                                          !Bin index
    integer :: i,j,k,l,m

    i=1;k=1;
    do
       do j=i,5
          r%spec_pair(:,k) = [prot%species(i),prot%species(j)]
          k = k+1
       end do
       i = i+1
       if (i > 5) exit
    end do
    n_pair = k-1                                            !Should be 15
    if (n_pair/=15) then
       write(*,*) 'n_pair is incorrect in Partial Dist Func'
       stop
    end if

    ! Building histogram with bin centered at input r_min and incrementing out to r_max
    dat_range = r_max-r_min
    r%r_min = r_min - 0.5_dp*dr
    n_bin = floor(dat_range/dr) + 1
    r%r_max = r_min + (n_bin-0.5)*dr
    r%n_bin = n_bin
    r%dr = dr
    allocate(r%dat(n_bin,n_pair))
    r%dat =0._dp

    allocate(mask1(size(prot%c_mask)),mask2(size(prot%c_mask)))
    ! Calculating and bining all internal distances for atoms in protein
    ! Loop maintains order of atom pairing
    do i=1,size(r%spec_pair,2)
       sym = r%spec_pair(:,i)
       call mask_spec_pair(mask1,mask2,sym,prot)
       if (sym(1) == sym(2)) then
          ! Avoids double counting for same species
          do j=1,size(mask1)
             do k=j,size(mask2)
                if(mask1(j) .and. mask2(k)) then
                   dist = norm2(prot%coordinates(:,j)-prot%coordinates(:,k))
                   if (dist>=r%r_max .or. dist<=r%r_min .or. dist ==0.) cycle
                   bin = ceiling((dist - r%r_min)/dr)
                   r%dat(bin,i) = r%dat(bin,i) + 1.0_dp
                end if
             end do
          end do
       else
          do j=1,size(mask1)
             do k=1,size(mask2)
                if(mask1(j) .and. mask2(k)) then
                   dist = norm2(prot%coordinates(:,j)-prot%coordinates(:,k))
                   if (dist>=r%r_max .or. dist<=r%r_min) cycle
                   bin = ceiling((dist - r%r_min)/dr)
                   r%dat(bin,i) = r%dat(bin,i) + 1.0_dp
                end if
             end do
          end do
       end if
    end do
  end function partial_rcf_
  ! -----------------------------------------------------------------------------

  function partial_rdf_(rcf,prot) result(r)
    ! FUNCITON: partial_rdf_(partial_rcf,protein)
    type(partial_rcf), intent(in) :: rcf                    !Partial RCF
    type(protein), intent(in) :: prot                       !Protein
    type(partial_rdf) :: r                                  !Returned partial RDF
    character(len=2) :: sym(2)                              !Atomic symbol pair
    integer :: i,j,m                                        !Dummy integers

    r%r_min = rcf%r_min
    r%r_max = rcf%r_max
    r%n_bin = rcf%n_bin
    r%dr = rcf%dr
    allocate(r%dat(r%n_bin,size(prot%species),size(prot%species)))

    r%dat = 0._dp
    do i=1, size(prot%species)
       do j=1, size(prot%species)
          sym(2) = prot%species(i)
          sym(1) = prot%species(j)
          r%spec_pair(:,j,i) = sym
          m  =search_spec_pair(sym)
          if(prot%spec_count(j)>0)then
             r%dat(:,j,i) = rcf%dat(:,m)/prot%spec_count(j)
          end if
       end do
    end do
  end function partial_rdf_
  ! -----------------------------------------------------------------------------

  function rho_(rdf,prot,option) result(r)
    ! FUNCTION rho_(partial_rdf,protein, [scattering factor option])
    ! Takes partial RDF (broken down by atom pairs) and creates Billinge
    ! rho(r) number density
    ! abscissa is dat(:,1) and ordinate is dat(:,2)
    ! OPTION for NEUTRON, XRD, NEUTRON+D
    type(partial_rdf), intent(in) :: rdf                    !Input radial dist function by atom pair
    type(protein), intent(in) :: prot                       !Input protein
    type(rho_r) :: r                                        !Return rho(r)
    character(len=*),intent(inout),optional :: option       !Option for experiment type
    real(dp), allocatable :: b(:,:)                         !Scattering factors
    real(dp), allocatable :: c(:)                           !Number concentration/fraction
    real(dp) :: b_avg                                       !Average scattering factor
    real(dp) :: b_i(2)                                      !Pertinant b
    real(dp) :: dr                                          !Range of bin
    integer :: N                                            !Total atom count
    integer :: i,j,k                                        !Dummy integers

    allocate(r%dat(rdf%n_bin,2),&
         b(size(rdf%dat,2),size(rdf%dat,3)),c(size(prot%spec_count)), &
         r%partials(rdf%n_bin,size(rdf%dat,2),size(rdf%dat,3)))

    ! Initializing for error
    b_i = 1.0_dp

    ! Creating combined b value of b1*b2
    if (present(option)) then
       call upper_case(option)
       select case(trim(option))
       case('XRD')
          do i=1,size(b,2)
             do j=1,size(b,1)
                do k=1,2
                   select case(rdf%spec_pair(k,j,i))
                   case('O ')
                      b_i(k) = 8.0_dp
                   case('C ')
                      b_i(k) = 6.0_dp
                   case('N ')
                      b_i(k) = 7.0_dp
                   case('H ')
                      b_i(k) = 1.0_dp
                   case('S ')
                      b_i(k) = 16.0_dp
                   end select
                end do
                b(j,i) = b_i(1)*b_i(2)
             end do
          end do
       case('NEUTRON')
          do i=1,size(b,1)
             do j=1,size(b,2)
                do k=1,2
                   select case(rdf%spec_pair(k,j,i))
                   case('O ')
                      b_i(k) = 5.803_dp
                   case('C ')
                      b_i(k) = 6.6460_dp
                   case('N ')
                      b_i(k) = 9.36_dp
                   case('H ')
                      b_i(k) = -3.7390_dp
                   case('S ')
                      b_i(k) = 2.847_dp
                   end select
                end do
                b(j,i) = b_i(1)*b_i(2)
             end do
          end do
       case('NEUTRON+D')
          do i=1,size(b,1)
             do j=1,size(b,2)
                do k=1,2
                   select case(rdf%spec_pair(k,j,i))
                   case('O ')
                      b_i(k) = 5.803_dp
                   case('C ')
                      b_i(k) = 6.6460_dp
                   case('N ')
                      b_i(k) = 9.36_dp
                   case('H ')
                      b_i(k) = 6.674_dp
                   case('S ')
                      b_i(k) = 2.847_dp
                   end select
                end do
                b(j,i) = b_i(1)*b_i(2)
             end do
          end do
       case default
          b=1.0_dp
       end select
    else
       b=1.0_dp
    end if

    ! Calculating number concentration and average b
    N=0
    do i=1,size(c)
       N = N+prot%spec_count(i)
    end do
    c = prot%spec_count/(N*1.0_dp)
    if (present(option)) then
       select case(trim(option))
       case('XRD')
          b_avg = c(1)*6.0_dp+c(2)*8.0_dp+c(3)*6.0_dp+c(4)+c(5)*16.0_dp
       case('NEUTRON')
          b_avg = c(1)*6.6460_dp+c(2)*5.803_dp+c(3)*9.36_dp-c(4)*3.7390_dp+c(5)*2.847_dp
       case('NEUTRON+D')
          b_avg = c(1)*6.6460_dp+c(2)*5.803_dp+c(3)*9.36_dp+c(4)*6.674_dp+c(5)*2.847_dp
       case default
          b_avg = 1.0_dp
       end select
    else
       b_avg = 1.0_dp
    end if

    ! Calculating the total number density function
    ! Sets r to the center of bin and rho(r) to correct sum
    r%dat = 0._dp
    do i=1,size(r%dat,1)
       r%dat(i,1) = rdf%r_min + (i-.5)*rdf%dr
    end do
    do i=1,size(b,2)
       do j=1,size(b,1)
          r%partials(:,j,i) = (b(j,i)/b_avg**2)*rdf%dat(:,j,i)
          r%dat(:,2) = r%dat(:,2) + c(j)*r%partials(:,j,i)
       end do
    end do

    r%b_ij = b
    r%b_avg = b_avg

    ! Finding average value of density
    b_avg =0._dp
    do i=1,size(r%dat,1)
       b_avg = b_avg + rdf%dr*r%dat(i,2)
    end do
    b_avg = b_avg/(r%dat(size(r%dat,1),1)-r%dat(1,1))

    r%rho0 = b_avg
  end function rho_
  
  
  function search_spec_pair(sym) result(r)
    ! FUNCTION: search_spec_pair(atom_symbol(2))
    character(len=2), intent(in) :: sym(2)                  !Atomic symbol pair to check
    integer ::  r
    select case(sym(1))
    case('C ')
       select case(sym(2))
       case('C ')
          r=1
       case('O ')
          r=2
       case('N ')
          r=3
       case('H ')
          r=4
       case('S ')
          r=5
       case default
          write(*,*)'Species search failure in Partial_Dist_Func.f90'
       end select
    case('O ')
       select case(sym(2))
       case('C ')
          r=2
       case('O ')
          r=6
       case('N ')
          r=7
       case('H ')
          r=8
       case('S ')
          r=9
       case default
          write(*,*)'Species search failure in Partial_Dist_Func.f90'
       end select
    case('N ')
       select case(sym(2))
       case('C ')
          r=3
       case('O ')
          r=7
       case('N ')
          r=10
       case('H ')
          r=11
       case('S ')
          r=12
       case default
          write(*,*)'Species search failure in Partial_Dist_Func.f90'
       end select
    case('H ')
       select case(sym(2))
       case('C ')
          r=4
       case('O ')
          r=8
       case('N ')
          r=11
       case('H ')
          r=13
       case('S ')
          r=14
       case default
          write(*,*)'Species search failure in Partial_Dist_Func.f90'
       end select
    case('S ')
       select case(sym(2))
       case('C ')
          r=5
       case('O ')
          r=9
       case('N ')
          r=12
       case('H ')
          r=14
       case('S ')
          r=15
       case default
          write(*,*)'Species search failure in Partial_Dist_Func.f90'
       end select
    case default
       write(*,*)'Species search failure in Partial_Dist_Func.f90'
    end select
  end function search_spec_pair
  ! --------------------------------------------------------------------------------------

  subroutine mask_spec_pair(mask1,mask2,sym,prot)
    ! SUBROUTINE: mask_spec_pair(mask1,mask2,atom_symbol(2),protein)
    logical,intent(inout) :: mask1(:)                         !Mask for first atom
    logical,intent(inout) :: mask2(:)                         !Mask for second atom
    character(len=2),intent(in) :: sym(2)                     !Atomic symbol pair to check
    type(protein),intent(in) :: prot                          !Protein

    if (sym(1) == 'C ') then
       mask1 = prot%c_mask
    else if (sym(1) == 'O ') then
       mask1 = prot%o_mask
    else if (sym(1) == 'N ') then
       mask1 = prot%n_mask
    else if (sym(1) == 'H ') then
       mask1 = prot%h_mask
    else if (sym(1) == 'S ') then
       mask1 = prot%s_mask
    else
       write(*,*) 'Species search failure in Partial_Dist_Func.f90 mask_spec_pair'
    end if
    if (sym(2) == 'C ') then
       mask2 = prot%c_mask
    else if (sym(2) == 'O ') then
       mask2 = prot%o_mask
    else if (sym(2) == 'N ') then
       mask2 = prot%n_mask
    else if (sym(2) == 'H ') then
       mask2 = prot%h_mask
    else if (sym(2) == 'S ') then
       mask2 = prot%s_mask
    else
       write(*,*) 'Species search failure in Partial_Dist_Func.f90 mask_spec_pair'
    end if
  end subroutine mask_spec_pair
  ! -----------------------------------------------------------------------------
end module pair_dist_funcs
