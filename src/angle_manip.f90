! LAST EDIT: PHil Maffettone 2017-05-01
! Updated for helix starting config setting. NEEDS more versitility than w=180.
module angle_manip
  use types
  use maths
  use random
  use class_proteins
  use cmd_io, only : init_config, mobile_range

  implicit none
  private
  public :: protein_bend,initial_config, omega_180

contains
  ! --------------------Public Contents-----------------------------------------
  ! SUBROUTINE: protein_bend(protein,residue_index,change_in_angle,angle_type,[info,errmsg])
  ! SUBROUTINE: initial_config(protein,trip_prby)
  ! SUBROUTINE: omega_180(peptide_list)
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  ! SUBROUTINE: melt(protein, log_2(melt_moves))
  ! ----------------------------------------------------------------------------

  subroutine protein_bend(prot,index,d_angle,option,info,errmsg)
    ! SUBROUTINE: protein_bend(protein,residue_index,change_in_angle,angle_type,[info,errmsg])
    ! Angle in radians, propogates torsional change throughout structure
    ! Currenly only works for changes to phi,psi,omega
    type(protein), intent(inout) :: prot                    !Input protein
    integer, intent(in) :: index                            !Index of residue change
    real(dp), intent(in) :: d_angle                         !Change in angle (radians)
    character(len=*), intent(in) :: option                  !Angle type
    integer,optional,intent(out) :: info                    !0:success, +:Faiure
    character(len=80),optional,intent(out) :: errmsg        !Message descrbing failures
    character(len=10) :: angle_type                         !Angle type
    real(dp),dimension(3) :: origin,axis,coords,vector      !Defined for changes
    real(dp) :: rot(3,3)                                    !Rotation matrix
    integer :: i_ca,i_n,i_c,i_o                             !Indecies of backbone atoms
    integer :: i,j                                          !Dummy integers

    if(present(info))info = 0
    if(present(errmsg))errmsg='SUCCESS'
    angle_type = option
    call upper_case(angle_type)

    i_ca = 0
    i_n = 0
    i_c = 0
    i_o = 0

    ! Case where groups (such as NH3) are added to end of residue chain
    ! Returns without doing anything
    if(prot%pep(index)%heteroatom) then
       if(present(info))info=1
       if(present(errmsg)) write(errmsg,"('Non-standard backbone at residue index ',i4)")index
       return
    end if

    ! Finding  internal residue backbone indicies
    do i=1,size(prot%pep(index)%backbone)
       if (get_label(prot%pep(index)%backbone(i)) == 'CA  ') then
          i_ca = i
       else if(get_label(prot%pep(index)%backbone(i)) == 'N   ') then
          i_n = i
       else if(get_label(prot%pep(index)%backbone(i)) == 'C   ') then
          i_c = i
       else if(get_label(prot%pep(index)%backbone(i)) == 'O   ') then
          i_o = i
       end if
    end do

    ! Returning error if all not found
    if(i_ca==0 .or. i_n==0 .or. i_c==0 .or. i_o==0) then
       if(present(info))info=2
       if(present(errmsg)) write(errmsg,"('Not all backbone members found at residue index ',i4)")index
       return
    end if

    ! Changing angle and setting origin and rotation axis
    select case(trim(angle_type))
    case('OMEGA')
       if(index == 1) then
          if(present(info))info=1
          if(present(errmsg)) write(errmsg,*)'Cannot change ',angle_type,' at index ',index
          return
       end if
       do i=1,size(prot%pep(index-1)%backbone)
          if(get_label(prot%pep(index-1)%backbone(i)) == 'C   ') then
             origin = get_coords(prot%pep(index-1)%backbone(i))
             axis = get_coords(prot%pep(index)%backbone(i_n))-origin
          end if
       end do
    case('PHI')
       if(index == 1)then
          if(present(info))info=1
          if(present(errmsg)) write(errmsg,*)'Cannot change ',angle_type,' at index ',index
          return
       end if
       origin = get_coords(prot%pep(index)%backbone(i_n))
       axis = get_coords(prot%pep(index)%backbone(i_ca))-origin
    case('PSI')
       if(index == prot%n_res)then
          if(present(info))info=1
          if(present(errmsg)) write(errmsg,*)'Cannot change ',angle_type,' at index ',index
          return
       end if
       origin = get_coords(prot%pep(index)%backbone(i_ca))
       axis = get_coords(prot%pep(index)%backbone(i_c))-origin
    end select
    call change_dihedral(prot%pep(index),d_angle,angle_type)

    call RotMatrix(d_angle,axis,rot)

    ! Propogating changes internal to residue
    ! Working left to right with exceptions for omega, then phi, then psi
    ! Exceptions include:
    ! No change to N
    ! Only omega changes alpha C and nitrogen H's
    ! Only omega and phi change C and alpha H
    ! PSI does not change side chain
    i=index
    do j=1,size(prot%pep(i)%backbone)
       if(j == i_n) cycle
       if(j == i_ca .and. .not. trim(angle_type) == 'OMEGA') cycle
       if(j == i_c .and. trim(angle_type) == 'PSI') cycle
       if(get_label(prot%pep(i)%backbone(j)) == 'HA  ' .and. trim(angle_type)=='PSI')cycle
       if(get_sym(prot%pep(i)%backbone(j)) =='H ' .and. .not. trim(angle_type)=='OMEGA')cycle
       vector = prot%pep(i)%backbone(j)%coords - origin
       vector = matmul(rot,vector)
       coords = origin+vector
       call change_coords(prot%pep(i)%backbone(j),coords,prot)
    end do

    if(size(prot%pep(i)%sidechain)>0 .and. .not. trim(angle_type) == 'PSI' ) then
       do j=1,size(prot%pep(i)%sidechain)
          vector = prot%pep(i)%sidechain(j)%coords - origin
          vector = matmul(rot,vector)
          coords = origin+vector
          call change_coords(prot%pep(i)%sidechain(j),coords,prot)
       end do
    end if


    ! Propogating changes rightward of residue
    do i=index+1,prot%n_res
       do j=1,size(prot%pep(i)%backbone)
          vector = prot%pep(i)%backbone(j)%coords - origin
          vector = matmul(rot,vector)
          coords = origin+vector
          call change_coords(prot%pep(i)%backbone(j),coords,prot)
       end do
       if(size(prot%pep(i)%sidechain)==0) cycle
       do j=1,size(prot%pep(i)%sidechain)
          vector = prot%pep(i)%sidechain(j)%coords - origin
          vector = matmul(rot,vector)
          coords = origin+vector
          call change_coords(prot%pep(i)%sidechain(j),coords,prot)
       end do
    end do

  end subroutine protein_bend
  ! -----------------------------------------------------------------------------

  subroutine initial_config(prot,trip_prby)
    ! Generates an initial configuration dependent on user request and Rama Prby
    ! Depends on init_config string from input_cmd
    ! Present options include: 'gaussian_all','random_all','gaussian', 'random', 'N/A'
    type(protein),intent(inout) :: prot                !Protein
    real(dp),intent(in) :: trip_prby(:,:,:)            !Triplet probability from routine
    real(dp) :: r                                      !random number
    real(dp) :: old_angle, new_angle                   !Old/New angles (degrees)
    real(dp) :: d_angle                                !Change in angle(radians)
    integer :: i,j,k                                   !Dummy integers
    integer :: info                                    !Info from subroutines
    character(len=10) :: angle_type                    !Angle type(Phi/Psi/Omega/)
    character(len=80) :: error                         !Error message


    ! First considers melting, then all other options
    if (trim(init_config(1:4)) == 'melt') then
       read(init_config(5:len_trim(init_config)),*,iostat=j) i
       if (j/= 0) then
          write(*,"('Bad value for melting magnitude in initial_config')")
          return
       end if
       call melt(prot, i)
       return
    end if
    
    select case(trim(init_config))
    case('gaussian_all')
       ! Samples a guassian distribution for a change in the angle (degrees)
       ! about the most probable angle given the triplet probabilities
       ! The present changes are 0.5*10.0 = +/-5.0 degrees
       angle_type = 'OMEGA'
       do i=2,prot%n_res
          call random_number(r)
          call rand_normal_dist(r)
          r = r-0.5_dp
          r = r*10.0_dp
          old_angle = get_dihedral(prot%pep(i),angle_type)
          new_angle = 180._dp + r
          d_angle = (new_angle-old_angle)*pi/180._dp
          if(d_angle>2*pi) d_angle = d_angle-2*pi
          call protein_bend(prot,i,d_angle,angle_type,info=info,errmsg=error)
          if(info/=0) then
             write(*,*)error
             stop 'Failure in initial_config for OMEGA'
          end if
       end do

       angle_type = 'PHI';j=1
       do i =2,prot%n_res
          k = maxloc(trip_prby(3,:,i),1)
          call random_number(r)
          call rand_normal_dist(r)
          r = r-0.5_dp
          r = r*10.0_dp
          old_angle = get_dihedral(prot%pep(i),angle_type)
          new_angle = trip_prby(j,k,i) + r
          d_angle = (new_angle - old_angle)*pi/180._dp
          if(d_angle>2*pi) d_angle = d_angle-2*pi
          call protein_bend(prot,i,d_angle,angle_type,info=info,errmsg=error)
          if(info/=0) then
             write(*,*)error
             stop 'Failure in initial_config for PHI'
          end if
       end do

       angle_type = 'PSI';j=2
       do i=1,prot%n_res-1
          k = maxloc(trip_prby(3,:,i),1)
          call random_number(r)
          call rand_normal_dist(r)
          r = r-0.5_dp
          r = r*10.0_dp
          old_angle = get_dihedral(prot%pep(i),angle_type)
          new_angle = trip_prby(j,k,i) + r
          d_angle = (new_angle - old_angle)*pi/180._dp
          if(d_angle>2*pi) d_angle = d_angle-2*pi
          call protein_bend(prot,i,d_angle,angle_type,info=info,errmsg=error)
          if(info/=0) then
             write(*,*)error
             stop 'Failure in initial_config for PSI'
          end if
       end do

    case('gaussian')
       ! Same as gaussian all for only phi and psi angles
       angle_type = 'PHI';j=1
       do i =2,prot%n_res
          k = maxloc(trip_prby(3,:,i),1)
          call random_number(r)
          call rand_normal_dist(r)
          r = r-0.5_dp
          r = r*10.0_dp
          old_angle = get_dihedral(prot%pep(i),angle_type)
          new_angle = trip_prby(j,k,i) + r
          d_angle = (new_angle - old_angle)*pi/180._dp
          if(d_angle>2*pi) d_angle = d_angle-2*pi
          call protein_bend(prot,i,d_angle,angle_type,info=info,errmsg=error)
          if(info/=0) then
             write(*,*)error
             stop 'Failure in initial_config for PHI'
          end if
       end do

       angle_type = 'PSI';j=2
       do i=1,prot%n_res-1
          k = maxloc(trip_prby(3,:,i),1)
          call random_number(r)
          call rand_normal_dist(r)
          r = r-0.5_dp
          r = r*10.0_dp
          old_angle = get_dihedral(prot%pep(i),angle_type)
          new_angle = trip_prby(j,k,i) + r
          d_angle = (new_angle - old_angle)*pi/180._dp
          if(d_angle>2*pi) d_angle = d_angle-2*pi
          call protein_bend(prot,i,d_angle,angle_type,info=info,errmsg=error)
          if(info/=0) then
             write(*,*)error
             stop 'Failure in initial_config for PSI'
          end if
       end do

    case('random_all')
       ! Randomizes all of the dihedral angles from [-180.0,180.0]
       angle_type = 'OMEGA'
       do i=2,prot%n_res
          call random_number(r)
          r = r-0.5_dp
          old_angle = get_dihedral(prot%pep(i),angle_type)
          new_angle = r*360._dp
          if (new_angle > 180) new_angle = new_angle - 360._dp
          if (new_angle < -180) new_angle = new_angle + 360._dp
          d_angle = (new_angle-old_angle)*pi/180._dp
          call protein_bend(prot,i,d_angle,angle_type,info=info)
          if(info/=0) stop 'Failure intial_config for OMEGA'
       end do

       angle_type = 'PHI';
       do i =2,prot%n_res
          call random_number(r)
          r = r-0.5_dp
          old_angle = get_dihedral(prot%pep(i),angle_type)
          new_angle = r*360._dp
          if (new_angle > 180) new_angle = new_angle - 360._dp
          if (new_angle < -180) new_angle = new_angle + 360._dp
          d_angle = (new_angle - old_angle)*pi/180._dp
          call protein_bend(prot,i,d_angle,angle_type,info=info)
          if(info/=0) stop 'Failure in initial_config for PHI'
       end do

       angle_type = 'PSI';
       do i=1,prot%n_res-1
          call random_number(r)
          r = r-0.5_dp
          old_angle = get_dihedral(prot%pep(i),angle_type)
          new_angle = r*360._dp
          if (new_angle > 180) new_angle = new_angle - 360._dp
          if (new_angle < -180) new_angle = new_angle + 360._dp
          d_angle = (new_angle - old_angle)*pi/180._dp
          call protein_bend(prot,i,d_angle,angle_type,info=info)
          if(info/=0) stop 'Failure in initial_config for PSI'
       end do


    case('random')
       ! Randomizes all of the dihedral angles, but not omega
       angle_type = 'PHI';
       do i =mobile_range(1),mobile_range(2)
          if (i==1)cycle
          call random_number(r)
          r = r-0.5_dp
          old_angle = get_dihedral(prot%pep(i),angle_type)
          new_angle = r*360._dp
          if (new_angle > 180) new_angle = new_angle - 360._dp
          if (new_angle < -180) new_angle = new_angle + 360._dp
          d_angle = (new_angle - old_angle)*pi/180._dp
          call protein_bend(prot,i,d_angle,angle_type,info=info)
          if(info/=0) stop 'Failure in initial_config for PHI'
       end do

       angle_type = 'PSI';
       do i=mobile_range(1),mobile_range(2)
          if(i==prot%n_res)cycle
          call random_number(r)
          r = r-0.5_dp
          old_angle = get_dihedral(prot%pep(i),angle_type)
          new_angle = r*360._dp
          if (new_angle > 180) new_angle = new_angle - 360._dp
          if (new_angle < -180) new_angle = new_angle + 360._dp
          d_angle = (new_angle - old_angle)*pi/180._dp
          call protein_bend(prot,i,d_angle,angle_type,info=info)
          if(info/=0) stop 'Failure in initial_config for PSI'
       end do

    case('helix')
       write(*,*)'HERE'
       angle_type = 'PHI';
       do i =mobile_range(1),mobile_range(2)
          if (i==1)cycle
          old_angle = get_dihedral(prot%pep(i),angle_type)
          new_angle = -60.0_dp
          d_angle = (new_angle - old_angle)*pi/180._dp
          call protein_bend(prot,i,d_angle,angle_type,info=info)
          if(info/=0) stop 'Failure in initial_config for PHI'
       end do

       angle_type = 'PSI';
       do i=mobile_range(1),mobile_range(2)
          if(i==prot%n_res)cycle
          old_angle = get_dihedral(prot%pep(i),angle_type)
          new_angle = -40._dp
          d_angle = (new_angle - old_angle)*pi/180._dp
          call protein_bend(prot,i,d_angle,angle_type,info=info)
          if(info/=0) stop 'Failure in initial_config for PSI'
       end do
    case('n/a')
       ! write(*,*) 'No changes desired from initial input structure
       return

    case default
       write(*,*) 'No or incorrect command received for initial_config. Initial structure will be as read.'

    end select
  end subroutine initial_config
  ! -----------------------------------------------------------------------------

  subroutine omega_180(prot)
    ! SUBROUTINE: omega_180(peptide_list)
    ! Sets all omega angles to 180
    ! Used particularly in reading sequence
    type(protein),intent(inout) :: prot                     !Polypeptide
    real(dp) :: old_angle, new_angle                        !Old/New angles (degrees)
    real(dp) :: d_angle                                     !Change in angle(radians)
    integer :: i,j,k                                        !Dummy integers
    integer :: info                                         !Info from subroutines
    character(len=10) :: angle_type                         !Angle type(Phi/Psi/Omega/)
    character(len=80) :: error                               !Error message

    angle_type = 'OMEGA'
    do i=2,prot%n_res
       old_angle = get_dihedral(prot%pep(i),angle_type)
       new_angle = 180._dp
       d_angle = (new_angle-old_angle)*pi/180._dp
       if(d_angle>2*pi) d_angle = d_angle-2*pi
       call protein_bend(prot,i,d_angle,angle_type,info=info,errmsg=error)
       if(info/=0) then
          write(*,*)error
          !stop 'Failure in initial_config for OMEGA'
       end if
    end do
  end subroutine omega_180
  ! -----------------------------------------------------------------------------

  subroutine melt(prot,mag)
    ! SUBROUTINE: melt(protein,log_2(melt_moves))
    ! Subroutine for melting a protein by an order of magnitude of moves.
    ! Generates and accepts 2^mag moves to start a MC
    type(protein), intent(inout) :: prot                    !Protein
    integer, intent(in) :: mag                              !Log_2(moves)
    integer :: moves                                        !Number melting moves
    integer :: i_move                                       !Residue integer
    integer :: pro_length                                   !Length of proteins
    real(dp) :: ran                                         !Random number
    real(dp) :: d_angle                                     !Change in angle (radians)
    character(len=3) :: option                              !Dihedral angle type
    integer :: i                                            !Dummy angle

    ! Initialization
    moves = 2**mag
    pro_length = prot%n_res
    do
       if(prot%pep(pro_length)%heteroatom) then
          pro_length = pro_length - 1
       else
          exit
       end if
    end do
    ! Sets mobile range to whole length if not or poorly set prior
    ! Duplicated in the main, in case melt is not used
    if (mobile_range(1)<= 0) mobile_range(1) = 1
    if (mobile_range(2)<=0 .or. mobile_range(2)>pro_length) mobile_range(2) = pro_length

    ! The melt
    do i=1,moves
       call random_number(ran)
       i_move = ceiling(ran * (mobile_range(2)-mobile_range(1)+1) + mobile_range(1)-1)
       call random_number(ran)
       if (ran > 0.5) then
          option = 'PHI'
       else
          option = 'PSI'
       end if
       call random_number(ran)
       d_angle = (ran - 0.5_dp)*360._dp*(pi/180._dp)
       call protein_bend(prot,i_move,d_angle,option)
    end do

  end subroutine melt
  ! -----------------------------------------------------------------------------
end module angle_manip
