!LAST EDIT: Phil Maffettone 2016-10-17
! Fixed Sampling
module ramachandran
  use types
  use maths
  use class_proteins

  implicit none

  private
  public :: Triplet_Prby, write_triplet,cumulative_prby,sample_cum, Uniform_Triplet_Prby

  type rama
     character(len=3) :: center, neigh
     real(dp) :: phi,psi,prby,log_p
  end type rama

contains

  ! --------------------Public Contents-----------------------------------------
  ! SUBROUTINE: Triplet_Prby(peptide_list, phi-psi-prby)
  ! SUBROUTINE: write_triplet(trip_prby,file_stem)
  ! FUNCTION: cumulative_prby(triplet_prby)
  ! FUNCTION: sample_cum(cumulative_prby,random,index)
  ! SUBROUTINE: Uniform_Triplet_Prby(peptide_list,phi-psi-prby)
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  ! SUBROUTINE: norm_P(probability_dist)
  ! ----------------------------------------------------------------------------

  subroutine Triplet_Prby(pep,trip_prby)
    ! SUBROUTINE: Triplet_Prby(phi-psi-prby,residue) are calculated
    ! The innermost labels the phi/psi combo and the respective probability.
    ! The second is all the bins. The 3rd dim runs over the residues of the protein.

    type(residue),intent(in) :: pep(:)                     !Protein
    real(dp), intent(out), allocatable :: trip_prby(:,:,:) !Triplet prby as defined above
    real(dp), allocatable :: prbyL(:,:)                    !Left hand relevant ln(P)
    real(dp), allocatable :: prbyR(:,:)                    !Right hand relevnt ln(P)
    real(dp), allocatable :: prbyALL(:)                    !Total ramachandran ln(P)
    type(rama), allocatable :: ramaL(:),ramaR(:)           !Objects for reading rama files
    character(len=3) :: left,center,right                  !Neighboring residue lables
    character(len=80) :: filename,line                     !Filename, and buffer line
    integer :: bin_size                                    !Angle increments in histogram
    integer :: n_bin                                       !Number of bins
    integer :: n_res                                       !Number of residues
    integer :: n_dist                                      !Number of distributions
    integer :: i,j,k,l,m                                   !Dummy integers
    integer :: unit                                        !File I/O unit
    logical :: file_exists

    bin_size = 5
    n_bin = (360/5)**2
    n_dist = 420                !420 distributions for 20*(20+1)
    n_res = size(pep)
    allocate(ramaL(n_bin*n_dist),ramaR(n_bin*n_dist),trip_prby(3,n_bin,n_res),&
         prbyL(3,n_bin),prbyR(3,n_bin),prbyALL(n_bin))
    trip_prby = 0._dp
    filename = '/Users/alggroup/Documents/Development/Common/Dunbrack/ramaTCBIGright.txt'
    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
       stop 'Ramachandran distributions have been removed from library'
    end if
    open(newunit=unit, file=filename, status='old')

    ! Working right first
    do i=1,n_bin*n_dist
       read(unit,*) ramaR(i)%center,line,ramaR(i)%neigh,ramaR(i)%phi,&
            ramaR(i)%psi,ramaR(i)%prby,ramaR(i)%log_p
    end do
    close(unit)

    filename = '/Users/alggroup/Documents/Development/Common/Dunbrack/ramaTCBIGleft.txt'
    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
       stop 'Ramachandran distributions have been removed from library'
    end if
    open(newunit=unit, file=filename, status='old')

    do i=1,n_bin*n_dist
       read(unit,*) ramaL(i)%center,line,ramaL(i)%neigh,ramaL(i)%phi,&
            ramaL(i)%psi,ramaL(i)%prby,ramaL(i)%log_p
    end do
    close(unit)

    ! Calculating triplet probabilties. The first is just the right hand probability
    ! The last is just the left prby, the others are left*right/all
    ! The log(p) is actually the -log(p)
    i=1
    prbyR = 0._dp
    prbyL = 0._dp
    prbyALL = 0._dp
    center = get_label(pep(i))
    right = get_label(pep(i+1))
    do j=1,20                   !Finding correct center label
       if (ramaR((j-1)*21*n_bin+1)%center == center) then
          do k=1,21             !Finding correct neighbor label
             if(ramaR((j-1)*n_bin*21+(k-1)*n_bin+1)%neigh == right) then
                m = (j-1)*n_bin*21+(k-1)*n_bin
                do l = 1,n_bin  !Saving ln(P) for correct neighbor set
                   prbyR(:,l) = [ramaR(m+l)%phi,ramaR(m+l)%psi,ramaR(m+l)%log_p]
                end do
             end if
          end do
       end if
    end do
    trip_prby(:,:,i) = prbyR
    trip_prby(3,:,i) = exp(-1._dp*trip_prby(3,:,i))

    do i=2,n_res-1
       if(pep(i+1)%heteroatom)exit                          !Only works for heteroatoms after chain
       prbyR = 0._dp
       prbyL = 0._dp
       prbyALL = 0._dp
       center = get_label(pep(i))
       right = get_label(pep(i+1))
       left = get_label(pep(i-1))
       do j=1,20
          if (ramaR((j-1)*21*n_bin+1)%center == center) then
             do k=1,21             !Finding correct neighbor label
                if(ramaR((j-1)*n_bin*21+(k-1)*n_bin+1)%neigh == right) then
                   m = (j-1)*n_bin*21+(k-1)*n_bin
                   do l = 1,n_bin  !Saving ln(P) for correct neighbor set
                      prbyR(:,l) = [ramaR(m+l)%phi,ramaR(m+l)%psi,ramaR(m+l)%log_p]
                   end do
                else if(ramaR((j-1)*n_bin*21+(k-1)*n_bin+1)%neigh == 'ALL') then
                   m = (j-1)*n_bin*21+(k-1)*n_bin
                   do l = 1,n_bin  !Saving ln(P) for total neighbor set
                      prbyALL(l) = ramaR(m+l)%log_p
                   end do
                end if
             end do
          end if
       end do
       do j=1,20
          if (ramaL((j-1)*21*n_bin+1)%center == center) then
             do k=1,21             !Finding correct neighbor label
                if(ramaL((j-1)*n_bin*21+(k-1)*n_bin+1)%neigh == left) then
                   m = (j-1)*n_bin*21+(k-1)*n_bin
                   do l = 1,n_bin  !Saving ln(P) for correct neighbor set
                      prbyL(:,l) = [ramaL(m+l)%phi,ramaL(m+l)%psi,ramaL(m+l)%log_p]
                   end do
                end if
             end do
          end if
       end do
       ! Check
       if (prbyR(1,1) /= prbyL(1,1)) then
          write(*,*) 'Left and Right Probability not alligned for ', get_label(pep(i))
          stop
       end if
       trip_prby(1:2,:,i) = prbyR(1:2,:)
       trip_prby(3,:,i) = exp(prbyALL(:)-prbyR(3,:)-prbyL(3,:))
    end do

    prbyR = 0._dp
    prbyL = 0._dp
    prbyALL = 0._dp
    center = get_label(pep(i))
    left = get_label(pep(i-1))
    do j=1,20
       if (ramaL((j-1)*21*n_bin+1)%center == center) then
          do k=1,21             !Finding correct neighbor label
             if(ramaL((j-1)*n_bin*21+(k-1)*n_bin+1)%neigh == left) then
                m = (j-1)*n_bin*21+(k-1)*n_bin
                do l = 1,n_bin  !Saving ln(P) for correct neighbor set
                   prbyL(:,l) = [ramaL(m+l)%phi,ramaL(m+l)%psi,ramaL(m+l)%log_p]
                end do
             end if
          end do
       end if
    end do
    trip_prby(:,:,i) = prbyL
    trip_prby(3,:,i) = exp(-1._dp*trip_prby(3,:,i))

    do i=1,n_res
       call norm_P(trip_prby(3,:,i))
    end do

  end subroutine Triplet_Prby
  ! -------------------------------------------------------------------------------

  subroutine Uniform_Triplet_Prby(pep,trip_prby)
    ! SUBROUTINE: Uniform_Triplet_Prby(residue,phi-psi-prby)
    ! Needs the inner dimensions to be set prior to use, but Triplet_Prby
    ! Creates a uniform probability to speed up error testing i/o. Should not be used in the implementation.
    ! The innermost labels the phi/psi combo and the respective probability.
    ! The second is all the bins. The 3rd dim runs over the residues of the protein.

    type(residue),intent(in) :: pep(:)                     !Protein
    real(dp), intent(out), allocatable :: trip_prby(:,:,:) !Triplet prby as defined above
    integer :: bin_size                                    !Angle increments in histogram
    integer :: n_bin                                       !Number of bins
    integer :: n_res                                       !Number of residues
    integer :: i                                           !Dummy integer
    bin_size = 5
    n_bin = (360/5)**2
    n_res = size(pep)
    allocate(trip_prby(3,n_bin,n_res))

    trip_prby(3,:,:) = 1._dp
    do i=1,n_res
       call norm_P(trip_prby(3,:,i))
    end do
  end subroutine Uniform_Triplet_Prby
  ! -----------------------------------------------------------------------------

  subroutine norm_P(p)
    ! Normalizes the probability distribution to sum to unity
    real(dp), intent(inout) :: p(:)
    real(dp) :: sum
    integer :: i
    sum=0._dp
    do i = 1,size(p)
       sum = sum+p(i)
    end do
    p=p/sum
  end subroutine norm_P
  ! -------------------------------------------------------------------------------

  subroutine write_triplet(trip_prby,stem)
    ! Writes all triplett probabilities into xyz files for each residue
    real(dp), intent(in) :: trip_prby(:,:,:)                !Triplet probabilities
    character(len=*),intent(in) :: stem                    !File stem
    character(len=80) :: filename                           !Filename
    character(len=3) :: fnum                                !Filenum
    integer :: i,j                                          !dummy integers
    integer :: unit                                         !file unit

    do i=1,size(trip_prby,3)
       write(fnum,"(I3.3)")i
       filename = trim(stem)//trim(fnum)//'.xyz'
       open(newunit=unit,file=filename,status='replace')
       do j=1,size(trip_prby,2)
          write(unit,"(2f8.2,ES15.5)")trip_prby(:,j,i)
       end do
       close(unit)
    end do
  end subroutine write_triplet
  !-------------------------------------------------------------------------------

  function cumulative_prby(trip_prby) result(r)
    ! FUNCTION: cumulative_prby(triplet_prby)
    ! Creates a cumulative probability distribution funciton from the triplets
    real(dp),intent(in) :: trip_prby(:,:,:)                 !Triplet Probabilities
    real(dp),allocatable :: r(:,:,:)                        !Cumulative Dist
    integer :: i,j                                          !Dumy integers

    allocate(r(size(trip_prby,1),size(trip_prby,2),size(trip_prby,3)))
    r=0._dp
    do i=1,size(r,3)                                        !Running over residues
       r(:,1,i) = trip_prby(:,1,i)
       do j=2,size(r,2)
          r(1:2,j,i) = trip_prby(1:2,j,i)
          r(3,j,i) = r(3,j-1,i)+trip_prby(3,j,i)
       end do
    end do
  end function cumulative_prby
  ! ------------------------------------------------------------------------------

  function sample_cum(cum,ran,index) result(r)
    ! FUNCTION: sample_cum(cumulative_prby,random,index)
    ! Searches the cumulative distribution for an appropriate return result
    real(dp), intent(in) :: cum(:,:,:)                      !Cumulative dist
    real(dp), intent(in) :: ran                             !Random number
    integer, intent(in) :: index                            !Residue index for change
    real(dp) :: r(2)                                        !Passed phi/psi combo
    real(dp) :: ran_loc                                     !Random fo
    integer :: i,j                                          !Dummy integer

    j=0
    do i = 1,size(cum,2)
       ! Stops at the first number where the random is less than the cum value
       ! By definition all other values will be less than the random
       if(ran <= cum(3,i,index)) then
          j = i
          exit
       end if
    end do
    if(j==0) then
       write(*,*)'Failure in sampling the cumulative distribution'
       write(*,*)'Failed in residue, ',index
       stop
    end if
    r = (cum(1:2,j,index))


    ! Casting the sampled angle in the range the dist offers
    ! The probabilities are given at the floor of a 5 degree range
    do i=1,2
       call random_number(ran_loc)
       r(i) = r(i) + 5.*ran_loc
    end do
  end function sample_cum
  !-------------------------------------------------------------------------------

end module ramachandran
