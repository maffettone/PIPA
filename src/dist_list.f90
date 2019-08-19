! LAST EDIT: Phil Maffettone 2016-04-08
! All distance lists default to using only alpha carbonsXSXS
! TODO: More comprehensive functions for dist mat (CC,CO,CN,etc) for incoperationg form factors
module dist_list

  use types
  use class_proteins
  use sorting

  implicit none

  private
  public dist_mat, sort_dl, chisq, change_dist_mat
contains
  ! --------------------Public Contents-----------------------------------------
  ! FUNCTION: dist_mat(protein, [flag])
  ! SUBROUTINE: change_dist_mat(distance_matrix,protein,residue_index,[flag])
  ! SUBROUTINE: sort_dl(distance_matrix,sorted_list)
  ! FUNCTION: chisq(sorted_experiment,sorted_simulation)
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  ! FUNCTION: bb_dist_mat(protein)
  ! FUNCTION: bbsc_dist_mat(protein)
  ! FUNCTION: ca_dist_mat(protein)
  ! FUNCTION: noh_dist_mat(prot)
  ! ----------------------------------------------------------------------------

  function dist_mat(prot,flag) result(r)
    ! FUNCTION: dist_mat(protein, flag)
    ! Directing funciton to for getting triangular distance matrix from protein
    ! Only used for the first utility of generating a previously unallocated matrix
    type(protein), intent(in) :: prot                       !Protein
    character(len=*), intent(in), optional :: flag          !Flag for dist inclusion
    real(dp), allocatable :: r(:,:)                         !Triangular matrix of distances


    if (present(flag)) then
       select case(trim(flag))
       case('ca')
          r = ca_dist_mat(prot)
       case('c_only')
          r = c_dist_mat(prot)
       case('n_only')
          r = n_dist_mat(prot)
       case('bb_generic')
          r = bb_dist_mat(prot)
       case('bb_noh')
          r = noh_dist_mat(prot)
       case('total_generic')
          r = bbsc_dist_mat(prot)
       case default
          r = ca_dist_mat(prot)
       end select
    else
       r = ca_dist_mat(prot)
    end if
  end function dist_mat
  ! -----------------------------------------------------------------------------

  function bb_dist_mat(prot) result(r)
    ! FUNCTION: bb_dist_mat(protein)
    ! Total backbone distance matrix irrespective of atom type
    type(protein), intent(in) :: prot                       !Protein
    real(dp), allocatable :: r(:,:)                         !Triangular matrix of distances
    integer :: i,j,k,l                                      !Dummy integers

    allocate(r(size(prot%coordinates,2),size(prot%coordinates,2)))
    r =0._dp

    k=1
    do i=1,size(prot%coordinates,2)
       l=1
       do j=1,i-1
          r(l,k) = norm2(prot%coordinates(:,j) - prot%coordinates(:,i))
          l=l+1
       end do
       k=k+1
    end do
  end function bb_dist_mat
  ! -----------------------------------------------------------------------------
  
  function noh_dist_mat(prot) result(r)
    ! FUNCTION: noh_dist_mat(protein)
    ! Total backbone distance matrix not including hydrogen (for x rays)
    type(protein), intent(in) :: prot                       !Protein
    real(dp), allocatable :: r(:,:)                         !Triangular matrix of distances
    integer :: i,j,k,l                                      !Dummy integers

    allocate(r(size(prot%coordinates,2),size(prot%coordinates,2)))
    r =0._dp

    k=1
    do i=1,size(prot%coordinates,2)
       l=1
       if(prot%h_mask(i)) cycle
       do j=1,i-1
          if(prot%h_mask(j)) cycle
          r(l,k) = norm2(prot%coordinates(:,j) - prot%coordinates(:,i))
          l=l+1
       end do
       k=k+1
    end do
  end function noh_dist_mat
  ! -----------------------------------------------------------------------------

  function bbsc_dist_mat(prot) result(r)
    ! FUNCTION: bbsc_dist_mat(protein)
    ! Total backbone and sidechain distance matrix irrespective of atom type
    type(protein), intent(in) :: prot                       !Protein
    real(dp), allocatable :: r(:,:)                         !Triangular matrix of distances

    write(*,*) 'Distance matrix function under development. Choose another.'
    stop
  end function bbsc_dist_mat
  ! -----------------------------------------------------------------------------

  function ca_dist_mat(prot) result(r)
    ! FUNCTION: ca_dist_mat(protein)
    ! Forms alpha carbon only distance matrix
    type(protein), intent(in) :: prot                       !Protein
    real(dp), allocatable :: r(:,:)                         !Triangular matrix of distances
    integer :: i,j,k,l                                      !Dummy integers

    k=0
    do i=1,size(prot%ca_mask)
       if(prot%ca_mask(i)) k=k+1
    end do
    allocate(r(k,k))

    r = 0._dp
    k=1;
    do i=1,size(prot%ca_mask)
       l=1;
       if(.not. prot%ca_mask(i)) cycle
       do j=1,i-1
          if(.not. prot%ca_mask(j)) cycle
          r(l,k) = norm2(prot%coordinates(:,j) - prot%coordinates(:,i))
          l=l+1
       end do
       k=k+1
    end do
  end function ca_dist_mat
  ! -----------------------------------------------------------------------------

  function c_dist_mat(prot) result(r)
    ! FUNCTION: c_dist_mat(protein)
    ! Forms the carbon only distance matrix
    type(protein), intent(in) :: prot                       !Protein
    real(dp), allocatable :: r(:,:)                         !Triangular matrix of distances
    integer :: i,j,k,l                                      !Dummy integers

    k=0
    do i=1,size(prot%c_mask)
       if(prot%c_mask(i)) k=k+1
    end do
    allocate(r(k,k))

    r = 0._dp
    k=1;
    do i=1,size(prot%c_mask)
       l=1;
       if(.not. prot%c_mask(i)) cycle
       do j=1,i-1
          if(.not. prot%c_mask(j)) cycle
          r(l,k) = norm2(prot%coordinates(:,j) - prot%coordinates(:,i))
          l=l+1
       end do
       k=k+1
    end do
  end function c_dist_mat
  ! -----------------------------------------------------------------------------

  function n_dist_mat(prot) result(r)
    ! FUNCTION: n_dist_mat(protein)
    ! Forms the nitrogen only distance matrix
    type(protein), intent(in) :: prot                       !Protein
    real(dp), allocatable :: r(:,:)                         !Triangular matrix of distances
    integer :: i,j,k,l                                      !Dummy integers

    k=0
    do i=1,size(prot%n_mask)
       if(prot%c_mask(i)) k=k+1
    end do
    allocate(r(k,k))

    r = 0._dp
    k=1;
    do i=1,size(prot%n_mask)
       l=1;
       if(.not. prot%n_mask(i)) cycle
       do j=1,i-1
          if(.not. prot%n_mask(j)) cycle
          r(l,k) = norm2(prot%coordinates(:,j) - prot%coordinates(:,i))
          l=l+1
       end do
       k=k+1
    end do
  end function n_dist_mat
  ! -----------------------------------------------------------------------------

  subroutine change_dist_mat(mat,prot,index,flag)
    ! SUBROUTINE: change_dist_mat(distance_matrix,protein,residue_index,[flag])
    ! Directing function to change in distance matrix from protein bend at
    ! residue index. Protein input is the protein after a change
    type(protein), intent(in) :: prot                       !Protein
    real(dp), intent(inout) :: mat(:,:)                     !Triangular matrix of distances
    integer, intent(in) :: index                            !Residue changed index
    character(len=*), intent(in), optional :: flag          !Flag for dist inclusion

    if (present(flag)) then
       select case(trim(flag))
       case('ca')
          call change_ca_dist(mat,prot,index)
       case('c_only')
          call change_c_dist(mat,prot,index)
       case('n_only')
          call change_n_dist(mat,prot,index)
       case('bb_noh')
          call change_noh_dist(mat,prot,index)
       case('bb_generic')
          call change_bb_dist(mat,prot,index)
       case default
          call change_ca_dist(mat,prot,index)
       end select
    else
       call change_ca_dist(mat,prot,index)
    end if
  end subroutine change_dist_mat
  ! -----------------------------------------------------------------------------

  
  subroutine change_bb_dist(mat,prot,index)
    ! SUBROUTINE: change_bb_dist(distance_matrix,protein,residue_index)
    type(protein), intent(in) :: prot                       !Protein
    real(dp), intent(inout) :: mat(:,:)                     !Triangular matrix of distances
    integer, intent(in) :: index                            !Residue changed index
    integer :: atom_index                                   !Index in the coordinate matrix
    integer :: i,j,k,l                                      !Dummy integers

    atom_index = prot%pep(index)%backbone(1)%idx
    k=atom_index
    do i = atom_index, size(prot%coordinates,2)
       l=1
       do j= 1,atom_index
          mat(l,k) = norm2(prot%coordinates(:,j) - prot%coordinates(:,i))
          l=l+1
       end do
       k=k+1
    end do
  end subroutine change_bb_dist
  ! -----------------------------------------------------------------------------

   subroutine change_noh_dist(mat,prot,index)
    ! SUBROUTINE: change_noh_dist(distance_matrix,protein,residue_index)
    type(protein), intent(in) :: prot                       !Protein
    real(dp), intent(inout) :: mat(:,:)                     !Triangular matrix of distances
    integer, intent(in) :: index                            !Residue changed index
    integer :: atom_index                                   !Index in the coordinate matrix
    integer :: i,j,k,l                                      !Dummy integers

    atom_index = prot%pep(index)%backbone(1)%idx
    k=atom_index
    do i = atom_index, size(prot%coordinates,2)
       l=1
       if(prot%h_mask(i))cycle
       do j= 1,atom_index
          if(prot%h_mask(j))cycle
          mat(l,k) = norm2(prot%coordinates(:,j) - prot%coordinates(:,i))
          l=l+1
       end do
       k=k+1
    end do
  end subroutine change_noh_dist
  ! -----------------------------------------------------------------------------
  
  subroutine change_ca_dist(mat,prot,index)
    ! SUBROUTINE: change_ca_dist(distance_matrix,protein,residue_index)
    type(protein), intent(in) :: prot                       !Protein
    real(dp), intent(inout) :: mat(:,:)                     !Triangular matrix of distances
    integer, intent(in) :: index                            !Residue changed index
    integer :: c_index                                      !Index in C coordinate matrix
    integer :: i,j,k,l                                      !Dummy integers
    integer :: true_count                                   !Count of mask prior to c_index

    c_index = 0
    do i = 1,size(prot%pep(index)%backbone)
       if(prot%pep(index)%backbone(i)%label == 'CA  ') then
          c_index = prot%pep(index)%backbone(i)%idx
       end if
    end do
    if (c_index == 0) then
       write(*,*) 'CA not found in changing CA correlation matrix'
       stop
    end if

    ! Deciding how to initialize k,l index for diagonal matrix
    true_count=0
    do i=1,c_index
       if(prot%ca_mask(i)) true_count = true_count +1
    end do

    k=true_count
    do i=c_index,size(prot%ca_mask)
       l=1
       if(.not. prot%ca_mask(i)) cycle
       do j=1,c_index
          if(.not. prot%ca_mask(j)) cycle
          mat(l,k) = norm2(prot%coordinates(:,j) - prot%coordinates(:,i))
          l=l+1
       end do
       k=k+1
    end do
  end subroutine change_ca_dist
  ! -----------------------------------------------------------------------------

  subroutine change_c_dist(mat,prot,index)
    ! SUBROUTINE: change_c_dist(distance_matrix,protein,residue_index)
    type(protein), intent(in) :: prot                       !Protein
    real(dp), intent(inout) :: mat(:,:)                     !Triangular matrix of distances
    integer, intent(in) :: index                            !Residue changed index
    integer :: c_index                                      !Index in C coordinate matrix
    integer :: i,j,k,l                                      !Dummy integers
    integer :: true_count                                   !Count of mask prior to c_index

    c_index = 0
    do i = 1,size(prot%pep(index)%backbone)
       if(prot%pep(index)%backbone(i)%sym == 'C ') then
          c_index = prot%pep(index)%backbone(i)%idx
       end if
    end do
    if (c_index == 0) then
       write(*,*) 'C not found in changing C correlation matrix'
       stop
    end if

    ! Deciding how to initialize k,l index for diagonal matrix
    true_count=0
    do i=1,c_index
       if(prot%c_mask(i)) true_count = true_count +1
    end do

    k=true_count
    do i=c_index,size(prot%c_mask)
       l=1
       if(.not. prot%c_mask(i)) cycle
       do j=1,c_index
          if(.not. prot%c_mask(j)) cycle
          mat(l,k) = norm2(prot%coordinates(:,j) - prot%coordinates(:,i))
          l=l+1
       end do
       k=k+1
    end do
  end subroutine change_c_dist
  ! -----------------------------------------------------------------------------

  subroutine change_n_dist(mat,prot,index)
    ! SUBROUTINE: change_n_dist(distance_matrix,protein,residue_index)
    type(protein), intent(in) :: prot                       !Protein
    real(dp), intent(inout) :: mat(:,:)                     !Triangular matrix of distances
    integer, intent(in) :: index                            !Residue changed index
    integer :: n_index                                      !Index in C coordinate matrix
    integer :: i,j,k,l                                      !Dummy integers
    integer :: true_count                                   !Count of mask prior to c_index

    n_index = 0
    do i = 1,size(prot%pep(index)%backbone)
       if(prot%pep(index)%backbone(i)%sym == 'N ') then
          n_index = prot%pep(index)%backbone(i)%idx
       end if
    end do
    if (n_index == 0) then
       write(*,*) 'N not found in changing N correlation matrix'
       stop
    end if

    ! Deciding how to initialize k,l index for diagonal matrix
    true_count=0
    do i=1,n_index
       if(prot%n_mask(i)) true_count = true_count +1
    end do

    k=true_count
    do i=n_index,size(prot%n_mask)
       l=1
       if(.not. prot%n_mask(i)) cycle
       do j=1,n_index
          if(.not. prot%c_mask(j)) cycle
          mat(l,k) = norm2(prot%coordinates(:,j) - prot%coordinates(:,i))
          l=l+1
       end do
       k=k+1
    end do
  end subroutine change_n_dist
  ! -----------------------------------------------------------------------------

  subroutine sort_dl(mat,r)
    ! SUBROUTINE: sort_dl(distance_matrix,sorted_list)
    ! Takes in the upper triangular distance matrix and converts
    ! it into a sorted distance list
    real(dp), intent(in) :: mat(:,:)                        !Uppr triangular matrix of distances
    real(dp), intent(inout) :: r(:)                         !Sorted distance list
    integer :: i,j,k                                        !Dummy integers

    k=1
    do i=1,size(mat,1)
       do j=1,i-1
          r(k) = mat(j,i)
          k=k+1
       end do
    end do
    call hpsort(r)
  end subroutine sort_dl
  ! -----------------------------------------------------------------------------

  function chisq(dat,sim) result(r)
    ! FUNCTION: chisq(sorted_experiment,sorted_simulation)
    real(dp), intent(in) :: dat(:)                          !Sorted dist list from experiment
    real(dp), intent(in) :: sim(:)                          !Sorted dist list from simulation
    real(dp) :: r                                           !Chi squared value
    integer :: i                                            !Dummy integer

    if (size(dat) /= size(sim)) then
       write(*,*) 'The number of distances in the data file does not match the simulation'
       write(*,*) size(dat),size(sim)
       stop
    end if
    
    r = 0._dp
    do i = 1,size(dat)
       r = r + ((dat(i) - sim(i))**2)/(dat(i))
    end do
  end function chisq
  ! -----------------------------------------------------------------------------


end module dist_list
