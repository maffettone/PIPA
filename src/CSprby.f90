! LAST EDIT: Phil Maffettone 2017-02-13
! Module for reading in CS probabilities from various programs, and
! incorperating into the triplet probabilities from the ramachandran distribution
module CSprby
  use types
  use maths
  use cmd_io

  implicit none

  private
  public :: incorpCS

contains

  ! --------------------Public Contents-----------------------------------------
  ! SUBROUTINE: incorpCS(triplet_probability, number_residues)
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  ! FUNCTION: dangle_read(directory,number_residues,triplet_probability)
  ! FUNCTION: dangle_old_read(directory,number_residues,triplet_probability)
  ! FUNCTION: talos_read(directory,number_residues,triplet_probability)
  ! ----------------------------------------------------------------------------

  subroutine incorpCS(trip_prby,n_res)
    ! SUBROUTINE: incorpCS(triplet_probability, number_residues)
    ! Selection procedure to change triplet probabilities according to chemical shifts

    real(dp), intent(inout) :: trip_prby(:,:,:)             !Triplet probability as defined in ramachandran.f90
    integer, intent(in) :: n_res                            !Number of residues
    real(dp), allocatable :: cs_prby(:,:)                   !Chemical shift probability
    integer :: i                                            !Dummy integer

    select case(cs_opt)
    case('dangle_old')
       cs_prby = dangle_old_read(cs_dir,n_res,trip_prby)
       if (NMR_weight >= 0.) then
          do i=1,n_res
             trip_prby(3,:,i) = trip_prby(3,:,i)*(cs_prby(:,i)**NMR_weight)
             trip_prby(3,:,i) = trip_prby(3,:,i)/sum(trip_prby(3,:,i))
          end do
       else
          do i=1,n_res
             trip_prby(3,:,i) = (trip_prby(3,:,i)**(abs(NMR_weight)))*cs_prby(:,i)
             trip_prby(3,:,i) = trip_prby(3,:,i)/sum(trip_prby(3,:,i))
          end do
       end if
    case('dangle')
       cs_prby = dangle_read(cs_dir,n_res,trip_prby)
       if (NMR_weight >= 0.) then
          do i=1,n_res
             trip_prby(3,:,i) = trip_prby(3,:,i)*(cs_prby(:,i)**NMR_weight)
             trip_prby(3,:,i) = trip_prby(3,:,i)/sum(trip_prby(3,:,i))
          end do
       else
          do i=1,n_res
             trip_prby(3,:,i) = (trip_prby(3,:,i)**(abs(NMR_weight)))*cs_prby(:,i)
             trip_prby(3,:,i) = trip_prby(3,:,i)/sum(trip_prby(3,:,i))
          end do
       end if
    case('talos')
       cs_prby = talos_read(cs_dir,n_res,trip_prby)
       if (NMR_weight >= 0.) then
          do i=1,n_res
             trip_prby(3,:,i) = trip_prby(3,:,i)*(cs_prby(:,i)**NMR_weight)
             trip_prby(3,:,i) = trip_prby(3,:,i)/sum(trip_prby(3,:,i))
          end do
       else
          do i=1,n_res
             trip_prby(3,:,i) = (trip_prby(3,:,i)**(abs(NMR_weight)))*cs_prby(:,i)
             trip_prby(3,:,i) = trip_prby(3,:,i)/sum(trip_prby(3,:,i))
          end do
       end if
    case('none')
       return
    case default
       return
    end select

  end subroutine incorpCS
  ! -----------------------------------------------------------------------------

  function dangle_read(dir,n_res,trip_prby) result(r)
    ! FUNCTION: dangle_read(directory,number_residues,triplet_probability)
    ! Reads set of Posterios_X.tsv files from Tim Steven's change in predictor.py
    ! From email on 24-11-2016
    character(len=*), intent(in) :: dir                     !Directory for .bay files
    integer, intent(in) :: n_res                            !Number of residues
    real(dp), intent(in) :: trip_prby(:,:,:)                !Previous triplet probability for size
    real(dp), allocatable :: r(:,:)                         !Return residue conformation probability
    integer :: phi, psi                                     !Dihedral angles
    real(dp) :: p                                           !Posterior probability of voxel
    integer :: info                                         !Error information. 0 for success.
    integer :: unit                                         !File unit
    integer :: i,j,k                                        !Dummy integers
    character(len=80) :: filename                           !Filename
    character(len=80) :: errmsg                             !Error Message

    info = 0
    ! Initialize the return probability to the same spacing as Ramachandran, and default to 1.0
    allocate(r(size(trip_prby,2), size(trip_prby,3)))
    r(:,:) = 1.0_dp

    do i=3,n_res-3
       r(:,i) = 0.0_dp                                      !All unavailable values are 0.
       write(filename, "(A10,I0,A4)") 'Posterior_',i,'.tsv'
       filename = trim(dir)//'/'//trim(filename)
       open(newunit=unit, file=trim(filename),status='old', iostat=info, iomsg=errmsg)
       if (info /= 0 )then
          info = 0
          write(*,*) 'DANGLE For Residue ',i, ': ', errmsg
          cycle
       end if
       do
          if(info<0)exit
          read(unit,*,iostat=info,iomsg=errmsg)phi,psi,p
          if(info<0)exit
          if(info>0) then
             write(*,*)'Error in reading DANGLE, Residue',i,errmsg
          end if
          do j=0,1
             do k=1,2
                r(72*(2*phi/10 + j) + 2*psi/10 + k,i) = p
             end do
          end do
       end do
    end do
    close(unit)
  end function dangle_read
  ! -----------------------------------------------------------------------------
  function dangle_old_read(dir, n_res,trip_prby) result(r)
    ! FUNCTION: dangle_old_read(directory,number_residues,triplet_probability)
    ! Reads set of .bay files for each available residue from WJF's copy of DANGLE
    ! Presently operates on the format which I assume from careful examination of
    ! Python source and matlab comparisons.
    character(len=*), intent(in) :: dir                     !Directory for .bay files
    integer, intent(in) :: n_res                            !Number of residues
    real(dp), intent(in) :: trip_prby(:,:,:)                !Previous triplet probability for size
    real(dp), allocatable :: r(:,:)                         !Return residue conformation probability
    real(dp) :: p                                           !Probability of a voxel
    integer :: info                                         !Error information. 0 for success.
    integer :: unit                                         !File unit
    integer :: i,j,k                                        !Dummy integers
    character(len=80) :: filename                           !Filename
    character(len=80) :: errmsg                             !Error Message

    info = 0
    ! Initialize the return probability to the same spacing as Ramachandran, and default to 1.0
    allocate(r(size(trip_prby,2), size(trip_prby,3)))
    r(:,:) = 1.0_dp

    do i=3,n_res-3
       write(filename, "(A4,I0,A4)") 'Res_',i,'.bay'
       filename = trim(dir)//'/'//trim(filename)
       open(newunit=unit, file=trim(filename),status='old', iostat=info, iomsg=errmsg)
       if (info /= 0 )then
          info = 0
          write(*,*) 'DANGLE For Residue ',i, ': ', errmsg
          cycle
       end if
       do j=1,36
          do k=1,36
             read(unit,*) p
             r(72*(2*(k-1)) + (72-(2*(j-1))) ,i) = p
             r(72*(2*(k-1)) + (72-(2*(j-1)+1)) ,i) = p
             r(72*(2*(k-1)+1) + (72-(2*(j-1))) ,i) = p
             r(72*(2*(k-1)+1)  + (72-(2*(j-1)+1)) ,i) = p
          end do
       end do
    end do
    close(unit)
  end function dangle_old_read
  ! -----------------------------------------------------------------------------

  function talos_read(dir,n_res,trip_prby) result(r)
    ! FUNCTION: talos_read(directory,number_residues,triplet_probability)
    ! Reads the file predABP.tab from a directory of TALOS-N output
    character(len=*), intent(in) :: dir                     !Directory for .bay files
    integer, intent(in) :: n_res                            !Number of residues
    real(dp), intent(in) :: trip_prby(:,:,:)                !Previous triplet probability for size
    real(dp), allocatable :: r(:,:)                         !Return residue conformation probability
    real(dp) :: dump                                        !Dump floats
    real(dp) :: list(324)                                   !List of
    integer :: info                                         !Error information. 0 for success.
    integer :: unit                                         !File unit
    integer :: i                                            !Read residue number
    integer :: j,k,l,m                                          !Dummy integers
    character :: a                                          !Dummy character
    character(len=80) :: filename                           !Filename
    character(len=80) :: line                               !Read line
    character(len=80) :: errmsg                             !Error Message

    info = 0
    ! Initialize the return probability to the same spacing as Ramachandran, and default to 1.0
    allocate(r(size(trip_prby,2), size(trip_prby,3)))
    r(:,:) = 1.0_dp
    filename = trim(dir)//'/predABP.tab'
    open(newunit=unit, file=trim(filename),status='old', iostat=info, iomsg=errmsg)
    if (info /= 0 )then
       info = 0
       write(*,*) 'No TALOS-N data found.'
       return
    end if

    do                                                      !Reading through eof
       if (info < 0) exit
       read(unit,'(A80)', iostat=info) line
       if (info < 0) exit
       line = adjustl(line)                                 !Ignoring leading blanks
       if (line(1:6) == 'FORMAT') then
          read(unit,*)

          ! Reading through long line of 'Q-scores'
          do
             read(unit,"(A1,I4,A1,A1,I3,I3,I4,F6.2,F8.3,323F6.3)",iostat=info,iomsg=errmsg)a,i,a,a,j,k,j,dump,list(:)
             if (info < 0) exit                             !EOF
             if (info > 0) then
                write(*,*)'Error in reading TALOS file at residue',i
                write(*,*)errmsg
                exit
             end if
             if (i>n_res) then
                write(*,*) 'TALOS-N file contains probability information for residues not in protein!'
                exit
             end if
             do j=1,18
                do k=1,18
                   do l=1,4
                      do m=0,3
                         r(72*(4*(j-1)+m) + 4*(k-1)+l,i) = list(18*(j-1)+k)
                      end do
                   end do
                end do
             end do
          end do

       end if
    end do
    close(unit)
  end function talos_read
  ! -----------------------------------------------------------------------------
end module CSprby
