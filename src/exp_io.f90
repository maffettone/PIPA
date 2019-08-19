! LAST EDIT: Phil Maffettone 2017_03_18
module exp_io

  use types
  use class_proteins

  implicit none

  private :: data, segment_data

  type data
     real(dp) :: x
     real(dp) :: y
     real(dp) :: e
     type(data),pointer :: next
  end type data

contains

  ! --------------------Public Contents-----------------------------------------
  ! SUBROUTINE: read_exp(protein,filename,G(r),[info_integer,error_message])
  ! SUBROUTINE: pdf2dl(g(r),distance_list)
  ! SUBROUTINE: read_p_r(filename,G(r),[info_integer,error_message])
  ! SUBROUTINE: write_xye(filename,xye_data,[info_integer,error_message])
  ! SUBROUTINE: read_xye(filename,G(r),[info_integer,error_message])
  ! SUBROUTINE: read_dist_list(filename,distance_list,[info_integer,error_message])
  ! SUBROUTINE: write_dist_list(filename, dist_list,[info_integer,error_message])
  ! SUBROUTINE: read_xy(filename,G(r),[info_integer,error_message])
  ! SUBROUTINE: read_csv(filename,G(r),[info_integer,error_message])
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  ! SUBROUTINE: segment_data(G(r))
  ! ----------------------------------------------------------------------------

  subroutine read_exp(prot,flag,filename,dl,info,errmsg)
    ! SUBROUTINE: read_exp(protein,flag,filename,dist_list,[info_integer,error_messafge])
    ! Choses appropriate read routine for experimental data
    type(protein),intent(in) :: prot                        !Protein
    character(len=*), intent(in), optional :: flag          !Flag for dist inclusion
    character(len=*),intent(in) :: filename                 !Filename
    integer,optional,intent(inout) :: info                  !Flag, 0 for success
    character(len=80),optional,intent(inout) :: errmsg      !Error message
    real(dp),allocatable,intent(out) :: dl(:)               !Distance list
    real(dp),allocatable:: g_r(:,:)                         !x,y,e data for appropriate func
    character(len=4) :: ext                                 !File extension
    character(len=80) :: errmsg_loc                         !Local error message
    integer  :: info_loc                                    !Local flag
    integer :: i,j                                          !Dummy integer


    info_loc = 0
    ! Work backwards from size of filename to period
    i = strlen(filename)
    do while (filename(i:i) /= '.')
       i = i-1
    end do

    if((strlen(filename)-i)>4) then
       if(present(info)) info = 1
       if(present(errmsg)) errmsg = 'No suitable length extension found. Filetype unknown.'
       return
    end if
    ext = filename((i+1):strlen(filename))

    ! Size for dist_list
    j=0;
    if (present(flag)) then
       select case(trim(flag))
       case('ca')
          do i=1,size(prot%ca_mask)
             if (prot%ca_mask(i)) j=j+1
          end do
       case('c_only')
          do i=1,size(prot%c_mask)
             if (prot%c_mask(i)) j=j+1
          end do
       case('n_only')
          do i=1,size(prot%n_mask)
             if (prot%n_mask(i)) j=j+1
          end do
       case('bb_noh')
          do i=1,size(prot%c_mask)
             if (prot%h_mask(i)) j=j+1
          end do
          j=size(prot%coordinates,2)-j
       case('bb_generic')
          j=size(prot%coordinates,2)
       case default
          do i=1,size(prot%ca_mask)
             if (prot%ca_mask(i)) j=j+1
          end do
       end select
    else
       do i=1,size(prot%ca_mask)
          if (prot%ca_mask(i)) j=j+1
       end do
    end if
    allocate(dl(j*(j-1)/2))

    ! Choosing apropriate read routine
    select case(trim(ext))
    case('xy')                                              !X-Y files
       call read_xy(filename,g_r,info_loc,errmsg_loc)
    case('xye')
       call read_xye(filename,g_r,info_loc,errmsg_loc)      !X-Y-Error files
    case('csv')
       call read_csv(filename,g_r,info_loc,errmsg_loc)      !CSV files
    case('out')
       call read_p_r(filename,g_r,info_loc,errmsg_loc)      !OUT from GNOM
    case default
       write(*,*)'Unable to parse file, ',filename
       write(*,*) 'Experimental data not read.'
    end select

    ! Fluffing all non-zero data to have a solid distance list
    j=size(g_r,1)
    do
       if (g_r(j,2) < 0.001) then
          j=j-1
       else
          exit
       end if
    end do
    do while (j<size(dl))
       write(*,*)'Caution: Spreading PDF data to generate enough distances'
       write(*,"('There are ',i10,' distances in the list.')")size(dl)
       call segment_data(g_r)
       j=size(g_r,1)
       do
          if (g_r(j,2) < 0.001) then
             j=j-1
          else
             exit
          end if
       end do
    end do
    call pdf2dl(g_r,dl)

    ! Sending along errors if they were provided
    if(present(info)) info = info_loc
    if(present(errmsg)) errmsg = errmsg_loc

  end subroutine read_exp
  ! ----------------------------------------------------------------------------

  subroutine pdf2dl(pdf,dl)
    ! SUBROUTINE: pdf2dl(g(r),distance_list)
    ! Converts from actual RDF data to distance list discretized by number of pairs
    ! This should be the RADIAL DISTRIBUTION FUNCION, such that the integral corresponds
    ! to the coordination number
    ! First creates a cumulative PDF which goes from 0 to 1, then samples distances for the list

    real(dp),intent(in) :: pdf(:,:)                         !Pair distribution function [:,3]
    real(dp),intent(inout) :: dl(:)                         !Distance list [n(n-1)/2]
    real(dp),allocatable :: cum(:,:)                        !Cumulative pdf
    real(dp) :: dy                                          !Change in y for dist list
    real(dp) :: search                                      !Number to search for
    integer :: i,j,k                                        !Dummy integers

    allocate(cum(size(pdf,1),2))
    cum=0._dp
    i=1
    cum(i,1) = pdf(i,1)
    cum(i,2) = pdf(i,2)
    do i=2,size(pdf,1)
       cum(i,1) = pdf(i,1)
       cum(i,2) = cum(i-1,2)+pdf(i,2)
    end do

    cum(:,2) = cum(:,2)/maxval(cum(:,2))
    dy = 1./size(dl)
    dl = 0.
    k=1
    do i=1,size(dl)
       do j=k,size(cum,1)
          search = i*dy
          if (search>1) search =1._dp
          if(cum(j,2) >= search) then
             k=j
             dl(i) = cum(j,1)
             exit
          end if
       end do
    end do
  end subroutine pdf2dl
  ! -----------------------------------------------------------------------------


  subroutine read_p_r(filename,r,info,errmsg)
    ! SUBROUTINE: read_p_r(filename,P(r),[info_integer,error_message])
    ! Reads in .out file from GNOM, to produce a P(r) function data

    character(len=*),intent(in) :: filename                 !Filename
    integer,optional,intent(out) :: info                    !Flag, 0 for success
    character(len=80),optional,intent(out) :: errmsg        !Error message
    integer :: info_loc                                     !Local flag
    character(len=80) :: errmsg_loc                         !Local error message
    character(len=80) :: line,scanfor                       !Buffer and search
    real(dp),allocatable,intent(out) :: r(:,:)              !x,y,e data
    type(data),pointer :: head                              !Pointer to head of linked list
    integer :: nvals = 0                                    !Number of data read
    type(data),pointer :: ptr                               !Temporary pointer
    type(data),pointer :: tail                              !pointer to tail of list
    real(dp) :: temp_x,temp_y,temp_e                        !Temporary variables
    integer :: i                                            !Dummy integer
    integer :: unit                                         !File unit
    
    info_loc = 0
    open(newunit=unit, file=trim(filename), status='old',IOSTAT = info_loc,IOMSG=errmsg_loc)
    if(info_loc /= 0) then
       if(present(info)) info=info_loc
       if(present(errmsg)) errmsg = errmsg_loc
       return
    end if

    ! Deallocating return data
    if(allocated(r)) deallocate(r)
    nullify(head,tail,ptr)

    ! Searching for P(r)
    scanfor = 'Distance distribution  function of particle'
    do
       read(unit,'(A80)')line
       if(trim(adjustl(line)) == scanfor) exit
    end do
    do i=1,4
       read(unit,*)
    end do

    do                                                      !Creating linked list
       read(unit,*,iostat = info_loc)temp_x,temp_y,temp_e
       if(info_loc /=0) exit
       nvals=nvals+1

       if(.not. associated(head))then                       !No values in list
          allocate(head,stat = info_loc)                    !Allocates new value
          tail => head                                      !Tail pts to new value
          nullify(tail%next)                                !Nullify ptr in new value
          tail%x = temp_x                                   !Store data
          tail%y = temp_y
          tail%e = temp_e
       else
          allocate(tail%next,stat = info_loc)               !Values already in list
          tail =>tail%next                                  !Tail points to new value
          nullify(tail%next)                                !Nullify ptr in new value
          tail%x = temp_x                                   !Store data
          tail%y = temp_y
          tail%e = temp_e
       end if
    end do

    close(unit)
    allocate(r(nvals,3))
    ptr =>head
    i=1
    do                                                      !Places data in return variable
       if(.not. associated(ptr))exit                        !Exits at null
       r(i,:) = [ptr%x,ptr%y,ptr%e]
       ptr=>ptr%next
       i=i+1
    end do

    do                                                      !Deallocates linked list
       ptr => head%next
       deallocate(head)
       if (.not. associated(ptr)) exit
       head => ptr
    end do
  end subroutine read_p_r
  ! ----------------------------------------------------------------------------

  subroutine read_xye(filename,r,info,errmsg)
    ! SUBROUTINE: read_p_r(filename,P(r),[info_integer,error_message])
    ! Reads in .xye file from to produce discrete function data with errors

    character(len=*),intent(in) :: filename                 !Filename
    integer,optional,intent(out) :: info                    !Flag, 0 for success
    character(len=80),optional,intent(out) :: errmsg        !Error message
    integer :: info_loc                                     !Local flag
    character(len=80) :: errmsg_loc                         !Local error message
    real(dp),allocatable,intent(out) :: r(:,:)              !x,y,e data
    type(data),pointer :: head                              !Pointer to head of linked list
    integer :: nvals                                        !Number of data read
    type(data),pointer :: ptr                               !Temporary pointer
    type(data),pointer :: tail                              !pointer to tail of list
    real(dp) :: temp_x,temp_y,temp_e                        !Temporary variables
    integer :: i                                            !Dummy integer
    integer :: unit                                         !File unit

    nvals = 0
    info_loc = 0
    ! Opening and returning for errors
    open(newunit=unit, file=trim(filename), status='old',IOSTAT = info_loc,IOMSG=errmsg_loc)
    if(info_loc /= 0) then
       if(present(info)) info=info_loc
       if(present(errmsg)) errmsg = errmsg_loc
       return
    end if

    ! Deallocating return data
    if(allocated(r)) deallocate(r)
    nullify(head,tail,ptr)

    ! Reading file to linked list then placing in return variable
    do                                                      !Creating linked list
       read(unit,*,iostat = info_loc)temp_x,temp_y,temp_e
       if(info_loc /=0) exit
       nvals=nvals+1

       if(.not. associated(head))then                       !No values in list
          allocate(head,stat = info_loc)                    !Allocates new value
          tail => head                                      !Tail pts to new value
          nullify(tail%next)                                !Nullify ptr in new value
          tail%x = temp_x                                   !Store data
          tail%y = temp_y
          tail%e = temp_e
       else
          allocate(tail%next,stat = info_loc)               !Values already in list
          tail =>tail%next                                  !Tail points to new value
          nullify(tail%next)                                !Nullify ptr in new value
          tail%x = temp_x                                   !Store data
          tail%y = temp_y
          tail%e = temp_e
       end if
    end do

    close(unit)
    allocate(r(nvals,3))
    ptr =>head
    i=1
    do                                                      !Places data in return variable
       if(.not. associated(ptr))exit                        !Exits at null
       r(i,:) = [ptr%x,ptr%y,ptr%e]
       ptr=>ptr%next
       i=i+1
    end do

    do                                                      !Deallocates linked list
       ptr => head%next
       deallocate(head)
       if (.not. associated(ptr)) exit
       head => ptr
    end do

  end subroutine read_xye
  ! ----------------------------------------------------------------------------


  subroutine read_xy(filename,r,info,errmsg)
    ! SUBROUTINE: read_p_r(filename,P(r),[info_integer,error_message])
    ! Reads in .xye file from to produce discrete function data with errors

    character(len=*),intent(in) :: filename                 !Filename
    integer,optional,intent(out) :: info                    !Flag, 0 for success
    character(len=80),optional,intent(out) :: errmsg        !Error message
    integer :: info_loc                                     !Local flag
    character(len=80) :: errmsg_loc                         !Local error message
    real(dp),allocatable,intent(out) :: r(:,:)              !x,y,e data
    type(data),pointer :: head                              !Pointer to head of linked list
    integer :: nvals                                        !Number of data read
    type(data),pointer :: ptr                               !Temporary pointer
    type(data),pointer :: tail                              !pointer to tail of list
    real(dp) :: temp_x,temp_y                               !Temporary variables
    integer :: i                                            !Dummy integer
    integer :: unit                                         !File unit

    nvals = 0
    info_loc = 0
    ! Opening and returning for errors
    open(newunit=unit, file=trim(filename), status='old',IOSTAT = info_loc,IOMSG=errmsg_loc)
    if(info_loc /= 0) then
       if(present(info)) info=info_loc
       if(present(errmsg)) errmsg = errmsg_loc
       return
    end if

    ! Deallocating return data
    if(allocated(r)) deallocate(r)
    nullify(head,tail,ptr)

    ! Reading file to linked list then placing in return variable
    do                                                      !Creating linked list
       read(unit,*,iostat = info_loc)temp_x,temp_y
       if(info_loc /=0) exit
       nvals=nvals+1

       if(.not. associated(head))then                       !No values in list
          allocate(head,stat = info_loc)                    !Allocates new value
          tail => head                                      !Tail pts to new value
          nullify(tail%next)                                !Nullify ptr in new value
          tail%x = temp_x                                   !Store data
          tail%y = temp_y
          tail%e = 0.
       else
          allocate(tail%next,stat = info_loc)               !Values already in list
          tail =>tail%next                                  !Tail points to new value
          nullify(tail%next)                                !Nullify ptr in new value
          tail%x = temp_x                                   !Store data
          tail%y = temp_y
          tail%e = 0.
       end if
    end do

    close(unit)
    allocate(r(nvals,3))
    ptr =>head
    i=1
    do                                                      !Places data in return variable
       if(.not. associated(ptr))exit                        !Exits at null
       r(i,:) = [ptr%x,ptr%y,ptr%e]
       ptr=>ptr%next
       i=i+1
    end do

    do                                                      !Deallocates linked list
       ptr => head%next
       deallocate(head)
       if (.not. associated(ptr)) exit
       head => ptr
    end do

  end subroutine read_xy
  ! ----------------------------------------------------------------------------

  subroutine read_csv(filename,r,info,errmsg)
    ! SUBROUTINE: read_csv(filename,G(r),[info_integer,error_message])
    ! Reads in .csv file from to produce discrete function data with errors
    ! WARNING: Mac OS doesn't always end lines with an eol deliminator. The i/o will
    ! fail if lines do not end with the correct deliminator.

    character(len=*),intent(in) :: filename                 !Filename
    integer,optional,intent(out) :: info                    !Flag, 0 for success
    character(len=80),optional,intent(out) :: errmsg        !Error message
    integer :: info_loc                                     !Local flag
    character(len=80) :: errmsg_loc                         !Local error message
    real(dp),allocatable,intent(out) :: r(:,:)              !x,y,e data
    type(data),pointer :: head                              !Pointer to head of linked list
    integer :: nvals                                        !Number of data read
    type(data),pointer :: ptr                               !Temporary pointer
    type(data),pointer :: tail                              !pointer to tail of list
    real(dp) :: temp_x,temp_y                               !Temporary variables
    integer :: i                                            !Dummy integer
    integer :: unit                                         !File unit

    nvals = 0
    info_loc = 0
    ! Opening and returning for errors
    open(newunit=unit, file=trim(filename), status='old',IOSTAT = info_loc,IOMSG=errmsg_loc)
    if(info_loc /= 0) then
       if(present(info)) info=info_loc
       if(present(errmsg)) errmsg = errmsg_loc
       return
    end if

    ! Deallocating return data
    if(allocated(r)) deallocate(r)
    nullify(head,tail,ptr)

    ! Reading file to linked list then placing in return variable
    do                                                      !Creating linked list
       read(unit,*,iostat = info_loc,IOMSG=errmsg_loc)temp_x,temp_y
       if(info_loc /=0) exit
       nvals=nvals+1

       if(.not. associated(head))then                       !No values in list
          allocate(head,stat = info_loc)                    !Allocates new value
          tail => head                                      !Tail pts to new value
          nullify(tail%next)                                !Nullify ptr in new value
          tail%x = temp_x                                   !Store data
          tail%y = temp_y
          tail%e = 0.
       else
          allocate(tail%next,stat = info_loc)               !Values already in list
          tail =>tail%next                                  !Tail points to new value
          nullify(tail%next)                                !Nullify ptr in new value
          tail%x = temp_x                                   !Store data
          tail%y = temp_y
          tail%e = 0.
       end if
    end do

    close(unit)
    allocate(r(nvals,3))
    ptr =>head
    i=1
    do                                                      !Places data in return variable
       if(.not. associated(ptr))exit                        !Exits at null
       r(i,:) = [ptr%x,ptr%y,ptr%e]
       ptr=>ptr%next
       i=i+1
    end do

    do                                                      !Deallocates linked list
       ptr => head%next
       deallocate(head)
       if (.not. associated(ptr)) exit
       head => ptr
    end do

  end subroutine read_csv
  ! -----------------------------------------------------------------------------§

  subroutine read_dist_list(filename,r,info,errmsg)
    ! SUBROUTINE: read_dist_list(filename,distance_list,[info_integer,error_message])
    ! Reads distance list from the filename instead of reading 2d pdf data

    character(len=*),intent(in) :: filename                 !Filename
    integer,optional,intent(out) :: info                    !Flag, 0 for success
    character(len=80),optional,intent(out) :: errmsg        !Error message
    integer :: info_loc                                     !Local flag
    character(len=80) :: errmsg_loc                         !Local error message
    real(dp),allocatable,intent(out) :: r(:)                !Distance list
    type(data),pointer :: head                              !Pointer to head of linked list
    integer :: nvals                                        !Number of data read
    type(data),pointer :: ptr                               !Temporary pointer
    type(data),pointer :: tail                              !pointer to tail of list
    real(dp) :: temp_x                                      !Temporary variable
    integer :: i                                            !Dummy integer
    integer :: unit                                         !File unit

    nvals = 0
    info_loc = 0
    ! Opening and returning for errors
    open(newunit=unit, file=trim(filename), status='old',IOSTAT = info_loc,IOMSG=errmsg_loc)
    if(info_loc /= 0) then
       if(present(info)) info=info_loc
       if(present(errmsg)) errmsg = errmsg_loc
       return
    end if

    ! Deallocating return data
    if(allocated(r)) deallocate(r)
    nullify(head,tail,ptr)

    ! Reading file to linked list then placing in return variable
    do                                                      !Creating linked list
       read(unit,*,iostat = info_loc)temp_x
       if(info_loc /=0) exit
       nvals=nvals+1

       if(.not. associated(head))then                       !No values in list
          allocate(head,stat = info_loc)                    !Allocates new value
          tail => head                                      !Tail pts to new value
          nullify(tail%next)                                !Nullify ptr in new value
          tail%x = temp_x                                   !Store data
       else
          allocate(tail%next,stat = info_loc)               !Values already in list
          tail =>tail%next                                  !Tail points to new value
          nullify(tail%next)                                !Nullify ptr in new value
          tail%x = temp_x                                   !Store data
       end if
    end do

    close(unit)
    allocate(r(nvals))
    ptr =>head
    i=1
    do                                                      !Places data in return variable
       if(.not. associated(ptr))exit                        !Exits at null
       r(i) = ptr%x
       ptr=>ptr%next
       i=i+1
    end do

    do                                                      !Deallocates linked list
       ptr => head%next
       deallocate(head)
       if (.not. associated(ptr)) exit
       head => ptr
    end do
  end subroutine read_dist_list
  ! -----------------------------------------------------------------------------
  
  subroutine write_xye(filename,dat,info,errmsg)
    ! SUBROUTINE: write_xye(filename, xye_data,[info_integer,error_message])
    ! Writes xye file from data for easier second use
    ! Also writes x,y file if need be
    character(len=*),intent(in) :: filename                 !Filename
    real(dp) :: dat(:,:)                                    !xye data
    integer,optional,intent(out) :: info                    !Flag, 0 for success
    character(len=80),optional,intent(out) :: errmsg        !Error message
    integer :: info_loc                                     !Local flag
    character(len=80) :: errmsg_loc                         !Local error message
    integer :: unit                                         !File unit
    integer :: i                                            !Dummy variable
    info_loc = 0
    open(newunit=unit, file=trim(filename),&
         status='replace',IOSTAT = info_loc,IOMSG=errmsg_loc)
    if(info_loc /= 0) then
       if(present(info)) info=info_loc
       if(present(errmsg)) errmsg = errmsg_loc
       return
    end if

    select case(size(dat,2))
    case(3)
       do i=1,size(dat,1)
          write(unit,'(3f12.5)') dat(i,:)
       end do
    case(2)
       do i=1,size(dat,1)
          write(unit,'(3f12.5)') dat(i,:),0._dp
       end do
    end select
    close(unit)
  end subroutine write_xye
  ! ----------------------------------------------------------------------------

  subroutine write_dist_list(filename,dat,info,errmsg)
    ! SUBROUTINE: write_dist_list(filename, dist_list,[info_integer,error_message])
    ! Writes dist list file from data for tracking

    character(len=*),intent(in) :: filename                 !Filename
    real(dp) :: dat(:)                                      !xye data
    integer,optional,intent(out) :: info                    !Flag, 0 for success
    character(len=80),optional,intent(out) :: errmsg        !Error message
    integer :: info_loc                                     !Local flag
    character(len=80) :: errmsg_loc                         !Local error message
    integer :: unit                                         !File unit
    integer :: i                                            !Dummy variable

    info_loc = 0
    open(newunit=unit, file=trim(filename),&
         status='replace',IOSTAT = info_loc,IOMSG=errmsg_loc)
    if(info_loc /= 0) then
       if(present(info)) info=info_loc
       if(present(errmsg)) errmsg = errmsg_loc
       return
    end if

    do i=1,size(dat)
       write(unit,'(3f12.5)')dat(i)
    end do
    close(unit)
  end subroutine write_dist_list
  ! -----------------------------------------------------------------------------

   subroutine write_thermo(filename,T,E,E2,tau,info,errmsg)
    ! SUBROUTINE: write_dist_list(filename, temperature,<E>, <E^2>,[info_integer,error_message])
    ! Writes thermodynamic data file from data for tracking

    character(len=*),intent(in) :: filename                 !Filename
    real(dp),intent(in) :: T,E,E2                           !xye data
    integer, intent(in) :: tau                              !autocorrelation time
    integer,optional,intent(out) :: info                    !Flag, 0 for success
    character(len=80),optional,intent(out) :: errmsg        !Error message
    integer :: info_loc                                     !Local flag
    character(len=80) :: errmsg_loc                         !Local error message
    logical, save :: FirstCall = .true.                     !Indicator for writing first time
    integer :: unit                                         !File unit

    if(FirstCall) then
       open(newunit=unit,file=trim(filename),status='replace',IOSTAT=info_loc,IOMSG=errmsg_loc)
       if(info_loc /= 0) then
          if(present(info)) info=info_loc
          if(present(errmsg)) errmsg = errmsg_loc
          return
       end if
       write(unit,"('Temperature',T25,'   <E>',T50,' C_v',T70, 'tau')")
       write(unit,"('----------------------------------------------------------------------------')")
       write(unit,"(ES12.3 T25 E12.3 T50 E12.3 T70 I5)") T, E, (E2-E**2)**2/T**2, tau
       FirstCall = .false.
    else
       open(newunit=unit,file=trim(filename),status='old',position='append',IOSTAT=info_loc,IOMSG=errmsg_loc)
       if(info_loc /= 0) then
          if(present(info)) info=info_loc
          if(present(errmsg)) errmsg = errmsg_loc
          return
       end if
       write(unit,"(ES12.3 T25 E12.3 T50 E12.3 T70 I5)") T, E, (E2-E**2)**2/T**2, tau
    end if
    close(unit)
  end subroutine write_thermo
  ! -----------------------------------------------------------------------------

  subroutine segment_data(g_r)
    ! SUBROUTINE: segment_data(G(r))
    ! Doubles the size of the function G(r) by discretizing the midpoints between data
    ! Uses a linear approximation
    real(dp), allocatable, intent(inout) :: g_r(:,:)        !x,y,e data of appropriate function
    real(dp), allocatable :: temp(:,:)                      !Temp holder for old G(r)
    real(dp) :: h                                           !Spacing of data
    real(dp) :: slope                                       !Slope of data
    integer :: n                                            !Number of rows in data
    integer :: m                                            !Number of columns in data
    integer :: i                                            !Dummy integer

    if (allocated(g_r)) then
       n = size(g_r,1)
       m = size(g_r,2)
       allocate(temp(n,m))
       temp = g_r
       deallocate(g_r)
    else
       stop 'Error in doubling data size in exp_io: Unallocated Call'
    end if

    allocate(g_r(2*n-1,m))

    do i=1,size(temp,1)-1
       g_r(2*(i-1)+1,1:2) = temp(i,1:2)
       h = temp(i+1,1)- temp(i,1)
       slope = (temp(i+1,2)- temp(i,2))/(temp(i+1,1)- temp(i,1))
       g_r(2*(i-1)+2,1) = temp(i,1) + h/2.
       g_r(2*(i-1)+2,2) = temp(i,2) + (h/2)*slope
       
       if (m>2) then
          g_r(2*(i-1)+1,3:size(temp,2)) = temp(i,3:size(temp,2))
          g_r(2*(i-1)+2,3:size(temp,2)) = temp(i,3:size(temp,2))
       end if
    end do
    g_r(2*(i-1)+1,1:2) = temp(i,1:2)
  end subroutine segment_data
  ! -----------------------------------------------------------------------------
  
end module exp_io
