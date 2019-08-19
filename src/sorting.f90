module sorting
  ! A routine for heap sorting a list of double precision numbers
  ! From NUMERICAL RECIPIES IN FORTRAN 77

  use types
  implicit none
  public
  
contains
  ! --------------------Public Contents-----------------------------------------
  ! SUBROUTINE: hpsort(real_array)
  ! ----------------------------------------------------------------------------
  ! --------------------Private Contents----------------------------------------
  ! 
  ! ----------------------------------------------------------------------------
  subroutine hpsort(ra)
    ! Sorts an array ra(1:n) into ascending numerical order using the
    ! heapsort algorithm
    real(dp), intent(inout) :: ra(:)                        !Array to sort
    integer :: n                                            !Size of array
    integer :: i, ir, j, l                                  !Dummy integers
    real(dp) :: rra                                         !Selection element

    n=size(ra)
    if (n .lt. 2) return
    ! The index l will be decremented from its initial value down to 1
    ! during the hiring (heap creation) phase. Once it reaches 1, the index ir
    !  will be decremented from its initial value down to 1 during the heap selection

    l = n/2 + 1
    ir = n
10  continue
    if (l .gt. 1) then                                      !Still in hiring phase
       l = l - 1
       rra = ra(l)
    else                                                    !In retiremen/promotion phase
       rra = ra(ir)                                         !Clear a space at the end of array
       ra(ir) = ra(1)                                       !Retire the top of the heap into it
       ir = ir - 1                                          !Decrease size of heap
       if (ir .eq. 1) then                                  !Done with last promotion
          ra(1) = rra                                       !The lowest value
          return
       end if
    end if

    ! In the hiring or promotion phase, we set up to sift down element rra
    ! to its proper level
    i = l
    j = l+1
20  if(j .le. ir) then
       if (j .lt. ir) then
          if(ra(j) .lt. ra(j+1)) j=j+1                      !Compare to better underling
       end if
       if (rra .lt. ra(j)) then                             !Demote rra
          ra(i) = ra(j)
          i = j
          j = j+j
       else                                              !This is rra's level. Set j to terminate sift-down
          j = ir+1
       end if
       goto 20
    end if
    ra(i) = rra                                             !Put rra into its slot
    goto 10
  end subroutine hpsort
end module sorting
