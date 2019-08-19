! LAST EDIT: Phil Maffettone 2016-12-13
module lists

  ! The following module is for creating doubly connected linked lists
  ! and the methods associted with them.
  ! Should be consistent with unlimited polymorphism
  ! A redundant dp_list is made for double precision variables, with easier implementation

  use types
  implicit none

  type list_node
     class(*), pointer :: item
     type(list_node), pointer :: next => null()
     type(list_node), pointer :: prev => null()
  end type list_node

  type list
     type(list_node), pointer, private :: first => null()
     type(list_node), pointer, private :: last => null()
     type(list_node), pointer, private :: iter => null()
     integer,private :: length = 0
   contains
     procedure :: finalize => list_finalize
     procedure :: add => list_add
     procedure :: iter_restart => list_iter_restart
     procedure :: iter_next => list_iter_next
     procedure :: for_each => list_for_each
     procedure :: len => list_len
     procedure :: push_front => list_push_front
     procedure :: push_back => list_push_back
     procedure :: pop_front => list_pop_front
     procedure :: pop_back => list_pop_back
     procedure :: display => list_display
  end type list

  type dp_node
     real(dp), pointer :: item
     type(dp_node), pointer :: next => null()
     type(dp_node), pointer :: prev => null()
  end type dp_node

  type dp_list
     type(dp_node), pointer, private :: first => null()
     type(dp_node), pointer, private :: last => null()
     type(dp_node), pointer, private :: iter => null()
     integer,private :: length = 0
   contains
     procedure :: finalize => dp_finalize
     procedure :: add => dp_add
     procedure :: iter_restart => dp_iter_restart
     procedure :: iter_next => dp_iter_next
     procedure :: len => dp_len
     procedure :: push_front => dp_push_front
     procedure :: push_back => dp_push_back
     procedure :: pop_front => dp_pop_front
     procedure :: pop_back => dp_pop_back
     procedure :: display => dp_list_display
  end type dp_list

  abstract interface
     subroutine foreach_sub(item)
       import
       class(*) :: item
     end subroutine foreach_sub
  end interface

contains

  subroutine list_add(self,item)
    class(list), intent(inout) :: self
    class(*), intent(in) :: item

    if (.not.associated(self%last)) then
       allocate(self%first)
       self%last => self%first
    else
       allocate(self%last%next)
       self%last%next%prev => self%last
       self%last => self%last%next
    end if

    allocate(self%last%item,source=item)
    self%length = self%length + 1
  end subroutine list_add
  ! -----------------------------------------------------------------------------

  recursive subroutine list_finalize(self)
    class(list), intent(inout) :: self
    type(list_node), pointer :: node, tmp

    node => self%first

    do while (associated(node))
       tmp => node
       node => node%next
       deallocate(tmp%item)
       deallocate(tmp)
    end do

    self%first => null()
    self%last => null()
    self%length = 0
  end subroutine list_finalize
  ! -----------------------------------------------------------------------------

  subroutine list_iter_restart(self)
    class(list), intent(inout) :: self
    self%iter => self%first
  end subroutine list_iter_restart
  ! -----------------------------------------------------------------------------

  subroutine list_iter_next(self, res)
    class(list), intent(inout) :: self
    class(*), pointer, intent(out) :: res
    if (associated(self%iter)) then
       res => self%iter%item
       self%iter => self%iter%next
    else
       res => null()
    end if
  end subroutine list_iter_next
  ! -----------------------------------------------------------------------------
  
  subroutine list_for_each(self,proc)
    class(list) :: self
    procedure(foreach_sub) :: proc
    type(list_node), pointer :: node

    node => self%first

    do while(associated(node))
       if(associated(node%item)) call proc(node%item)
       node => node%next
    end do
  end subroutine list_for_each
  ! -----------------------------------------------------------------------------

  pure integer function list_len(self)
    class(list),intent(in) :: self
    list_len = self%length
  end function list_len
  ! -----------------------------------------------------------------------------

  subroutine list_push_back(self,item)
    class(list), intent(inout) :: self
    class(*), pointer, intent(in) :: item

    call self%add(item)
  end subroutine list_push_back
  ! -----------------------------------------------------------------------------

  subroutine list_push_front(self,item)
    class(list), intent(inout) :: self
    class(*), pointer, intent(in) :: item
    type(list_node), pointer :: node

     if (associated(self%first)) then
       node => self%first
       allocate(node%prev)
       self%first => node%prev
       self%first%next => node
       allocate(self%first%item,source=item)
    else
       allocate(self%first)
       self%last => self%first
       allocate(self%last%item,source=item)
    end if
    self%length = self%length + 1

  end subroutine list_push_front
  ! -----------------------------------------------------------------------------

  subroutine list_pop_front(self)
    class(list), intent(inout) :: self
    if (associated(self%first)) then
       if (associated(self%first%next)) then
          self%first => self%first%next
          deallocate(self%first%prev)
          self%first%prev => null()
       else
          deallocate(self%first)
          self%first => null()
          self%last => null()
       end if
    end if
  end subroutine list_pop_front
  ! -----------------------------------------------------------------------------

  subroutine list_pop_back(self)
    class(list), intent(inout) :: self

    if (associated(self%last)) then
       self%last => self%last%prev
       deallocate(self%last%next)
       self%last%next => null()
    end if
  end subroutine list_pop_back
  ! -----------------------------------------------------------------------------

  subroutine list_display(self)
    class(list), intent(inout) :: self
    type(list_node), pointer :: node
    class(*), pointer :: assoc 

    node => self%first
    do while (associated(node))
       assoc => node%item
       select type(assoc)
       type is(real(dp))
          write(*,*) assoc
       type is(real(qp))
          write(*,*) assoc
       type is(integer)
          write(*,*) assoc
       type is (character(len=*))
          write(*,*) assoc
       class default
          return
       end select
       node => node%next
    end do
  end subroutine list_display
  ! -----------------------------------------------------------------------------


  ! Begining Double Precision routines
  recursive subroutine dp_finalize(self)
    class(dp_list), intent(inout) :: self
    type(dp_node), pointer :: node, tmp

    node => self%first

    do while (associated(node))
       tmp => node
       node => node%next
       deallocate(tmp%item)
       deallocate(tmp)
    end do

    self%first => null()
    self%last => null()
    self%length = 0
  end subroutine dp_finalize
  ! -----------------------------------------------------------------------------

  subroutine dp_add(self,item)
    class(dp_list), intent(inout) :: self
    real(dp), intent(in) :: item

    if (.not.associated(self%last)) then
       allocate(self%first)
       self%last => self%first
    else
       allocate(self%last%next)
       self%last%next%prev => self%last
       self%last => self%last%next
    end if

    allocate(self%last%item,source=item)
    self%length = self%length + 1
  end subroutine dp_add
  ! -----------------------------------------------------------------------------

  subroutine dp_iter_restart(self)
    class(dp_list), intent(inout) :: self
    self%iter => self%first
  end subroutine dp_iter_restart
  ! -----------------------------------------------------------------------------

  subroutine dp_iter_next(self, res)
    class(dp_list), intent(inout) :: self
    real(dp), intent(out) :: res
    if (associated(self%iter)) then
       res = self%iter%item
       self%iter => self%iter%next
    else
       res = 0.
    end if
  end subroutine dp_iter_next
  ! -----------------------------------------------------------------------------

  pure integer function dp_len(self)
    class(dp_list),intent(in) :: self
    dp_len = self%length
  end function dp_len

  subroutine dp_push_back(self,item)
    class(dp_list), intent(inout) :: self
    real(dp), intent(in) :: item

    call self%add(item)
  end subroutine dp_push_back
  ! -----------------------------------------------------------------------------

  subroutine dp_push_front(self,item)
    class(dp_list), intent(inout) :: self
    real(dp), intent(in) :: item
    type(dp_node), pointer :: node

    if (associated(self%first)) then
       node => self%first
       allocate(node%prev)
       self%first => node%prev
       self%first%next => node
       allocate(self%first%item,source=item)
    else
       allocate(self%first)
       self%last => self%first
       allocate(self%last%item,source=item)
    end if
    self%length = self%length + 1
  end subroutine dp_push_front
  ! -----------------------------------------------------------------------------

  subroutine dp_pop_front(self)
    class(dp_list), intent(inout) :: self
    if (associated(self%first)) then
       if (associated(self%first%next)) then
          self%first => self%first%next
          deallocate(self%first%prev)
          self%first%prev => null()
       else
          deallocate(self%first)
          self%first => null()
          self%last => null()
       end if
       self%length = self%length - 1
    end if
  end subroutine dp_pop_front
  ! -----------------------------------------------------------------------------

  subroutine dp_pop_back(self)
    class(dp_list), intent(inout) :: self

    if (associated(self%last)) then
       self%last => self%last%prev
       deallocate(self%last%next)
       self%last%next => null()
       self%length = self%length - 1
    end if
  end subroutine dp_pop_back
  ! -----------------------------------------------------------------------------

  subroutine dp_list_display(self)
    class(dp_list),intent(in) :: self
    type(dp_node), pointer :: node

    node =>self%first
    do while (associated(node))
       write(*,*) node%item
       node => node%next
    end do
  end subroutine dp_list_display
  ! -----------------------------------------------------------------------------
end module lists
