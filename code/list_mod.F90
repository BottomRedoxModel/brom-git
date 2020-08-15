!-----------------------------------------------------------------------
! SPBM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the SPBM distribution.
!-----------------------------------------------------------------------
! Original author(s): Shamil Yakubov
!-----------------------------------------------------------------------

module list_mod
  use item_mod

  implicit none
  private
  public list

  type,abstract:: list
    private
    class(item), pointer:: first_item => null()
    class(item), pointer:: current_item => null()
  contains
    procedure,non_overridable:: add_item
    procedure,non_overridable:: get_item
    procedure,non_overridable:: set_item
    procedure,non_overridable:: next
    procedure,non_overridable:: moreitems
    procedure,non_overridable:: reset
    procedure,non_overridable:: count_list
    procedure,non_overridable:: delete_list
  end type
contains
  subroutine add_item(self, var)
    class(list),intent(inout):: self
    class(*)   ,intent(in)   :: var
    class(item),pointer:: new_item

    if (.not.associated(self%first_item)) then
      self%first_item => item(var, null())
    else
      new_item => item(var,self%first_item)
      self%first_item => new_item
    end if
  end subroutine

  function get_item(self)
    class(list),intent(in):: self
    class(*),pointer:: get_item

    get_item => self%current_item%get_item()
  end function

  subroutine set_item(self,new_var)
    class(list),intent(inout):: self
    class(*),intent(in):: new_var

    call self%current_item%set_item(new_var)
  end subroutine

  subroutine next(self)
    class(list),intent(inout) :: self

    self%current_item => self%current_item%next_item()
  end subroutine

  function moreitems(self)
    class(list),intent(in) :: self
    logical moreitems

    moreitems = associated(self%current_item)
  end function

  subroutine reset(self)
    class(list),intent(inout) :: self

    self%current_item => self%first_item
  end subroutine

  integer function count_list(self)
    class(list),intent(inout):: self

    count_list = 0
    call self%reset()
    do
      count_list = count_list+1
      call self%next()
      if (.not.self%moreitems()) exit
    end do
  end function

  subroutine delete_list(self)
    class(list),intent(inout):: self

    call self%reset()
    do
      if (.not.associated(self%first_item)) exit
      call self%next()
      call delete(self%first_item)
      self%first_item=>self%current_item
    end do
  end subroutine
end module
