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

module item_mod
  implicit none
  private
  public:: item,delete

  type item
    private
    class(*),pointer:: var => null()
    type(item),pointer:: next => null()
  contains
    procedure,non_overridable:: next_item
    procedure,non_overridable:: get_item
    procedure,non_overridable:: set_item
  end type

  interface item
    module procedure item_constructor
  end interface
contains
  function item_constructor(var,next)
    class(*),intent(in):: var
    class(item),pointer,intent(in):: next
    class(item),pointer:: item_constructor

    allocate(item_constructor)
    allocate(item_constructor%var,source=var)
    item_constructor%next => next
  end function

  function next_item(self)
    class(item),intent(in):: self
    class(item),pointer:: next_item

    next_item => self%next
  end function

  function get_item(self)
    class(item),intent(in):: self
    class(*),pointer:: get_item

    get_item => self%var
  end function

  subroutine set_item(self,new_var)
    class(item),intent(inout):: self
    class(*),intent(in):: new_var

    deallocate(self%var)
    allocate(self%var,source=new_var)
  end subroutine

  !non type-bound procedure because of pointer
  subroutine delete(any_item)
    class(item),pointer,intent(inout):: any_item

    deallocate(any_item)
  end subroutine
end module
