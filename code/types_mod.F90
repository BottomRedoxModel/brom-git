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

module types_mod
  use fabm_types,only: rk
  use fabm_driver
  use list_mod

  implicit none
  type,abstract:: variable
    character(len=64):: name  = ''
    character(len=64):: long_name  = ''
    character(len=64):: units = ''
  contains
    procedure,non_overridable:: inverse
    procedure,non_overridable:: print_value
    procedure,non_overridable:: print_name
    procedure,non_overridable:: variable_1d_size
  end type

  type,extends(variable):: alone_variable
    real(rk):: value = 0._rk
  end type

  type,extends(variable):: variable_1d
    real(rk),allocatable,dimension(:):: value
  end type

  type,extends(variable):: variable_2d
    real(rk),allocatable,dimension(:,:):: value
  end type

  type,abstract,extends(list):: list_variables
  contains
    procedure:: get_var_by_number
    procedure:: get_var
    procedure:: set_var
    procedure:: inv_var
    procedure:: get_1st_dim_length
    procedure:: print_var
    procedure:: print_list_variables
    procedure:: get_array
    procedure:: get_column
    procedure:: get_value
  end type
contains
  subroutine inverse(self)
    class(variable),intent(inout):: self
    integer:: temp,temp2,i

    select type(self)
    class is(alone_variable)
      return
    class is(variable_1d)
      self%value = self%value(size(self%value):1:-1)
    class is(variable_2d)
      temp = size(self%value,2)
      temp2 = size(self%value,1)
      do i=1,temp
        self%value(:,i) = &
        self%value(temp2:1:-1,i)
      end do
    class default
      call fatal_error("Inverse value","Wrong type")
    end select
  end subroutine

  subroutine print_value(self)
    class(variable),intent(in):: self

    select type(self)
    class is(alone_variable)
      write(*,*) self%value
    class is(variable_1d)
      write(*,'(f15.13)') self%value
    class is(variable_2d)
      write(*,'(f10.5)') self%value(:,1)
    class default
      call fatal_error("Print value","Wrong type")
    end select
  end subroutine

  subroutine print_name(self)
    class(variable),intent(in):: self

    write(*,*) self%name
  end subroutine
  !
  integer function variable_1d_size(self)
    class(variable),intent(in):: self

    select type(self)
    class is(alone_variable)
      variable_1d_size = 1
    class is(variable_1d)
      variable_1d_size = size(self%value,1)
    class is(variable_2d)
      variable_1d_size = size(self%value,1)
    class default
      call fatal_error("Variable size function","Wrong type")
    end select
  end function

  subroutine get_var_by_number(self,number,get_variable)
    class(list_variables),target,intent(in):: self
    integer                     ,intent(in):: number
    class(variable),allocatable ,intent(out):: get_variable
    class(*)             ,pointer:: curr
    class(list_variables),pointer:: temporary
    integer i

    !to make self intent(in)
    i = 1
    temporary => self
    call temporary%reset()
    do
      curr=>temporary%get_item()
      select type(curr)
      class is(variable)
        if (i==number) then
          allocate(get_variable,source=curr)
          return
        end if
      end select
      call temporary%next()
      if (.not.temporary%moreitems()) then
        call fatal_error("Getting variables",&
                         "can't find variable")
      end if
      i = i+1
    end do
  end subroutine

  subroutine get_var(self,inname,get_variable)
    class(list_variables),target,intent(in):: self
    character(len=*)            ,intent(in):: inname
    class(variable),allocatable ,intent(out):: get_variable
    class(*)             ,pointer:: curr
    class(list_variables),pointer:: temporary

    !to make self intent(in)
    temporary => self
    call temporary%reset()
    do
      curr=>temporary%get_item()
      select type(curr)
      class is(variable)
        if (trim(curr%name)==trim(inname)) then
          allocate(get_variable,source=curr)
          return
        end if
      end select
      call temporary%next()
      if (.not.temporary%moreitems()) then
        call fatal_error("Getting variables",&
                         "can't find "//inname//&
                         " variable")
      end if
    end do
  end subroutine

  subroutine set_var(self,inname,new_var)
    class(list_variables),intent(inout):: self
    character(len=*)     ,intent(in)   :: inname
    class(variable),allocatable,intent(in):: new_var
    class(*),pointer:: curr

    call self%reset()
    do
      curr=>self%get_item()
      select type(curr)
      class is(variable)
        if (trim(curr%name)==trim(inname)) then
          call self%set_item(new_var)
          return
        end if
      end select
      call self%next()
      if (.not.self%moreitems()) then
        call fatal_error("Setting variables",&
                         "can't find '"//inname//&
                         "' variable")
      end if
    end do
  end subroutine

  subroutine inv_var(self,inname)
    class(list_variables),intent(inout):: self
    character(len=*)     ,intent(in)   :: inname
    class(*),pointer:: curr

    call self%reset()
    do
      curr=>self%get_item()
      select type(curr)
      class is(variable)
        if (trim(curr%name)==trim(inname)) then
          call curr%inverse()
          return
        end if
      end select
      call self%next()
      if (.not.self%moreitems()) then
        call fatal_error("Inverting variables",&
                         "can't find '"//inname//&
                         "' variable")
      end if
    end do
  end subroutine

  integer function get_1st_dim_length(self,inname)
    class(list_variables),intent(in):: self
    character(len=*)     ,intent(in):: inname
    class(variable),allocatable:: get_variable

    call self%get_var(inname,get_variable)
    select type(get_variable)
    class is(alone_variable)
      get_1st_dim_length=1
    class is(variable_1d)
      get_1st_dim_length=size(get_variable%value,1)
    class is(variable_2d)
      get_1st_dim_length=size(get_variable%value,1)
    end select
  end function

  subroutine print_var(self,inname)
    class(list_variables),intent(in):: self
    character(len=*)     ,intent(in):: inname
    class(variable),allocatable:: var

    call self%get_var(inname,var)
    call var%print_name()
    call var%print_value()
  end subroutine

  subroutine print_list_variables(self,message)
    class(list_variables),intent(inout):: self
    character(len=*)     ,intent(in)   :: message
    class(*),pointer:: curr
    logical first

    call self%reset()
    first = self%moreitems()
    write(*,*) message
    if (.not.first) then
      write(*,*) 'Empty'
    else
      do
        curr=>self%get_item()
        select type(curr)
        class is(variable)
          call curr%print_name()
        end select
        call self%next()
        if (.not.self%moreitems()) exit
      end do
    end if
  end subroutine

  function get_array(self,inname)
    real(rk),allocatable:: get_array(:,:)
    class(list_variables),intent(in):: self
    character(len=*)     ,intent(in):: inname

    class(variable),allocatable:: get_variable

    call self%get_var(inname,get_variable)
    select type(get_variable)
    class is(variable_2d)
      allocate(get_array,source=get_variable%value)
    class default
      call fatal_error("Getting column failed",&
                       "Wrong variable")
    end select
  end function

  function get_column(self,inname,column)
    real(rk),allocatable:: get_column(:)
    class(list_variables),intent(in):: self
    character(len=*)     ,intent(in):: inname
    integer,intent(in),optional     :: column

    class(variable),allocatable:: get_variable

    call self%get_var(inname,get_variable)
    select type(get_variable)
    class is(variable_1d)
      !get_column = get_variable%value(:)
      allocate(get_column,source=get_variable%value(:))
    class is(variable_2d)
      if (.not.present(column)) call fatal_error(&
                       "Getting column failed",&
                       "Column should be present")
      !get_column = get_variable%value(:,column)
      allocate(get_column,source=get_variable%value(:,column))
    class default
      call fatal_error("Getting column failed",&
                       "Wrong variable")
    end select
  end function

  function get_value(self,inname,id)
    real(rk),allocatable:: get_value
    class(list_variables),intent(in):: self
    character(len=*)     ,intent(in):: inname
    integer,optional     ,intent(in):: id

    class(variable),allocatable:: get_variable

    call self%get_var(inname,get_variable)
    select type(get_variable)
    class is(alone_variable)
      !get_value = get_variable%value
      allocate(get_value,source=get_variable%value)
    class is(variable_1d)
      if (.not.present(id)) call fatal_error(&
                       "Getting value failed",&
                       "id should be present")
      allocate(get_value,source=get_variable%value(id))
    class default
      call fatal_error("Getting value failed",&
                       "Wrong variable")
    end select
  end function
end module
