!-----------------------------------------------------------------------
! IPBM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the IPBM distribution.
!-----------------------------------------------------------------------
! Original author(s): Shamil Yakubov
!-----------------------------------------------------------------------

module yaml_mod
  use list_mod
  use yaml_types
  use yaml,yaml_parse=>parse,yaml_error_length=>error_length

  implicit none
  private

  type:: yaml_ipbm_type
    character(len=64):: name
    class(*):: value
  end type
  
  type,extends(list):: type_yaml_input
  contains
    private
    procedure:: add_input
  end type

  interface get_ipbm_parameter
    module procedure get_ipbm_parameter_text
    module procedure get_ipbm_parameter_value
  end interface

  public read_ipbm_configuration, get_ipbm_parameter

  class(type_yaml_input):: ipbm_parameter

contains
  !
  !
  !
  subroutine add_input(self, var)
    class(type_yaml_input):: self
    class(variable)     :: var
    class(*),allocatable:: temp

    allocate(temp,source=var)
    call self%add_item(temp)
  end subroutine
  !
  !
  !
  function get_ipbm_parameter_text(name) result(text)
    character(len=*),  intent(in):: name
    character(len=64),intent(out):: text
  
    class(*),pointer:: curr
    
    call ipbm_parameter%reset()
    do
      curr => ipbm_parameter%get_item()
      select type(curr)
      type is(yaml_ipbm_type)
        if (trim(curr%name) == trim(name)) then
          result = curr
          return
        end if
      end select
      call self%next()
      if (.not.self%moreitems()) exit
    end do
    
  end function get_ipbm_parameter_text
  !
  !
  !
  function get_ipbm_parameter_value(name) result(value)
    character(len=*), intent(in):: name
    real(rk)        ,intent(out):: value
  
  end function get_ipbm_parameter_value
  !
  !
  !
  subroutine read_ipbm_configuration()
  !adapted from the FABM code, see at fabm.net: fabm_config.F90

    class (type_node),pointer        :: node
    character(len=yaml_error_length) :: yaml_error
    character(len=256)               :: path_eff

    path_eff = 'ipbm.yaml'
    ! Parse YAML file.
    node => yaml_parse(trim(path_eff),1,yaml_error)
    if (yaml_error/='') call fatal_error('read ipbm configuration: ',trim(yaml_error))
    if (.not.associated(node)) call fatal_error('read ipbm configuration: ', &
    'No configuration information found in '//trim(path_eff)//'.')
    !call node%dump(output_unit,0)

    select type (node)
      class is (type_dictionary)
        call parse(node)
      class is (type_node)
        call fatal_error('read_ipbm_configuration', trim(path_eff)//' must contain a dictionary &
      &at the root (non-indented) level, not a single value. Are you missing a trailing colon?')
    end select
  contains
    !
    !
    !
    subroutine parse(mapping)
      class (type_dictionary), intent(in) :: mapping
      
      class (type_node), pointer          :: node
      type  (type_key_value_pair), pointer:: pair
      type  (yaml_ipbm_type)              :: ipbm_type
        
      logical             :: success
      real(rk)            :: realvalue
      character(len=64)   :: initname

      node => mapping%get('initialization')
      if (.not.associated(node)) &
        call fatal_error('read_ipbm_configuration', 'No "initialization" dictionary found at root level.')
      select type (node)
      class is (type_dictionary)
        pair => node%first
      class is (type_node)
        nullify(pair)
        call fatal_error('read_ipbm_configuration',trim(node%path)// &
            ' must be a dictionary (initialization:)')
      end select

      do while (associated(pair))
        select type (value=>pair%value)
        class is (type_scalar)
          initname = trim(pair%key)
          ipbm_type%name = initname
          realvalue = value%to_real(default=real(0,real_kind),success=success)
          if (.not.success) then !Assume the value is a string
            !realvalue = huge(1.0_rk)
            ipbm_type%value = value%string
            call ipbm_parameter%add_input(ipbm_type)
          else
            !outname = repeat('z',64)
            ipbm_type%value = realvalue
            call ipbm_parameter%add_input(realvalue)
          end if
        class is (type_null)
          call fatal_error('parse_initialization',trim(value%path)//' must be set to a real number, not to null.')
        class is (type_dictionary)
          call fatal_error('parse_initialization',trim(value%path)//' must be set to a real number, not to a dictionary.')
        end select
      end do
    end subroutine parse
  end subroutine read_ipbm_configuration
  
end module
