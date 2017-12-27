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
  use fabm_driver
  use fabm_types,only: rk
  use list_mod
  use yaml_types
  use yaml,yaml_parse=>parse,yaml_error_length=>error_length

  implicit none
  private

  type:: yaml_ipbm_type
    character(len=64):: name
    real(rk)         :: real_value
    character(len=64):: char_value
  end type

  type,extends(list):: type_yaml_input
  end type

  public read_ipbm_configuration
  public get_ipbm_real_parameter
  public get_ipbm_string_parameter

  type(type_yaml_input):: ipbm_parameters

contains
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
        call fatal_error('read_ipbm_configuration', trim(path_eff)//' &
          must contain a dictionary at the root (non-indented) level, &
          not a single value. Are you missing a trailing colon?')
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
            ipbm_type%char_value = value%string
          else
            ipbm_type%real_value = realvalue
          end if
          call ipbm_parameters%add_item(ipbm_type)
        class is (type_null)
          call fatal_error('parse_initialization',trim(value%path)//' must be set to a real number, not to null.')
        class is (type_dictionary)
          call fatal_error('parse_initialization',trim(value%path)//' must be set to a real number, not to a dictionary.')
        end select
        pair => pair%next
      end do
    end subroutine parse
  end subroutine read_ipbm_configuration

  real(rk) function get_ipbm_real_parameter(name)
    character(len=*):: name
    class(*),pointer:: curr

    call ipbm_parameters%reset()
    do
      curr=>ipbm_parameters%get_item()
      select type(curr)
      class is(yaml_ipbm_type)
        if (trim(curr%name)==trim(name)) then
          get_ipbm_real_parameter = curr%real_value
          return
        end if
      end select
      call ipbm_parameters%next()
      if (.not.ipbm_parameters%moreitems()) then
        call fatal_error("Getting variables",&
                         "can't find "//name//&
                         " variable")
      end if
    end do
  end function get_ipbm_real_parameter

  character(len=64) function get_ipbm_string_parameter(name)
    character(len=*):: name
    class(*),pointer:: curr

    call ipbm_parameters%reset()
    do
      curr=>ipbm_parameters%get_item()
      select type(curr)
      class is(yaml_ipbm_type)
        if (trim(curr%name)==trim(name)) then
          get_ipbm_string_parameter = curr%char_value
          return
        end if
      end select
      call ipbm_parameters%next()
      if (.not.ipbm_parameters%moreitems()) then
        call fatal_error("Getting variables",&
                         "can't find "//name//&
                         " variable")
      end if
    end do
  end function get_ipbm_string_parameter
end module
