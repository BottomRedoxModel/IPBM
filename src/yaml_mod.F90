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

#include "../include/spbm.h"

module yaml_mod
  use fabm_driver
  use fabm_types,only: rk
  use list_mod
  use yaml_types
  use yaml,yaml_parse=>parse,yaml_error_length=>error_length

  implicit none
  private

  type:: yaml_spbm_type
    character(len=64):: name
    real(rk)         :: real_value
    character(len=64):: char_value
  end type

  type,extends(list):: type_yaml_input
  end type

  !interface get_spbm_parameter
  !  module procedure get_spbm_real_parameter
  !  module procedure get_spbm_string_parameter
  !end interface

  public read_spbm_configuration
  public get_spbm_real_parameter
  public get_spbm_string_parameter

  type(type_yaml_input):: spbm_parameters

contains
  !
  !
  !
  subroutine read_spbm_configuration()
  !adapted from the FABM code, see at fabm.net: fabm_config.F90

    class (type_node),pointer        :: node
    character(len=yaml_error_length) :: yaml_error
    character(len=256)               :: path_eff

    path_eff = _SPBM_FILE_NAME_
    ! Parse YAML file.
    node => yaml_parse(trim(path_eff),1,yaml_error)
    if (yaml_error/='') call fatal_error('read spbm configuration: ',trim(yaml_error))
    if (.not.associated(node)) call fatal_error('read spbm configuration: ', &
    'No configuration information found in '//trim(path_eff)//'.')
    !call node%dump(output_unit,0)

    select type (node)
      class is (type_dictionary)
        call parse(node)
      class is (type_node)
        call fatal_error('read_spbm_configuration', trim(path_eff)//' &
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
      type  (yaml_spbm_type)              :: spbm_type

      logical             :: success
      real(rk)            :: realvalue
      character(len=64)   :: initname

      node => mapping%get('initialization')
      if (.not.associated(node)) &
        call fatal_error('read_spbm_configuration', 'No "initialization" dictionary found at root level.')
      select type (node)
      class is (type_dictionary)
        pair => node%first
      class is (type_node)
        nullify(pair)
        call fatal_error('read_spbm_configuration',trim(node%path)// &
            ' must be a dictionary (initialization:)')
      end select

      do while (associated(pair))
        select type (value=>pair%value)
        class is (type_scalar)
          initname = trim(pair%key)
          spbm_type%name = initname
          realvalue = value%to_real(default=real(0,real_kind),success=success)
          if (.not.success) then !Assume the value is a string
            spbm_type%char_value = value%string
          else
            spbm_type%real_value = realvalue
          end if
          call spbm_parameters%add_item(spbm_type)
        class is (type_null)
          call fatal_error('parse_initialization',trim(value%path)//' must be set to a real number, not to null.')
        class is (type_dictionary)
          call fatal_error('parse_initialization',trim(value%path)//' must be set to a real number, not to a dictionary.')
        end select
        pair => pair%next
      end do
    end subroutine parse
  end subroutine read_spbm_configuration

  real(rk) function get_spbm_real_parameter(name)
    character(len=*):: name
    class(*),pointer:: curr

    call spbm_parameters%reset()
    do
      curr=>spbm_parameters%get_item()
      select type(curr)
      class is(yaml_spbm_type)
        if (trim(curr%name)==trim(name)) then
          get_spbm_real_parameter = curr%real_value
          return
        end if
      end select
      call spbm_parameters%next()
      if (.not.spbm_parameters%moreitems()) then
        call fatal_error("Getting variables",&
                         "can't find "//name//&
                         " variable")
      end if
    end do
  end function get_spbm_real_parameter

  character(len=64) function get_spbm_string_parameter(name)
    character(len=*):: name
    class(*),pointer:: curr

    call spbm_parameters%reset()
    do
      curr=>spbm_parameters%get_item()
      select type(curr)
      class is(yaml_spbm_type)
        if (trim(curr%name)==trim(name)) then
          get_spbm_string_parameter = curr%char_value
          return
        end if
      end select
      call spbm_parameters%next()
      if (.not.spbm_parameters%moreitems()) then
        call fatal_error("Getting variables",&
                         "can't find "//name//&
                         " variable")
      end if
    end do
  end function get_spbm_string_parameter
end module
