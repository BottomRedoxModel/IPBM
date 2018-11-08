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
#include "../include/parameters.h"

module variables_mod
  use types_mod
  use input_mod
  use ice_mod
  use fabm_driver
  use yaml_mod !to use a function mentioned in spbm.h

  implicit none
  !NaN value
  !REAL(rk), PARAMETER :: D_QNAN = &
  !          TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_rk)
  real(rk) D_QNAN

  type,extends(list_variables):: spbm_standard_variables
    type(ice) type_ice
  contains
    private
    procedure:: initialize=>initialize_standard_variables
    procedure:: add_var=>add_standard_var
    procedure:: add_ice_thickness
    procedure:: add_grid_on_faces
    procedure:: add_grid_on_centers
    procedure:: add_day_number
    procedure:: add_layer_thicknesses
    procedure:: add_constant_in_sed
    procedure:: add_porosity
    procedure:: add_diffusivity
    procedure,public:: first_day
  end type

  type,extends(variable_1d):: spbm_state_variable
    logical  is_solid
    logical  is_gas
    integer  use_bound_up
    integer  use_bound_low
    real(rk) bound_up
    real(rk) bound_low
    real(rk) density
    real(rk),allocatable,dimension(:):: sinking_velocity
    real(rk),allocatable,dimension(:):: fickian_fluxes
  contains
    procedure:: set_spbm_state_variable
    procedure:: print_state_variable
  end type

  interface spbm_standard_variables
    module procedure spbm_standard_variables_constructor
  end interface
contains
  function spbm_standard_variables_constructor()
    type(spbm_standard_variables):: spbm_standard_variables_constructor

    !NaN
    D_QNAN = 0._rk
    D_QNAN = D_QNAN / D_QNAN

    call spbm_standard_variables_constructor%initialize()
  end function spbm_standard_variables_constructor
  !
  !Initialize standard variables list
  !
  subroutine initialize_standard_variables(self)
    class(spbm_standard_variables),intent(inout):: self
    type(type_input):: kara_input

    !open 'input_file_name'nc from spbm.yaml, make list with all variables
    kara_input = type_input(_FILE_NAME_)
    !horizontal variables
    call self%add_var(kara_input,_OCEAN_TIME_)
    call self%add_var(kara_input,_SHORTWAVE_RADIATION_)
    call self%add_day_number("day_number")
    !ice variables
    call self%add_ice_thickness(kara_input)
    call self%add_var(kara_input,_SNOW_THICKNESS_)
    call self%add_var(kara_input,_ICE_SURFACE_TEMPERATURE_)
    self%type_ice = ice(self%get_1st_dim_length("day_number"),&
      self%get_column(_ICE_THICKNESS_))
    !vertical variables
    call self%add_grid_on_faces(kara_input,&
      self%type_ice%get_number_of_layers())
    call self%add_grid_on_centers("middle_layer_depths","dz",&
                                  "number_of_layers","air_ice_indexes")
    call self%add_layer_thicknesses("layer_thicknesses")
    !2d variables
    !Add variables which are constants in sediments
    call self%add_constant_in_sed(kara_input,_TEMPERATURE_)
    call self%add_constant_in_sed(kara_input,_SALINITY_)
    call self%add_constant_in_sed(kara_input,_RHO_)
    !other variables
    call self%add_porosity()
    call self%add_diffusivity(kara_input)
    call self%print_list_variables('Allocated spbm_standard_variables:')
    !call self%print_var("porosity")
    !call self%print_var("porosity_on_interfaces")
    !delete unneeded list
    call kara_input%delete_list()
  end subroutine initialize_standard_variables

  subroutine add_standard_var(self,name_input,inname)
    class(spbm_standard_variables),intent(inout):: self
    type(type_input)              ,intent(in)   :: name_input
    character(len=*)              ,intent(in)   :: inname
    class(variable),allocatable:: var
    class(*),allocatable:: temp

    call name_input%get_var(inname,var)
    !memory allocation problems occur without it
    allocate(temp,source=var)
    call self%add_item(temp)
  end subroutine add_standard_var
  !
  !Add discretized ice_thicknesses
  !
  subroutine add_ice_thickness(self,name_input)
    class(spbm_standard_variables),intent(inout):: self
    type(type_input)              ,intent(in)   :: name_input
    class(variable),allocatable:: var
    type(variable_1d) new_var_1d
    real(rk),dimension(:),allocatable:: value_1d

    call name_input%get_var(_ICE_THICKNESS_,var)
    select type(var)
    class is(variable_1d)
      allocate(value_1d(size(var%value,1)))
      value_1d = var%value-mod(var%value,_ICE_LAYERS_RESOLUTION_)
      new_var_1d = variable_1d(var%name,var%long_name,var%units,value_1d)
      call self%add_item(new_var_1d)
    end select
  end subroutine add_ice_thickness
  !
  !Adds bbl and sediments to depths of layers faces
  !
  subroutine add_grid_on_faces(self,name_input,ice_layers)
    class(spbm_standard_variables),intent(inout):: self
    type(type_input),intent(in):: name_input
    integer         ,intent(in):: ice_layers

    class(variable)        ,allocatable:: var

    real(rk),dimension(:)  ,allocatable:: ice_thickness
    real(rk),dimension(:,:),allocatable:: value_2d

    type(alone_variable) new_var
    type(variable_2d) new_var_2d

    integer length,bbl_count,sediments_count
    integer ice_water,water_bbl,bbl_sediments,total_boundaries
    integer time
    integer i

    real(rk):: width_bbl
    real(rk):: resolution_bbl
    real(rk):: width_sediments
    real(rk):: resolution_sediments

    width_bbl = _WIDTH_BBL_
    resolution_bbl = _RESOLUTION_BBL_
    width_sediments = _WIDTH_SEDIMENTS_
    resolution_sediments = _RESOLUTION_SEDIMENTS_

    call name_input%get_var(_DEPTH_ON_BOUNDARY_,var)
    select type(var)
    class is(variable_1d)
      length=size(var%value,1)
    end select
    time = self%get_1st_dim_length("day_number")
    allocate(ice_thickness(time))
    ice_thickness = self%get_column(_ICE_THICKNESS_)

    bbl_count = width_bbl/resolution_bbl
    sediments_count = width_sediments/resolution_sediments
    allocate(value_2d(ice_layers+length+bbl_count+sediments_count,time))
    water_bbl = 1+bbl_count+sediments_count
    bbl_sediments = 1+sediments_count
    ice_water = length+bbl_count+sediments_count
    total_boundaries = ice_layers+length+bbl_count+sediments_count

    !adding indexes of inner boundaries
    new_var = alone_variable("ice_water_index","ice_water_index",'Layer',ice_water)
    call self%add_item(new_var)
    new_var = alone_variable("water_bbl_index","water_bbl_index",'Layer',water_bbl)
    call self%add_item(new_var)
    new_var = alone_variable("bbl_sediments_index","bbl_sediments_index",'Layer',bbl_sediments)
    call self%add_item(new_var)
    new_var = alone_variable("number_of_boundaries","number_of_boundaries",'Layer',total_boundaries)
    call self%add_item(new_var)

    select type(var)
    class is(variable_1d)
      value_2d = 0._rk
      forall (i=1:time) value_2d(water_bbl:ice_water,i) = var%value
      value_2d(ice_water+1:,:) = self%type_ice%do_grid(ice_thickness)
      value_2d(bbl_sediments,:) = var%value(1)
      do i = bbl_sediments+1,water_bbl
        value_2d(i,:) = value_2d(i-1,:)-resolution_bbl
      end do
      if (value_2d(water_bbl,1)<=value_2d(water_bbl+1,1)) then
        call fatal_error("BBL configurating","Wrong _BBL_WIDTH_")
      end if
      value_2d(1,:) = value_2d(bbl_sediments,:)+width_sediments
      do i = 2,(bbl_sediments-1)
        value_2d(i,:) = value_2d(i-1,:)-resolution_sediments
      end do
      new_var_2d = variable_2d(var%name,var%long_name,var%units,value_2d)
      call self%add_item(new_var_2d)
    class default
      call fatal_error("Adding layers","Wrong type")
    end select
  end subroutine add_grid_on_faces

  subroutine add_grid_on_centers(self,inname,inname_increments,&
                                 number_of_layers,air_ice_indexes)
    class(spbm_standard_variables),intent(inout):: self

    character(len=*),intent(in):: inname
    character(len=*),intent(in):: inname_increments
    character(len=*),intent(in):: number_of_layers
    character(len=*),intent(in):: air_ice_indexes

    class(variable)        ,allocatable:: var
    real(rk),dimension(:,:),allocatable:: mid_depths
    !increments between layer midpoints
    real(rk),dimension(:,:),allocatable:: dz

    type(alone_variable):: new_var
    type(variable_1d):: new_var_1d
    type(variable_2d):: new_var_2d

    integer i,j,length,time
    integer ice_water_i
    integer,dimension(:),allocatable:: active_layers

    time = self%get_1st_dim_length("day_number")
    allocate(active_layers(time))
    ice_water_i = self%get_value("ice_water_index")
    active_layers = ice_water_i+self%type_ice%get_active_layers()
    new_var_1d = variable_1d(air_ice_indexes,air_ice_indexes,'Layers',active_layers)
    call self%add_item(new_var_1d)

    length = self%get_value("number_of_boundaries")-1._rk
    new_var = alone_variable(number_of_layers,number_of_layers,'Layer',length)
    call self%add_item(new_var)

    call self%get_var(_DEPTH_ON_BOUNDARY_,var)
    select type(var)
    class is(variable_2d)
      allocate(mid_depths(length,time))
      mid_depths = D_QNAN
      do j = 1,time
        do i = 1,(active_layers(j)-1)
          mid_depths(i,j) = abs((var%value(i+1,j)+&
                          var%value(i,j))/2._rk)
        end do
      end do
      new_var_2d = variable_2d(inname,inname,var%units,mid_depths)
      call self%add_item(new_var_2d)
      !for example: for 4 boundaries exists 3 layers
      !and 3 midpoints of these layers and 2 dz
      allocate(dz(length-1,time))
      dz = D_QNAN
      do j = 1,time
        do i = 2,(active_layers(j)-1)
          dz(i-1,j) = abs(mid_depths(i,j)-mid_depths(i-1,j))
        end do
      end do
      new_var_2d = variable_2d(inname_increments,inname_increments,var%units,dz)
      call self%add_item(new_var_2d)
    class default
      call fatal_error("Add grid on centers",&
        "Wrong type")
    end select
  end subroutine add_grid_on_centers

  subroutine add_day_number(self,inname)
    class(spbm_standard_variables),intent(inout):: self
    character(len=*)              ,intent(in)   :: inname
    class(variable),allocatable:: var

    call self%get_var(_OCEAN_TIME_,var)
    select type(var)
    class is(variable_1d)
      var%value = var%value/86400._rk
      var%name = inname
      call self%add_item(var)
    class default
      call fatal_error("Add day number",&
        "Wrong type")
    end select
  end subroutine add_day_number

  integer function first_day(self)
    class(spbm_standard_variables),intent(in):: self
    class(variable),allocatable:: var

    call self%get_var("day_number",var)
    select type(var)
    class is(variable_1d)
      first_day = int(var%value(1))
    end select
  end function first_day

  subroutine add_layer_thicknesses(self,inname)
    class(spbm_standard_variables),intent(inout):: self
    character(len=*)              ,intent(in)   :: inname
    class(variable)        ,allocatable:: var
    integer ,dimension(:)  ,allocatable:: air_ice_indexes
    real(rk),dimension(:,:),allocatable:: value_2d
    type(variable_2d):: new_var
    integer i,j,length,time

    call self%get_var(_DEPTH_ON_BOUNDARY_,var)
    length = self%get_value("number_of_layers")
    time = self%get_1st_dim_length("day_number")
    allocate(value_2d(length,time))
    value_2d = D_QNAN
    allocate(air_ice_indexes(time))
    air_ice_indexes = int(self%get_column("air_ice_indexes"))

    select type(var)
    class is(variable_2d)
      !forall(i = 1:length)&
      !  value_2d(i,:) = abs(var%value(i+1,:)-var%value(i,:))
      do j = 1,time
        do i = 1,(air_ice_indexes(j)-1)
          value_2d(i,j) = abs(var%value(i+1,j)-&
                          var%value(i,j))
        end do
      end do
      new_var = variable_2d(inname,inname,var%units,value_2d)
      call self%add_item(new_var)
    class default
      call fatal_error("Add layer_thicknesses",&
        "Wrong type")
    end select
  end subroutine add_layer_thicknesses

  subroutine add_constant_in_sed(self,name_input,inname)
    class(spbm_standard_variables),intent(inout):: self
    type(type_input)              ,intent(in)   :: name_input
    character(len=*)              ,intent(in)   :: inname
    class(variable),allocatable        :: var
    real(rk),dimension(:,:),allocatable:: value_2d
    real(rk),dimension(:,:),allocatable:: salinity
    real(rk),dimension(:)  ,allocatable:: air_temp,ice_thickness
    type(variable_2d):: new_var
    integer i,length,time
    integer water_bbl_index,ice_water_index
    character(64):: temperature_name, salinity_name, density_name

    temperature_name = trim(_TEMPERATURE_)
    salinity_name = trim(_SALINITY_)
    density_name = trim(_RHO_)

    length          = self%get_value("number_of_layers")
    water_bbl_index = self%get_value("water_bbl_index")
    ice_water_index = self%get_value("ice_water_index")
    time            = self%get_1st_dim_length("day_number")
    allocate(value_2d(length,time))
    value_2d = 0._rk

    call name_input%get_var(inname,var)
    select type(var)
    class is(variable_2d)
      value_2d(water_bbl_index:ice_water_index-1,:time) = var%value
      if (inname == temperature_name) then
        allocate(air_temp(time))
        allocate(ice_thickness(time))
        air_temp      = self%get_column(_ICE_SURFACE_TEMPERATURE_)
        ice_thickness = self%get_column(_ICE_THICKNESS_)
        value_2d(ice_water_index:,:time) = self%type_ice%do_ice_temperature(&
          air_temp,value_2d(ice_water_index-1,:time),ice_thickness)
      else if (inname == salinity_name) then
        value_2d(ice_water_index:,:time) = self%type_ice%do_ice_brine_salinity(&
          value_2d(ice_water_index-1,:time))
      else if (inname == density_name) then
        allocate(salinity(length,time))
        salinity = self%get_array(salinity_name)
        value_2d(ice_water_index:,:time) = self%type_ice%do_ice_brine_density(&
          salinity(ice_water_index:,:time))/1000._rk
      else
        forall (i = 1:time)&
          value_2d(ice_water_index:,i) = &
          value_2d(ice_water_index-1,i)
      end if
      forall (i = 1:time)&
        value_2d(:water_bbl_index-1,i) =&
        value_2d(water_bbl_index,i)
      new_var = variable_2d(inname,var%long_name,var%units,value_2d)
      call self%add_item(new_var)
    class default
      call fatal_error("Adding constant in sediments variable",&
                       "Wrong type")
    end select
  end subroutine add_constant_in_sed
  !
  !adopted from Phil Wallhead (PW)
  !Adds porosity,tortuosity and porosity factors, PW:
  !"These allow us to use a single equation to model diffusivity
  !updates in the water column and sediments, for both solutes and solids:
  !dC/dt = d/dz(pF2*kzti*d/dz(pF1*C))
  !where C has units [mass per unit total volume (water+sediments)]"
  !
  subroutine add_porosity(self)
    !
    ! equations for solutes, PW:
    ! dC/dt = d/dz(kzti*dC/dz)            in the water column
    ! dC/dt = d/dz(phi*kzti*d/dz(C/phi))  in the sediments
    ! equations for solids, PW:
    ! dC/dt = d/dz(kzti*dC/dz)                   in the water column
    ! dC/dt = d/dz((1-phi)*kzti*d/dz(C/(1-phi))) in the sediments
    !
    class(spbm_standard_variables),intent(inout):: self

    real(rk),dimension(:)  ,allocatable:: w_b,u_b
    real(rk),dimension(:,:),allocatable:: porosity
    real(rk),dimension(:,:),allocatable:: tortuosity
    real(rk),dimension(:,:),allocatable:: porosity_factor
    type(alone_variable):: dphidz_SWI
    type(variable_1d):: var_1d
    type(variable_2d):: new_var

    integer ice_water_index,swi_index,length,time,i
    real(rk) max_porosity,min_porosity,porosity_decay

    integer ,dimension(:)  ,allocatable:: air_ice_indexes
    real(rk),dimension(:)  ,allocatable:: swi_depth
    real(rk),dimension(:,:),allocatable:: depth_center
    real(rk),dimension(:,:),allocatable:: depth_boundary

    max_porosity   = _MAX_POROSITY_
    min_porosity   = _MIN_POROSITY_
    porosity_decay = _POROSITY_DECAY_

    ice_water_index = self%get_value("ice_water_index")
    swi_index       = self%get_value("bbl_sediments_index")
    length          = self%get_value("number_of_layers")
    time            = self%get_1st_dim_length("day_number")
    allocate(air_ice_indexes(time))
    air_ice_indexes = int(self%get_column("air_ice_indexes"))
    allocate(depth_center(length, time))
    depth_center = self%get_array("middle_layer_depths")
    !allocate(depth_center,source=self%get_array("middle_layer_depths"))
    allocate(depth_boundary(length+1, time))
    depth_boundary = self%get_array(_DEPTH_ON_BOUNDARY_)
    !allocate(depth_boundary,source=self%get_array(_DEPTH_ON_BOUNDARY_))
    allocate(swi_depth,source = depth_boundary(swi_index,:))

    !Indices of layer interfaces in the sediments
    var_1d = variable_1d("k_sed","Indices of layer interfaces",&
                         "",(/(i,i=1,swi_index-1)/))
    call self%add_item(var_1d)
    !Indices of layer interfaces in the sediments including the SWI
    var_1d = variable_1d("k_sed1","Indices of layer interfaces with SWI",&
                         "",(/(i,i=1,swi_index)/))
    call self%add_item(var_1d)

    !alone variable
    !dphi/dz on SWI for sedimentation velocity caclulations
    dphidz_SWI = alone_variable("dphidz_SWI",&
                 "dphi/dz on SWI for sedimentation velocity caclulations",&
                 "m^-1",1.0_rk*(max_porosity-min_porosity)/porosity_decay)
    call self%add_item(dphidz_SWI)

    !for layers
    allocate(porosity(length,time))
    porosity = 1._rk
    do i = 1,time
      porosity(1:swi_index-1,i) = min_porosity+(&
        max_porosity-min_porosity)*exp(-1._rk*(&
        depth_center(1:swi_index-1,i)-swi_depth(i))/&
        porosity_decay)
    end do
    porosity(ice_water_index:,:) = self%type_ice%do_brine_relative_volume(&
          .true.,self%get_column(_ICE_THICKNESS_))
    new_var = variable_2d("porosity","porosity",'-',porosity)
    call self%add_item(new_var)

    !porosity factor 1 for solutes
    allocate(porosity_factor(length,time))
    porosity_factor = 1._rk/porosity
    !PW:
    !Factor to convert [mass per unit total volume]
    !to [mass per unit volume pore water] for solutes in sediments
    new_var = variable_2d("porosity_factor_solutes_1",&
                          "Total vol. to pore water vol.",&
                          '-',porosity_factor)
    call self%add_item(new_var)

    !porosity factor 1 for solids
    porosity_factor = 1._rk
    do i = 1,time
      porosity_factor(air_ice_indexes(i):,i) = D_QNAN
    end do
    porosity_factor(1:swi_index-1,:) = 1._rk/&
      (1._rk-porosity(1:swi_index-1,:))
    !PW:
    !Factor to convert [mass per unit total volume]
    !to [mass per unit volume solids] for solids in sediments
    new_var = variable_2d("porosity_factor_solids_1",&
                          "Total vol. to vol. solids",&
                          '-',porosity_factor)
    call self%add_item(new_var)

    deallocate(porosity)
    deallocate(porosity_factor)

    !for boundaries
    allocate(porosity(length+1,time))

    !background vertical advective velocities of particulates and
    !solutes on layer interfaces in the sediments (w_b, u_b)
    !these assume steady state compaction and neglect reaction terms
    allocate(w_b(swi_index))
    allocate(u_b(swi_index))

    porosity = 1._rk
    do i = 1,self%get_1st_dim_length("day_number")
      porosity(1:swi_index,i) = min_porosity+(&
        max_porosity-min_porosity)*exp(-1._rk*(&
        depth_boundary(1:swi_index,i)-swi_depth(i))/&
        porosity_decay)
      if (ice_water_index==air_ice_indexes(i)) then
        porosity(ice_water_index,i) = 1._rk
      else
        porosity(ice_water_index,i) = 0.5_rk
      end if
    end do
    porosity(ice_water_index+1:,:) = &
          self%type_ice%do_brine_relative_volume(&
          .false.,self%get_column(_ICE_THICKNESS_))
    new_var = variable_2d("porosity_on_interfaces","porosity_on_interfaces",'-',porosity)
    call self%add_item(new_var)
    w_b = 0._rk
    u_b = 0._rk
    w_b(1:swi_index) = -((1._rk-min_porosity)/&
                       (1._rk-porosity(1:swi_index,1)))*_BURIAL_VELOCITY_
    u_b(1:swi_index) = -(min_porosity)/&
                       (porosity(1:swi_index,1))*_BURIAL_VELOCITY_
    var_1d = variable_1d("w_b","Particulates velocities","m s^-1",w_b)
    call self%add_item(var_1d)
    var_1d = variable_1d("u_b","Solutes velocities","m s^-1",u_b)
    call self%add_item(var_1d)

    !porosity factor 2 for solutes
    allocate(porosity_factor(length+1,time))
    porosity_factor = porosity
    !PW:
    !Porosity-related area restriction factor for fluxes
    !across layer interfaces
    new_var = variable_2d("porosity_factor_solutes_2",&
                          "Area restriction factor - solutes",&
                          '-',porosity_factor)
    call self%add_item(new_var)

    !porosity factor 2 for solids
    porosity_factor = 1._rk
    do i = 1,time
      porosity_factor(air_ice_indexes(i)+1:,i) = D_QNAN
    end do
    porosity_factor(1:swi_index,:) = &
      1._rk-porosity(1:swi_index,:)
    !PW:
    !Porosity-related area restriction factor for fluxes
    !across layer interfaces
    new_var = variable_2d("porosity_factor_solids_2",&
                          "Area restriction factor - solids",&
                          '-',porosity_factor)
    call self%add_item(new_var)

    !tortuosity on layer interfaces
    !Boudreau 1996, eq. 4.120
    allocate(tortuosity(length+1,time))
    tortuosity = sqrt(1._rk-2._rk*log(porosity))
    new_var = variable_2d("tortuosity_on_interfaces",&
                          "Tortousity on layer interfaces",&
                          '-',tortuosity)
    call self%add_item(new_var)
  end subroutine add_porosity

  subroutine add_diffusivity(self,name_input)
    class(spbm_standard_variables),intent(inout):: self
    type(type_input),intent(in):: name_input

    class(variable)        ,allocatable:: var

    real(rk),dimension(:,:),allocatable:: eddy_kz
    real(rk),dimension(:,:),allocatable:: value_2d
    real(rk),dimension(:,:),allocatable:: tortuosity
    real(rk),dimension(:,:),allocatable:: porosity
    real(rk),dimension(:,:),allocatable:: depth_boundary
    integer ,dimension(:)  ,allocatable:: air_ice_indexes

    type(variable_2d):: new_var_2d

    real    z_conv
    integer i,j,time
    integer number_of_boundaries
    integer water_bbl_index,ice_water_index
    integer bbl_sediments_index

    number_of_boundaries &
      = self%get_value("number_of_boundaries")
    water_bbl_index = self%get_value("water_bbl_index")
    ice_water_index = self%get_value("ice_water_index")
    bbl_sediments_index &
      = self%get_value("bbl_sediments_index")
    time = self%get_1st_dim_length("day_number")
    allocate(air_ice_indexes(time))
    air_ice_indexes = int(self%get_column("air_ice_indexes"))
    allocate(eddy_kz(number_of_boundaries,time))
    eddy_kz = 0._rk

    !add eddy diffusivity
    allocate(depth_boundary(number_of_boundaries,time))
    depth_boundary = self%get_array(_DEPTH_ON_BOUNDARY_)
    call name_input%get_var(_TURBULENCE_,var)
    select type(var)
    class is(variable_2d)
      !write(*,*) size(var%value, 1)
      !write(*,*) size(var%value, 2)
      eddy_kz(water_bbl_index:ice_water_index,:) = var%value

      !linear interpolation in bbl
      do i = bbl_sediments_index+1,water_bbl_index
        eddy_kz(i,:) = &
        !(i-bbl_sediments_index)*eddy_kz(water_bbl_index+1,:)/(&
        !water_bbl_index+1-bbl_sediments_index)
        eddy_kz(bbl_sediments_index,:)+&
        ((eddy_kz(water_bbl_index+1,:)-eddy_kz(bbl_sediments_index,:))/&
        (depth_boundary(water_bbl_index+1,:)-&
         depth_boundary(bbl_sediments_index,:)))*&
        (depth_boundary(i,:)-depth_boundary(bbl_sediments_index,:))
      end do

      do i = 1,time
        eddy_kz(air_ice_indexes(i)+1:,i) = D_QNAN
      end do

      !to start not from the bottom most interface
      eddy_kz(2:bbl_sediments_index,:) = _DISPERSION_COEFFICIENT_

      new_var_2d = variable_2d(_TURBULENCE_,&
                               'Eddy diffusivity',var%units,eddy_kz)
      call self%add_item(new_var_2d)
    class default
      call fatal_error("Adding turbulence",&
        "Wrong type")
    end select

    !add molecular diffusivity
    allocate(tortuosity(number_of_boundaries,time))
    tortuosity = self%get_array("tortuosity_on_interfaces")
    allocate(value_2d(number_of_boundaries,time))
    value_2d = 0._rk
    value_2d(bbl_sediments_index+1:number_of_boundaries,:) &
      = _INFINITE_DILLUTION_MOLECULAR_DIFFUSIVITY_
    value_2d(2:bbl_sediments_index,:) = &
      _RELATIVE_DYNAMIC_VISCOSITY_*&
      _INFINITE_DILLUTION_MOLECULAR_DIFFUSIVITY_/&
      tortuosity(2:bbl_sediments_index,:)**2
    do i = 1,time
      value_2d(ice_water_index:air_ice_indexes(i),i) = 0._rk
      value_2d(air_ice_indexes(i)+1:,i) = D_QNAN
    end do
    new_var_2d = variable_2d("molecular_diffusivity",&
                             "Molecular diffusivity",var%units,value_2d)
    call self%add_item(new_var_2d)

    !add bioturbation diffusivity
    value_2d = 0._rk
    do i = 1,bbl_sediments_index
    !accuracy problem
      where (depth_boundary(i,:)-depth_boundary(bbl_sediments_index,:)&
          <_MIXED_LAYER_DEPTH_)
        value_2d(i,:) = _MAX_BIOTURBATION_DIFFUSIVITY_
      else where
        value_2d(i,:) = _MAX_BIOTURBATION_DIFFUSIVITY_*exp(-1._rk*(&
          depth_boundary(i,:)-depth_boundary(bbl_sediments_index,:)-&
          _MIXED_LAYER_DEPTH_)/_DECAY_BIOTURBATION_SCALE_)
      end where
    end do
    value_2d(1,:) = 0._rk
    do i = 1,time
      value_2d(air_ice_indexes(i)+1:,i) = D_QNAN
    end do
    new_var_2d = variable_2d("bioturbation_diffusivity",&
                             "Bioturbation diffusivity",var%units,value_2d)
    call self%add_item(new_var_2d)

    !add brine channel gravity drainage diffusivity
    allocate(porosity(number_of_boundaries,time))
    value_2d = 0._rk
    porosity = self%get_array("porosity_on_interfaces")
    do i = 1,time
      z_conv = 0._rk
      !if no ice, go to the next day
      if (air_ice_indexes(i)/=ice_water_index) then
        do j = ice_water_index,air_ice_indexes(i)
          if (porosity(j,i)>_REQUIRED_VOLUME_) then
            z_conv = z_conv+_ICE_LAYERS_RESOLUTION_
          else
            exit
          end if
        end do
        value_2d(ice_water_index:j-1,i) = &
          z_conv*_GRAVITY_DRAINAGE_
      end if
    end do
    new_var_2d = variable_2d("ice_gravity_drainage",&
                             "Ice gravity drainage diffusivity",&
                              var%units,value_2d)
    call self%add_item(new_var_2d)
  end subroutine add_diffusivity

  pure subroutine set_spbm_state_variable(self,is_solid,&
      is_gas,use_bound_up,use_bound_low,bound_up,&
      bound_low,density,sinking_velocity,value,layer)
    class(spbm_state_variable),intent(inout):: self
    logical,optional          ,intent(in)   :: is_solid
    logical,optional          ,intent(in)   :: is_gas
    integer,optional          ,intent(in)   :: use_bound_up
    integer,optional          ,intent(in)   :: use_bound_low
    real(rk),optional         ,intent(in)   :: bound_up
    real(rk),optional         ,intent(in)   :: bound_low
    real(rk),optional         ,intent(in)   :: density
    real(rk),allocatable,dimension(:),optional,intent(in)&
                                            :: sinking_velocity
    real(rk),optional         ,intent(in)   :: value
    integer ,optional         ,intent(in)   :: layer

    if(present(is_solid)) self%is_solid = is_solid
    if(present(is_gas)) self%is_gas = is_gas
    if(present(use_bound_up)) self%use_bound_up = use_bound_up
    if(present(use_bound_low)) self%use_bound_low = use_bound_low
    if(present(bound_up)) self%bound_up = bound_up
    if(present(bound_low)) self%bound_low = bound_low
    if(present(density)) self%density = density
    if(present(sinking_velocity)) self%sinking_velocity = sinking_velocity
    if(present(value).and.present(layer)) self%value(layer) = value
  end subroutine set_spbm_state_variable

  subroutine print_state_variable(self)
    class(spbm_state_variable),intent(in):: self

    integer i

    do i = 1,72!size(self%value,1)
      write(*,'(f9.7,2x,i3)') self%value(i),i
    end do
    write(*,*) self%name
    !write(*,*) "is_solid:",self%is_solid
    !write(*,*) "density:",self%density
    !write(*,'(f9.3)') (/ self%value(size(self%value,1):1:-1) /)
    !write(*,'(x,a,2x,i2,f9.6)') 'up',self%use_bound_up,self%bound_up
    !write(*,'(x,a,i2,f9.6)') 'down',self%use_bound_low,self%bound_low
    !write(*,*) 'sinking',self%sinking_velocity
    _LINE_
  end subroutine print_state_variable
end module
