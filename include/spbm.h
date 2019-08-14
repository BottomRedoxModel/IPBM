#define _LINE_      write(*,*) repeat('*',79)
#define _PAUSE_     read(*,*)

#define _SPBM_FILE_NAME_ 'spbm.yaml'

#define _FABM_YAML_NAME_			get_spbm_string_parameter('fabm_yaml_name')
#define _FILE_NAME_                 get_spbm_string_parameter('input_file_name')
#define _INITIALIZATION_SINCE_YEAR_ get_spbm_real_parameter('initialization_since_year') 
#define _DEPTH_ON_BOUNDARY_         get_spbm_string_parameter('depth_on_boundary')
#define _OCEAN_TIME_                get_spbm_string_parameter('time')
#define _TEMPERATURE_               get_spbm_string_parameter('temperature')
#define _SALINITY_                  get_spbm_string_parameter('salinity')
#define _RHO_                       get_spbm_string_parameter('density')
#define _TURBULENCE_                get_spbm_string_parameter('eddy_diffusivity')
#define _ICE_THICKNESS_             get_spbm_string_parameter('ice_thickness')
#define _SNOW_THICKNESS_            get_spbm_string_parameter('snow_thickness')
#define _ICE_SURFACE_TEMPERATURE_   get_spbm_string_parameter('ice_surface_temperature')
#define _SHORTWAVE_RADIATION_       get_spbm_string_parameter('shortwave_radiation')

#define _RE_DAY_					get_spbm_real_parameter('repeat_day') 
#define _RE_YEAR_					get_spbm_real_parameter('repeat_year')
#define _SECONDS_PER_CIRCLE_        get_spbm_real_parameter('seconds_per_circle')
#define _IS_SHORTWAVE_RADIATION_IS_PAR_ get_spbm_real_parameter('is_par')

#define _REQUIRED_VOLUME_           get_spbm_real_parameter('required_volume')
#define _IALGAE_VELOCITY_           get_spbm_real_parameter('ice_algae_velocity')

#define _FILE_NAME_ICE_             'ice.nc'
#define _FILE_NAME_WATER_           'water.nc'
#define _FILE_NAME_SEDIMENTS_       'sediments.nc'

#define _BACKGROUND_ATTENUATION_    get_spbm_real_parameter('background_attenuation')
#define _SILT_ATTENUATION_          get_spbm_real_parameter('silt_attenuation')
#define _SILT_CONCENTRATION_        get_spbm_real_parameter('silt_concentration')

#define _LATITUDE_                  get_spbm_real_parameter('latitude')

#define _WIDTH_BBL_                 get_spbm_real_parameter('width_bbl')
#define _RESOLUTION_BBL_            get_spbm_real_parameter('resolution_bbl')
#define _DISPERSION_COEFFICIENT_    get_spbm_real_parameter('dispersion_coefficient')

#define _WIDTH_SEDIMENTS_           get_spbm_real_parameter('width_sediments')
#define _RESOLUTION_SEDIMENTS_      get_spbm_real_parameter('resolution_sediments')
#define _MAX_POROSITY_              get_spbm_real_parameter('max_porosity')
#define _MIN_POROSITY_              get_spbm_real_parameter('min_porosity')
#define _BURIAL_VELOCITY_           get_spbm_real_parameter('burial_velocity')

#define _ICE_LAYERS_RESOLUTION_     get_spbm_real_parameter('ice_layers_resolution')

#define _Alk_                       get_spbm_string_parameter('Alk')
#define _DIC_                       get_spbm_string_parameter('DIC')

#define _PO4_                       get_spbm_string_parameter('PO4')
#define _NH4_                       get_spbm_string_parameter('NH4')
#define _NO3_                       get_spbm_string_parameter('NO3')
#define _SO4_                       get_spbm_string_parameter('SO4')
#define _Si_                        get_spbm_string_parameter('Si')

#define _O2_                        get_spbm_string_parameter('O2')
#define _CH4_                       get_spbm_string_parameter('CH4')

#define _DOM_                       get_spbm_string_parameter('DOM')
#define _DOC_                       get_spbm_string_parameter('DOC')

#define _CaCO3_                     get_spbm_string_parameter('CaCO3')
#define _S0_                        get_spbm_string_parameter('S0')
#define _Fe3_                       get_spbm_string_parameter('Fe3')
#define _FeCO3_                     get_spbm_string_parameter('FeCO3')
#define _FeS_                       get_spbm_string_parameter('FeS')
#define _FeS2_                      get_spbm_string_parameter('FeS2')
#define _Mn4_                       get_spbm_string_parameter('Mn4')
#define _MnCO3_                     get_spbm_string_parameter('MnCO3')
#define _MnS_                       get_spbm_string_parameter('MnS')
#define _Sipart_                    get_spbm_string_parameter('Sipart')
#define _Phy_                       get_spbm_string_parameter('Phy')
#define _Het_                       get_spbm_string_parameter('Het')
#define _POM_                       get_spbm_string_parameter('POM')
#define _POC_                       get_spbm_string_parameter('POC')
#define _SmallPOM_                  get_spbm_string_parameter('SmallPOM')
#define _MediumPOM_                 get_spbm_string_parameter('MediumPOM')
#define _LargePOM_                  get_spbm_string_parameter('LargePOM')
#define _iDiatoms_                  get_spbm_string_parameter('iDiatoms')
#define _Diatoms_                   get_spbm_string_parameter('Diatoms')
#define _NanoPhy_                   get_spbm_string_parameter('NanoPhy')
#define _PicoPhy_                   get_spbm_string_parameter('PicoPhy')
#define _MicroPhy_                  get_spbm_string_parameter('MicroPhy')
#define _MesoZoo_                   get_spbm_string_parameter('MesoZoo')
#define _MicroZoo_                  get_spbm_string_parameter('MicroZoo')
#define _NanoFlag_                  get_spbm_string_parameter('NanoFlag')

#define _RELAXATION_FILE_NAME_      get_spbm_string_parameter('relaxation_file_name')
#define _RELAXATION_PARAMETER_      get_spbm_real_parameter('relaxation_parameter')
#define _DX_ get_spbm_real_parameter('dx')

#define _DIC_rel_                   get_spbm_string_parameter('DICrel')
#define _PO4_rel_                   get_spbm_string_parameter('PO4rel')
#define _NH4_rel_                   get_spbm_string_parameter('NH4rel')
#define _NO3_rel_                   get_spbm_string_parameter('NO3rel')
#define _SO4_rel_                   get_spbm_string_parameter('SO4rel')
#define _Si_rel_                    get_spbm_string_parameter('Sirel')
#define _O2_rel_                    get_spbm_string_parameter('O2rel')
#define _CH4_rel_                   get_spbm_string_parameter('CH4rel')

#define _CH4_flux_                  get_spbm_string_parameter('CH4flux')
#define _DOM_flux_                  get_spbm_string_parameter('DOMflux')
#define _DOC_flux_                  get_spbm_string_parameter('DOCflux')
#define _POM_flux_                  get_spbm_string_parameter('POMflux')
#define _POC_flux_                  get_spbm_string_parameter('POCflux')
#define _DIC_flux_                  get_spbm_string_parameter('DICflux')
#define _Alk_flux_                  get_spbm_string_parameter('Alkflux')
