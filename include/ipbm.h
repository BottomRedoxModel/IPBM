#define _FABM_YAML_NAME_			      get_ipbm_string_parameter('fabm_yaml_name')
#define _FILE_NAME_                 get_ipbm_string_parameter('input_file_name')
#define _INITIALIZATION_SINCE_YEAR_ get_ipbm_real_parameter('initialization_since_year') 
#define _DEPTH_ON_BOUNDARY_         get_ipbm_string_parameter('depth_on_boundary')
#define _OCEAN_TIME_                get_ipbm_string_parameter('time')
#define _TEMPERATURE_               get_ipbm_string_parameter('temperature')
#define _SALINITY_                  get_ipbm_string_parameter('salinity')
#define _RHO_                       get_ipbm_string_parameter('density')
#define _TURBULENCE_                get_ipbm_string_parameter('eddy_diffusivity')
#define _ICE_THICKNESS_             get_ipbm_string_parameter('ice_thickness')
#define _SNOW_THICKNESS_            get_ipbm_string_parameter('snow_thickness')
#define _ICE_SURFACE_TEMPERATURE_   get_ipbm_string_parameter('ice_surface_temperature')

#define _FILE_NAME_ICE_             'ice.nc'
#define _FILE_NAME_WATER_           'water.nc'
#define _FILE_NAME_SEDIMENTS_       'sediments.nc'

#define _SECONDS_PER_CIRCLE_        300

#define _LONGITUDE_                 get_ipbm_real_parameter(longitude) 77.8_rk
#define _LATITUDE_                  get_ipbm_real_parameter(latitude) 77.3_rk

#define _WIDTH_BBL_                 get_ipbm_real_parameter(width_bbl) 0.5_rk
#define _RESOLUTION_BBL_            get_ipbm_real_parameter(resolution_bbl) 0.1_rk

#define _WIDTH_SEDIMENTS_           get_ipbm_real_parameter(width_sediments) 0.1_rk
#define _RESOLUTION_SEDIMENTS_      get_ipbm_real_parameter(resolution_sediments) 0.02_rk

#define _ICE_LAYERS_RESOLUTION_     get_ipbm_real_parameter(ice_layers_resolution) 0.06_rk

#define _O2_                        'O2_o'
#define _CaCO3_                     'L2_c'

#define _Alk_                       'B_C_Alk'
#define _DIC_                       'B_C_DIC'
#define _PO4_                       'B_NUT_PO4'
#define _NO3_                       'B_NUT_NO3'
#define _Si_                        'B_NUT_Si'
#define _NH4_                       'B_NUT_NH4'
#define _H2S_                       'B_S_H2S'
#define _CH4_                       'B_CH4_CH4'
#define _S0_                        'B_S_S0'
#define _Fe3_                       'B_Fe_Fe3'
#define _FeCO3_                     'B_Fe_FeCO3'
#define _FeS_                       'B_Fe_FeS'
#define _FeS2_                      'B_Fe_FeS2'
#define _Mn4_                       'B_Mn_Mn4'
#define _MnCO3_                     'B_Mn_MnCO3'
#define _MnS_                       'B_Mn_MnS'
#define _Sipart_                    'B_Si_Sipart'
#define _Phy_                       'B_BIO_Phy'
#define _PON_                       'B_BIO_PON'
