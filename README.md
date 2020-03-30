# A 1-Dimensional Sympagic-Pelagic-Benthic transport model, (SPBM)
## About
Coupled simulation of ice, water column, and sediment biogeochemistry.

## Supported compilers:
* recent gfortran compiler (part of GCC)
* Intel Fortran Compiler version 12.1 or higher

## How to use
This manual is referred to the SPBM v0.2.
At first you must have compliant compiler, [Git], [CMake], and [NetCDF] Fortran library compiled with the same Fortran compiler as used for compiling SPBM.
For the VisualStudio solution under Windows pre-compiled NetCDF libraries are provided.
In case of you do not use any biogeochemical model you do not need to specify a corresponding environment variable.

Then:

## Linux(bash shell):
1. Download all required programs and switch to the spbm_v0.2 tag branch everywhere except [ERSEM]:

  download SPBM;

  ```
  git clone https://github.com/BottomRedoxModel/SPBM.git
  cd SPBM
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```
  
  download configuration files and data for SPBM;
  
  ```
  git clone https://github.com/BottomRedoxModel/SPBM_data.git
  cd SPBM_data
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```
   
  download [FABM] from this repository, this version has SPBM registered as the host program already;

  ```
  git clone https://github.com/BottomRedoxModel/fabm.git
  cd fabm
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```

  download biogeochemistry model BROM;
  
  ```
  git clone https://github.com/BottomRedoxModel/brom_niva_module
  cd brom_niva_module
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```

  download [ERSEM] from the official cite.

2. Add FABMDIR, ERSEMDIR, BROMDIR, NetCDF_ROOT, NetCDF_modules environment variables.

  For example you can add to `~/.bashrc` current lines:

  ```
  export FABMDIR='/path/to/fabm'
  export ERSEMDIR='/path/to/ersem'
  export BROMDIR='/path/to/brom_niva_module'
  export NetCDF_ROOT='/path/to/NetCDF/bin'
  export NetCDF_modules='/path/to/NetCDF/modules'
  ```

  NetCDF_ROOT should point to the nf-config location.
  NetCDF_modules should point to netcdf module `*.mod` files location.  
  For some OS these variables point to the similar location.
  Somtimes, in case of installing NetCDF libraries into the standard directory it is not
  necessary to specify their location.
  Reload .bashrc `$ source ~/.bashrc`

3. Make a build.

  Enter SPBM folder and execute `$ bash build_release.sh`, it will make a build
  directory, copy there all necessary data files, and launch cmake to make a build.

4. Compile the code.

  From build folder execute `$ make`

5. Run SPBM.

  SPBM program needs to be launched:
  1. The data file with forcing (for example `ROMS_Laptev_Sea.nc`);
  2. Two configuration files: SPBM configuration file itself - spbm.yaml and the second [FABM] configuration file specified in spbm.yaml (for example fabm_ersem.yaml).
  For the [ERSEM] test case rename `spbm_ersem.yaml` to `spbm.yaml`.
  Respectively rename `spbm_brom.yaml` to `spbm.yaml` for the BROM test case.
  Then execute `$ ./SPBM`

## Windows 10, 8:

1. Download all required programs and switch to the spbm_v0.2 tag branch everywhere except [ERSEM]:

  Right-click in Windows Explorer within the directory where you want to place the SPBM directory, and choose "Git Bash Here", or use PowerShell program.
  At first download SPBM;

  ```
  git clone https://github.com/BottomRedoxModel/SPBM.git
  cd SPBM
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```
  
  download all necessary configuration files and data for SPBM;
  
  ```
  git clone https://github.com/BottomRedoxModel/SPBM_data.git
  cd SPBM_data
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```
   
  then download [FABM] from this repository, this version has SPBM registered as the host program already;

  ```
  git clone https://github.com/BottomRedoxModel/fabm.git
  cd fabm
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```
    
  download BROM - biogeochemistry model;
  
  ```
  git clone https://github.com/BottomRedoxModel/brom_niva_module
  cd brom_niva_module
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```

  download [ERSEM] from the official cite.
  
2. Add SPBMDIR environment variable (only if you are going to use pre-compiled NetCDF libraries).

  * In Search, search for and then select: Environment variables or something similar
  * In the user variables specify the name **SPBMDIR** and the value **path:\to\SPBM**

3. Make a build.

  * Start "CMake"
  * Browse the **Where is the source code** to the **path:\to\SPBM\src**
  * Browse the **Where to build the binaries** - e.g. **path:\to\SPBM\build**
  * Click the **Configure** button.
  Select a build system generator, if you use Intel Visual Fortran with Visual Studio integration and want to use NetCDF libraries that come with SPBM please select a 32-bit generator.
  * Now all configuration variables for the build system are listed and you can change them according to your preferences.
  You need set **FABM\_BASE** variable to the directory where you have downloaded [FABM].
  Then click the **Configure** button again.
  Specify **FABM\_ERSEM\_BASE** as well.
  Then click the **Configure** button again.
  Select **Advanced** option and specify **-DFABM\_NIVA\_BASE** to **path/to/brom_niva_module** also.
  * Click the **Configure** button until no new (red-coloured) configuration variables appear, then press **Generate** button.

4. Compile the code.

  After generating the build system, you should build the software.
  You can do either by opening Visual Studio and choosing **Build All** (after opening **path:\to\SPBM\build\SPBM.sln**, right click on SPBM in **Solution Explorer** and select **Set as StartUp Project**) or typing **make** if using a build system based on makefiles.

5. Run SPBM.

  SPBM program needs to be launched:
  1. The data file with forcing (for example `ROMS_Laptev_Sea.nc`);
  2. Two configuration files: SPBM configuration file itself - `spbm.yaml` and the second [FABM] configuration file specified in `spbm.yaml` (for example `fabm_ersem.yaml`).
  All these files are available in the data folder (`SPBM_data`).
  For the [ERSEM] test case rename `spbm_ersem.yaml` to `spbm.yaml`.
  Respectively rename `spbm_brom.yaml` to `spbm.yaml` for the BROM test case.
  Now you have `SPBM.exe` file in your `path:\to\SPBM\build\Debug(Release)` directory.
  You should place all files into the one folder (for example `SPBM.exe`, `spbm.yaml`, `fabm_ersem.yaml`, `ROMS_Laptev_Sea.nc`) and launch `SPBM.exe`.

## Results visualisation

For the SPBM output files visualisation you can use the python [script] or any other software recognizing [NetCDF] files.

[Git]:https://git-scm.com/downloads
[FABM]:http://fabm.net
[CMake]:https://cmake.org/
[NetCDF]:http://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html
[ERSEM]:https://gitlab.ecosystem-modelling.pml.ac.uk/stable/ERSEM/tree/master
[script]:https://github.com/lisapro/ice_brom_pic
