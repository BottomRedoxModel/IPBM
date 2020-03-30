# A 1-Dimensional Sympagic-Pelagic-Benthic transport model, (SPBM)
## About
Coupled simulation of ice, water column, and sediment biogeochemistry.

## Supported compilers:
* recent gfortran compiler (part of GCC)
* Intel Fortran Compiler version 12.1 or higher

## How to use
At first you must have compliant compiler, [Git], [CMake], and [NetCDF] Fortran library compiled with the same Fortran compiler as used for compiling SPBM.
For the VisualStudio solution under Windows pre-compiled NetCDF libraries are provided.
In case of you do not use any biogeochemical model you do not need to specify a corresponding environment variable.

Then:

## Linux(bash shell):
1. Download all required programs:

  download SPBM, configuration files and data for SPBM;

  ```
  git clone https://github.com/BottomRedoxModel/SPBM.git
  git clone https://github.com/BottomRedoxModel/SPBM_data.git
  ```
   
  download [FABM] from the following repository, that version has SPBM registered as the host program already;

  ```
  git clone https://github.com/BottomRedoxModel/fabm.git
  cd fabm
  git checkout spbm
  cd ..
  ```

  download biogeochemistry model BROM;
  
  ```
  git clone https://github.com/BottomRedoxModel/brom_niva_module
  cd brom_niva_module
  git checkout dev-sham
  cd ..
  ```

2. Add FABMDIR, BROMDIR, NetCDF_ROOT, NetCDF_modules environment variables.

  For example you can add to `~/.bashrc` current lines:

  ```
  export FABMDIR='/path/to/fabm'
  export BROMDIR='/path/to/brom_niva_module'
  export NetCDF_ROOT='/path/to/NetCDF/bin'
  export NetCDF_modules='/path/to/NetCDF/modules'
  ```

  NetCDF_ROOT should point to the nf-config location.
  NetCDF_modules should point to netcdf module `*.mod` files location.  
  For some OS these variables point to the similar location.
  Sometimes, in case of installing NetCDF libraries into the standard directory it is not
  necessary to specify their location.
  Reload .bashrc `$ source ~/.bashrc`

3. Make a build.

  Enter SPBM folder and execute `$./prepare_build.sh north`, it will make a build
  directory with all necessary data files.

4. Compile the code.

  From build folder execute `$ make`

5. Run SPBM: `$ ./SPBM`.

  SPBM program needs to be launched:
  1. The data file with forcing (for example `wadden_sea_out.nc`);
  2. Two configuration files: SPBM configuration file itself - spbm.yaml and the second [FABM] configuration file specified in spbm.yaml (for example fabm_brom.yaml).

## Windows 10, 8:
1. Download all required programs: 

  Right-click in Windows Explorer within the directory where you want to place the SPBM directory, and choose "Git Bash Here", or use PowerShell program.
  Then input the commands from the previous section.
  
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
  Select **Advanced** option and specify **-DFABM\_NIVA\_BASE** to **path/to/brom_niva_module** also.
  * Click the **Configure** button until no new (red-coloured) configuration variables appear, then press **Generate** button.

4. Compile the code.

  After generating the build system, you should build the software.
  You can do either by opening Visual Studio and choosing **Build All** (after opening **path:\to\SPBM\build\SPBM.sln**, right click on SPBM in **Solution Explorer** and select **Set as StartUp Project**) or typing **make** if using a build system based on makefiles.

5. Run SPBM.

## Results visualisation

For the SPBM output files visualisation you can use the python [script] or any other software recognizing [NetCDF] files.

[Git]:https://git-scm.com/downloads
[FABM]:http://fabm.net
[CMake]:https://cmake.org/
[NetCDF]:http://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html
[script]:https://github.com/lisapro/ice_brom_pic
