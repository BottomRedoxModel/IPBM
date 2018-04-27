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

Then:

## Linux(bash shell):
1. Download all required programs into the one folder:

  download SPBM

  ```
  git clone https://github.com/BottomRedoxModel/SPBM.git
  cd SPBM
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```
  
  download configuration files and data for SPBM
  
  ```
  git clone https://github.com/BottomRedoxModel/SPBM_data.git
  cd SPBM_data
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```
   
  download [FABM] from this repository

  ```
  git clone https://github.com/BottomRedoxModel/fabm.git
  cd fabm
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```
  
  This [FABM] version has SPBM registered as the host program already.
  
  download biogeochemistry model BROM
  
  ```
  git clone https://github.com/BottomRedoxModel/brom_niva_module
  cd brom_niva_module
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```

  download [ERSEM] from the official cite.

2. Add FABMDIR and NetCDF_ROOT environment variables

  For example you can add to `~/.bashrc` current lines:

  ```
  export FABMDIR='/path/to/FABM'
  export NetCDF_ROOT='/path/to/NetCDF/bin'
  ```
  
  Reload .bashrc `$ source ~/.bashrc`

3. Make a build 

  Enter SPBM folder and execute `$ bash build_release.sh` - it will make release for the test case using BROM only.

4. Compile the code

  From build folder execute `$ make`

5. Run SPBM

  From build folder execute `$ ./SPBM`

## Windows 10, 8:

1. Download all required programs 

  Right-click in Windows Explorer within the directory where you want to place the SPBM directory, and choose "Git Bash Here", or use PowerShell program.
  At first download SPBM

  ```
  git clone https://github.com/BottomRedoxModel/SPBM.git
  cd SPBM
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```
  
  download all nessesary configuration files and data for SPBM
  
  ```
  git clone https://github.com/BottomRedoxModel/SPBM_data.git
  cd SPBM_data
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```
   
  then download [FABM] from this repository

  ```
  git clone https://github.com/BottomRedoxModel/fabm.git
  cd fabm
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```
  
  This [FABM] version has SPBM registered as the host program already.
  
  download BROM - biogeochemistry model
  
  ```
  git clone https://github.com/BottomRedoxModel/brom_niva_module
  cd brom_niva_module
  git checkout tags/spbm_v0.2 -b spbm_v0.2
  cd ..
  ```

  download [ERSEM] from the official cite.
  
2. Add SPBMDIR environment variable (only if you are going to use pre-compiled NetCDF libraries)

  * In Search, search for and then select: Environment variables or something similar
  * In the user variables specify the name **SPBMDIR** and the value **path:\to\SPBM**

3. Make a build

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
  Select **Advanced** option and specify **-DFABM\_NIVA\_BASE** to `path/to/brom_niva_module` also.
  * Click the **Configure** button until no new (red-coloured) configuration variables appear, then press **Generate** button.

4. Compile the code

  After generating the build system, you should build the software.
  You can do either by opening Visual Studio and choosing **Build All** (after opening **path:\to\SPBM\build\SPBM.sln**, right click on SPBM in **Solution Explorer** and select **Set as StartUp Project**) or typing **make** if using a build system based on makefiles.

5. Run SPBM

  Now you have **SPBM.exe** file in your `path:\to\SPBM\build\Debug(Release)` directory.

## Results visualisation

For the SPBM output files visualisation you can use the python [script].

[Git]:https://git-scm.com/downloads
[FABM]:http://fabm.net
[CMake]:https://cmake.org/
[NetCDF]:http://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html
[ERSEM]:https://gitlab.ecosystem-modelling.pml.ac.uk/stable/ERSEM/tree/master
[script]:https://github.com/lisapro/ice_brom_pic
