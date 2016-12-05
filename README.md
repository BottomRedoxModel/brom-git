# BROM
## About
Bottom RedOx Model (BROM): a coupled benthic-pelagic model for simulation of water and sediment biogeochemistry.

## Supported compilers:
* gfortran 4.7 or higher (part of GCC)
* Intel Fortran Compiler version 12.1 or higher

## How to use
At first you must download [FABM] and do all prerequisites it needs (you should have compliant compiler, [Git], [CMake], and [NetCDF] Fortran library compiled with the same Fortran compiler as used for compiling BROM & FABM. For the VisualStudio solution under Windows pre-compiled NetCDF libraries are provided.) Then:

## Linux(bash shell):
1. Download BROM

  `$ git clone https://github.com/e-yakushev/brom-git.git brom-git`

2. Add FABMDIR and NetCDF_ROOT environment variables

  For example you can add to `~/.bashrc` current lines:

  ```
  export FABMDIR='/path/to/FABM'
  export NetCDF_ROOT='/path/to/NetCDF/bin'
  ```

  You may also need in case of compiling NetCDF libraries from source not to standart folder add 2 more lines:

  ```
  export NCDIR='/path/to/netcdf'
  export LD_LIBRARY_PATH=$NCDIR/lib:$LD_LIBRARY_PATH
  ```

  Don't forget reload .bashrc `$ source ~/.bashrc`

3. Make a build 

  Enter brom-git folder and execute `$ bash build_brom.sh`

  or

  You can greate build folder manually: creat it, copy there files from `/data` folder, then execute command:

  `$ cmake path/to/BROM/code -DFABM_BASE=$FABMDIR`

  from it.

4. Compile the code

  From build folder execute `$ make`

5. Run BROM

  From build folder execute `$ ./brom`

## Windows:

1. Download BROM

  Right-click in Windows Explorer within the directory where you want to place the BROM direcrory, and choose "Git Bash Here". In the terminal that appears type:

  `$ git clone https://github.com/e-yakushev/brom-git.git brom-git`

  if using other software, use URL `https://github.com/e-yakushev/brom-git.git` and recommended directory is `..\brom-git`

2. Add BROMDIR environment variable

  Windows 10 and Windows 8

  * In Search, search for and then select: System (Control Panel)
  * Click the **Advanced system settings** link
  * Click **Environment Variables**. In the section **System Variables** click **New**
  * In the **New System Variable** specify the name **BROMDIR** and the value **path:\to\brom-git**

3. Make a build

  * Start "CMake"
  * Browse the **Where is the source code** to the **path:\to\brom-git\code**
  * Browse the **Where to build the binaries** - e.g. **path:\to\brom-git\build**
  * Click the **Configure** button. Select a build system generator, if you use Intel Visual Fortran with Visual Studio integration and want to use NetCDF libraries that come with BROM please select a 32-bit generator.
  * Now all configuration variables for the build system are listed and you can change them according to your preferences. You need set FABM_BASE variable to the directory where you have downloaded [FABM]
  * Click the **Configure** button until no new (red-coloured) configuration variables appear, then press **Generate** button.

4. Compile the code

  After generating the build system, you should build the software. You can do either by opening Visual Studio and choosing **Build All**, or typing **make** if using a build system based on makefiles.

5. Run BROM

  Now you have **brom.exe** file in your `path:\to\brom-git\build\Debug(Release)` directory. It needs 4 files as input data: brom.yaml, fabm.yaml, nns_annual.nc, start.dat. You can find it in `..\brom-git\data` folder. In case of running BROM under Visual Studio remember to specify the working directory (`..\brom-git\data`).

[Git]:https://git-scm.com/downloads
[FABM]:http://fabm.net
[CMake]:https://cmake.org/
[NetCDF]:http://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html
