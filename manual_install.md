## Manual Installation:
If running ns_install() does not work, please email the author (who wants to know what went wrong).  You can also try the following manual installation instructions.

Below, PATH_NEUTRAL_SURFACE is a generic stand-in for the path to the folder containing this README.md file.
PATH_TOPOBARIC_SURFACE is a generic stand-in for the path to the topobaric-surface subfolder in PATH_NEUTRAL_SURFACE.

If you will use a bash shell, start by entering the following commands,
modified to use the path to the location of this file, in a bash shell: 
```
  PATH_NEUTRAL_SURFACE=~/work/projects-gfd/neutral-surfaces
  PATH_TOPOBARIC_SURFACE=$PATH_NEUTRAL_SURFACE/topobaric-surface
```

### Step 1: Check MATLAB version
  In MATLAB, run the following command:
  ```
  >> version 
  ```
  Ensure that the output is something greater than or equal to than 9.1 (R2016b).
  If not, you must update MATLAB to use Topobaric Surface. 

### Step 2: Install recon.jar 
In MATLAB, run
```
>> version('-java')
```
The first numbers below indicate MATLAB's Java Runtime Environment (JRE) version.

If this is 1.7 or 1.8:

Copy (or move) recon_1.7.jar or recon_1.8.jar to recon.jar, located in $PATH_TOPOBARIC_SURFACE/src/recon/build/  
For example, run the following in bash (changing 1.8 to 1.7 as appropriate):
```
cp $PATH_TOPOBARIC_SURFACE/src/recon/build/recon_1.8 $PATH_TOPOBARIC_SURFACE/src/recon/build/recon.jar
```
Continue to Step 3.

Otherwise:
You can compile ReCon yourself, or contact the author requesting an update.

To compile ReCon yourself,
- Install Java Development Kit (JDK) and Apache Ant:
  - On Ubuntu:
    sudo apt install default-jdk ant
  - On MacOS:
    - Install (the most recent version of) the JDK following instructions from
      https://www.oracle.com/technetwork/java/javase/downloads/index.html
    - Install homebrew following instructions from 
      https://brew.sh/
    - Install Apache Ant by the following command in Terminal:
      brew install ant
  - On Windows:
    - Install (the most recent version of) the JDK following instructions from
      https://www.oracle.com/technetwork/java/javase/downloads/index.html
    - Install Apache Ant, following instructions from
      https://ant.apache.org/manual/install.html
- Edit the text file PATH_TOPOBARIC_SURFACE/src/recon/build.xml and replace each instance of "1.8" in the <javac ...> lines to specify MATLAB's JRE version.
    Alternatively, run the following commands in bash, replacing ??? below with MATLAB's JRE version -- with the period escaped, e.g. 1.8 is entered as 1\.8 
    ```
    sed -i "s/1\.8/???/g" $PATH_TOPOBARIC_SURFACE/src/recon/build.xml
    ```
- Compile ReCon: run the following command in a shell
  ```
  ant -buildfile $PATH_TOPOBARIC_SURFACE/src/recon/build.xml
  ```
- Confirm the existence of the file $PATH_TOPOBARIC_SURFACE/src/recon/build/recon.jar

### Step 3: Add ReCon to MATLAB's javaclasspath
In MATLAB, run
```
>> prefdir
```
to determine MATLAB's preferences directory. 

Create a text file in this directory called javaclasspath.txt if it does not already exist.

Add the following line to javaclasspath.txt
```
$PATH_TOPOBARIC_SURFACE/src/recon/build/recon.jar
```
Alternatively, run the following commands in a shell, replacing ~/Documents/MATLAB with the output of MATLAB's prefdir:
```
PREFDIR=~/Documents/MATLAB
echo $PATH_TOPOBARIC_SURFACE/src/recon/build/recon.jar" >> $PREFDIR/javaclasspath.txt
```
  

### Step 4: Increase memory that MATLAB allocates to Java
In MATLAB, navigate to Preferences -> General -> Java HEAP Memory

Move the slider far to the right.

e.g. the default value is 512 MB, and you increase it to 2048 MB.

### Step 5: Setup MEX with C
In MATLAB, run
```
>> mex('-setup', 'C')
```
If that fails, ensure you have a C compiler installed on your system.
- On Ubuntu: sudo apt install gcc
- On MacOS: Install Xcode from the App Store (to get clang)
- On Windows: Install gcc via cygwin following instructions (and links therein) at
   https://gcc.gnu.org/install/binaries.html
   
Then try the above mex command in MATLAB again (perhaps restarting MATLAB prior).

### Step 6: Add Topobaric Surface to MATLAB's search path
  In MATLAB, run
  ```
  >> userpath
  ```
  The output is your user path. 
  Create a text file in this directory called startup.m if it does not already exist.
  Add the following lines to the end of startup.m, replacing PATH_NEUTRAL_SURFACE with
  the path to the neutral-surfaces directory. 
  ```
  run('PATH_NEUTRAL_SURFACE' filesep 'ns_add_to_path.m')
  ```

### Step 7: Restart MATLAB. Confirm that the output of
```
>> javaclasspath
```
has recon.jar in the last line. 

### Step 8: Download and link data 
This step is optional, only necessary to run the scripts in $PATH_TOPOBARIC_SURFACE/run

Download the ECCO2 and OCCA files by running the following in a bash shell, replacing the paths assigned to $ECCO2 and $OCCA with your own choices:
```
ECCO2=~/work/data/ECCO2/latlon
wget -P $ECCO2/SALT/ ftp://ecco.jpl.nasa.gov/ECCO2/cube92_latlon_quart_90S90N/SALT.nc/SALT.1440x720x50.20021223.nc
wget -P $ECCO2/THETA/ ftp://ecco.jpl.nasa.gov/ECCO2/cube92_latlon_quart_90S90N/THETA.nc/THETA.1440x720x50.20021223.nc 
wget -P $ECCO2/UVEL/ ftp://ecco.jpl.nasa.gov/ECCO2/cube92_latlon_quart_90S90N/UVEL.nc/UVEL.1440x720x50.20021223.nc 
wget -P $ECCO2/VVEL/ ftp://ecco.jpl.nasa.gov/ECCO2/cube92_latlon_quart_90S90N/VVEL.nc/VVEL.1440x720x50.20021223.nc 
wget -P $ECCO2/SSH/ ftp://ecco.jpl.nasa.gov/ECCO2/cube92_latlon_quart_90S90N/SSH.nc/SSH.1440x720.20021222.nc
wget -P $ECCO2/SSH/ ftp://ecco.jpl.nasa.gov/ECCO2/cube92_latlon_quart_90S90N/SSH.nc/SSH.1440x720.20021223.nc
wget -P $ECCO2/SSH/ ftp://ecco.jpl.nasa.gov/ECCO2/cube92_latlon_quart_90S90N/SSH.nc/SSH.1440x720.20021224.nc
wget -P $ECCO2/GAMMA/ https://ndownloader.figshare.com/files/14536058
mv $ECCO2/GAMMA/14536058 $ECCO2/GAMMA/ GAMMA.1440x720x50.20021223.mat
wget -P $ECCO2/omega_v1.1gjs https://ndownloader.figshare.com/files/14536061
mv $ECCO2/omega_v1.1gjs/14536061 $ECCO2/omega_v1.1gjs/omega.1440x720.20021223.from_SIGMA1000_through_(180,0,1000).Boussinesq.mat
wget -P $ECCO2/omega_v1.1gjs https://ndownloader.figshare.com/files/14536064
mv $ECCO2/omega_v1.1gjs/14536064 $ECCO2/omega_v1.1gjs/omega.1440x720.20021223.from_SIGMA2000_through_(180,0,2000).Boussinesq.mat

OCCA=~/work/data/OCCA
wget -P $OCCA/ ftp://mit.ecco-group.org/ecco_for_las/OCCA_1x1_v2/2004-6/annual/DDsalt.0406annclim.nc
wget -P $OCCA/ ftp://mit.ecco-group.org/ecco_for_las/OCCA_1x1_v2/2004-6/annual/DDtheta.0406annclim.nc
wget -P $OCCA/ ftp://mit.ecco-group.org/ecco_for_las/OCCA_1x1_v2/2004-6/annual/DDphihyd.0406annclim.nc
wget -P $OCCA/ ftp://mit.ecco-group.org/ecco_for_las/OCCA_1x1_v2/2004-6/annual/DDetan.0406annclim.nc
wget -P $OCCA/omega_v1.1gjs https://ndownloader.figshare.com/files/14536133
mv $OCCA/omega_v1.1gjs/14536133 $OCCA/omega_v1.1gjs/    omega.0406annclim.from_SIGMA1000_through_(180,0,1000).nonBoussinesq.mat

printf $ECCO2 > $PATH_NEUTRAL_SURFACE/lib/dat/PATH_ECCO2.txt
printf $OCCA > $PATH_NEUTRAL_SURFACE/lib/dat/PATH_OCCA.txt
```

Alternatively, manually download the data through a web browser.  Go to ftp://ecco.jpl.nasa.gov/ECCO2/cube92_latlon_quart_90S90N/ and to ftp://mit.ecco-group.org/ecco_for_las/OCCA_1x1_v2/2004-6/annual/ 

Download the files named in the wget commands above, maintaining the directory structure as on the ftp site. 
That is, if ~/work/data/ECCO2/latlon/ is your base ECCO2 data directory, then the SALT data file should be in 
  ~/work/data/ECCO2/latlon/SALT.nc/SALT.1440x720x50.20021223.nc
and similarly for others. 
If ~/work/data/OCCA/ is your base OCCA data directory, then the SALT data file should be in
  ~/work/data/OCCA/DDsalt.0406annclim.nc
and similarly for others. 

Then, download files from the following URLs to your local computer, editing local paths as necessary:
- https://ndownloader.figshare.com/files/14536058 to ~/work/data/ECCO2/latlon/GAMMA/GAMMA.1440x720x50.20021223.mat
- https://ndownloader.figshare.com/files/14536061 to ~/work/data/ECCO2/latlon/omega_v1.1gjs/omega.1440x720.20021223.from_SIGMA1000_through_(180,0,1000).Boussinesq.mat
- https://ndownloader.figshare.com/files/14536064 to ~/work/data/ECCO2/latlon/omega_v1.1gjs/omega.1440x720.20021223.from_SIGMA2000_through_(180,0,2000).Boussinesq.mat
- https://ndownloader.figshare.com/files/14536133 to ~/work/data/OCCA/omega_v1.1gjs/omega.0406annclim.from_SIGMA1000_through_(180,0,1000).nonBoussinesq.mat

Then create a text file $PATH_NEUTRAL_SURFACE/lib/dat/PATH_ECCO2.txt containing, as single line of text, the base ECCO2 data directory.
Also create a text file $PATH_NEUTRAL_SURFACE/lib/dat/PATH_OCCA.txt containing, as single line of text, the base OCCA data directory.

## Copyright:
This file is part of Neutral Surfaces.

Copyright (C) 2019  Geoff Stanley

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Author(s) : Geoff Stanley

Email     : g.stanley@unsw.edu.au

Email     : geoffstanley@gmail.com

Version   : 2.0.0
