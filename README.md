# Topobaric Surface
## Software to compute topobaric surfaces and geostrophic streamfunctions

Topobaric surfaces are approximately neutral surfaces that are highly accurate, fast to
compute, and possess an exact geostrophic streamfunction. Underlying topobaric surfaces
is a multivalued functional relationship between the in-situ density and the pressure on
a neutral surface. The single-valued branches of this function are valid on regions that
are determined by the Reeb graph of the pressure on the surface. This software generates
that Reeb graph from any other approximately neutral surface, then empirically fits the 
density to the pressure using simple functions, then updates the surface such that the
density on the surface matches that given by the these simple functions, repeating this
procedure iteratively until the surface converges. 

The exact geostrophic streamfunction on a hypothetical neutral surface exists but is 
ill-defined. It is approximated by the topobaric geostrophic streamfunction, which 
is well-defined and services any approximately neutral surface. It, too, is built upon
a multivalued function relationship between in-situ density and pressure on a neutral
surface. 

A simpler variant of the topobaric geostrophic streamfunction is the orthobaric Montgomery
potential, which inverts the problem and fits the pressure as a single-valued function of
the density on a surface.


## References:
Stanley, G. J. (2019a) Neutral surface topology. Ocean Modelling, accepted.

Stanley, G. J. (2019b) The exact geostrophic streamfunction for neutral surfaces. Ocean Modelling, submitted.



## Requirements:
MATLAB 2016b or higher (tested on 2017b and 2018b) with
 - the Image Processing Toolbox (for bwconncomp), and
 - the Optimization Toolbox (for lsqlin)


## Contents:
- ./COPYING         - GNU general public license
- ./COPYING.lesser  - GNU lesser general public license
- ./fex/            - software from the MATLAB file exchange
- ./lib/            - extra libraries for analysing surfaces, and computing geostrophic streamfunctions
- ./lib/eos/        - alternative versions of the equation of state (JMD1995, TEOS-10)
- ./run/            - example scripts
- ./src/            - source code for topobaric surfaces
- ./src/private/    - internal functions for topobaric surfaces
- ./src/recon/      - modified ReCon software
- ./README.md       - this file

In particular,
- ./src/topobaric_surface.m     - create topobaric surfaces (Stanley 2019a)
- ./src/topobaric_geostrf.m     - create topobaric geostrophic streamfunction (Stanley 2019b)
- ./lib/orthobaric_montgomery.m - create orthobaric Montgomery potentials (Stanley 2019b)
- ./lib/zhanghogg92.m           - create Zhang and Hogg (1992) geostrophic streamfunctions
- ./lib/mcdougallklocker10.m    - create McDougall and Klocker (2010) geostrophic streamfunctions
- ./lib/cunningham00.m          - create Cunningham (2000) geostrophic streamfunctions
- ./lib/isopycnal.m             - create potential density surfaces
- ./lib/deltasurf.m             - create specific volume anomaly or in-situ density anomaly surfaces


## Installation:
Run the following commands in MATLAB, replacing 
  ~/work/dphil/projects/topobaric_surface/
with the path to this README.md file
```
>> cd('~/work/dphil/projects/topobaric_surface/run')
>> topobaric_surface_install();
```

## Manual Installation:
If the above installation function does not work, try the following manual installation instructions.

Note: the Windows instructions have not been properly tested. 
Please contact the author if there are any problems (with Windows or otherwise).

Below, PATH_TOPOBARIC_SURFACE is a generic stand-in for the path to this file.

If you will use a bash shell, start by entering the following command,
modified to use the path to the location of this file, in a bash shell: 
```
  PATH_TOPOBARIC_SURFACE=~/work/dphil/projects/topobaric_surface
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
Note, this step is optional, but highly recommended for speed of execution.

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
  Add the following lines to the end of startup.m:
  ```
  % Add Topobaric Surface
  tmp_STARTUP_PATH = pwd();
  cd(['$PATH_TOPOBARIC_SURFACE' filesep 'src']);
  topobaric_surface_add_to_path();
  cd(tmp_STARTUP_PATH);
  clearvars tmp_STARTUP_PATH
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

printf PATH_TOPOBARIC_SURFACE:$PATH_TOPOBARIC_SURFACE: > $PATH_TOPOBARIC_SURFACE/run/paths.txt
printf ECCO2:$ECCO2: >> $PATH_TOPOBARIC_SURFACE/run/paths.txt
printf OCCA:$OCCA: >> $PATH_TOPOBARIC_SURFACE/run/paths.txt
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

Then create a text file $PATH_TOPOBARIC_SURFACE/run/paths.txt containing the following text
```
PATH_TOPOBARIC_SURFACE~/work/dphil/projects/topobaric_surface:ECCO2:~/work/data/ECCO2/latlon/:OCCA:~/work/data/OCCA/:
```
replacing the local paths as appropriate.


## Usage:
./run/examples.m gives examples to create a topobaric surface, a topobaric geostrophic streamfunction, the orthobaric Montgomery potential, and other geostrophic streamfunctions from the literature.

./run/run_ECCO2.m generates most figures in “Neutral surface topology” and “The exact geostrophic streamfunction for neutral surfaces”.

./run/run_OCCA.m generates Figure 3 in “Neutral surface topology”.

./run/pitch.m generates Figure B.7 in “Neutral surface topology”.

Note, these were originally run on MATLAB 2017b and 2018b on a Mac. 
Running on Linux is known to produce cosmetic colour and font differences.



## MATLAB software:
Additional MATLAB software from the MATLAB File Exchange is provided in ./fex/

* binsrchn  by  Geoff Stanley
    https://www.mathworks.com/matlabcentral/fileexchange/70108
* bisectguess  by  Geoff Stanley
    https://www.mathworks.com/matlabcentral/fileexchange/69710
* bfs, scomponents in the GAIMC toolbox  by  David F. Gleich
    https://www.mathworks.com/matlabcentral/fileexchange/24134
* catstruct  by  Jos van der Geest
    https://www.mathworks.com/matlabcentral/fileexchange/7842
* CC2periodic  by  Geoff Stanley
    https://www.mathworks.com/matlabcentral/fileexchange/66079
* columncalculus  by  Geoff Stanley
    https://www.mathworks.com/matlabcentral/fileexchange/69713
* % distinguishable_colors  by   Tim Holy
    https://www.mathworks.com/matlabcentral/fileexchange/29702
* % export_fig  by  Yair Altman
    https://www.mathworks.com/matlabcentral/fileexchange/23629
* PriorityQueue  by  Andrew Woodward
    https://www.mathworks.com/matlabcentral/fileexchange/69142
* splinefit, ppint  by  Jonas Lundgren
    https://www.mathworks.com/matlabcentral/fileexchange/13812
* % subaxis  by  Aslak Grinsted
    https://www.mathworks.com/matlabcentral/fileexchange/3696
* % textborder  by  Joao Henriques
    https://www.mathworks.com/matlabcentral/fileexchange/27383
* % wprctile  by  Durga Lal Shrestha
    https://www.mathworks.com/matlabcentral/fileexchange/16920

% optional: not required to compute topobaric surfaces, but used by scripts in ./run/

Note, a small modification has been made to line 25 of GAIMC’s bfs.m, which originally read:
  if ~exist('target','var') || isempty(full), target=0; end
It now reads:
  if nargin < 3 || isempty(target), target = 0; end


## ReCon
Topobaric Surface uses ReCon to compute the Reeb Graph.

ReCon is available at
http://vgl.serc.iisc.ernet.in/software/software.php?pid=003

For more information on ReCon, see
Doraiswamy, H. & Natarajan, V. Computing Reeb Graphs as a Union of Contour Trees. IEEE Transactions on Visualization and Computer Graphics 19, 249–262 (2013).

The ReCon code included with Topobaric Surface has been modified to work with double precision and to allow the input and output to be from memory rather than files on the hard disk.

Note that ReCon requires Java 1.5 or higher, though it must be built for the version of the Java Runtime Environment that is packaged with MATLAB.  


## Copyright:
Copyright 2019 Geoff Stanley

This file is part of Topobaric Surface.

Topobaric Surface is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

Topobaric Surface is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Topobaric Surface.  If not, see
<https://www.gnu.org/licenses/>.

Author(s) : Geoff Stanley

Email     : g.stanley@unsw.edu.au 

Email     : geoffstanley@gmail.com

Version   : 1.0
