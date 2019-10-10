Recon - A fast algorithm to compute Reeb graphs
===============================================

Recon is a library to compute Reeb graphs. It implements the algorithm outlined in the paper 

Harish Doraiswamy and Vijay Natarajan, "Computing Reeb Graphs as a Union of Contour Trees", IEEE Transactions on Visualization and Computer Graphics, 19(2), 2013, 249--262.

Recon is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

This folder contains the following files / folders:
1. src folder
2. loaders.xml
3. computeReebGraph.sh
4. input.properties
5. Eclipse project files, in case you want to import project into Eclipse
6. build.xml
7. license
8. copyright
9. and this Readme

Requirements
------------
This program requires **Java 1.5** or higher.
To compile this code using the given build.xml Apache Ant is required. (tested on Apache Ant 1.8.0)

Compiling the code
------------------
Compile the code by running the "ant" command
This will create recon.jar file in the build folder.

Usage
----- 
run the computeReebGraph.sh after setting the following in the input.properties file
* loader - Specifies the input mesh loader type. It should be one of the loaders mentioned in loaders.xml. Currently OFF, TET and SIM are supported
* inputFile
* inputFunction - Input function should be 0, or i, where i is the co-ordinate index (to be used for the height function along the ith axis)
* output - The output file (optional) to store the Reeb graph.

Input
------
The library currently supports the following three formats for the input mesh:

OFF
***
1. Optional first line containing "OFF"
2. Next line specifies the no. of vertices (nv) followed by the number of triangles (nt) (space seperated)
3. The next nv lines contains
   x y z [f]
   where x, y & z specify the co-ordinates of the vertex and f specifies the function value. (If the input type is not f, then the function value is optional)
4. the next nt lines has 
   [3] v1 v2 v3 
   where v1, v2 and v3 are the vertex indices of the vertices that form the triangles (the 3 is optional)


TET
***
1. First line specifies the no. of vertices (nv) followed by the number of tetrahedrons (nt) (space seperated)
2. The next nv lines contains
   x y z [f]
   where x, y & z specify the co-ordinates of the vertex and f specifies the function value. (If the input type is not 0, then the function value is optional)
3. the next nt lines has 
   v1 v2 v3 v4
   where v1, v2, v3 and v4 are the vertex indices of the vertices that form the tetrahedron.


SIM
***
1. First line specifies the dimension (d) of the input
2. The next line specifies the no. of vertices (nv) followed by the number of simplices (ns) (space seperated)
3. The next nv lines contains
   c1 c2 ... cd [f]
   where ci specifies the ith co-ordinate of the vertex and f specifies the function value.
4. The next ns lines has
   (l + 1) v1 v2 ... v{l+1}
   where l is the dimension of the simplex and vi is the index of the ith vertex of the simplex.



Extending code to support other file formats
---------------------------------------------

1. You need to write a loader that extend the interface vgl.iisc.external.loader.MeshLoader
2. Add this loader source in the src/meshloader folder in the appropriate package
3. Register this new loader by providing its name and class information in the loaders.xml file present in the root folder.
4. You can now use this loader by providing its registered name in the input.properties file.
See the javadoc for more details regarding the Loader interface. 


In case of any errors
---------------------
1. If you get a OutOfMemory (java heap) exception, try increasing the memory allocated to the jvm in the run.sh file (there is a -Xmx parameter) and run again.
2. Any other error/exception, kindly let me know of the error. It would be great if you can provide
   the stack trace of the exception (if it is an exception that has occured) along with the input data.



