<html>
<h2>2 Site Model</h2>
<p> The files contained in the directory are for the simulation of the tight binding hubbard model with interactions for an ensemble of 2 site systems. The program will calculate the density of states (DOS), generalized inverse participaion ratio (GIPR) and the filling of the DOS.</p>
<h3>Compilling and Running</h3>
<h4>Basic</h4>
<p> Download all the files to a directory or use your forked repository (information on how to do this in README in AH-onsite directory). If you plan on modifying the files and wanting the update the files on GitHub we recommend forking. </p>
<p>In the 2site folder type the command "<em>make</em>". This will compile the executable. To get remove auxiliary type "<em>make clean</em>". The executable is called <em>main.e</em> and it can be run using the command "<em>./main.e</em>". </p>
<p>To modify the input variables open main.f90 in a text editor. In the section labeled input parameters (line 25-33) modify the variables. The description of each variable can be found in the comments beside it.
<h4>Advanced</h4>
<p>The compiller specified in the makefile is gfortran if you would like to use another compiller open the makefile in a text editor and set the variable CC to the compiller of your choice. To modify the optimization of the code or add additional flags for debugging or other purposes add them to the CFLAGS variable in the makefile. Optimization -O0 is lowest and -O3 is highest. </p>
<p>The only library that needs to be linked is the LAPACK library (linear algebra package). If not already owned it can be downloaded <a href="http://www.netlib.org/lapack/#_lapack_version_3_5_0">here</a>.</p>
<h3>Data and Graphs</h3>
Previous data and graphs that have been using this code can be found in the <em>data</em> directory. The graphs are stored as .agr files and can be opened using xmgrace and some as also saved in pdf format. The files have names based on the parameters of the simulation that created them. </p>
<p>On the top of the data files more information regarding the simulation can be found as well as the filling of the DOS.</p>
</html>
