This repository holds the materials needed for our hands-on GraphBLAS tutorial.  

Clone the repo and test that everything works by going into the src directory

     cd src

And then build our basic graph program

     make
 
You can run this by typing

     ./BuildGraph.exe
     
We include in this repo pre-built versions of the SuiteSparse library for 
OSX and for Linux.   I'm sorry if you use Windows.   The creators of this tutorial
do not have Windows systems so we could not build/text the tutorial on Windows.

You can go to the SuiteSparse web site 

       http://faculty.cse.tamu.edu/davis/suitesparse.html

And download the source and build for you won system if our prebuild OSX and Linux
versions do not cover your own system.  Reverse engineering our makefile for a 
windows or other non-Linux like system should be straightforward.  The software
dependencies are by design quite simple.  Just build the code, link the suiteSparse
libries, and reference the suiteSparse and our own tutorial-include-file.






