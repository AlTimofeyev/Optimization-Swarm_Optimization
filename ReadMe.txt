******************************************************************************************************************
Author:	Al Timofeyev
Date:	May 19, 2019
Desc:	This is how to compile and execute the code.
	This is an implementation of 3 Swarm Algorithms: Particle Swarm Optimization, Firefly Algorithm, and
	Harmony Search. Each of these algorithms are modular and easily importable to other projects. The
	main.cpp file provided with the source code is an example of how these algorithms can be used/implemented;
	either through file reading or inputing the parameters manually in the code.
******************************************************************************************************************

**********************************************************************************
NOTES TO PROFESSOR (if any) ARE AT THE VERY BOTTOM OF THIS README!!
**********************************************************************************


**********************************************************************************
------------------------------- CODING ENVIRONMENT -------------------------------
Windows 10
CLion version 2019.1.3
CMake version 3.14.3 (bundled with CLion)
cygwin version 3.0.4
cygwin GDB version 8.1.1
gcc version 7.4.0 
g++ version 7.4.0

******
NOTE:
1)	CLion generated a CMakeLists.txt file included with the source code.
	cmake_minimum_required(VERSION 3.13)
2)	The program was written in C++.
******
**********************************************************************************


**********************************************************************************
--------------------------- SETUP CONFIGURATION FILES ----------------------------
----------------------------------------------------------------------------------
---- Structure of Particle Swarm Optimization configuration file

<Dimensions>			---- The number of dimensions.
<Population Size>		---- The size of the population.
<Max Number of Iterations>	---- The maximum number of iterations.
<k>				---- The K Dampening Factor used to dampen the velocities.
<c1>				---- The c1 Scaling Factor (used to bring particle closer to personal best).
<c2>				---- The c2 Scaling Factor (used to bring particle closer to global best).
<list of parameters>		---- The list of function IDs and upper and lower bounds.
Use only a comma (,) delimiter, no spaces between values.

-- Example:
30		---- Number of Dimensions
500		---- Population size
500		---- Max number of iterations
0.9		---- k
0.3		---- c1
0.7		---- c2
1,-500,40	---- Function ID 1, with -500 to 40 min/max bounds.
5,-32,100	---- Function ID 5, with -32 to 100 min/max bounds.
8,0,pi		---- Function ID 8, with 0 to pi min/max bounds.


----------------------------------------------------------------------------------
----------------------------------------------------------------------------------
---- Structure of Firefly Algorithm configuration file

<Dimensions>			---- The number of dimensions.
<Population Size>		---- The size of the population.
<Max Number of Iterations>	---- The maximum number of iterations.
<alpha>				---- The alpha scaling factor.
<beta>				---- The beta scaling factor.
<gamma>				---- The gamma scaling factor.
<list of parameters>		---- The list of function IDs and upper and lower bounds.
Use only a comma (,) delimiter, no spaces between values.

-- Example:
30		---- Number of Dimensions
500		---- Population size
500		---- Max number of iterations
0.1		---- alpha
0.2		---- beta
1.0		---- gamma
1,-500,40	---- Function ID 1, with -500 to 40 min/max bounds.
5,-32,100	---- Function ID 5, with -32 to 100 min/max bounds.
8,0,pi		---- Function ID 8, with 0 to pi min/max bounds.

----------------------------------------------------------------------------------
----------------------------------------------------------------------------------
---- Structure of Harmony Search configuration file

<Dimensions>			---- The number of dimensions.
<Population Size>		---- The size of the population.
<Max Number of Iterations>	---- The maximum number of iterations.
<HMCR>				---- The Harmony Memory Consideration Rate.
<PAR>				---- The Pitch Adjustment Rate.
<bandwidth>			---- The bandwidth range (pitch scaling factor).
<list of parameters>		---- The list of function IDs and upper and lower bounds.
Use only a comma (,) delimiter, no spaces between values.

-- Example:
30		---- Number of Dimensions
500		---- Population size
500		---- Max number of iterations
0.9		---- HMCR
0.4		---- PAR
0.2		---- bandwidth
1,-500,40	---- Function ID 1, with -500 to 40 min/max bounds.
5,-32,100	---- Function ID 5, with -32 to 100 min/max bounds.
8,0,pi		---- Function ID 8, with 0 to pi min/max bounds.

******
NOTE:
1)	Depending on which IDE you are running, configuration files should be either
	in the same folder as source code or in build folder.
2)	Configuration files can be passed as command line parameters or use the default
	configuration file (just alter the default files as needed).
******
**********************************************************************************


**********************************************************************************
------------------------------ COMPILE AND EXECUTE -------------------------------
---- To compile for an IDE project.
To Compile:
You could use CMake to compile CMakeLists.txt file that's included with source code.

To Execute:
run main.cpp


---- I'm assuming it could also be compiled and run from command line:
To Compile:
g++ -o main main.cpp

To Execute:
./main				---- Default psConfig.txt, faConfig.txt, or hsConfig is used as initialization file.
./main frogconfig.txt		---- frogconfig.txt is initialization file example.
./main configFile2.txt		---- configFile2.txt is initialization file example.
./main blabla.txt		---- blabla.txt is initialization file example.

**** NOTE:
1)	The code in main needs to be slightly modified depending on which of the
	swarm algorithms you are planning to use.
**********************************************************************************


**********************************************************************************
------------------------------ NOTES TO PROFESSOR --------------------------------
1)	The Swarm Algorithm tester functions in main.cpp have been set to only
	print and save the analysis of each algorithm. If you'd like to see/save
	the actual results, please uncomment code in their respected tester functions.

2)	Analysis of population stagnation will only be performed on the first
	benchmark function f(1) of the Particle Swarm (if it's included in the report).
**********************************************************************************