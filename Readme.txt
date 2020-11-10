Thank you for downloading this file. 

Here we describe the procedure of executing the slow update exact stochastic simulation algorithm (SUESSA),
slow update exact sorting stochastic simulation algorithm (SUESSSA) and block searched stochastic simulation algorithm (BlSSSA). 

The description of the necessary files are given below.

xyz.txt : This file contains the stoichiometry matrix of each network. 
k.txt : This file contains the reaction rates. 
conc.txt : This file contains the population of each species. 
outputfile.csv : This file contains the output of the simulation in terms of the population of the species. 



1. Colloidal aggregation model with 100 species and 5000 reactions.
2. B cell receptor signaling network with 1122 species and 24388 reactions.
3. FceRI signaling network with 380 species and 3862 reactions. 
4. Linear chain model with 100 species and 99 reactions. 
5. Kinetic model of 1,3-Butadiene Oxidation. This network contains 92 species and 613 reactions. The 
database given in "http://ignis.usc.edu/Mechanisms/USC-Mech%20II/USC_Mech%20II.htm" contains 95 species and 615 reactions for this networks. Therefore, this network will be executed with 95 no. of species and 615 reactions. 


#Java must have been installed in the system before simulating the algorithms.


#The procedure for executing the colloidal aggregation model in SUESSA. 

Step1: Open the folder named code. Then, open the folder named Colloidal aggregation model.
Step2: Copy the files named xyz.txt, conc.txt and k.txt and paste them in any directory, say bin for example, the directory bin may have the path c:\user\bin
Step3: Copy the file named SUESSA.jar to the directory bin having path c:\user\bin
Step4: In command line, move to the path c:\user\bin
Step5: Execute the command given by
c:\user\bin> java -jar SUESSA.jar
Enter the values of the number of species and reactions 
100
5000
Reading File from Java code
Input the fluctuation interval for species populations. For example 0.1 for 10%
0.1
Enter total number of simulation steps
10000000
If you want to write the output in file then, press 1 and if not then, press 0
1
Enter the number of simulation steps you want to be written in the file. Steps must be less than total no. of simulation steps 
100
 Search time 6047ms
 Update time 1055ms
 Total simulation time 8855ms




#The procedure for executing the B cell receptor signaling network in SUESSSA.

Step1: Open the folder named code. Then, open the folder named B cell receptor signaling network.
Step2: Copy the files named xyz.txt, conc.txt and k.txt and paste them in any directory, say bin for example, the directory bin may have the path c:\user\bin
Step3: Copy the file named SUESSSA.jar to the directory bin having path c:\user\bin
Step4: In command line, move to the path c:\user\bin
Step5: Execute the command given by
c:\user\bin> java -jar SUESSSA.jar
Enter the values of the number of species and reactions 
1122
24388
Reading File from Java code
Input the fluctuation interval for species populations. For example 0.1 for 10%
0.1
Enter total number of simulation steps
10000000
If you want to write the output in file then, press 1 and if not then, press 0
0
 Search time 4645ms
 Update time 695ms
 Total simulation time 7394ms


#The procedure for executing the Kinetic model of 1,3-Butadiene Oxidation in BlSSSA.

Step1: Open the folder named code. Then, open the folder named kineticmodelofbutadineoxidation.
Step2: Copy the files named xyz.txt, conc.txt and k.txt and paste them in any directory, say bin for example, the directory bin may have the path c:\user\bin
Step3: Copy the file named BlSSSA.jar to the directory bin having path c:\user\bin
Step4: In command line, move to the path c:\user\bin
Step5: Execute the command given by
c:\user\bin> java -jar BlSSSA.jar
Enter the values of the number of species and reactions 
95
615
Reading File from Java code
Input the fluctuation interval for species populations. For example 0.1 for 10%
0.1
Enter total number of simulation steps
10000000
If you want to write the output in file then, press 1 and if not then, press 0
0
 Intervals = 9

 Search time 1638ms
 Update time 528ms
 Total simulation time 3433ms

Thank you. 

