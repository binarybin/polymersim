# System requirements:
1.The program is run in macOS High Sierra and Catalina, with most recent 10.15.3 Beta version of macOS Catalina
2.Check github repo. Each branch stands for:
* master->polymer and polymer interaction
* rubisco->polymer and rigid 4x2 interaction
* square->polymer and rigid 2x2 square interaction



# Installation guide:
use Makefile to compile polymersim

run the file with appropriate parameters and reference tasklist files

tasklist.txt specifies the parameters used in the simulated annealing process

for example: ./polymersim 262 0 125 4 4 8 50 50 tasklist.txt scan empty. Check main.cpp to see the meaning of each input parameter. 

The expected output of each step is a file describes the occupation status of each point on the lattice model. Then we use analysis program to compute the cluster size and etc. 

Run time could vary by computing power and the input parameters. For example, more polymers and bigger space would result in longer running time. Usually it takes 10+ hours to run a rigid polymer simulation and less time to run flexible polymer simulations.

Link to the 3D simulation repo: https://github.com/BenjaminWeiner/magic-numbers


