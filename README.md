## Cite this paper
If you find this code useful or use it in your work, please cite the following paper:
Dominik Richard, Michael Morin, Claude-Guy Quimper (2024), Local Search for Coverage Path Planning with Imperfect Extended Detection (submitted for review to ICRA 2025).

## Contact info
dominik.richard.1@ulaval.ca

## Dataset
We provide our benchmarked dataset in the directory data/testdata.

## Config file
The input of the executables are toml config files tructued as:
insert here example

## Dependencies
Insert here the link to concorde solver.
Insert here the link to qsopt.
Insert here the link to cplex.

## How to run
Compile using cmake with the CMakeList.txt provided in abuild directory.
./build/initialization_algorithm init_config.toml
./build/local_search_algorithm ls_config.toml
