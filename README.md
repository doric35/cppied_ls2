## Cite this paper
If you find this code useful or use it in your work, please cite the following paper:\
Dominik Richard, Michael Morin, Claude-Guy Quimper (2024), Local Search for Coverage Path Planning with Imperfect Extended Detection\
(submitted for review to ICRA 2025).

## Contact info
dominik.richard.1@ulaval.ca

## Dataset
We provide our configuration dataset in the directory data/config.\
We provide our benchmarked dataset in the directory data/test.\
The name of each subdirectory is an instance id in the format: s{n}{n}_{s,r}{id1}_lrc{0,1}{mod_id}{r}\
n: instance size; {s,r}: structured(s) or noisy(r); id1: instance id unique relative to size and type;\
id2: Lateral range curve model; mod_id; Lateral range curve model approximation id; r: lateral range;

## Config file
The input of the executables are toml config files structured as:\
\[problem_parameters\] \
problem_file = "./data/test/instance_id/cppied_problem.txt" # Seabed encoding source\
solution_file = "/data/test/instance_id/cppied_solution.txt" # Solution encoding target\
coverage_file = "/data/test/instance_id/cppied_coverage.txt" # Coverage encoding target\
pod_file = "/data/test/instance_id/cppied_pod.txt" # Lateral range curve approximation source\
required_coverage = 0.9 # w_{min}\
tolerance = 0.01 # For numerical stability\
id = "instance_id" # Instance id\
\
\[algorithm_parameters\]\
solver_temp_files_dump = "./" # Where to store intermediate files\
log_directory = "./" # Where to store the run log\
solver_executable = "./External/concorde/LINKERN/linkern" # TSP solver executable; Use dps_NN for nearest neighbor; dps_NN3Opt for Nearest neighbor with a 3opt iteration;\
save_solution = true # Either to store solution in targets\
seed = 0 # Random number generator seed (Mersenne Twister 19937)\
uuid = "icra1" # Run unique id\
n_perturbations = 15 # k_p\
n_improvements = 10 # I\
max_time = 600 # Max allowed time in seconds\ 
neighborhood_width = 2 # W\
n_neighbors = 9 # k_i\
epsilon = 0.5 # See the paper for algorithm parameters details.\

## Dependencies
[Concorde solver](https://www.math.uwaterloo.ca/tsp/concorde.html)\
[Qsopt](https://www.math.uwaterloo.ca/~bico/qsopt/)\
[Cplex](https://www.ibm.com/products/ilog-cplex-optimization-studio)\

## How to run
Compile using cmake with the CMakeList.txt (include -O3 option) in "build" directory.\
Once the config files are formated, run:\
mkdir build\
cmake ..\
make\
cd ..\
./build/initialization_algorithm init_config.toml\
./build/local_search_algorithm ls_config.toml
