#include "../local_search_algorithms/ls2.hpp"

int main(int argc, char *argv[]){
    if (argc != 2)
        throw std::runtime_error("Need a configuration file.");
    auto config = toml::parse_file(argv[1]);
    std::string problem_path = config["problem_parameters"]["problem_file"].value_or("");
    std::string pod_path = config["problem_parameters"]["pod_file"].value_or("");
    std::string coverage_path = config["problem_parameters"]["coverage_file"].value_or("");
    std::string solution_path = config["problem_parameters"]["solution_file"].value_or("");
    std::string opt_path = config["problem_parameters"]["opt_file"].value_or("");
    std::string solver_path = config["algorithm_parameters"]["solver_executable"].value_or("");
    std::string tsp_path_prefix = config["algorithm_parameters"]["solver_temp_files_dump"].value_or("../temp_files/");
    std::string log_path_prefix = config["algorithm_parameters"]["log_directory"].value_or("../log/");
    std::string uuid = config["algorithm_parameters"]["uuid"].value_or("");
    std::string id = config["problem_parameters"]["id"].value_or("");
    bool save_solution = config["algorithm_parameters"]["save_solution"].value_or(true);
    std::filesystem::path tsp_base = tsp_path_prefix;
    tsp_base = tsp_base / "tsp_instance.tsp"; std::string tsp_base_str = tsp_base.string();
    std::vector<std::vector<unsigned short>> problem;
    utils::read_matrix<unsigned short>(problem_path, problem);
    std::vector<std::vector<double>> coverage;
    utils::read_matrix<double>(coverage_path, coverage);
    std::vector<std::vector<double>> pod;
    utils::read_matrix<double>(pod_path, pod);
    std::vector<csp::active_node> solution;
    utils::read_struct<csp::active_node>(solution_path, solution, 3);
    csp::active_node start_point = solution[0];
    local_search::ls2 optimizer(config["algorithm_parameters"]["n_perturbations"].value_or(10),
                                config["algorithm_parameters"]["n_improvements"].value_or(300),
                                config["algorithm_parameters"]["max_time"].value_or(3600),
                                config["algorithm_parameters"]["neighbourhood_width"].value_or(1),
                                config["algorithm_parameters"]["n_neighbors"].value_or(3600),
                                coverage,
                                pod,
                                problem,
                                solution,
                                config["problem_parameters"]["required_coverage"].value_or(0.9),
                                config["problem_parameters"]["tolerance"].value_or(0.001),
                                pod[0].size(),
                                start_point,
                                solver_path,
                                tsp_base_str,
                                uuid,
                                pod[0].size(),
                                1.0,
                                config["algorithm_parameters"]["seed"].value_or(0));
    std::cout << "Improving solution with LS2 heuristic" << std::endl;
    int err = utils::solve(optimizer, log_path_prefix, uuid);
    if (err !=0) throw std::runtime_error("Could not complete local search, return code " + std::to_string(err));
    std::cout << "Solution Cost: " << optimizer.get_lexicographic_cost() << " found in " <<
                optimizer.get_resolution_time<std::chrono::milliseconds>() << "milliseconds" << std::endl;
    if (save_solution) {
        utils::write_struct(solution_path, optimizer.get_solution(), 3);
        utils::write_matrix<double>(coverage_path, optimizer.get_coverage());
    }

    optimizer.info[0].emplace_back("id"); optimizer.info[0].emplace_back("uuid");
    optimizer.info[1].push_back(id); optimizer.info[1].push_back(uuid);
    std::filesystem::path base_path = log_path_prefix;
    std::filesystem::path log_name = base_path / "ls2.csv";
    utils::write_info(log_name.string(), optimizer.info);
    log_name = base_path / "ls2_solver.csv";
    utils::write_info(log_name.string(), optimizer.solver_info);

    return 0;
}


