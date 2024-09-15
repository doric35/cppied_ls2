#include "../initialization_algorithms/dp_sweeper.hpp"

int main(int argc, char *argv[]){
    if (argc != 2){
        throw std::runtime_error("Need at least configuration file.");
    }
    auto config = toml::parse_file(argv[1]);
    std::string problem_path = config["problem_parameters"]["problem_file"].value_or("");
    std::string pod_path = config["problem_parameters"]["pod_file"].value_or("");
    std::string coverage_path = config["problem_parameters"]["coverage_file"].value_or("");
    std::string solution_path = config["problem_parameters"]["solution_file"].value_or("");
    std::string tsp_path_prefix = config["algorithm_parameters"]["solver_temp_files_dump"].value_or("../temp_files/");
    std::string uuid = config["algorithm_parameters"]["uuid"].value_or("");
    std::filesystem::path tsp_base = tsp_path_prefix;
    tsp_base = tsp_base / "tsp_instance.tsp"; std::string tsp_base_str = tsp_base.string();
    std::string id = config["problem_parameters"]["id"].value_or("");
    std::string log_path_prefix = config["algorithm_parameters"]["log_directory"].value_or("../log/");
    std::string solver_path = config["algorithm_parameters"]["solver_executable"].value_or("");
    std::vector<std::vector<unsigned short>> problem;
    utils::read_matrix<unsigned short>(problem_path, problem);
    std::vector<std::vector<double>> coverage;
    coverage = std::vector<std::vector<double>>(problem.size(), std::vector<double>(problem[0].size(), 0.0f));
    std::vector<std::vector<double>> pod;
    utils::read_matrix<double>(pod_path, pod);
    std::vector<csp::active_node> solution {{0, (unsigned short)pod[0].size(), 2}};
    csp::active_node start_point = solution[0];
    init::dp_sweeper initializer(coverage,
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
    std::cout << "Initializing CPPIED Solution with DpSweeper" << std::endl;
    int err = utils::solve(initializer, log_path_prefix, uuid);
    if (err !=0) throw std::runtime_error("Could not complete initialization, return code " + std::to_string(err));
    std::cout << "Initial Solution Cost: " << initializer.get_lexicographic_cost() << " found in " <<
                    initializer.get_resolution_time<std::chrono::milliseconds>() << "milliseconds" << std::endl;
    utils::write_struct<csp::active_node>(solution_path, initializer.get_solution(), 3);
    utils::write_matrix<double>(coverage_path, initializer.get_coverage());
    initializer.info[0].emplace_back("id"); initializer.info[0].emplace_back("uuid");
    initializer.info[1].push_back(id); initializer.info[1].push_back(uuid);
    std::filesystem::path base_path = log_path_prefix;
    std::filesystem::path log_name = "dp_sweeper.csv";
    log_name = log_path_prefix / log_name;
    utils::write_info(log_name.string(), initializer.info);
    std::cout<< "done" << std::endl;
    return 0;
}
