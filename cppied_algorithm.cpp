#include "cppied_algorithm.hpp"
void cppied_algorithm::init(){
    pad();
    set_lexicographic_cost();
    for (auto &i: mPod) {
        for (auto &j: i)
            if (j > mRequiredCoverage + cFLAT_EPSILON) j = mRequiredCoverage + cFLAT_EPSILON;
    }
}

void cppied_algorithm::pad() {
    csp::active_node pad_node {mPadSize, mPadSize, 0};
    mGridDimensions.first += 2*mPadSize; mGridDimensions.second += 2*mPadSize;
    for (auto& node:mSolution) node += pad_node;
    mStartNode += pad_node;
    mSeabed=utils::pad_matrix(mSeabed, mPadSize, (unsigned short )0);
    mCoverage=utils::pad_matrix(mCoverage, mPadSize, mPadValue);
}

void cppied_algorithm::unpad() {
    csp::active_node pad_node {mPadSize, mPadSize, 0};
    mGridDimensions.first -= 2*mPadSize; mGridDimensions.second -= 2*mPadSize;
    for (auto& node:mSolution) node-=pad_node;
    mStartNode -= pad_node;
    utils::unpad_matrix(mSeabed, mPadSize);
    utils::unpad_matrix(mCoverage, mPadSize);
}

cppied_algorithm::~cppied_algorithm()=default;

void cppied_algorithm::reset_coverage()
{
    utils::unpad_matrix(mCoverage, mPadSize);
    std::fill(mCoverage.begin(), mCoverage.end(),
              std::vector<double>(mCoverage[0].size(), 0.0));
    mCoverage = utils::pad_matrix(mCoverage, mPadSize, mPadValue);
    for (auto node_it= mSolution.begin(); node_it!=mSolution.end(); node_it++){
        set_insertion_coverage(node_it, std::next(node_it));
    }
}

void cppied_algorithm::solve_tsp(const std::optional<std::string>& pSolverCommand)
{
    write_tsp_instance();
    int err;
    if (pSolverCommand.has_value()) err = std::system(pSolverCommand.value().c_str());
    else err = std::system((mSolverPath + " -o " + mTspPathPrefix + ".sol " + mTspPathPrefix + ".tsp").c_str());
    if (err) throw std::runtime_error("Could not solve the tsp instance " + mTspPathPrefix);
    read_tsp_solution();
}

void cppied_algorithm::push_shortest_path(const csp::active_node& pDestination)
{
    csp::active_node current = mSolution.back();
    csp::active_node next{}, next_best{};
    unsigned best_cost, next_cost;
    do {
        best_cost=cINF_EDGE;
        for(auto& neighbor:cNEIGHBORS[current.dir]){
            next = current + neighbor;
            if (valid_node(next)) {
                next_cost = get_distance(next, pDestination, true) + (current.dir != next.dir) + cLVALUE;
                if (next_cost < best_cost) {next_best = next; best_cost = next_cost;}
            }
        }
        mSolution.push_back(next_best); current = mSolution.back();
    } while (current.x != pDestination.x || current.y != pDestination.y || current.dir != pDestination.dir);
}

std::string cppied_algorithm::get_tsp_options(){
    std::string tsp_options;
    if (mSolverPath.find("linkern") != std::string::npos){
        mTspSolType = "tour";
        tsp_options = "-s " + std::to_string(mSeed);
    } else if (mSolverPath.find("edgegen")!= std::string::npos)
        throw std::runtime_error("Cppied solver does not support edgegen yet.\n Must be linkern, concorde or dps_NN for initialization.\n");
    else if (mSolverPath.find("concorde") != std::string::npos) {
        mTspSolType = "tour";
        tsp_options = "-f -s " + std::to_string(mSeed);
    } else throw std::runtime_error("Cppied solver does not support the tsp solver: " +
                                    mSolverPath + " must be linkern, edgegen or concorde.\n");
    return  tsp_options;
}
