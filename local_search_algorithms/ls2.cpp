#include "ls2.hpp"

/*****************************Initialization Methods*****************************/
void local_search::ls2::generate_nbh()
{
    short min_dx, max_dx, min_dy, max_dy;
    max_dy = (short)std::min(mNeighbourhoodWidth+1, 2*mLateralRange);
    min_dy = (short)std::max(-mNeighbourhoodWidth, -2*mLateralRange+1);
    for (short dy = min_dy; dy<max_dy; dy++){
        mLateralNbh.push_back({0, dy, 2});
        mLateralNbh.push_back({0, dy, 3});
    }
    max_dx = (short)std::min(mNeighbourhoodWidth+1, mLateralRange+1);
    min_dx = (short)std::max(-mNeighbourhoodWidth+1, -mLateralRange+1);
    max_dy = (short)std::min(mLateralRange,mNeighbourhoodWidth);
    min_dy = (short)-std::min(mLateralRange,mNeighbourhoodWidth);
    for (short dy = min_dy; dy<max_dy; dy++){
        for (short dx = min_dx; dx<max_dx; dx++){
            mLateralNbh.push_back({dx, dy, 0});
            mLateralNbh.push_back({dx, dy, 1});
        }
    }

    max_dx = (short)std::min(mNeighbourhoodWidth+1, 2*mLateralRange);
    min_dx = (short)std::max(-mNeighbourhoodWidth, -2*mLateralRange+1);
    for (short dx = min_dx; dx<max_dx; dx++){
        mVerticalNbh.push_back({dx, 0, 0});
        mVerticalNbh.push_back({dx, 0, 1});
    }
    max_dx = (short)std::min(mLateralRange, mNeighbourhoodWidth);
    min_dx = (short)-std::min(mLateralRange, mNeighbourhoodWidth);
    max_dy = (short)std::min(mLateralRange+1,mNeighbourhoodWidth+1);
    min_dy = (short)std::max(-mLateralRange+1,-mNeighbourhoodWidth+1);
    for (short dy = min_dy; dy<max_dy; dy++){
        for (short dx = min_dx; dx<max_dx; dx++){
            mVerticalNbh.push_back({dx, dy, 2});
            mVerticalNbh.push_back({dx, dy, 3});
        }
    }
}

void local_search::ls2::pad() {
    cppied_algorithm::pad();
    mMutableCoverage = mCoverage;
    mMutableSolution.clear();
    std::copy(mSolution.begin(), mSolution.end(), std::back_inserter(mMutableSolution));
    mXCoordHorizontalRNodes = std::uniform_int_distribution<int>(mPadSize, (int)mCoverage[0].size() - (1 + mPadSize));
    mYCoordHorizontalRNodes = std::uniform_int_distribution<int>(mPadSize, (int)mCoverage.size() - mPadSize);
    mXCoordVerticalRNodes = std::uniform_int_distribution<int>(mPadSize, (int)mCoverage[0].size() - mPadSize);
    mYCoordVerticalRNodes = std::uniform_int_distribution<int>(mPadSize, (int)mCoverage.size()-(1 +mPadSize));
}

void local_search::ls2::unpad() {
    cppied_algorithm::unpad();
    mMutableCoverage = mCoverage;
    mMutableSolution.clear();
    std::copy(mSolution.begin(), mSolution.end(), std::back_inserter(mMutableSolution));
    mXCoordHorizontalRNodes = std::uniform_int_distribution<int>(0, (int)mCoverage[0].size() - 1);
    mYCoordHorizontalRNodes = std::uniform_int_distribution<int>(0, (int)mCoverage.size());
    mXCoordVerticalRNodes = std::uniform_int_distribution<int>(0, (int)mCoverage[0].size());
    mYCoordVerticalRNodes = std::uniform_int_distribution<int>(0, (int)mCoverage.size()-1);
}

void local_search::ls2::init() {
    cppied_algorithm::init();
}

/*****************************Solver Methods*****************************/
void local_search::ls2::d_solve(){
    mStart = std::chrono::system_clock::now();
    bool optimal_improved;
    std::chrono::time_point<std::chrono::system_clock> solver_start, solver_end;
    unsigned n_bad_iter=0, iter = 0;
    auto [n_moves, n_turns] = get_bi_cost();
    //initial upper bound
    std::vector<csp::active_node> s_upper_bound=mSolution;
    solver_info.push_back({std::to_string(iter++), std::to_string(mSolutionCost), std::to_string(n_moves),
                           std::to_string(n_turns),
                           std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(mEnd - mStart).count())});
    do {
        optimal_improved = solve_perturb_and_optimize(iter);
        //A solution may be optimized in the solve and perturb wrt heuristics but wrt upper bound
        //This check leads to more time under perturbation
        optimal_improved = optimal_improved && !equal_solutions(s_upper_bound);
        if (optimal_improved){
            //save upper bound before recomputing a feasible path
            s_upper_bound = mSolution;
            n_bad_iter=0;
            optimize_path();
            mEnd = std::chrono::system_clock::now();
            std::tie(n_moves, n_turns) = get_bi_cost();
            solver_info.push_back({std::to_string(iter++), std::to_string(mSolutionCost), std::to_string(n_moves),
                                   std::to_string(n_turns),
                                   std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(mEnd - mStart).count())});
        } else {n_bad_iter++;iter++;}

    } while(std::chrono::duration_cast<std::chrono::seconds>(mEnd - mStart).count() < mMaxTime);// & n_bad_iter < mMaxBadIter);
    build_feasible();
}

void local_search::ls2::build_feasible(){
    std::vector<csp::active_node> solution_cp = mSolution;
    mSolution.clear();
    mSolution.push_back(solution_cp[0]);
    for (std::size_t i = 1; i< solution_cp.size(); i++)
        push_shortest_path(solution_cp[i]);
    //solution cost should not have changed
    reset_coverage();
}

bool local_search::ls2::equal_solutions(const std::vector<csp::active_node> &pSolution) const {
    if (pSolution.size() != mSolution.size()) return false;
    for (std::size_t i =0; i< mSolution.size(); i++)
        if (mSolution[i] != pSolution[i]) return false;
    return true;
}

void local_search::ls2::terminate() {
    unpad();
    //Info: Solver, cost, n_moves, n_turns, Total time;
    std::string solver, cost, n_moves, n_turns, total_time;
    solver = mSolverPath.substr(mSolverPath.find_last_of('/')+1);
    auto [ui_n_moves, ui_n_turns] = get_bi_cost();
    info.push_back({"solver", "cost", "n_moves", "n_turns", "total_time"});
    info.push_back({solver, std::to_string(get_lexicographic_cost()), std::to_string(ui_n_moves), std::to_string(ui_n_turns),
                    std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(mEnd-mStart).count())});
}

void local_search::ls2::optimize_path() {
    solve_lkh();
    reset_coverage();
    mMutableSolution.clear();
    std::copy(mSolution.begin(), mSolution.end(), std::back_inserter(mMutableSolution));
    mMutableCoverage = mCoverage;
    set_lexicographic_cost();
}

void local_search::ls2::set_lexicographic_cost() {
    cppied_algorithm::set_lexicographic_cost();
    mLocalCost = mSolutionCost;
}

bool local_search::ls2::restart(bool pPerturbate, bool pHeuristicImprovement) {
    bool improved = false;
    if (mLocalCost < mSolutionCost || pHeuristicImprovement){
        //update upper bound
        mSolutionCost = mLocalCost;
        improved = true;
        mSolution.clear();
        std::copy(mMutableSolution.begin(), mMutableSolution.end(), std::back_inserter(mSolution));
        mCoverage = mMutableCoverage;
    }
    else{
        //reset to upper bound
        mMutableSolution.clear();
        mLocalCost = mSolutionCost;
        std::copy(mSolution.begin(), mSolution.end(), std::back_inserter(mMutableSolution));
        mMutableCoverage = mCoverage;
    }
    if (pPerturbate)
        perturbation_procedure();
    return improved;
}

bool local_search::ls2::solve_perturb_and_optimize(std::size_t pIter){
    bool best_improved, improved, imp_improved, local_improved, two_opt_improved;
    best_improved= false;
    //std::chrono::time_point<std::chrono::system_clock> solver_start, solver_end;
    for (int i=0; i<mNImprovements; i++){
        imp_improved = false;
        do
        {
            improved = improvement_procedure();
            imp_improved = imp_improved || improved;
        } while (improved);
        mEnd = std::chrono::system_clock::now();
        two_opt_improved=two_opt_procedure();
        if (mLocalCost < mSolutionCost){
            mEnd = std::chrono::system_clock::now();
            solver_info.push_back({std::to_string(pIter), std::to_string(std::min(mSolutionCost, mLocalCost)), "",
                                   "",
                                   std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(mEnd - mStart).count())});
        }
        local_improved = restart(i<mNImprovements-1,
                                 (mLocalCost==mSolutionCost && (imp_improved || two_opt_improved)));
        best_improved = local_improved || best_improved;
    }
    return best_improved;
}

template <class ForwardIter>
unsigned local_search::ls2::set_mutable_extraction_coverage(ForwardIter begin, ForwardIter end)
{
    unsigned short i, j, n_unsat=0;
    for (; begin != end; begin++){
        if ((*begin).dir>1) {
            j = (*begin).x;
            for (i = (*begin).y-mLateralRange; i<(*begin).y + mLateralRange; i++){
                mMutableCoverage[i][j] -= map_pod(i, j, {.x = j, .y = (*begin).y, .dir=2});
                mMutableCoverage[i][j] /= (1.0 - map_pod(i, j, {.x = j, .y = (*begin).y, .dir=2}));
                n_unsat += mMutableCoverage[i][j] < (mRequiredCoverage - mTolerance);
            }
        } else {
            i = (*begin).y;
            for (j = (*begin).x-mLateralRange; j < (*begin).x+mLateralRange; ++j) {
                mMutableCoverage[i][j] -= map_pod(i, j,{.x = (*begin).x, .y = i, .dir=0});
                mMutableCoverage[i][j] /= (1.0 - map_pod(i, j, {.x = (*begin).x, .y = i, .dir=0}));
                n_unsat += mMutableCoverage[i][j] < (mRequiredCoverage - mTolerance);
            }
        }
    }
    return n_unsat;
}

void local_search::ls2::sequence_satisfy(const csp::active_node& n1, const csp::active_node& n2,
                                         const csp::active_node& ref, bool* response) {
    int i,j; double p_ij1;
    response[0]= true; response[1]= true; response[2]= true;
    if (ref.dir >1) {
        j = ref.x;
        for (i = ref.y - mLateralRange; i < ref.y + mLateralRange; i++){
            if (mMutableCoverage[i][j] < mRequiredCoverage - mTolerance) {
                p_ij1 = mMutableCoverage[i][j] + (1.0 - mMutableCoverage[i][j]) * map_pod(i, j, n1);
                response[0] = response[0] && node_cover(i, j, n1) &&
                              (p_ij1 > mRequiredCoverage - mTolerance);
                response[1] = response[1] && node_cover(i, j, n2) &&
                              (mMutableCoverage[i][j] + (1.0 - mMutableCoverage[i][j]) * map_pod(i, j, n2) >
                               mRequiredCoverage - mTolerance);
                response[2] = response[0] || response[1] ||
                              (node_cover(i, j, n1) && node_cover(i, j, n2) &&
                               (p_ij1 + (1.0 - p_ij1) * map_pod(i, j, n2) > mRequiredCoverage - mTolerance));
                if (!response[2]) break;
            }
        }
    } else{
        i = ref.y;
        for (j = ref.x - mLateralRange; j < ref.x + mLateralRange; j++){
            if (mMutableCoverage[i][j] < mRequiredCoverage - mTolerance) {
                p_ij1 = mMutableCoverage[i][j] + (1.0 - mMutableCoverage[i][j]) * map_pod(i, j, n1);
                response[0] = response[0] && node_cover(i, j, n1) &&
                        (p_ij1 > mRequiredCoverage - mTolerance);
                response[1] = response[1] && node_cover(i, j, n2) &&
                        (mMutableCoverage[i][j] + (1.0 - mMutableCoverage[i][j]) * map_pod(i, j, n2) >
                               mRequiredCoverage - mTolerance);
                response[2] = response[0] || response[1] ||
                        (node_cover(i, j, n1) && node_cover(i, j, n2) &&
                        (p_ij1 + (1.0 - p_ij1) * map_pod(i, j, n2) > mRequiredCoverage - mTolerance));
                if (!response[2]) break;
            }
        }
    }
}

bool local_search::ls2::improvement_procedure(){
    bool skip_iteration, best_improve= false, local_improve;
    unsigned n_uncovered, extraction_gain, count=0;
    double random_real, heuristic_weight;
    auto node_it = std::next(mMutableSolution.begin());
    do {
        n_uncovered = set_mutable_extraction_coverage(node_it, std::next(node_it));
        extraction_gain = get_distance(*std::prev(node_it), *node_it, false) +
                          get_distance(*node_it, *std::next(node_it), false) -
                          get_distance(*std::prev(node_it), *std::next(node_it), false);
        if (n_uncovered==0){
            node_it = mMutableSolution.erase(node_it);
            mLocalCost -= extraction_gain;
            best_improve = true; continue;
        }
        random_real = mRealGenerator(mRng);
        heuristic_weight = (((double)(n_uncovered))/((double)(2*mLateralRange)+cEPSILON));
        skip_iteration = ( random_real > heuristic_weight );
        if (extraction_gain==0 || skip_iteration) {
            set_mutable_insertion_coverage(node_it, std::next(node_it));
            node_it++; continue;
        }
        auto next = std::next(node_it);
        local_improve = local_nbh_search_procedure(node_it, std::next(node_it));
        best_improve = local_improve || best_improve;
        node_it = next;
    } while (node_it != std::prev(mMutableSolution.end()) && node_it!=mMutableSolution.end());
    return best_improve;
}

template <class ForwardIter>
bool local_search::ls2::local_nbh_search_procedure(ForwardIter begin, ForwardIter end)
{
    return random_nbh_search_procedure(begin, end);
}

template <class ForwardIter>
bool local_search::ls2::random_nbh_search_procedure(ForwardIter begin, ForwardIter end){
    csp::active_node node1{}, node2{};
    std::vector<csp::active_node_neighbors>& nbh = this -> mVerticalNbh; if ((*begin).dir > 1) nbh = this -> mLateralNbh;
    std::shuffle(nbh.begin(), nbh.end(), mRng); if (nbh.size()%2==1) nbh.pop_back();
    std::vector<ForwardIter> insertion_positions;
    std::vector<csp::active_node> init_sequence(begin, end); csp::active_node init_node = *begin;
    int heuristic_cost_bound = 4*cLVALUE*mLateralRange + 4;
    unsigned extraction_gain = get_distance(*std::prev(begin), *begin, false) +
                            get_distance(*begin, *std::next(begin), false) -
                            get_distance(*std::prev(begin), *std::next(begin), false);
    auto init_seq_pos = mMutableSolution.erase(begin, end);
    filter_insertion_positions(init_node, insertion_positions, extraction_gain+heuristic_cost_bound);
    for (int i=0; i < std::min((int)nbh.size()-2, mMaxSubsequences); i += 2) {
        bool nodes_satisfy[3]{true, true, true};
        node1 = {static_cast<unsigned short>((init_sequence.front()).x + nbh[i].x),
                 static_cast<unsigned short>((init_sequence.front()).y + nbh[i].y),
                 nbh[i].dir};
        node2 = {static_cast<unsigned short>((init_sequence.back()).x + nbh[i + 1].x),
                 static_cast<unsigned short>((init_sequence.back()).y + nbh[i + 1].y),
                 nbh[i + 1].dir};
        if (valid_node(node1) && valid_node(node2)) {
            sequence_satisfy(node1, node2, init_node, nodes_satisfy);
            if (best_insertion_search(node1, node2, insertion_positions, nodes_satisfy, extraction_gain)) return true;
        }
    }
    mMutableSolution.insert(init_seq_pos, init_sequence.begin(), init_sequence.end());
    set_mutable_insertion_coverage(init_sequence.begin(), init_sequence.end());
    return false;
}

template<class ForwardIter>
bool local_search::ls2::best_insertion_search(csp::active_node& n1, csp::active_node& n2,
                                              std::vector<ForwardIter>& pInsertionPositions,
                                              const bool* pNodesSatisfy, unsigned pMaxCost) {
    unsigned insertion_cost, second_insertion_cost;
    if (pNodesSatisfy[2]) {
        if (pNodesSatisfy[0]) {
            insertion_cost = best_insertion_search(n1, pInsertionPositions, pMaxCost);
            if (insertion_cost < pMaxCost) {
                set_mutable_insertion_coverage(n1);
                mLocalCost -= (pMaxCost - insertion_cost);
                return true; }
        } else if (pNodesSatisfy[1]) {
            insertion_cost = best_insertion_search(n2, pInsertionPositions, pMaxCost);
            if (insertion_cost < pMaxCost){
                set_mutable_insertion_coverage(n2);
                mLocalCost -= (pMaxCost - insertion_cost);
                return true; }
        }
        else {
            insertion_cost = best_insertion_search(n1, pInsertionPositions, pMaxCost);
            if (insertion_cost < pMaxCost) {
                second_insertion_cost = best_insertion_search(n2, pInsertionPositions, pMaxCost - insertion_cost);
                if (second_insertion_cost < pMaxCost - insertion_cost) {
                    set_mutable_insertion_coverage(n1);
                    set_mutable_insertion_coverage(n2);
                    mLocalCost -= (pMaxCost - (insertion_cost+ second_insertion_cost));
                    return true; }
                else {mMutableSolution.erase(pInsertionPositions.back()); pInsertionPositions.pop_back();}
            }
        }
    }
    return false;
}

void local_search::ls2::set_mutable_insertion_coverage(csp::active_node &pNode) {
    if (pNode.dir>1) {
        for (int i = pNode.y-mLateralRange; i<pNode.y + mLateralRange; i++){
            mMutableCoverage[i][pNode.x] += (1.0 - mMutableCoverage[i][pNode.x]) * map_pod(i, pNode.x, pNode);
        }
    } else {
        for (int j = pNode.x-mLateralRange; j < pNode.x+mLateralRange; ++j) {
            mMutableCoverage[pNode.y][j] += (1.0 - mMutableCoverage[pNode.y][j]) * map_pod(pNode.y, j,pNode);
        }
    }
}

template<class ForwardIter>
void local_search::ls2::set_mutable_insertion_coverage(ForwardIter begin, ForwardIter end){
    unsigned short i, j;
    for (; begin != end; begin++){
        if ((*begin).dir>1) {
            j = (*begin).x;
            for (i = (*begin).y-mLateralRange; i<(*begin).y + mLateralRange; i++){
                mMutableCoverage[i][j] += (1.0 - mMutableCoverage[i][j]) * map_pod(i, j, {.x = j, .y = (*begin).y, .dir=2});
            }
        } else {
            i = (*begin).y;
            for (j = (*begin).x-mLateralRange; j < (*begin).x+mLateralRange; ++j) {
                mMutableCoverage[i][j] += (1.0 - mMutableCoverage[i][j]) * map_pod(i, j,{.x = (*begin).x, .y = i, .dir=0});
            }
        }
    }
}

template <class ForwardIter>
unsigned local_search::ls2::best_insertion_search(csp::active_node& n1, std::vector<ForwardIter>& pInsertionPositions,
                                                  unsigned max_cost){
    unsigned insertion_cost;
    if (!valid_node(n1)) return max_cost;
    for (auto pos: pInsertionPositions){
        insertion_cost = get_distance(*pos, n1, false) +
                         get_distance(n1, *std::next(pos), false) -
                         get_distance(*pos, *std::next(pos), false);
        if (insertion_cost < max_cost){
            auto it = mMutableSolution.insert(std::next(pos), n1);
            pInsertionPositions.push_back(it);
            return insertion_cost;
        }
    }
    return max_cost;
}

void local_search::ls2::perturbation_procedure()
{
    unsigned best_insertion_cost, local_insertion_cost;
    csp::active_node r_node{};
    for (int i =0; i< mNPerturbations; i++){
        best_insertion_cost = cINF_EDGE;
        r_node = random_node_gen();
        auto insertion_pos = mMutableSolution.begin();
        for (auto pos = mMutableSolution.begin(); pos!= std::prev(mMutableSolution.end()); pos++){
            local_insertion_cost = get_distance(*pos, r_node, false) +
                         get_distance(r_node, *std::next(pos), false) -
                         get_distance(*pos, *std::next(pos), false);
            if (local_insertion_cost < best_insertion_cost){
                best_insertion_cost = local_insertion_cost;
                insertion_pos = pos;
                if (local_insertion_cost==0){
                    break;
                }
            }
        }
        mLocalCost += best_insertion_cost;
        auto new_pos = mMutableSolution.insert(std::next(insertion_pos), r_node);
        set_mutable_insertion_coverage(new_pos, std::next(new_pos));
    }
}

csp::active_node local_search::ls2::random_node_gen()
{
    unsigned short x, y, direction;
    direction = mDirectionRNodes(mRng);
    if (direction > 1){
        x = mXCoordHorizontalRNodes(mRng);
        y = mYCoordHorizontalRNodes(mRng);
    }
    else{
        x = mXCoordVerticalRNodes(mRng);
        y = mYCoordVerticalRNodes(mRng);
    }
    return {x, y, direction};
}

bool local_search::ls2::two_opt_procedure(){
    bool improved= false, opt_improved= false;
    unsigned current_cost, swap_cost;
    do
    {
        improved = false;
        for (auto it = mMutableSolution.begin(); it != std::prev(mMutableSolution.end(),3); it++){
            for (auto it2 = std::next(it, 2); it2 != std::prev(mMutableSolution.end(),2); it2++){
                csp::active_node swaped_node_2 = {(*it2).x, (*it2).y, cSWAP[(*it2).dir]};
                csp::active_node swaped_node_ = {(*std::next(it)).x, (*std::next(it)).y,cSWAP[(*std::next(it)).dir]};
                current_cost = get_distance(*it, *std::next(it), false) +
                               get_distance(*it2, *std::next(it2), false);
                swap_cost = get_distance(*it, swaped_node_2, false) +
                            get_distance(swaped_node_, *std::next(it2), false);
                if (swap_cost < current_cost){
                    swap_directions(std::next(it), std::next(it2));
                    std::reverse(std::next(it), std::next(it2));
                    improved = true;
                    mLocalCost -= (current_cost - swap_cost);
                }
            }
        }
        opt_improved = opt_improved || improved;
    } while (improved);

    return opt_improved;
}

template<class ForwardIter>
void local_search::ls2::filter_insertion_positions(const csp::active_node& ref, std::vector<ForwardIter>& container, unsigned pUpperbound) {
    unsigned cost;
    for (auto node_it = mMutableSolution.begin(); node_it != std::prev(mMutableSolution.end(), 2); node_it++) {
        cost = (abs((*node_it).x - ref.x) + abs((*node_it).y - ref.y)) +
               (abs((*std::next(node_it)).x - ref.x) + abs((*std::next(node_it)).y - ref.y)) -
               (abs((*node_it).x - (*std::next(node_it)).x) + abs((*node_it).y - (*std::next(node_it)).y));
        if (cost*cLVALUE <= pUpperbound) container.push_back(node_it);
    }
}

/*****************************TSP Methods*****************************/
void local_search::ls2::solve_lkh(){
    std::string tsp_options = get_tsp_options() + " -y " + mTspPathPrefix + ".tour";
    std::string solver_command = mSolverPath + " " + tsp_options + " -o " + mTspPathPrefix + ".sol " + mTspPathPrefix;
    cppied_algorithm::solve_tsp(solver_command);
}

void local_search::ls2::write_tsp_instance()
{
    write_tsp_initial_tour();
    std::ofstream o_tsp_instance;
    o_tsp_instance.open(mTspPathPrefix, std::ios::out);
    if (o_tsp_instance.is_open()){
        unsigned dimension = (3*mSolution.size());
        o_tsp_instance << "NAME : " << mUniqueId << std::endl;
        o_tsp_instance << "TYPE : TSP" << std::endl;
        o_tsp_instance << "DIMENSION : " << dimension << std::endl;
        o_tsp_instance << "EDGE_WEIGHT_TYPE : EXPLICIT" << std::endl;
        o_tsp_instance << "EDGE_WEIGHT_FORMAT : FULL_MATRIX" << std::endl;
        o_tsp_instance << "EDGE_WEIGHT_SECTION" << std::endl;
        write_tsp_matrix(o_tsp_instance);
        o_tsp_instance << "EOF" << std::endl;
        o_tsp_instance.close();
    }
    else {
        throw std::runtime_error("Error: Unable to open file");
    }
    o_tsp_instance.close();
}

void local_search::ls2::write_tsp_matrix(std::ofstream& ostream) {
    write_tsp_start_edges(ostream);
    for (std::size_t i=1; i< mSolution.size(); i++) {
        write_tsp_node_edges(ostream, mSolution[i], i);
        ostream << cINF_EDGE << " ";
        for (int j = 1; j < i; j++) {
            ostream << cINF_EDGE << " " << cINF_EDGE << " " << cINF_EDGE << " ";
        }
        ostream << 0 << " " << cINF_EDGE << " " << 0 << " ";
        for (int j = i + 1; j < mSolution.size(); j++) {
            ostream << cINF_EDGE << " " << cINF_EDGE << " " << cINF_EDGE << " ";
        }
        ostream << cINF_EDGE << " " << cINF_EDGE << std::endl;

        write_tsp_node_edges(ostream, {mSolution[i].x, mSolution[i].y, cSWAP[mSolution[i].dir]}, i);
    }
    write_tsp_dummy_edges(ostream);
}

void local_search::ls2::write_tsp_start_edges(std::ofstream &ostream) {
    ostream << cINF_EDGE << " ";
    for (int i=1; i< mSolution.size(); i++){
        ostream << get_distance(mStartNode, mSolution[i], false) << " "
                << cINF_EDGE << " "
                << get_distance(mStartNode, {mSolution[i].x, mSolution[i].y, cSWAP[mSolution[i].dir]}, false) << " ";
    }
    ostream << cINF_EDGE << " " << 0 << std::endl;
}

void local_search::ls2::write_tsp_node_edges(std::ofstream &ostream, const csp::active_node &pNode, const std::size_t pPosition) {
    csp::active_node reversed_node = {pNode.x, pNode.y, cSWAP[pNode.dir]};
    ostream << get_distance(mStartNode, pNode, false) << " ";
    for (int j = 1; j < pPosition; j++) {
        ostream << get_distance(reversed_node, mSolution[j], false) << " "
                << cINF_EDGE << " "
                << get_distance(reversed_node, {mSolution[j].x, mSolution[j].y, cSWAP[mSolution[j].dir]}, false) << " ";
    }
    ostream << cINF_EDGE << " " << 0 << " " << cINF_EDGE << " ";
    for (int j = pPosition + 1; j < mSolution.size(); j++) {
        ostream << get_distance(reversed_node, mSolution[j], false) << " "
                << cINF_EDGE << " "
                << get_distance(reversed_node, {mSolution[j].x, mSolution[j].y, cSWAP[mSolution[j].dir]}, false) << " ";
    }
    ostream << 0 << " " << cINF_EDGE << std::endl;
}

void local_search::ls2::write_tsp_dummy_edges(std::ofstream& ostream){
    //Dummy node to close loop
    ostream << cINF_EDGE << " ";
    for (int i=1; i<mSolution.size(); i++)
        ostream << 0 << " " << cINF_EDGE << " " << 0 << " ";
    ostream << cINF_EDGE << " " << 0 << " " << std::endl;
    ostream << 0 << " ";
    for (int i=1; i<mSolution.size(); i++)
        ostream << cINF_EDGE << " " << cINF_EDGE << " " << cINF_EDGE << " ";
    ostream << 0 << " " << cINF_EDGE << " " << std::endl;
}

void local_search::ls2::write_tsp_initial_tour() {
    std::ofstream o_tsp_solution;
    o_tsp_solution.open(mTspPathPrefix + ".tour", std::ios::out);
    if (o_tsp_solution.is_open()){
        unsigned dimension = (3*mSolution.size());
        o_tsp_solution << dimension << std::endl;
        for (std::size_t i=0; i< dimension-1; i++) o_tsp_solution << i << ' ';
        o_tsp_solution << dimension - 1 << std::endl << "EOF" << std::endl;
    } else throw std::runtime_error("Error: Unable to open file");
    o_tsp_solution.close();
}

void local_search::ls2::read_tsp_solution() {
    std::ifstream i_tsp_solution; bool reverse_sol; unsigned offset, dimension;
    i_tsp_solution.open(mTspPathPrefix + ".sol", std::ios::in);
    if (i_tsp_solution.is_open()) std::tie(reverse_sol, offset, dimension)= read_tsp_solution(i_tsp_solution);
    else throw std::runtime_error("Error: Unable to open tsp (lkh) sol file");
    i_tsp_solution.close();
    std::vector<csp::active_node> solution_copie = mSolution; mSolution.clear();
    if (reverse_sol){
        offset = dimension - (offset+1);
        std::reverse(mTspSolution.begin(), mTspSolution.end());
    }
    if (mTspSolution[0] != 0)
        std::rotate(mTspSolution.begin(), mTspSolution. begin()+offset, mTspSolution.end());
    assertm(mTspSolution[0]==0 && mTspSolution[dimension-2]==dimension-2 && mTspSolution[dimension-1]==dimension-1 , "Good Swap operations");
    mTspSolution.pop_back(); mTspSolution.pop_back();
    tsp_to_solution(solution_copie);
}

std::tuple<bool, std::size_t, std::size_t> local_search::ls2::read_tsp_solution(std::ifstream& i_tsp_solution){
    unsigned dimension, n_cols, source, dest, cost, offset; bool reverse_sol;
    i_tsp_solution >> dimension >> n_cols;
    assertm(dimension==n_cols, "Error while reading tsp sol file.\n");
    mTspSolution.clear(); mTspSolution.reserve(dimension);
    for(int i=0; i<dimension; i++) {
        i_tsp_solution >> source >> dest >> cost;
        mTspSolution.push_back(source);
        if (mTspSolution[i] == 0) { reverse_sol = (dest==dimension-1); offset=i;}
    }
    return std::make_tuple(reverse_sol, offset, dimension);
}

void local_search::ls2::tsp_to_solution(const std::vector<csp::active_node>& pSolution) {
    mSolution.clear();
    assertm(pSolution[0]==mStartNode, "The first node is different from the starting point.\n");
    mSolution.reserve(pSolution.size());mSolution.push_back(pSolution[0]);
    for (std::size_t i = 1; i<mTspSolution.size(); i+=3){
        if (mTspSolution[i]%3==1) push_shortest_path(pSolution[(mTspSolution[i]+2)/3]);
        else if (mTspSolution[i]%3 ==0)
            push_shortest_path({pSolution[mTspSolution[i]/3].x, pSolution[mTspSolution[i]/3].y, cSWAP[pSolution[mTspSolution[i]/3].dir]});
        else throw std::runtime_error("Entered a dummy edge in the TSP.\n");
    }
}

/*****************************************************************************/







