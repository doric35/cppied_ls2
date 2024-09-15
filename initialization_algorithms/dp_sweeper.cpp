#include "dp_sweeper.hpp"

void init::dp_sweeper::d_solve() {
    mStart = std::chrono::system_clock::now();
    generate_feasible_segments();
    mTspStart = std::chrono::system_clock::now();
    if (mSolverPath.find("dps_NN") == std::string::npos) {
        std::string tsp_options = get_tsp_options();
        solve_tsp(mSolverPath + " " + tsp_options + " -o " + mTspPathPrefix + ".sol " + mTspPathPrefix);
    } else nn3opt(mSolverPath.find("3Opt") != std::string::npos);
    mTspEnd = std::chrono::system_clock::now();
    mCoverage = std::vector<std::vector<double>>(mGridDimensions.first, std::vector<double>(mGridDimensions.second, 0.0));
    cppied_algorithm::set_insertion_coverage(mSolution.begin(), mSolution.end());
    mEnd = std::chrono::system_clock::now();
}

void init::dp_sweeper::terminate() {
    unpad();
    set_lexicographic_cost();
    std::string solver = mSolverPath.substr(mSolverPath.find_last_of('/')+1);
    auto [ui_n_moves, ui_n_turns] = get_bi_cost();
    info.push_back({"solver", "cost", "n_moves", "n_turns", "dp_time(ms)", "tsp_time(ms)", "total_time(ms)"});
    info.push_back({solver, std::to_string(get_lexicographic_cost()), std::to_string(ui_n_moves), std::to_string(ui_n_turns),
                    std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(mTspStart-mStart).count()),
                    std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(mTspEnd - mTspStart).count()),
                    std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(mEnd-mStart).count())});
}

void init::dp_sweeper::generate_feasible_segments() {
    std::size_t h_cardinality, v_cardinality;
    do {
        dp_h_forward();
        dp_v_forward();
        h_cardinality = mHSegments.size(); v_cardinality = mVSegments.size(); unsigned short prev_covered = mNcovered;
        if (mHdpGain.back().first > mVdpGain.back().first + cEPSILON){
            dp_h_backward();
            for (std::size_t i = h_cardinality; i<mHSegments.size(); i++)
                set_segment_insertion_coverage(mHSegments[i].first,mHSegments[i].second);
        } else if (mHdpGain.back().first + cEPSILON < mVdpGain.back().first){
            dp_v_backward();
            for (std::size_t i = v_cardinality; i<mVSegments.size(); i++)
                set_segment_insertion_coverage(mVSegments[i].first,mVSegments[i].second);
        } else
            genererate_segments_from_cardinality(h_cardinality, v_cardinality);
        if (mNcovered==prev_covered && h_cardinality == mHSegments.size() && v_cardinality == mVSegments.size())
            mPenalty/=(double)(2*mLateralRange);
    } while (mNcovered < mNcells);
}

void init::dp_sweeper::genererate_segments_from_cardinality(unsigned short prev_h_cardinality, unsigned short prev_v_cardinality) {
    dp_h_backward();
    dp_v_backward();
    if (mHSegments.size()-prev_h_cardinality <= mVSegments.size() -prev_v_cardinality){
        for (std::size_t i = prev_h_cardinality; i<mHSegments.size(); i++)
            set_segment_insertion_coverage(mHSegments[i].first,mHSegments[i].second);
        mVSegments.resize(prev_v_cardinality);
    } else{
        for (std::size_t i = prev_v_cardinality; i<mVSegments.size(); i++)
            set_segment_insertion_coverage(mVSegments[i].first,mVSegments[i].second);
        mHSegments.resize(prev_h_cardinality);
    }
}

std::pair<double, std::pair<unsigned short , unsigned short>> init::dp_sweeper::h_kadane(unsigned short line_num) {
    double local_gain, max_current=0.0, max_global=0.0;
    unsigned short start=mPadSize, end=mGridDimensions.second - mPadSize - 1, temps_start=mPadSize;
    for (unsigned short j=mPadSize; j<mGridDimensions.second - mPadSize; j++){
        local_gain=0.0;
        for (unsigned short i = line_num - mLateralRange; i < line_num + mLateralRange; i++)
            local_gain += (mRequiredCoverage > mCoverage[i][j]) ?
                    std::min(mRequiredCoverage - mCoverage[i][j], (1.0 - mCoverage[i][j])*
                                                                    map_pod(i, j,{j, line_num, 2})) :-mPenalty;
        if (local_gain <= cEPSILON) {
            if (max_current > max_global) {
                start = temps_start;
                end = j - 1;
                max_global = max_current;
            }
            temps_start = j + 1;
            max_current = 0.0;
        } else max_current += local_gain;
    }
    if (end == mGridDimensions.second - mPadSize - 1 && start == mPadSize) max_global=max_current;
    return std::make_pair(max_global, std::make_pair(start, end));
}

std::pair<double, std::pair<unsigned short , unsigned short >> init::dp_sweeper::v_kadane(unsigned short col_num) {
    double local_gain, max_current=0.0, max_global=0.0;
    unsigned short start=mPadSize, end=mGridDimensions.second - mPadSize - 1, temps_start=mPadSize;
    for (unsigned short i=mPadSize; i<mGridDimensions.first - mPadSize; i++){
        local_gain=0.0;
        for (unsigned short j = col_num-mLateralRange; j < col_num + mLateralRange; j++)
            local_gain += (mRequiredCoverage > mCoverage[i][j]) ?
                          std::min(mRequiredCoverage - mCoverage[i][j], (1.0 - mCoverage[i][j])*
                                                                        map_pod(i, j,{col_num, i, 0})) :-mPenalty;
        if (local_gain<=cEPSILON){
            if (max_current > max_global + cEPSILON){
                start = temps_start; end = i-1; max_global = max_current;
            }
            temps_start = i+1; max_current = 0.0;
        } else max_current += local_gain;
    }
    if (end == mGridDimensions.first - mPadSize - 1 && start == mPadSize) max_global=max_current;
    return std::make_pair(max_global, std::make_pair(start, end));
}

void init::dp_sweeper::dp_h_forward() {
    std::pair<double, std::pair<unsigned short , unsigned short >> max_subarray;
    mHdpGain[0] = h_kadane(mPadSize);
    for (unsigned short i=mPadSize+1; i < static_cast<unsigned short>(std::min(mPadSize + 2*mLateralRange,mGridDimensions.first-mPadSize+1)); i++)
        mHdpGain[i-mPadSize] = std::max(mHdpGain[i-(mPadSize+1)],h_kadane(i));
    for (unsigned short i=static_cast<unsigned short>(std::min(mPadSize + 2*mLateralRange,mGridDimensions.first-mPadSize+1)); i< static_cast<unsigned short>(mGridDimensions.first-mPadSize+1); i++){
        max_subarray = h_kadane(i);
        mHdpGain[i-mPadSize] = std::make_pair(std::max(mHdpGain[i-mPadSize-2*mLateralRange].first + max_subarray.first,
                                                            mHdpGain[i-mPadSize-1].first),
                                                   max_subarray.second);
    }
}

void init::dp_sweeper::dp_v_forward() {
    std::pair<double, std::pair<unsigned short , unsigned short >> max_subarray;
    mVdpGain[0] = v_kadane(mPadSize);
    for (unsigned short j=mPadSize+1; j < static_cast<unsigned short>(std::min(mPadSize + 2*mLateralRange,mGridDimensions.second-mPadSize+1)); j++)
        mVdpGain[j-mPadSize] = std::max(mVdpGain[j-(mPadSize+1)],v_kadane(j));
    for (unsigned short j=static_cast<unsigned short>(std::min(mPadSize + 2*mLateralRange,mGridDimensions.second-mPadSize+1));j< static_cast<unsigned short>(mGridDimensions.second-mPadSize+1); j++){
        max_subarray = v_kadane(j);
        mVdpGain[j-mPadSize] = std::make_pair(std::max(mVdpGain[j-(mPadSize+2*mLateralRange)].first + max_subarray.first,
                                                   mVdpGain[j-(mPadSize+1)].first),
                                              max_subarray.second);
    }
}

void init::dp_sweeper::dp_h_backward() {
    auto i = static_cast<short>(mGridDimensions.first - mPadSize);
    double diff;
    while (i > mPadSize){
        diff = mHdpGain[i-mPadSize].first - mHdpGain[i-(1+mPadSize)].first;
        if (diff > cEPSILON){
            mHSegments.push_back({{mHdpGain[i-mPadSize].second.first, (unsigned short) i, 2},
                                  {mHdpGain[i-mPadSize].second.second, (unsigned short) i, 3}});
            i-=2*mLateralRange;
        }
        else if (mHdpGain[i-mPadSize].first <= cEPSILON) i=static_cast<short>(mPadSize-1);
        else i-=1;
    }
    if (i==mPadSize && mHdpGain[0].first>cEPSILON)
        mHSegments.push_back({{mHdpGain[0].second.first, (unsigned short) i, 2},
                              {mHdpGain[0].second.second, (unsigned short) i, 3}});
}

void init::dp_sweeper::dp_v_backward() {
    auto j = static_cast<short>(mGridDimensions.second - mPadSize);
    double diff;
    while (j > mPadSize){
        diff = mVdpGain[j-mPadSize].first - mVdpGain[j-(1+mPadSize)].first;
        if (diff > cEPSILON){
            mVSegments.push_back({{(unsigned short)j, mVdpGain[j-mPadSize].second.first, 1},
                                  {(unsigned short)j, mVdpGain[j-mPadSize].second.second, 0}});
            j-=2*mLateralRange;
        }
        else if (mVdpGain[j-mPadSize].first <= cEPSILON) j=static_cast<short>(mPadSize-1);
        else j-=1;
    }
    if (j==mPadSize && mVdpGain[0].first > cEPSILON)
        mVSegments.push_back({{(unsigned short)j, mVdpGain[0].second.first, 1},
                              {(unsigned short)j, mVdpGain[0].second.second, 0}});
}

void init::dp_sweeper::set_segment_insertion_coverage(const csp::active_node& pFirst,
                                              const csp::active_node& pLast) {
    double p_ij;
    if (pFirst.dir>1) {
        for (unsigned short i = pFirst.y-mLateralRange; i<pFirst.y + mLateralRange; i++){
            for (unsigned short j = pFirst.x; j<=pLast.x; j++){
                p_ij = mCoverage[i][j];
                mCoverage[i][j] += (1.0 - mCoverage[i][j])* map_pod(i, j, {.x = j, .y = pFirst.y, .dir=2});
                mNcovered += (p_ij <(mRequiredCoverage-mTolerance)) & (mCoverage[i][j]>= (mRequiredCoverage-mTolerance));
            }
        }
    } else {
        for (unsigned short i = pFirst.y; i <= pLast.y; i++){
            for (unsigned short j = pFirst.x-mLateralRange; j < pFirst.x+mLateralRange; ++j){
                p_ij = mCoverage[i][j];
                mCoverage[i][j] += (1.0 - mCoverage[i][j])* map_pod(i, j, {.x = pFirst.x, .y = i, .dir=0});
                mNcovered += (p_ij <(mRequiredCoverage-mTolerance)) & (mCoverage[i][j]>= (mRequiredCoverage-mTolerance));
            }
        }
    }
}

/******************************************TSP Computation*******************************************/

void init::dp_sweeper::nn3opt(bool pLocalSearch) {
    //Make data struct from HDp and VDp
    std::vector<std::pair<csp::active_node, csp::active_node>> segment_tour = {{mSolution[0], mSolution[0].swap()}};
    segment_tour.insert(segment_tour.end(), mHSegments.begin(), mHSegments.end());
    segment_tour.insert(segment_tour.end(), mVSegments.begin(), mVSegments.end());
    for (std::size_t i=0; i<segment_tour.size()-2; i++) move_nearest_segment(segment_tour, i);

    if (!pLocalSearch) {
        std::size_t i=1;
        if (get_distance(segment_tour.front().second.swap(), segment_tour[1].first, true)==0)
            push_shortest_path(segment_tour[i++].second.swap());
        for (; i< segment_tour.size(); i++) {
            push_shortest_path(segment_tour[i].first);
            push_shortest_path(segment_tour[i].second.swap());
        }
    } else {dps3opt(segment_tour);}
}

void init::dp_sweeper::dps3opt(std::vector<std::pair<csp::active_node,csp::active_node>>& segment_path) {
    unsigned i_cost;
    std::size_t dimension = segment_path.size();
    for (std::size_t i=1; i<dimension-1; i++) {
        for (std::size_t j=i+2; j< dimension; j++){
            i_cost = dps2opt_move(segment_path,  i, j);
            for (std::size_t k=j+2; k<dimension; k++) {
                i_cost = dps3opt_move(segment_path, i, j, k, i_cost);
            }
        }
    }
    dps3opt_reconstruct_path(segment_path);
}

void init::dp_sweeper::dps3opt_reconstruct_path(std::vector<std::pair<csp::active_node, csp::active_node>>& segment_path) {
    //Rearrange segment path
    if (get_distance(mSolution[0], segment_path[1].first, true)==0)
        push_shortest_path(segment_path[1].second.swap());
    else push_shortest_path(segment_path[1].first);
    for (std::size_t i=2; i< segment_path.size(); i++){
        push_shortest_path(segment_path[i].first);
        push_shortest_path(segment_path[i].second.swap());
    }

}

unsigned init::dp_sweeper::dps2opt_move(std::vector<std::pair<csp::active_node, csp::active_node>>& segment_path,
                                         std::size_t i, std::size_t j){
    //(1.0, 1.1) (2.0, 2.1) INIT
    unsigned best_cost = get_distance(segment_path[i].second.swap(),segment_path[i+1].first, false) +
                        get_distance(segment_path[j].second.swap(),segment_path[j+1].first, false);
    //(1.0, 2.0) (1.1, 2.1)
    unsigned current_cost = get_distance(segment_path[i].second.swap(), segment_path[j].second, false) +
                         get_distance(segment_path[i+1].first.swap(), segment_path[j+1].first, false);
    if (current_cost < best_cost){
        opt_flip(std::next(segment_path.begin(), i+1), std::next(segment_path.begin(), j+1));
        best_cost = current_cost;
    }
    return best_cost;
}

unsigned init::dp_sweeper::dps3opt_move(std::vector<std::pair<csp::active_node,csp::active_node>>& segment_path,
                                    std::size_t i, std::size_t j, std::size_t k, unsigned i_cost){
    //(1.0, 1.1) (2.0, 2.1) (3.0, 3.1) INIT
    unsigned new_i_cost, current_cost;
    unsigned best_cost = i_cost + get_distance(segment_path[k].second.swap(),segment_path[k+1].first, false);
    //(1.0, 2.0) (1.1, 3.0) (2.1, 3.1) OK
    new_i_cost = get_distance(segment_path[i].second.swap(), segment_path[j].second, false) +
                   + get_distance(segment_path[i+1].first.swap(),segment_path[k].first, false);
    current_cost = new_i_cost + get_distance(segment_path[j+1].first.swap(),segment_path[k+1].first, false);
    if (current_cost < best_cost){
        opt_flip(std::next(segment_path.begin(), i+1), std::next(segment_path.begin(), j+1));
        opt_flip(std::next(segment_path.begin(), j+1), std::next(segment_path.begin(), k+1));
        return new_i_cost;
    }
    //(1.0, 3.0) (2.1, 1.1) (2.0, 3.1) OK
    new_i_cost = get_distance(segment_path[i].second.swap(),
                                segment_path[k].second, false)
                   +get_distance(segment_path[j+1].first.swap(),
                                 segment_path[i+1].first, false);
    current_cost = new_i_cost + get_distance(segment_path[j].second.swap(),
                                 segment_path[k+1].first, false);
    if (current_cost < best_cost){
        opt_flip(std::next(segment_path.begin(), i+1), std::next(segment_path.begin(), k+1));
        j = (i+1) + (k-(j+1));
        opt_flip(std::next(segment_path.begin(), j+1), std::next(segment_path.begin(), k+1));
        return new_i_cost;
    }
    //(1.0, 2.1) (3.0, 2.0) (1.1, 3.1) OK
    new_i_cost = get_distance(segment_path[i].second.swap(),
                                segment_path[j+1].first, false)
                   +get_distance(segment_path[k].second.swap(),
                                 segment_path[j].second, false);
    current_cost = new_i_cost + get_distance(segment_path[i+1].first.swap(),
                                 segment_path[k+1].first, false);
    if (current_cost < best_cost){
        opt_flip(std::next(segment_path.begin(), j+1), std::next(segment_path.begin(), k+1));
        opt_flip(std::next(segment_path.begin(), i+1), std::next(segment_path.begin(), k+1));
        return new_i_cost;
    }
    //(1.0, 2.1) (3.0, 1.1) (2.0, 3.1) OK
    new_i_cost = get_distance(segment_path[i].second.swap(),
                                segment_path[j+1].first, false)
                   +get_distance(segment_path[k].second.swap(),
                                 segment_path[i+1].first, false);
    current_cost = new_i_cost + get_distance(segment_path[j].second.swap(),
                                 segment_path[k+1].first, false);
    if (current_cost < best_cost){
        opt_flip(std::next(segment_path.begin(), i+1), std::next(segment_path.begin(), j+1));
        opt_flip(std::next(segment_path.begin(), j+1), std::next(segment_path.begin(), k+1));
        opt_flip(std::next(segment_path.begin(), i+1), std::next(segment_path.begin(), k+1));
        return new_i_cost;
    } else return i_cost;
}

void init::dp_sweeper::opt_flip(std::vector<std::pair<csp::active_node, csp::active_node>>::iterator first,
                                std::vector<std::pair<csp::active_node, csp::active_node>>::iterator last) {
    std::reverse(first, last);
    for (; first != last; first++){
        std::swap(first->first, first->second);
    }
}

void init::dp_sweeper::move_nearest_segment(std::vector<std::pair<csp::active_node, csp::active_node>> &segment_path,
                                            std::size_t pos) {
    unsigned best_cost = cINF_EDGE;
    std::size_t nearest_nb=pos+1;
    bool accept_zero = pos==0, flip= false;
    std::pair<csp::active_node, csp::active_node> temp;
    unsigned upper_bound = cLVALUE - (accept_zero*cLVALUE);
    for (std::size_t next_pos=pos+1; next_pos<segment_path.size(); next_pos++) {
        unsigned cost1 = get_distance(segment_path[pos].second.swap(), segment_path[next_pos].first, accept_zero);
        unsigned cost2 = get_distance(segment_path[pos].second.swap(), segment_path[next_pos].second, accept_zero);
        if (cost1 < best_cost) {nearest_nb = next_pos; best_cost = cost1; flip=false;}
        if (cost2 < best_cost) {nearest_nb = next_pos; best_cost = cost2; flip=true;}
        if (best_cost<= upper_bound) break;
    }
    if (flip) {
        temp = segment_path[pos+1];
        segment_path[pos+1] = {segment_path[nearest_nb].second, segment_path[nearest_nb].first};
        segment_path[nearest_nb] = temp;
    }
    else {
        temp = segment_path[pos+1];
        segment_path[pos+1] = segment_path[nearest_nb];
        segment_path[nearest_nb] = temp;
    }
}

void init::dp_sweeper::write_tsp_instance() {
    std::ofstream o_tsp_instance;
    o_tsp_instance.open(mTspPathPrefix, std::ios::out);
    if (o_tsp_instance.is_open()){
        unsigned dimension = mHSegments.size() * 3 + mVSegments.size() * 3 + 3;
        o_tsp_instance << "NAME : "<< mUniqueId << std::endl;
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

void init::dp_sweeper::write_tsp_matrix(std::ofstream &ostream) {
    write_tsp_start(ostream);
    write_h_tsp(ostream);
    write_v_tsp(ostream);
}

void init::dp_sweeper::write_tsp_start(std::ofstream &ostream) {
    ostream << cINF_EDGE << " " << 0 << " " << cINF_EDGE << " ";
    for (auto & mHSegment : mHSegments){
        ostream << get_distance(mStartNode, mHSegment.first, true)
                << " " << cINF_EDGE << " "
                << get_distance(mStartNode, mHSegment.second, true)
                << " ";
    }
    for (auto & mVSegment : mVSegments){
        ostream << get_distance(mStartNode, mVSegment.first, true)
                << " " << cINF_EDGE << " "
                << get_distance(mStartNode, mVSegment.second, true)
                << " ";
    }
    ostream << std::endl;
    ostream << 0 << " " << cINF_EDGE << " " << 0 << " ";
    for (std::size_t j = 0; j < mHSegments.size(); j++)
        ostream << cINF_EDGE << " " << cINF_EDGE << " " << cINF_EDGE << " ";

    for (std::size_t j = 0; j < mVSegments.size(); j++)
        ostream << cINF_EDGE << " " << cINF_EDGE << " " << cINF_EDGE << " ";

    ostream << std::endl;
    ostream << cINF_EDGE << " " << 0 << " " << cINF_EDGE << " ";
    for (std::size_t j=0; j< mHSegments.size(); j++){
        ostream << 0 << " " << cINF_EDGE << " " << 0 << " ";
    }
    for (std::size_t j=0; j< mVSegments.size(); j++){
        ostream << 0 << " " << cINF_EDGE << " " << 0 << " ";
    }
    ostream << std::endl;
}

void init::dp_sweeper::write_h_tsp(std::ofstream &ostream) {
    for (std::size_t i=0; i< mHSegments.size(); i++){
        // start
        ostream << get_distance(mStartNode, mHSegments[i].first, true)
                << " " << cINF_EDGE << " " << 0 << " ";
        // source
        for (std::size_t j = 0; j < i; j++)
            ostream << get_distance({mHSegments[i].first.x, mHSegments[i].first.y, cSWAP[mHSegments[i].first.dir]},
                                    mHSegments[j].first, false) << " " << cINF_EDGE << " "
                    << get_distance({mHSegments[i].first.x, mHSegments[i].first.y, cSWAP[mHSegments[i].first.dir]},
                                    mHSegments[j].second, false) << " ";
        ostream << cINF_EDGE << " " << 0 << " " << cINF_EDGE << " ";
        for (std::size_t j = i+1; j < mHSegments.size(); j++)
            ostream << get_distance({mHSegments[i].first.x, mHSegments[i].first.y, cSWAP[mHSegments[i].first.dir]},
                                    mHSegments[j].first, false) << " " << cINF_EDGE << " "
                    << get_distance({mHSegments[i].first.x, mHSegments[i].first.y, cSWAP[mHSegments[i].first.dir]},
                                    mHSegments[j].second, false) << " ";
        for (auto &segment: mVSegments)
            ostream << get_distance({mHSegments[i].first.x, mHSegments[i].first.y, cSWAP[mHSegments[i].first.dir]},
                                    segment.first, false) << " " << cINF_EDGE << " "
                    << get_distance({mHSegments[i].first.x, mHSegments[i].first.y, cSWAP[mHSegments[i].first.dir]},
                                    segment.second, false) << " ";
        ostream << std::endl;
        // dummy
        ostream << cINF_EDGE << " " << cINF_EDGE << " " << cINF_EDGE << " ";
        for (std::size_t j = 0; j < mHSegments.size(); j++)
            ostream << ((i == j) ? 0 : cINF_EDGE) << " " << cINF_EDGE << " " << ((i == j) ? 0 : cINF_EDGE) << " ";
        for (std::size_t j = 0; j < mVSegments.size(); j++)
            ostream << cINF_EDGE << " " << cINF_EDGE << " " << cINF_EDGE << " ";
        ostream << std::endl;
        //target
        ostream << get_distance(mStartNode, mHSegments[i].second, true)
                << " " << cINF_EDGE << " " << 0 << " ";
        for (std::size_t j = 0; j < i; j++)
            ostream << get_distance({mHSegments[i].second.x, mHSegments[i].second.y, cSWAP[mHSegments[i].second.dir]},
                                    mHSegments[j].first, false) << " " << cINF_EDGE << " "
                    << get_distance({mHSegments[i].second.x, mHSegments[i].second.y, cSWAP[mHSegments[i].second.dir]},
                                    mHSegments[j].second, false) << " ";
        ostream << cINF_EDGE << " " << 0 << " " << cINF_EDGE << " ";
        for (std::size_t j=i+1; j< mHSegments.size(); j++)
            ostream << get_distance({mHSegments[i].second.x, mHSegments[i].second.y, cSWAP[mHSegments[i].second.dir]},
                                    mHSegments[j].first, false) << " " << cINF_EDGE << " "
                    << get_distance({mHSegments[i].second.x, mHSegments[i].second.y, cSWAP[mHSegments[i].second.dir]},
                                    mHSegments[j].second, false) << " ";
        for (auto &segment: mVSegments)
            ostream << get_distance({mHSegments[i].second.x, mHSegments[i].second.y, cSWAP[mHSegments[i].second.dir]},
                                    segment.first, false) << " " << cINF_EDGE << " "
                    << get_distance({mHSegments[i].second.x, mHSegments[i].second.y, cSWAP[mHSegments[i].second.dir]},
                                    segment.second, false) << " ";
        ostream << std::endl;
    }
}

void init::dp_sweeper::write_v_tsp(std::ofstream &ostream) {
    for (std::size_t i=0; i< mVSegments.size(); i++){
        // start
        ostream << get_distance(mStartNode, mVSegments[i].first, true)
                << " " << cINF_EDGE << " " << 0 << " ";
        // source
        for (auto &segment: mHSegments)
            ostream << get_distance({mVSegments[i].first.x, mVSegments[i].first.y, cSWAP[mVSegments[i].first.dir]},
                                    segment.first, false) << " " << cINF_EDGE << " "
                    << get_distance({mVSegments[i].first.x, mVSegments[i].first.y, cSWAP[mVSegments[i].first.dir]},
                                    segment.second, false) << " ";
        for (std::size_t j = 0; j < i; j++)
            ostream << get_distance({mVSegments[i].first.x, mVSegments[i].first.y, cSWAP[mVSegments[i].first.dir]},
                                    mVSegments[j].first, false) << " " << cINF_EDGE << " "
                    << get_distance({mVSegments[i].first.x, mVSegments[i].first.y, cSWAP[mVSegments[i].first.dir]},
                                    mVSegments[j].second, false) << " ";
        ostream << cINF_EDGE << " " << 0 << " " << cINF_EDGE << " ";
        for (std::size_t j = i+1; j < mVSegments.size(); j++)
            ostream << get_distance({mVSegments[i].first.x, mVSegments[i].first.y, cSWAP[mVSegments[i].first.dir]},
                                    mVSegments[j].first, false) << " " << cINF_EDGE << " "
                    << get_distance({mVSegments[i].first.x, mVSegments[i].first.y, cSWAP[mVSegments[i].first.dir]},
                                    mVSegments[j].second, false) << " ";
        ostream << std::endl;
        // dummy
        ostream << cINF_EDGE << " " << cINF_EDGE << " " << cINF_EDGE << " ";
        for (std::size_t j = 0; j < mHSegments.size(); j++)
            ostream << cINF_EDGE << " " << cINF_EDGE << " " << cINF_EDGE << " ";
        for (std::size_t j = 0; j < mVSegments.size(); j++)
            ostream << ((i == j) ? 0 : cINF_EDGE) << " " << cINF_EDGE << " " << ((i == j) ? 0 : cINF_EDGE) << " ";
        ostream << std::endl;
        //target
        ostream << get_distance(mStartNode, mVSegments[i].second, true)
                << " " << cINF_EDGE << " " << 0 << " ";
        for (auto &segment: mHSegments)
            ostream << get_distance({mVSegments[i].second.x, mVSegments[i].second.y, cSWAP[mVSegments[i].second.dir]},
                                    segment.first, false) << " " << cINF_EDGE << " "
                    << get_distance({mVSegments[i].second.x, mVSegments[i].second.y, cSWAP[mVSegments[i].second.dir]},
                                    segment.second, false) << " ";
        for (std::size_t j = 0; j < i; j++)
            ostream << get_distance({mVSegments[i].second.x, mVSegments[i].second.y, cSWAP[mVSegments[i].second.dir]},
                                    mVSegments[j].first, false) << " " << cINF_EDGE << " "
                    << get_distance({mVSegments[i].second.x, mVSegments[i].second.y, cSWAP[mVSegments[i].second.dir]},
                                    mVSegments[j].second, false) << " ";
        ostream << cINF_EDGE << " " << 0 << " " << cINF_EDGE << " ";
        for (std::size_t j=i+1; j< mVSegments.size(); j++)
            ostream << get_distance({mVSegments[i].second.x, mVSegments[i].second.y, cSWAP[mVSegments[i].second.dir]},
                                    mVSegments[j].first, false) << " " << cINF_EDGE << " "
                    << get_distance({mVSegments[i].second.x, mVSegments[i].second.y, cSWAP[mVSegments[i].second.dir]},
                                    mVSegments[j].second, false) << " ";
        ostream << std::endl;
    }
}

void init::dp_sweeper::read_tsp_solution() {
    std::ifstream i_tsp_solution;
    i_tsp_solution.open(mTspPathPrefix + ".sol", std::ios::in);
    std::size_t dimension = mHSegments.size() * 3 + mVSegments.size() * 3 + 3, n_lines, n_cols;
    mSolution.clear(); mTspSolution.clear(); mTspSolution.reserve(dimension);
    i_tsp_solution >> n_lines >> n_cols;
    assertm( n_cols==n_lines & n_lines==dimension, "Good number of nodes");
    if (i_tsp_solution.is_open()) {
        if (mTspSolType == "tour") read_tsp_tour(i_tsp_solution, dimension);
        else if (mTspSolType == "edges") read_tsp_edges(i_tsp_solution, dimension);
        else throw std::runtime_error("Unrecognized tsp solution type.\n");
    } else throw std::runtime_error("Error: Unable to open tsp sol file. \n");
    i_tsp_solution.close();
    mSolution.push_back(mStartNode);
    assertm(mTspSolution[0]==0 & mTspSolution[1]!=1 & mTspSolution[dimension-1]==1 & mTspSolution[dimension-2]==2, "Good Swap operations");
    mTspSolution.pop_back(); mTspSolution.pop_back();
    tsp_to_solution();
}

void init::dp_sweeper::read_tsp_tour(std::ifstream& pIfStream, std::size_t dimension) {
    unsigned offset, source, dest, cost; bool reverse_sol;
    for(std::size_t i=0; i<dimension; i++) {
        pIfStream >> source >> dest >> cost;
        mTspSolution.push_back(source);
        if (mTspSolution[i] == 0) { reverse_sol = (dest==1); offset = i;}
    }
    if (reverse_sol) {
        offset = dimension - (offset + 1);
        std::reverse(mTspSolution.begin(), mTspSolution.end());
    }
    if (mTspSolution[0] != 0)
        std::rotate(mTspSolution.begin(), mTspSolution. begin()+offset, mTspSolution.end());
}

void init::dp_sweeper::read_tsp_edges(std::ifstream& pIfStream, std::size_t dimension) {
    unsigned source, dest, cost;
    std::vector<std::vector<std::pair<unsigned, unsigned>>> edges(dimension);
    for(std::size_t i=0; i< dimension; i++){
        pIfStream >> source >> dest >> cost;
        edges[source].emplace_back(source, dest);
        edges[dest].emplace_back(source, dest);
    }
    source = 1;
    dest=0;
    for(std::size_t i=0; i< dimension; i++){
        mTspSolution.push_back(dest);
        auto edge1 = edges[dest][0], edge2 = edges[dest][1];
        if ((edge1.first == source && edge1.second == dest) || (edge1.first == dest && edge1.second == source)){
            if (edge2.first == dest){source=dest; dest=edge2.second;}
            else{source=dest; dest=edge2.first;}
        }
        else if (edge1.first==dest){source=dest; dest=edge1.second;}
        else{source=dest; dest=edge1.first;}
    }
    assertm(dest==0, "Error, loop was not completed properly.\n");
}

void init::dp_sweeper::tsp_to_solution() {
    push_tsp_first_segment();
    for (std::size_t i = 4; i < mTspSolution.size(); i+=3){
        if (mTspSolution[i] >= 3 + 3*mHSegments.size()) push_tsp_vertical_segment(i);
        else push_tsp_horizontal_segment(i);
    }
}

void init::dp_sweeper::push_tsp_first_segment() {
    csp::active_node node_s{}, node_t{};
    if (mTspSolution[1] >= 3 + 3*mHSegments.size()){
        if (mTspSolution[1] % 3 == 0) {
            node_s = mVSegments[mTspSolution[1] / 3 - (1 + mHSegments.size())].first;
            node_t = mVSegments[mTspSolution[1] / 3 - (1 + mHSegments.size())].second;
        }
        else {
            node_s = mVSegments[mTspSolution[1] / 3 - (1 + mHSegments.size())].second;
            node_t = mVSegments[mTspSolution[1] / 3 - (1 + mHSegments.size())].first;
        }
    } else{
        if (mTspSolution[1] % 3 == 0) {
            node_s = mHSegments[mTspSolution[1]/3 - 1].first;
            node_t = mHSegments[mTspSolution[1]/3 - 1].second;
        }
        else {
            node_s = mHSegments[mTspSolution[1]/3 - 1].second;
            node_t = mHSegments[mTspSolution[1]/3 - 1].first;
        }
    }
    node_t = {node_t.x, node_t.y, cSWAP[node_t.dir]};
    if (mStartNode.x != node_s.x | mStartNode.y != node_s.y | mStartNode.dir != node_s.dir) push_shortest_path(node_s);
    if (node_t.x != node_s.x | node_t.y != node_s.y) push_shortest_path(node_t);
}

void init::dp_sweeper::push_tsp_vertical_segment(std::size_t i) {
    std::pair<csp::active_node, csp::active_node> segment = mVSegments[mTspSolution[i]/3 - (1 + mHSegments.size())];
    if (mTspSolution[i] % 3 == 0) {
        push_shortest_path(segment.first);
        if (segment.first.y != segment.second.y)
            push_shortest_path({segment.second.x, segment.second.y, cSWAP[segment.second.dir]});
    }
    else {
        push_shortest_path(segment.second);
        if (segment.first.y != segment.second.y)
            push_shortest_path({segment.first.x, segment.first.y, cSWAP[segment.first.dir]});
    }
}

void init::dp_sweeper::push_tsp_horizontal_segment(std::size_t i) {
    std::pair<csp::active_node, csp::active_node> segment = mHSegments[mTspSolution[i]/3 - 1];
    if (mTspSolution[i] % 3 == 0) {
        push_shortest_path(segment.first);
        if (segment.first.x != segment.second.x)
            push_shortest_path({segment.second.x, segment.second.y, cSWAP[segment.second.dir]});
    }
    else {
        push_shortest_path(segment.second);
        if (segment.first.x != segment.second.x)
            push_shortest_path({segment.first.x, segment.first.y, cSWAP[segment.first.dir]});
    }
}
