#ifndef CPPIED_ALGORITHM_HPP
#define CPPIED_ALGORITHM_HPP

// Include any necessary headers
#include "utils/utils.hpp"

// Define any necessary macros or constants
class cppied_algorithm
{
protected:
    /* data */
    std::vector<std::vector<double>> mCoverage;
    std::vector<std::vector<double>> mPod;
    std::vector<std::vector<unsigned short>> mSeabed;
    std::vector<csp::active_node> mSolution;
    double mRequiredCoverage;
    double mTolerance;
    unsigned short mLateralRange;
    unsigned short mPadSize;
    double mPadValue;
    csp::active_node mStartNode;
    std::pair<unsigned short, unsigned short> mGridDimensions;
    unsigned mSolutionCost;
    /*Algorithm Utils*/
    std::string mSolverPath, mTspPathPrefix, mUniqueId, mTspSolType;
    unsigned mSeed;
    std::chrono::time_point<std::chrono::system_clock> mStart, mEnd;

    virtual void d_solve()=0;
    virtual void init();
    virtual void pad();
    virtual void unpad();
    virtual void terminate() = 0;
public:
    std::vector<std::vector<std::string>> info;
    cppied_algorithm(std::vector<std::vector<double>>& pCoverage,
                     std::vector<std::vector<double>>& pPod,
                     std::vector<std::vector<unsigned short>>& pSeabed,
                     std::vector<csp::active_node>& pSolution,
                        double pRequiredCoverage,
                        double pTolerance,
                        unsigned pLateralRange,
                        csp::active_node& pStartNode,
                        std::string& pSolverPath,
                        std::string& pTspPathPrefix,
                        std::string& pUniqueId,
                        unsigned pPadSize,
                        double pPadValue,
                        unsigned pSeed) : mCoverage(pCoverage), mPod(pPod), mSeabed(pSeabed), mSolution(pSolution),
                        mRequiredCoverage(pRequiredCoverage), mTolerance(pTolerance), mLateralRange(pLateralRange),
                        mStartNode(pStartNode), mSolverPath(pSolverPath), mSolutionCost(0), mTspPathPrefix(pTspPathPrefix),
                        mUniqueId(pUniqueId), mPadSize(std::max(pPadSize, pLateralRange)), mPadValue(pPadValue), info(),
                        mGridDimensions(std::make_pair(pCoverage[0].size(), pCoverage.size())), mSeed(pSeed), mTspSolType("tour"){}
    ~cppied_algorithm();

    /*Solver Methods*/
    void solve(){
        init();
        d_solve();
        terminate();
    }
    virtual void write_tsp_instance() = 0;
    virtual void read_tsp_solution() = 0;
    void solve_tsp(const std::optional<std::string>& pSolverCommand=std::nullopt);
    template<class ForwardIter> unsigned unfeasible_lexicographic_cost(ForwardIter begin, ForwardIter end){
        unsigned cost=0;
        for (; begin!=std::prev(end); begin++)
            cost += get_distance(*begin, *std::next(begin), false);
        return cost;
    }
    template<class ForwardIter> unsigned feasible_lexicographic_cost(ForwardIter begin, ForwardIter end){
        unsigned cost = 0;
        for (; begin!=std::prev(end); begin++)
            cost += cLVALUE + ((*begin).dir != (*std::next(begin)).dir);
        return cost;
    }
    template<class ForwardIter> void swap_directions(ForwardIter begin, ForwardIter end){
        for (; begin != end; begin++)
            (*begin).dir = cSWAP[(*begin).dir];
    }
    template<typename T> friend void utils::unpad_matrix(std::vector<std::vector<T>>& pMatrix, int pPadSize);
    void reset_coverage();
    template<class ForwardIter> void set_insertion_coverage(ForwardIter begin, ForwardIter end){
        unsigned short i, j;
        for (; begin != end; begin++){
            if ((*begin).dir>1) {
                j = (*begin).x;
                for (i = (*begin).y-mLateralRange; i<(*begin).y + mLateralRange; i++)
                    mCoverage[i][j] += (1.0f - mCoverage[i][j]) * map_pod(i, j, {.x = j, .y = (*begin).y, .dir=2});
            } else {
                i = (*begin).y;
                for (j = (*begin).x-mLateralRange; j < (*begin).x+mLateralRange; ++j)
                    mCoverage[i][j] += (1.0f - mCoverage[i][j]) * map_pod(i, j,{.x = (*begin).x, .y = i, .dir=0});
            }
        }
    }
    template<class ForwardIter> void set_extraction_coverage(ForwardIter begin, ForwardIter end){
        unsigned short i,j;
        for (; begin!=end; begin++) {
            if (begin.dir > 1) {
                j = begin.x;
                for (i = begin.y - mLateralRange; i < begin.y + mLateralRange; i++) {
                    mCoverage[i][j] -= map_pod(i, j, {.x = j, .y = (*begin).y, .dir=2});
                    mCoverage[i][j] /= (1.0f - map_pod(i, j, {.x = j, .y = (*begin).y, .dir=2}));
                }
            } else {
                i = (*begin).y;
                for (j = (*begin).x - mLateralRange; j < (*begin).x + mLateralRange; ++j) {
                    mCoverage[i][j] -= map_pod(i, j, {.x = (*begin).x, .y = i, .dir=0});
                    mCoverage[i][j] /= (1.0f - map_pod(i, j, {.x = (*begin).x, .y = i, .dir=0}));
                }
            }
        }
    }
    virtual void set_lexicographic_cost(){mSolutionCost = feasible_lexicographic_cost(mSolution.begin(), mSolution.end());}
    [[nodiscard]] unsigned get_lexicographic_cost() const{return mSolutionCost;}
    [[nodiscard]] std::pair<unsigned, unsigned> get_bi_cost() const{
        unsigned n_turns = mSolutionCost - ((mSolution.size()-1)*cLVALUE);
        return std::make_pair((mSolutionCost-n_turns)/cLVALUE, n_turns);
    }
    unsigned cost(std::vector<csp::active_node>& pSolution) {
        unsigned n_turns=0;
        for (int i=0; i<pSolution.size()-1; i++){
            n_turns += (pSolution[i].dir != pSolution[i+1].dir);
        }
        return (cLVALUE*(pSolution.size()-1)) + n_turns;
    }
    [[nodiscard]] const std::vector<csp::active_node>& get_solution() const{return mSolution;}
    [[nodiscard]] const std::vector<std::vector<double>>& get_coverage() const{return mCoverage;}
    void push_shortest_path(const csp::active_node& pDestination);
    [[nodiscard]] inline double map_pod(const unsigned short line, const unsigned short column,
                                       const csp::active_node& agent_position) const{
        unsigned short distance = (unsigned short)std::max(std::floor(std::abs(agent_position.x-0.5-column)),
                                     std::floor(std::abs(agent_position.y-0.5-line)));
        return mPod[mSeabed[line][column]][distance];
    }
    [[nodiscard]] inline bool node_cover(const unsigned short line, const unsigned short column,
                           const csp::active_node& agent_position) const{
        if (agent_position.dir >1){
            return agent_position.x==column && std::floor(std::abs(agent_position.y-0.5-line)) < mLateralRange;
        } else{
            return agent_position.y==line && std::floor(std::abs(agent_position.x-0.5-column)) < mLateralRange;
        }
    }
    [[nodiscard]] inline unsigned short get_distance(const csp::active_node& pSource, const csp::active_node& pDestination,
                                 bool pAcceptZeroDistance) const
    {
        if (pAcceptZeroDistance && (pSource.dir==pDestination.dir) && (pSource.y == pDestination.y) && (pSource.x == pDestination.x))
            return 0;

        if (pSource.dir == 2){
            if (pDestination.dir == 2){
                if (pSource.y !=  pDestination.y){
                    if (pSource.x<pDestination.x) return (2*cLVALUE + 2) + ((abs(pSource.y-pDestination.y)-1)*cLVALUE) + ((pDestination.x-pSource.x-1)*cLVALUE);
                    else return (4*cLVALUE + 4) + ((pSource.x-pDestination.x)*cLVALUE)+ (abs(abs(pDestination.y-pSource.y)-2)*cLVALUE);}

                else{
                    if (pSource.x<pDestination.x) return (pDestination.x-pSource.x)*cLVALUE;
                    else return (4*cLVALUE + 4) + ((pSource.x-pDestination.x)*cLVALUE);}
            }
            else if (pDestination.dir == 0){
                if (pSource.y < pDestination.y){
                    if (pSource.x<pDestination.x) return (3*cLVALUE + 3) + ((pDestination.y-pSource.y)*cLVALUE) + (abs(pDestination.x-pSource.x-2)*cLVALUE);
                    else return (3*cLVALUE + 3) + ((pDestination.y-pSource.y)*cLVALUE) + ((pSource.x-pDestination.x)*cLVALUE);}
                else if (pSource.y > pDestination.y){
                    if (pSource.x<pDestination.x) return cLVALUE + 1 + ((pSource.y-pDestination.y-1)*cLVALUE) + ((pDestination.x-pSource.x-1)*cLVALUE);

                    else if  (pSource.x==pDestination.x) return (pSource.y< mGridDimensions.first-mPadSize or pSource.y-pDestination.y>1) ? 3*cLVALUE + 3 + (abs(pSource.y-pDestination.y-2)*cLVALUE) : 6*cLVALUE + 5;
                    else return (pSource.y< mGridDimensions.first-mPadSize or pDestination.y<(mGridDimensions.first-1-mPadSize)) ? (3*cLVALUE + 3) + (abs(pSource.y-pDestination.y-2)*cLVALUE) + ((pSource.x-pDestination.x)*cLVALUE) : 5*cLVALUE + 5 + (pSource.x - pDestination.x-1)*cLVALUE;
                }
                else{
                    if (pSource.x<pDestination.x) return (pSource.x<mGridDimensions.second-1-mPadSize) ? (3*cLVALUE + 3) + (abs(pDestination.x-pSource.x-2)*cLVALUE) : 6*cLVALUE + 5;
                    else return (3*cLVALUE + 3) + ((pSource.x-pDestination.x)*cLVALUE);
                }
            }
            else if (pDestination.dir == 1){
                if (pSource.y > pDestination.y){
                    if (pSource.x>pDestination.x) return (3*cLVALUE + 3) + ((pSource.y-pDestination.y-1)*cLVALUE) + ((pSource.x-pDestination.x)*cLVALUE);
                    else if  (pSource.x<pDestination.x) return ((pSource.x < mGridDimensions.second-1-mPadSize)) ? (3*cLVALUE + 3) + ((pSource.y-pDestination.y-1)*cLVALUE) + (abs(pDestination.x-pSource.x-2)*cLVALUE) : 5*cLVALUE + 5 + (abs(pSource.y-pDestination.y-2)*cLVALUE);
                    else return (3*cLVALUE + 3) + ((pSource.y-pDestination.y-1)*cLVALUE);
                }
                else if (pSource.y < pDestination.y){
                    if (pSource.x<pDestination.x) return cLVALUE + 1 + ((pDestination.x-pSource.x-1)*cLVALUE) + ((pDestination.y-pSource.y)*cLVALUE);
                    else if (pDestination.x<pSource.x) return (3*cLVALUE + 3) + ((pSource.x-pDestination.x)*cLVALUE) + ((pDestination.y-pSource.y-1)*cLVALUE);
                    else return (3*cLVALUE + 3) + ((pDestination.y-pSource.y-1)*cLVALUE);
                }
                else{
                    if (pSource.x<pDestination.x) return cLVALUE + 1 + ((pDestination.x-pSource.x-1)*cLVALUE);
                    else if (pDestination.x<pSource.x) return (3*cLVALUE + 3) + ((pSource.x-pDestination.x)*cLVALUE) + cLVALUE;
                    else return (pSource.y==0) ? 6*cLVALUE + 5 : 4*cLVALUE + 3;
                }
            }
            else{
                if (pSource.y !=  pDestination.y) return (2*cLVALUE + 2) + (abs(pDestination.x-pSource.x)*cLVALUE) + (abs((abs(pDestination.y-pSource.y)-1))*cLVALUE);
                else return (pSource.x < mGridDimensions.second-1-mPadSize or pSource.x!=pDestination.x) ? (4*cLVALUE + 4) + (abs(abs(pDestination.x-pSource.x)-1)*cLVALUE) : 7*cLVALUE + 6;
            }
        }
        else if (pSource.dir== 0){
            if (pDestination.dir == 2){
                if (pSource.y < pDestination.y){
                    if (pSource.x<pDestination.x) return (3*cLVALUE + 3) + ((pDestination.y-pSource.y-1)*cLVALUE) + ((pDestination.x-pSource.x-1)*cLVALUE);
                    else if (pDestination.x<pSource.x) return (3*cLVALUE + 3) +  ((pDestination.y-pSource.y-1)*cLVALUE) + ((pSource.x-pDestination.x-1)*cLVALUE);
                    else return (pSource.x>0) ? 4*cLVALUE + 3 + ((pDestination.y-pSource.y-1)*cLVALUE) : 5*cLVALUE + 5 + (abs(pDestination.y-pSource.y-2)*cLVALUE);
                }
                else if (pDestination.y < pSource.y){
                    if (pSource.x < pDestination.x) return cLVALUE + 1 + ((pDestination.x-pSource.x)*cLVALUE) + ((pSource.y-pDestination.y)*cLVALUE);
                    else if (pDestination.x<pSource.x) return (3*cLVALUE + 3) + ((pSource.x-pDestination.x-1)*cLVALUE) + ((pSource.y-pDestination.y-1)*cLVALUE);
                    else return cLVALUE + 1 + ((pSource.y-pDestination.y)*cLVALUE);
                }
                else{
                    if (pDestination.x<pSource.x) return (pSource.y>0) ? (3*cLVALUE + 3) + ((pSource.x-pDestination.x-1)*cLVALUE) + cLVALUE : 5*cLVALUE + 5 + (abs(pSource.x-pDestination.x-2)*cLVALUE);
                    else return cLVALUE + 1 + ((pDestination.x-pSource.x)*cLVALUE);
                }
            }
            else if (pDestination.dir == 0){
                if (pSource.y < pDestination.y){
                    if (pSource.x>pDestination.x) return (4*cLVALUE + 4) + (abs(pSource.x-pDestination.x-2)*cLVALUE) + ((pDestination.y-pSource.y)*cLVALUE);
                    else if (pSource.x<pDestination.x) return (4*cLVALUE + 4) + (abs(pDestination.x-pSource.x-2)*cLVALUE) + ((pDestination.y-pSource.y)*cLVALUE);
                    else return (4*cLVALUE + 4) + ((pDestination.y-pSource.y)*cLVALUE);
                }
                else if (pSource.y > pDestination.y){
                    if (pSource.x<pDestination.x) return (2*cLVALUE + 2) + ((pDestination.x-pSource.x-1)*cLVALUE) + ((pSource.y-pDestination.y-1)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (2*cLVALUE + 2) + ((pSource.x-pDestination.x-1)*cLVALUE) + ((pSource.y-pDestination.y-1)*cLVALUE);
                    else return (pSource.y-pDestination.y)*cLVALUE;
                }
                else{
                    if (pSource.x<pDestination.x) return (4*cLVALUE + 4) + (abs(pDestination.x-pSource.x-2)*cLVALUE);
                    else return (4*cLVALUE + 4) + (abs(pSource.x-pDestination.x-2)*cLVALUE);
                }
            }
            else if (pDestination.dir == 1){
                if (pSource.y < pDestination.y){
                    if (pSource.x<pDestination.x) return (2*cLVALUE + 2) + ((pDestination.x-pSource.x-1)*cLVALUE) + ((pDestination.y-pSource.y)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (2*cLVALUE + 2) + ((pSource.x-pDestination.x-1)*cLVALUE) + ((pDestination.y-pSource.y)*cLVALUE);
                    else return (4*cLVALUE + 4) + ((pDestination.y-pSource.y-1)*cLVALUE);
                }
                else if (pSource.y > pDestination.y){
                    if (pSource.x<pDestination.x) return (2*cLVALUE + 2) + ((pSource.y-pDestination.y)*cLVALUE) + ((pDestination.x-pSource.x-1)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (2*cLVALUE + 2) + ((pSource.y-pDestination.y)*cLVALUE) + ((pSource.x-pDestination.x-1)*cLVALUE);
                    else return (4*cLVALUE + 4) + ((pSource.y-pDestination.y-1)*cLVALUE);
                }
                else{
                    if (pSource.x<pDestination.x) return (2*cLVALUE + 2) + ((pDestination.x-pSource.x-1)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (2*cLVALUE + 2) + ((pSource.x-pDestination.x-1)*cLVALUE);
                    else return (pSource.y == 0) ? 7*cLVALUE + 6 : 5*cLVALUE +4;
                }
            }
            else {
                if (pSource.y < pDestination.y){
                    if (pSource.x<pDestination.x) return (3*cLVALUE + 3) + ((pDestination.y-pSource.y-1)*cLVALUE) + ((pDestination.x-pSource.x)*cLVALUE);
                    else if (pDestination.x<pSource.x) return (pSource.x < mGridDimensions.second-mPadSize or (pSource.x-pDestination.x>1)) ? (3*cLVALUE + 3) + (abs(pSource.x-pDestination.x-2)*cLVALUE) + ((pDestination.y-pSource.y-1)*cLVALUE) : 5*cLVALUE + 5 + (abs(pDestination.y-pSource.y-2)*cLVALUE);
                    else return (3*cLVALUE + 3) + ((pDestination.y-pSource.y-1)*cLVALUE);
                }
                else if (pDestination.y<pSource.y){
                    if (pSource.x<pDestination.x) return (3*cLVALUE + 3) + ((pSource.y-pDestination.y-1)*cLVALUE) + ((pDestination.x-pSource.x)*cLVALUE);
                    else if (pDestination.x<pSource.x) return cLVALUE + 1+ (abs(pSource.x-pDestination.x-1)*cLVALUE) + ((pSource.y-pDestination.y)*cLVALUE);
                    else return (3*cLVALUE + 3) + ((pSource.y-pDestination.y-1)*cLVALUE);
                }
                else{
                    if (pDestination.x>pSource.x) return (3*cLVALUE + 3) + ((pDestination.x-pSource.x)*cLVALUE) + cLVALUE;
                    else if (pSource.x>pDestination.x) return cLVALUE + 1 + ((pSource.x-pDestination.x-1)*cLVALUE);
                    else return (pSource.y==0) ? 6*cLVALUE + 5 : 4*cLVALUE + 3;
                }
            }
        }
        else if (pSource.dir== 1){
            if (pDestination.dir == 2){
                if (pSource.y < pDestination.y){
                    if (pSource.x<pDestination.x) return cLVALUE + 1 + ((pDestination.x-pSource.x)*cLVALUE) + ((pDestination.y-pSource.y-1)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (pSource.y < mGridDimensions.first-1-mPadSize) ? (3*cLVALUE + 3) + ((pSource.x-pDestination.x-1)*cLVALUE) + (abs(pDestination.y-pSource.y-2)*cLVALUE) : 5*cLVALUE + 5 + (abs(pSource.x-pDestination.x-2)*cLVALUE);
                    else return cLVALUE + 1 + ((pDestination.y-pSource.y-1)*cLVALUE);
                }
                else if (pSource.y > pDestination.y){
                    if (pSource.x<pDestination.x) return (3*cLVALUE + 3) + ((pDestination.x-pSource.x-1)*cLVALUE) + ((pSource.y-pDestination.y)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (3*cLVALUE + 3) + ((pSource.x-pDestination.x-1)*cLVALUE) + ((pSource.y-pDestination.y)*cLVALUE);
                    else return (pSource.x>0) ? 4*cLVALUE + 3 + ((pSource.y-pDestination.y)*cLVALUE) : 5*cLVALUE + 5 + ((pSource.y-pDestination.y-1)*cLVALUE);
                }
                else{
                    if (pSource.x<pDestination.x) return (3*cLVALUE + 3) + ((pDestination.x-pSource.x-1)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (3*cLVALUE + 3) + ((pSource.x-pDestination.x-1)*cLVALUE);
                    else return (pSource.x==0) ? 6*cLVALUE + 5 : 4*cLVALUE + 3;
                }
            }
            else if (pDestination.dir == 0){
                if (pSource.y < pDestination.y){
                    if (pSource.x<pDestination.x) return (2*cLVALUE + 2) + ((pDestination.y-pSource.y)*cLVALUE) + ((pDestination.x-pSource.x-1)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (2*cLVALUE + 2) + ((pDestination.y-pSource.y)*cLVALUE) + ((pSource.x-pDestination.x-1)*cLVALUE);
                    else return (4*cLVALUE + 4) + ((pDestination.y-pSource.y-1)*cLVALUE);
                }
                else if (pSource.y > pDestination.y){
                    if (pSource.x<pDestination.x) return (2*cLVALUE + 2) + ((pSource.y-pDestination.y)*cLVALUE) + ((pDestination.x-pSource.x-1)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (2*cLVALUE + 2) + ((pSource.y-pDestination.y)*cLVALUE) + ((pSource.x-pDestination.x-1)*cLVALUE);
                    else return (4*cLVALUE + 4) + ((pSource.y-pDestination.y-1)*cLVALUE);
                }
                else{
                    if (pSource.x<pDestination.x) return (2*cLVALUE + 2) + ((pDestination.x-pSource.x-1)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (2*cLVALUE + 2) + ((pSource.x-pDestination.x-1)*cLVALUE);
                    else return (pSource.y == mGridDimensions.first-1-mPadSize) ? 7*cLVALUE + 6 : 5*cLVALUE +4;
                }
            }
            else if (pDestination.dir == 1){
                if (pSource.y < pDestination.y){
                    if (pSource.x<pDestination.x) return (2*cLVALUE + 2) + ((pDestination.y-pSource.y-1)*cLVALUE)+((pDestination.x-pSource.x-1)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (2*cLVALUE + 2) + ((pDestination.y-pSource.y-1)*cLVALUE)+((pSource.x-pDestination.x-1)*cLVALUE);
                    else return (pDestination.y-pSource.y)*cLVALUE;
                }
                else if (pSource.y > pDestination.y){
                    if (pSource.x<pDestination.x) return (4*cLVALUE + 4) + (abs(pDestination.x-pSource.x-2)*cLVALUE) + ((pSource.y-pDestination.y)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (4*cLVALUE + 4) +(abs(pSource.x-pDestination.x-2)*cLVALUE) + ((pSource.y-pDestination.y)*cLVALUE);
                    else return (4*cLVALUE + 4) + ((pSource.y-pDestination.y)*cLVALUE);
                }
                else return (pSource.x<pDestination.x) ? (4*cLVALUE + 4) + (abs(pDestination.x-pSource.x-2)*cLVALUE) : (4*cLVALUE + 4) + (abs(pSource.x-pDestination.x-2)*cLVALUE);
            }
            else {
                if (pSource.y < pDestination.y){
                    if (pSource.x<pDestination.x) return (pSource.y < mGridDimensions.first-1-mPadSize) ? (3*cLVALUE + 3) + ((pDestination.x-pSource.x)*cLVALUE) + (abs(pDestination.y-pSource.y-2)*cLVALUE) : 5*cLVALUE + 5 + (pDestination.x-pSource.x-1)*cLVALUE;
                    else if (pSource.x>pDestination.x) return cLVALUE + 1 + ((pDestination.y-pSource.y-1)*cLVALUE) + ((pSource.x-pDestination.x-1)*cLVALUE);
                    else return (pSource.y < mGridDimensions.first-1-mPadSize) ? (3*cLVALUE + 3) + (abs(pDestination.y-pSource.y-2)*cLVALUE) : 6*cLVALUE + 5;
                }
                else if (pSource.y > pDestination.y){
                    if (pSource.x<pDestination.x) return (3*cLVALUE + 3) + ((pDestination.x-pSource.x)*cLVALUE) + ((pSource.y-pDestination.y)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (3*cLVALUE + 3) + (abs(pSource.x-pDestination.x-2)*cLVALUE) + ((pSource.y-pDestination.y)*cLVALUE);
                    else return 3*cLVALUE + 3 + ((pSource.y-pDestination.y)*cLVALUE);
                }
                else{
                    if (pSource.x<pDestination.x) return (3*cLVALUE + 3) + ((pDestination.x-pSource.x)*cLVALUE);
                    else if (pSource.x>pDestination.x) return (pDestination.x < mGridDimensions.second-1-mPadSize) ? (3*cLVALUE + 3) + (abs(pSource.x-pDestination.x-2)*cLVALUE) : 6*cLVALUE + 5;
                    else return 3*cLVALUE + 3;
                }
            }
        }
        else {
            if (pDestination.dir == 2){
                if (pSource.y !=  pDestination.y) return (2*cLVALUE + 2) + (abs(pDestination.x-pSource.x)*cLVALUE) + ((abs(pDestination.y-pSource.y)-1)*cLVALUE);
                else{
                    if (pSource.x !=0 or pSource.x!=pDestination.x) return (4*cLVALUE + 4) + (abs(abs(pDestination.x-pSource.x)-1)*cLVALUE);
                    else return (pSource.x==0) ? 7*cLVALUE + 5 : 5*cLVALUE +4;
                }
            }
            else if (pDestination.dir == 0){
                if (pSource.y < pDestination.y){
                    if (pSource.x>pDestination.x) return (3*cLVALUE + 3) + ((pDestination.y-pSource.y)*cLVALUE) + ((pSource.x-pDestination.x-1)*cLVALUE);
                    else if (pDestination.x>pSource.x) return (3*cLVALUE + 3) + ((pDestination.y-pSource.y)*cLVALUE) + ((pDestination.x-pSource.x-1)*cLVALUE);
                    else return (pSource.x>0) ? (3*cLVALUE + 3) + ((pDestination.y-pSource.y)*cLVALUE) + cLVALUE : 5*cLVALUE + 5 + ((pDestination.y-pSource.y-1)*cLVALUE);
                }
                else if (pSource.y > pDestination.y){
                    if (pSource.x>=pDestination.x) return cLVALUE + 1 + ((pSource.y-pDestination.y-1)*cLVALUE) + ((pSource.x-pDestination.x)*cLVALUE);
                    else return ((pSource.y-pDestination.y) > 1) ? (3*cLVALUE + 3) + ((pSource.y-pDestination.y-2)*cLVALUE) + ((pDestination.x-pSource.x-1)*cLVALUE) : (pSource.y < mGridDimensions.first-mPadSize) ? 4*cLVALUE + 3 + ((pDestination.x-pSource.x-1)*cLVALUE)
                                                                                                                                                        : ((pDestination.x-pSource.x) > 1) ? 5*cLVALUE + 5 + ((pDestination.x-pSource.x-2)*cLVALUE) : 6*cLVALUE + 5;
                }
                else{
                    if (pSource.x>pDestination.x) return (3*cLVALUE + 3) + ((pSource.x-pDestination.x-1)*cLVALUE);
                    else if (pDestination.x > pSource.x) return (3*cLVALUE + 3) + ((pDestination.x-pSource.x-1)*cLVALUE);
                    else return (pSource.x!=0) ? 4*cLVALUE + 3 : 6*cLVALUE + 5;
                }
            }
            else if (pDestination.dir == 1){
                if (pSource.y > pDestination.y){
                    if (pSource.x>pDestination.x) return (3*cLVALUE + 3) + ((pSource.y-pDestination.y-1)*cLVALUE) + ((pSource.x-pDestination.x-1)*cLVALUE);
                    else if (pDestination.x>pSource.x) return (3*cLVALUE + 3) + ((pDestination.x-pSource.x-1)*cLVALUE) + ((pSource.y-pDestination.y-1)*cLVALUE);
                    else return (pSource.x!=0) ? (3*cLVALUE + 3) + ((pSource.y-pDestination.y-1)*cLVALUE) + cLVALUE : 5*cLVALUE + 5 + (abs(pSource.y-pDestination.y-2)*cLVALUE);
                }
                else if (pSource.y < pDestination.y){
                    if (pSource.x>pDestination.x) return cLVALUE + 1 + ((pSource.x-pDestination.x)*cLVALUE) + ((pDestination.y-pSource.y)*cLVALUE);
                    else if (pDestination.x>pSource.x) return (3*cLVALUE + 3) +((pDestination.x-pSource.x-1)*cLVALUE) + ((pDestination.y-pSource.y-1)*cLVALUE);
                    else return cLVALUE + 1 + ((pDestination.y-pSource.y)*cLVALUE);
                }
                else{
                    if (pSource.x>pDestination.x) return cLVALUE + 1 + ((pSource.x-pDestination.x)*cLVALUE);
                    else if (pDestination.x>pSource.x) return (pSource.y>0) ? 4*cLVALUE + 3 + ((pDestination.x-pSource.x-1)*cLVALUE) : 5*cLVALUE + 5 + (abs(pDestination.x-pSource.x-2)*cLVALUE);
                    else return cLVALUE + 1;
                }
            }
            else{
                if (pSource.y !=  pDestination.y){
                    if (pSource.x>pDestination.x) return (2*cLVALUE + 2) + ((pSource.x-pDestination.x-1)*cLVALUE)+ ((abs(pSource.y-pDestination.y)-1)*cLVALUE);
                    else if (pDestination.x>pSource.x) return (4*cLVALUE + 4) + (abs(abs(pDestination.y-pSource.y)-2)*cLVALUE) + ((pDestination.x-pSource.x)*cLVALUE);
                    else return (4*cLVALUE + 4) + (abs(abs(pDestination.y-pSource.y)-2)*cLVALUE);
                }
                else{
                    if (pSource.x>pDestination.x) return (pSource.x-pDestination.x)*cLVALUE;
                    else return (4*cLVALUE + 4) + ((pDestination.x-pSource.x)*cLVALUE);
                }
            }
        }
    }
    [[nodiscard]] bool inline valid_node(csp::active_node pNode) const{
        if (pNode.dir>1)
            return (pNode.x>=mPadSize & pNode.x < mGridDimensions.second-mPadSize &
                    pNode.y>=mPadSize & pNode.y <= mGridDimensions.first-mPadSize);
        else
            return (pNode.x>=mPadSize & pNode.x <= mGridDimensions.second-mPadSize &
                    pNode.y>=mPadSize & pNode.y < mGridDimensions.first-mPadSize);
    }
    template<typename T>
    [[nodiscard]] long get_resolution_time() const {return std::chrono::duration_cast<T>(mEnd - mStart).count();}
    std::string get_tsp_options();

    /*Constants*/
    const unsigned short cINF_EDGE = 65535;
    const unsigned short cLVALUE = 10;
    const double cFLAT_EPSILON = 0.01;
    const std::vector<std::vector<csp::active_node_neighbors>> cNEIGHBORS = {{{0, -1, 0}, {-1,0, 3}, {0, 0, 2}},
                                         {{0, 1, 1}, {-1, 1, 3}, {0, 1, 2}},
                                         {{1, -1, 0}, {1, 0, 1}, {1, 0, 2}},
                                         {{0, -1, 0}, {0, 0, 1}, {-1, 0, 3}}};
    const unsigned short cSWAP[4] {1, 0, 3, 2};
};
#endif // CPPIED_ALGORITHM_HPP
