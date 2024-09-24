#ifndef LS2_HPP
#define LS2_HPP

// Include any necessary headers
#include "../cppied_algorithm.hpp"

// Define any necessary classes or functions
namespace local_search{
    class ls2 : public cppied_algorithm
    {
    private:
        /* data */
        std::list<csp::active_node> mMutableSolution;
        std::vector<std::vector<double>> mMutableCoverage;
        unsigned mLocalCost;

        /*Algorithm parameters*/
        int mNPerturbations, mNImprovements, mMaxTime, mMaxSubsequences;
        unsigned short mNeighbourhoodWidth;
        std::vector<unsigned> mTspSolution;
        std::vector<csp::active_node_neighbors> mLateralNbh;
        std::vector<csp::active_node_neighbors> mVerticalNbh;

        /*Algorithm Utils*/
        std::mt19937 mRng;
        std::uniform_int_distribution<int> mDirectionRNodes;
        std::uniform_int_distribution<int> mXCoordHorizontalRNodes;
        std::uniform_int_distribution<int> mYCoordHorizontalRNodes;
        std::uniform_int_distribution<int> mXCoordVerticalRNodes;
        std::uniform_int_distribution<int> mYCoordVerticalRNodes;
        std::uniform_real_distribution<double> mRealGenerator;
    protected:
        void d_solve() override;
        void init() override;
        /*TSP utils*/
        void write_tsp_initial_tour();
        void read_tsp_solution() override;
        std::tuple<bool, std::size_t, std::size_t> read_tsp_solution(std::ifstream& ostream);
        void write_tsp_instance() override;
        void write_tsp_matrix(std::ofstream& ostream);
        void write_tsp_start_edges(std::ofstream& ostream);
        void write_tsp_dummy_edges(std::ofstream& ostream);
        void write_tsp_node_edges(std::ofstream& ostream, const csp::active_node& pNode, std::size_t pPosition);
        void tsp_to_solution(const std::vector<csp::active_node>& pSolution);
        void build_feasible();
    public:
        std::vector<std::vector<std::string>> ls_info;
        std::vector<std::vector<std::string>> twoopt_info;
        std::vector<std::vector<std::string>> solver_info;
        /*Constructor Methods*/
        void generate_nbh();
        template<typename... Args>
        ls2(int pNPerturbations,
            int pNImprovements, int pMaxTime,
            int pNeighbourhoodWidth, int pMaxSubsequences,
            Args&&... pCppiedArgs) :
                cppied_algorithm(std::forward<Args>(pCppiedArgs)...),
                mNPerturbations(pNPerturbations),
                mNImprovements(pNImprovements), mMaxTime(pMaxTime), mMutableSolution(), mMaxSubsequences(pMaxSubsequences), mNeighbourhoodWidth(pNeighbourhoodWidth),
                mLateralNbh(), mVerticalNbh(), mTspSolution(), mLocalCost(0)
        {
            ls_info.push_back({"iteration", "cost", "iteration_time"});
            twoopt_info.push_back({"iteration", "cost", "iteration_time"});
            solver_info.push_back({"iteration", "cost", "n_moves", "n_turns", "iteration_time"});
            mLocalCost = mSolutionCost;
            mRng.seed(mSeed);
            mMutableCoverage = mCoverage;
            mTspSolution.reserve(mSolution.size()+1);
            mMutableSolution.clear();
            std::copy(mSolution.begin(), mSolution.end(), std::back_inserter(mMutableSolution));
            mXCoordHorizontalRNodes = std::uniform_int_distribution<int>(0, (int)mCoverage[0].size() - 1);
            mYCoordHorizontalRNodes = std::uniform_int_distribution<int>(0, (int)mCoverage.size());
            mXCoordVerticalRNodes = std::uniform_int_distribution<int>(0, (int)mCoverage[0].size());
            mYCoordVerticalRNodes = std::uniform_int_distribution<int>(0, (int)mCoverage.size()-1);
            mDirectionRNodes = std::uniform_int_distribution<int>(0, 3);
            mRealGenerator = std::uniform_real_distribution(0.0,1.0);
            generate_nbh();
        }
        ~ls2()=default;

        void pad() override;
        void unpad() override;
        void set_lexicographic_cost() override;

        /*Solver methods*/
        void terminate() override;
        void solve_lkh();
        void optimize_path();
        bool solve_perturb_and_optimize(std::size_t pIter);
        bool improvement_procedure();
        void perturbation_procedure();
        bool two_opt_procedure();
        bool restart(bool pPerturbate, bool pHeuristicImprovement);
        template<class ForwardIter> unsigned set_mutable_extraction_coverage(ForwardIter begin, ForwardIter end);
        void set_mutable_insertion_coverage(csp::active_node& pNode);
        template<class ForwardIter> void set_mutable_insertion_coverage(ForwardIter begin, ForwardIter end);
        template<class ForwardIter> bool local_nbh_search_procedure(ForwardIter begin,
                                                                        ForwardIter end);
        void sequence_satisfy(const csp::active_node& n1, const csp::active_node& n2, const csp::active_node& ref, bool* response);
        template<class ForwardIter> bool random_nbh_search_procedure(ForwardIter begin, ForwardIter end);
        template<class ForwardIter>
        bool best_insertion_search(csp::active_node& n1, csp::active_node& n2, std::vector<ForwardIter>& pInsertionPositions,
                                   const bool* pNodesSatisfy, unsigned pMaxCost);
        template<class ForwardIter>
        unsigned best_insertion_search(csp::active_node& n1, std::vector<ForwardIter>& pInsertionPositions, unsigned max_cost);

        template<class ForwardIter>
        void filter_insertion_positions(const csp::active_node& ref, std::vector<ForwardIter>& container, unsigned pUpperbound);
        bool equal_solutions(const std::vector<csp::active_node>& pSolution) const;
        csp::active_node random_node_gen();

        /*Constants*/
        const double cEPSILON = 0.1;
};
} // namespace local_search

#endif // LS2_HPP
