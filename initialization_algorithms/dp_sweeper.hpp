#include "../cppied_algorithm.hpp"

#ifndef DP_SWEEPER_HPP
#define DP_SWEEPER_HPP

namespace init{
    class dp_sweeper: public cppied_algorithm{
    private:
        std::vector<std::pair<csp::active_node, csp::active_node>> mHSegments, mVSegments;
        std::vector<std::pair<double, std::pair<unsigned short, unsigned short>>> mHdpGain, mVdpGain;
        unsigned mNcovered, mNcells;
        double mPenalty;
        std::vector<unsigned> mTspSolution;
        std::chrono::time_point<std::chrono::system_clock> mTspStart, mTspEnd;
    protected:
        void d_solve() override;
        void terminate() override;
        static void opt_flip(std::vector<std::pair<csp::active_node, csp::active_node>>::iterator first,
                      std::vector<std::pair<csp::active_node, csp::active_node>>::iterator last);
    public:
        template<typename... Args>
        dp_sweeper(Args&&... pCppiedArgs) :
                cppied_algorithm(std::forward<Args>(pCppiedArgs)...), mHSegments(), mVSegments(), mTspSolution(), mHdpGain(),
                mVdpGain(), mPenalty(mTolerance/2), mNcovered(0)
        {
                    assertm(mSolution.size()==1, "Solver only support empty solution at initialization.");
                    mCoverage = std::vector<std::vector<double>>(mGridDimensions.first,
                                                                std::vector<double>(mGridDimensions.second, 0.0f));
                    mNcells = mGridDimensions.first * mGridDimensions.second;
                    mPenalty = mTolerance;
                    mHdpGain.resize(mGridDimensions.first + 1);
                    mVdpGain.resize(mGridDimensions.second + 1);
        }
        //TSP functions
        void write_tsp_instance() override;
        void read_tsp_solution() override;
        void read_tsp_tour(std::ifstream& pIfStream, std::size_t dimension);
        void read_tsp_edges(std::ifstream& pIfStream, std::size_t dimension);
        void write_tsp_matrix(std::ofstream& ostream);
        void write_tsp_start(std::ofstream& ostream);
        void write_h_tsp(std::ofstream& ostream);
        void write_v_tsp(std::ofstream& ostream);
        void tsp_to_solution();
        void push_tsp_first_segment();
        void push_tsp_horizontal_segment(std::size_t i);
        void push_tsp_vertical_segment(std::size_t i);

        //NN3Opt methods
        void nn3opt(bool pLocalSearch);
        void dps3opt(std::vector<std::pair<csp::active_node,csp::active_node>>& segment_path);
        unsigned dps3opt_move(std::vector<std::pair<csp::active_node,csp::active_node>>& segment_path,
                              std::size_t i, std::size_t j, std::size_t k, unsigned i_cost);
        unsigned dps2opt_move(std::vector<std::pair<csp::active_node, csp::active_node>>& segment_path,
                              std::size_t i, std::size_t j);
        void dps3opt_reconstruct_path(std::vector<std::pair<csp::active_node, csp::active_node>>& segment_path);

        void move_nearest_segment(std::vector<std::pair<csp::active_node,csp::active_node>>& pSegmentSet, std::size_t pos);

        //DP Functions
        void generate_feasible_segments();
        void genererate_segments_from_cardinality(unsigned short prev_h_cardinality, unsigned short prev_v_cardinality);
        void dp_v_forward();
        void dp_h_forward();
        void dp_v_backward();
        void dp_h_backward();
        std::pair<double, std::pair<unsigned short, unsigned short>> v_kadane(unsigned short col_num);
        std::pair<double, std::pair<unsigned short, unsigned short>> h_kadane(unsigned short line_num);
        void set_segment_insertion_coverage(const csp::active_node& pFirst, const csp::active_node& pLast);

        const double cEPSILON = 1e-6f;
    };
}

#endif //ALGORITHMS_DP_SWEEPER_HPP
