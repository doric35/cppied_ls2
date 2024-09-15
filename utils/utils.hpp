#include <vector>
#include <iostream>
#include <ostream>
#include <istream>
#include <iomanip>
#include <fstream>
#include <tuple>
#include <map>
#include <exception>
#include <stdexcept>
#include <chrono>
#include <cstdint>
#include <algorithm>
#include <functional>
#include <list>
#include <iterator>
#include <cassert>
#include <optional>
#include <random>
#include <cmath>
#include <filesystem>
#include <cassert>
#include <fcntl.h>
#include "toml.hpp"

#ifndef assertm
#define assertm(exp, msg) assert(((void)msg, exp))
#endif

namespace tsp{
    template<typename T>
    struct node{
        T encoding; bool reversed;
        unsigned id, rank;
    };

    template<typename T>
    struct segment_list{
        unsigned id, rank;
        bool reversed;
        std::list<node<T>> segment_nodes;
    };

    template<typename T>
    struct segment_vector{
        unsigned id,rank;
        bool reversed;
        std::vector<node<T>> segment_nodes;
    };

    template<typename T>
    struct tour_vector{
        unsigned dimension;
        std::vector<segment_vector<T>> tour;
    };

    template<typename T>
    struct tour_list{
        unsigned dimension;
        std::list<segment_list<T>> tour;
    };
}

namespace csp{
    const unsigned short cSWAP[4] {1, 0, 3, 2};
    struct active_node_neighbors{
        short x, y;
        unsigned short dir;
    };
    struct active_node{
        unsigned short x, y, dir;
        [[nodiscard]] active_node swap() const{return {x, y, cSWAP[dir]};}
        active_node operator+(const active_node& node) const{return {static_cast<unsigned short>(x+node.x), static_cast<unsigned short>(y+node.y), dir};}
        void operator+=(const active_node& node){x += node.x; y+=node.y;}
        void operator+=(const active_node_neighbors& dnode){x +=dnode.x; y+=dnode.y;dir = dnode.dir;}
        active_node operator+(const active_node_neighbors& dnode) const{return {static_cast<unsigned short>(x+dnode.x), static_cast<unsigned short>(y+dnode.y), dnode.dir};}
        bool operator==(const active_node& node) const{return x==node.x && y==node.y && dir==node.dir;}
        bool operator!=(const active_node& node) const{return x!=node.x || y!=node.y || dir!=node.dir;}
        active_node operator-(const active_node& node) const{return {static_cast<unsigned short>(x-node.x), static_cast<unsigned short>(y-node.y), dir};}
        void operator-=(const active_node& node){x -= node.x; y-=node.y;}
        void operator-=(const active_node_neighbors& dnode){x -=dnode.x; y-=dnode.y;dir = dnode.dir;}
        active_node operator-(const active_node_neighbors& dnode) const{return {static_cast<unsigned short>(x-dnode.x), static_cast<unsigned short>(y-dnode.y), dnode.dir};}
        //TODO: reverse operation on a node (i.e. reverse direction)
        friend std::istream &operator>>(std::istream& is, active_node& node){
            is >> node.x >> node.y >> node.dir;
            return is;
        }
        friend std::ostream &operator<<(std::ostream& os, const active_node& node){
            os << node.x << " " << node.y << " " << node.dir;
            return os;
        }
    };
}

namespace utils{
    template <typename T>
    class has_solve{
        //SFINAE to check for the presence of a method named 'solve'
        template <typename U>
        static auto test(int) -> decltype(std::declval<U>().solve(), std::true_type{});

        // Fallback for types that do not have a 'solve' method
        template <typename>
        static std::false_type test(...);

        public:
            static constexpr bool value = decltype(test<T>(0))::value;
    };

    template <typename T>
    int solve(T& pSolver, std::string& pErrPath, std::string& uuid){
        static_assert(has_solve<T>::value, "Solver must have a solve() method. \n");
        try {
            pSolver.solve();
        } catch (const std::exception& e) {
            std::ofstream err_stream(pErrPath + "/err/hypersearch_errors.txt", std::ios::app);
            if (err_stream) err_stream << uuid << "," << e.what() << std::endl;
            else std::cerr << uuid << "," << e.what() << std::endl;
            err_stream.close();
            return 1;
        } catch (...) {
            std::ofstream err_stream(pErrPath + "/err/hypersearch_errors.txt", std::ios::app);
            if (err_stream) err_stream << uuid << "," << "Caught an unknown exception!" << std::endl;
            else std::cerr << uuid << "," << "Caught an unknown exception!" << std::endl;
            err_stream.close();
            return 1;
        }
        return 0;
    }

    inline void string_strip(const char pChar, std::string& pString){
        if (pString[0] == pChar) pString.erase(pString.begin());
        if (pString[pString.size()-1] == pChar) pString.erase(std::prev(pString.end()));
    }

    template <typename T>
    void read_matrix(const std::string& pPath, std::vector<std::vector<T>>& pMatrix){
        std::ifstream mat_stream; std::ostringstream err;
        mat_stream.open(pPath, std::ios::in);
        int m, n; T value; char hash; std::string eof;
        mat_stream >> hash >> n >> m;
        if (!mat_stream.is_open()) {
            err << "Could not open matrix file " << pPath << std::endl;
            throw std::runtime_error(err.str());
        }
        if (hash != '#') {
            err << "Invalid header marker " << pPath << std::endl;
            throw std::runtime_error(err.str());
        }
        pMatrix.resize(n, std::vector<T>(m));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) { mat_stream >> pMatrix[i][j]; }
        }
        std::getline(mat_stream, eof);
        if (std::getline(mat_stream, eof)) {
            if (eof != "# EOF" && eof != "# EOF\n" && eof !="# EOF \n") {
                err << "Invalid EOF marker " << eof << " " << pPath << std::endl;
                throw std::runtime_error(err.str());
            }
        }
        else throw std::runtime_error("Failed to read EOF marker \n");
        mat_stream.close();
    }

    template <typename T>
    void write_matrix(const std::string& pPath, const std::vector<std::vector<T>>& pMatrix){
        std::ofstream mat_stream; std::ostringstream err;
        mat_stream.open(pPath, std::ios::out);
        if (!mat_stream.is_open()) {
            err << "Could not open matrix file " << pPath << std::endl;
            throw std::runtime_error(err.str());
        }
        int n = pMatrix.size(); int m = n > 0 ? pMatrix[0].size() : 0;
        mat_stream << "# " << n << " " << m << std::endl;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m-1; ++j)
                mat_stream << pMatrix[i][j] << " ";
            mat_stream << pMatrix[i][m-1] << std::endl;
        }

        // Write the EOF marker
        mat_stream << "# EOF" << std::endl;
        mat_stream.close();
    }

    inline void write_info(const std::string& pPath, const std::vector<std::vector<std::string>>& pInfo,
                           std::ios_base::openmode pMode=std::ios::out){
        std::ofstream csv_file;
        csv_file.open(pPath, pMode);
        if (csv_file.is_open()){
            for (auto data_row: pInfo){
                for (int i=0; i< data_row.size()-1; i++)
                    csv_file << data_row[i] << ",";
                csv_file << data_row.back() << std::endl;
            }
        } else throw std::runtime_error("Could not open the solver info log file " + pPath + " \n");
        csv_file.close();
    }

    inline std::vector<std::vector<std::string>> read_info(const std::string& pPath){
        std::ifstream csv_file;
        csv_file.open(pPath, std::ios::in);
        std::vector<std::vector<std::string>> info{{}};
        std::string line, field;
        if (csv_file.is_open()){
            while (std::getline(csv_file, line)){
                std::stringstream line_stream(line);
                while (std::getline(line_stream, field, ',')) {
                    string_strip(' ', field);
                    info.back().push_back(field);
                }
                info.emplace_back();
            }
            info.pop_back();
        } else throw std::runtime_error("Could not open the csv info file.");
        csv_file.close();
        return info;
    }

    template<typename T>
    void read_struct(const std::string& pPath, std::vector<T>& pStruct, unsigned short pNFields){
        std::ifstream struct_stream; std::ostringstream err; std::string eof;
        struct_stream.open(pPath, std::ios::out);
        if (!struct_stream.is_open()){
            err << "Could not open struct file " << pPath << std::endl;
            throw std::runtime_error(err.str());
        }
        char hash; size_t n, m;
        struct_stream >> hash >> n >> m;
        if (m != pNFields){
            err << "Wrong number of fields in the struct file" << pPath << std::endl;
            throw std::runtime_error(err.str());
        }
        T structure; pStruct.reserve(n);
        for (int i = 0; i < n; ++i) {
            struct_stream >> structure; pStruct.push_back(structure);
        }
        std::getline(struct_stream, eof);
        if (std::getline(struct_stream, eof)) {
            if (eof != "# EOF" && eof != "# EOF\n" && eof !="# EOF \n") {
                err << "Invalid EOF marker " << eof << " " << pPath << std::endl;
                throw std::runtime_error(err.str());
            }
        }
        else throw std::runtime_error("Failed to read EOF marker \n");
        struct_stream.close();
    }

    template<typename T>
    void write_struct(const std::string& pPath, const std::vector<T>& pStruct, unsigned short pNFields){
        std::ofstream struct_stream; std::ostringstream err;
        struct_stream.open(pPath, std::ios::out);
        if (!struct_stream.is_open()){
            err << "Could not open struct file " << pPath << std::endl;
            throw std::runtime_error(err.str());
        }
        struct_stream << "# " << pStruct.size() << " " << pNFields << std::endl;
        for (auto& structure: pStruct) struct_stream << structure << std::endl;
        struct_stream << "# EOF" << std::endl;
        struct_stream.close();
    }

    template<typename T>
    inline bool compare_opt_solution(const std::string& pNewPath,
                                     const std::string& pOptPath,
                                     T& pSolver){
        std::vector<csp::active_node> opt_solution; bool best;
        read_struct<csp::active_node>(pOptPath, opt_solution, 3);
        auto current_cost = pSolver.get_lexicographic_cost();
        auto opt_cost = pSolver.cost(opt_solution);
        best = opt_cost > current_cost;
        if (best) write_struct<csp::active_node>(pNewPath, pSolver.get_solution(), 3);
        return best;
    }

    inline bool update_pareto_front_2d(const std::pair<int, int> pBiSolutionCost, const std::string& pInstanceParetoLog,
                                 const std::string& pCost1Id, const std::string& pCost2Id, const std::string& uuid){
        std::vector<std::vector<std::string>> info = read_info(pInstanceParetoLog);
        std::size_t field1, field2;
        std::vector<std::string> place_holder(info[0].size());
        place_holder[0] = uuid;
        bool solution_dominates = true;
        for (std::size_t i=0; i<info[0].size(); i++){
            if (info[0][i] == pCost1Id) field1 = i;
            else if (info[0][i] == pCost2Id) field2 = i;
            else continue;
        }
        place_holder[field1] = std::to_string(pBiSolutionCost.first);
        place_holder[field2] = std::to_string(pBiSolutionCost.second);
        for (int i =1; i< info.size(); i++){
            if (std::stoi(info[i][field1]) > pBiSolutionCost.first && std::stoi(info[i][field2]) > pBiSolutionCost.second)
                info.erase(info.begin() + i);
            else if (std::stoi(info[i][field1]) <= pBiSolutionCost.first && std::stoi(info[i][field2]) <= pBiSolutionCost.second)
                solution_dominates = false;
        }
        //Partial dominance over all elements or complete dominance over one.
        if (solution_dominates) {
            info.push_back(place_holder);
            write_info(pInstanceParetoLog, info);
        }
        return solution_dominates;
    }

    template <typename T>
    std::vector<std::vector<T>> pad_matrix(const std::vector<std::vector<T>>& pMatrix, const int pPadSize, const T pPadValue){
        std::vector<std::vector<T>> padded_matrix(2*pPadSize + pMatrix.size(), std::vector<T>(2*pPadSize + pMatrix[0].size(), pPadValue));
        for (std::size_t i = pPadSize; i< pPadSize + pMatrix.size(); i++){
            for (std::size_t j= pPadSize; j<pPadSize + pMatrix[0].size(); j++)
                padded_matrix[i][j] = pMatrix[i-pPadSize][j-pPadSize];
        }
        return padded_matrix;
    }

    template<typename T>
    void unpad_matrix(std::vector<std::vector<T>>& pMatrix, const int pPadSize){
        std::vector<std::vector<T>> unpadded_matrix(pMatrix.size()-2*pPadSize, std::vector<T>(pMatrix[0].size()-2*pPadSize));
        for (size_t i = 0; i< pMatrix.size()-2*pPadSize; i++){
            for (size_t j= 0; j<pMatrix[0].size()-2*pPadSize; j++)
                unpadded_matrix[i][j] = pMatrix[i+pPadSize][j+pPadSize];
        }
        pMatrix = unpadded_matrix;
    }
}
