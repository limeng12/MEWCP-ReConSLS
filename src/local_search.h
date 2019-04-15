#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H
#include <limits.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <math.h>
#include <unistd.h>
#include <string.h>
using namespace std;

//#define format_konect
#define format_dimacs

//#define construct_select_with_step
#define select_with_step

#define init_with_degree
//#define init_with_step

#define greedy_add_construct
//#define no_greedy_add_construct

#define is_simplify

extern unsigned long long * vertex_step;
extern long* add_candidate_index;
extern long* drop_candidate_index;

//clock_t start
extern long vertex_count;
extern long edge_count;

extern int init_num;
extern int last_improved_no;
extern long long add_num;
extern long long swap_num;
extern long long drop_num;
//int degree_weight = 0;

struct weighted_edge {
    long vertex;
    long weight;
};
struct is_add_nei_element {
    long weight;
    bool is_nei;
};

struct best_element {
    long vertex_num;
    long index_num;
};
// The struct to maintain the remaining vertices
struct Remaining_vertex {
    vector<long> vertex;
    vector<vector<long>::size_type> index;

    vector<long>::iterator begin() {
        return vertex.begin();
    }
    vector<long>::iterator end() {
        return vertex.end();
    }
    void init(vector<long>::size_type vertex_size) {
        vertex.reserve(vertex_size);
        index.resize(vertex_size + 1);
        for (vector<long>::size_type i = 1; i <= vertex_size; ++i) {
            vertex.push_back((long)i);
            index[i] = i - 1;
        }
    }
    void remove(long v) {
        index[*vertex.rbegin()] = index[v];
        vertex[index[v]] = *vertex.rbegin();
        vertex.pop_back();
    }

    vector<long>::size_type size() {
        return vertex.size();
    }

    bool empty() {
        return vertex.empty();
    }
};
extern Remaining_vertex remaining_vertex;
// Build
extern vector<vector<weighted_edge>> vertex_adjacency_list;
extern long* vertex_edge_weight;
extern long* vertex_degree;

// Construction
extern bool is_new_graph;
extern vector<long> candidates;
extern long* cand_neighbor_weight;
extern bool* is_in_candidates;
extern is_add_nei_element* is_addv_neighbor;
extern int ttt;

extern vector<long> start_vertices;
extern long untest_pointer;
extern long head_p, tail_p;

// Reduction
extern long* vertex_upper_bound;

#define min_num(x, y) (x <= y ? x : y)
//Select vertex
extern best_element* best_add_array;
extern best_element* best_swap_array;
extern best_element* best_drop_array;

// Solutions
extern vector<long> solution;
extern vector<long> best_solution;
extern long solution_weight;
extern long Waim;
extern long best_solution_weight;
extern double best_solution_time;
extern long best_solution_try;

// Main loop
extern bool better_since_cons;
extern long simp_count;

// BMS
extern long start_bms_count;
extern long min_bms_count;
extern long max_bms_count;
extern long real_bms_count;
extern long add_bms_count;
extern long swap_bms_count;

// Definition of global variables
extern long tries;
extern long max_tries;
extern int cutoff_time;
extern int t;


//local search
struct swap_element {
    long v1;
    long v2;
    long score;
};
struct is_neighbor_element
{
    long weight;
    bool is_in_neighbor;
};
struct vertex_connect_clique
{
    int number;
    long weight;
    /* data */
};
struct swap_number {
    long v1;
    long v2;
};
// add 
extern vector<long> add_candidate_set;
extern long* add_vertex_score;
// drop
extern long* drop_vertex_score;
// swap
extern vector<swap_element> swap_candidate_set;

extern bool* is_in_solution;
extern vertex_connect_clique* vertex_num_connect_to_clique; 
extern bool* tmp_is_v1;
extern is_neighbor_element* is_in_neighbor;
extern bool* tmp_v2_in_swap;

extern bool* conf_change;

long soft_upper_bound(long);

void print_solution();

bool comp(const swap_element &a, const swap_element &b);// 

void update_best_solution();

long find_v1_from_solution(long v2);

// index: add_v's index in add_candidate_set
void add(long add_v, vector<long>::size_type index, bool change_conf, char ccc);

//index: del_v's index in solution
void drop(long del_v, vector<long>::size_type index);

void init_add();

void initGreedyConstruct();

void local_search();
#endif
