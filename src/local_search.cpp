#include "local_search.h"
const int MY_RAND_MAX_INT = 10000000;
const int prob1 =  0.0*100000;  //random add 0.1%
const int prob2 =  0.0*100000;  //random swap 0.1%
const int prob3 =  50*100000;  //random drop 50%
unsigned long long total_step = 0;
unsigned long long * vertex_step = NULL;
long* add_candidate_index = NULL;
long* drop_candidate_index = NULL;
long start_vertex = 0;
long vertex_count = 0;
long edge_count = 0;
int init_num = 0;
long continuous_noimprove = 0;
long continuous_drop = 0;
int last_improved_no = 0;

long long add_num = 0;
long long swap_num = 0;
long long drop_num = 0;

// Build
vector<vector<weighted_edge>> vertex_adjacency_list;
long* vertex_edge_weight = NULL;
long* vertex_degree = NULL;

// Construction
bool is_new_graph = false;
vector<long> candidates;
long* cand_neighbor_weight;
bool* is_in_candidates;
is_add_nei_element* is_addv_neighbor;
int ttt = 2;

vector<long> start_vertices;
long untest_pointer = 0;
long head_p, tail_p;

//Select vertex
best_element* best_add_array = NULL;
best_element* best_swap_array = NULL;
best_element* best_drop_array = NULL;

// Solutions
vector<long> solution;
vector<long> best_solution;
long solution_weight;
long Waim = INT_MAX;
long best_solution_weight=0;
double best_solution_time;
long best_solution_try;

// Main loop
bool better_since_cons = true;
long simp_count=0;

// BMS
long start_bms_count = 1;
long min_bms_count = 8;
long max_bms_count = 128;
long real_bms_count;
long add_bms_count = 128;
long swap_bms_count = 128;

// Definition of global variables
long tries;
long max_tries = 2000000000;
int cutoff_time = 100;
int t=2;

// add 
vector<long> add_candidate_set;
long* add_vertex_score = NULL;
// drop
long* drop_vertex_score = NULL;
// swap
vector<swap_element> swap_candidate_set;

bool* is_in_solution = NULL;
vertex_connect_clique* vertex_num_connect_to_clique = NULL; 
bool* tmp_is_v1 = NULL;
is_neighbor_element* is_in_neighbor = NULL;
bool* tmp_v2_in_swap = NULL;

bool* conf_change = NULL;

Remaining_vertex remaining_vertex;

void print_solution() {
    for (int i=0; i<solution.size(); i++) {
        cout << solution[i] << " ";
    }
    cout<<endl;
}

bool comp(const swap_element &a, const swap_element &b) {
    return a.v2 > b.v2;
}

void update_best_solution() {
    best_solution.clear();
    for (auto v : solution) {
        best_solution.push_back(v);
    }
}

long find_v1_from_solution(long v2) {
    for (long i=0; i<solution.size(); i++) {
        tmp_is_v1[solution[i]] = 1;
    }
    for (auto u : vertex_adjacency_list[v2]) {
        if (tmp_is_v1[u.vertex]) {
            tmp_is_v1[u.vertex] = 0;
        }
    }
    for (long i=0; i<solution.size(); i++) {
        if(tmp_is_v1[solution[i]]) {
            tmp_is_v1[solution[i]] = 0;
            return solution[i];
        }
    }
    return 0;
}

void add(long add_v, vector<long>::size_type index, bool change_conf, char ccc) {
    //is_addv_neighbor.clear();
    vertex_step[add_v] = total_step;
    swap_element se;
    solution.push_back(add_v);
    drop_candidate_index[add_v] = solution.size()-1;
    conf_change[add_v] = 0;
    solution_weight += add_vertex_score[add_v];
    if (solution_weight > best_solution_weight) {
        better_since_cons = true;
        last_improved_no = init_num;
        best_solution_weight = solution_weight;
        cout << "o " << best_solution_weight << endl;//<< " " << solution.size() << " " << init_num << " " << ccc << " " << start_bms_count<< " " << real_bms_count << endl;
    
        update_best_solution();
        if (best_solution_weight == Waim)
        {
            cout << "finished!" << endl;
            exit(0);
        }
        continuous_noimprove = 0;
    }
    is_in_solution[add_v] = 1;

    add_candidate_set[index]=*(add_candidate_set.end() - 1);
    add_candidate_index[*(add_candidate_set.end()-1)] = add_candidate_index[add_v];
    add_candidate_set.pop_back();
    add_vertex_score[add_v] = 0;

    is_in_candidates[add_v] = false;
    if (change_conf) {
        for (auto u : vertex_adjacency_list[add_v]) {
            conf_change[u.vertex] = 1;
            vertex_num_connect_to_clique[u.vertex].number ++;
            vertex_num_connect_to_clique[u.vertex].weight += u.weight;
            is_in_neighbor[u.vertex].is_in_neighbor = 1;
            is_in_neighbor[u.vertex].weight = u.weight;
            if (is_in_candidates[u.vertex]) {
                cand_neighbor_weight[u.vertex] -= u.weight;
                is_addv_neighbor[u.vertex].is_nei = true;
                is_addv_neighbor[u.vertex].weight = u.weight;
            }
            if (1 == is_in_solution[u.vertex]) {
                drop_vertex_score[u.vertex] -= u.weight;
            }
        }
    }
    else { 
        for (auto u : vertex_adjacency_list[add_v]) {
            vertex_num_connect_to_clique[u.vertex].number ++;
            vertex_num_connect_to_clique[u.vertex].weight += u.weight;
            is_in_neighbor[u.vertex].is_in_neighbor = 1;
            is_in_neighbor[u.vertex].weight = u.weight;
            if (is_in_candidates[u.vertex]) {
                cand_neighbor_weight[u.vertex] -= u.weight;
                is_addv_neighbor[u.vertex].is_nei = true;
                is_addv_neighbor[u.vertex].weight = u.weight;
            }
            if (1 == is_in_solution[u.vertex]) {
                drop_vertex_score[u.vertex] -= u.weight;
            }
        }
    }
    drop_vertex_score[add_v] = 0 - vertex_num_connect_to_clique[add_v].weight;

    // remove swap_candidate_set
    for (vector<swap_element>::size_type i = 0; i < swap_candidate_set.size();) {
        long v2 = swap_candidate_set[i].v2;
        long v1;
        if (0 == is_in_neighbor[v2].is_in_neighbor) {
            swap_candidate_set[i] = *(swap_candidate_set.end() - 1);
            swap_candidate_set.pop_back();
            tmp_v2_in_swap[v2] = false;
        }
        else {
            v1 = swap_candidate_set[i].v1;
            swap_candidate_set[i].score += (is_in_neighbor[v2].weight - is_in_neighbor[v1].weight);
            i++;
        }
    }

    if (2 == solution.size()) {
        for (auto u : vertex_adjacency_list[add_v]) {
            // add swap_candidate_set; if the number of the new clique is 2, then 
            if (!is_in_candidates[u.vertex] && !is_in_solution[u.vertex]) {
                se.v1 = solution[0];
                se.v2 = u.vertex;
                se.score = u.weight - solution_weight;
                swap_candidate_set.push_back(se);
                tmp_v2_in_swap[se.v2] = true;
            }
        }
    }
    // remove add_candidate_set; add swap_candidate_set
    for (vector<long>::size_type i = 0; i < add_candidate_set.size(); ) {
        long cur_v = add_candidate_set[i];

        if (!is_addv_neighbor[cur_v].is_nei){
            se.v1 = add_v;
            se.v2 = cur_v;
            //update
            
            is_in_candidates[cur_v] = false; 
            add_candidate_set[i] = *(add_candidate_set.end() - 1);
            add_candidate_index[add_candidate_set[i]] = add_candidate_index[cur_v];
            add_candidate_set.pop_back();
            if (ccc == 'c') {
                cand_neighbor_weight[cur_v] = 0;
                for (auto v : vertex_adjacency_list[cur_v]) {
                    if (is_in_candidates[v.vertex]) {
                        cand_neighbor_weight[v.vertex] -= v.weight;
                    }
                }
            }
            add_vertex_score[cur_v] = 0;
            swap_candidate_set.push_back(se);
            tmp_v2_in_swap[se.v2] = true;
            swap_candidate_set[swap_candidate_set.size()-1].score = vertex_num_connect_to_clique[cur_v].weight - vertex_num_connect_to_clique[add_v].weight;
        }
        else {
            // modify add score
            add_vertex_score[cur_v] += is_addv_neighbor[cur_v].weight; 
            is_addv_neighbor[cur_v].is_nei = false; // erase the value here, making this array clean after the step.
            i++;
        }
    } 

    for (auto u : vertex_adjacency_list[add_v]) {
        is_in_neighbor[u.vertex].is_in_neighbor = 0;
    }
}

void drop(long del_v, vector<long>::size_type index) {
    //is_addv_neighbor.clear();
    vertex_step[del_v] = total_step;
    solution[index] = *(solution.end() - 1);
    solution.pop_back();
    drop_candidate_index[solution[index]] = drop_candidate_index[del_v];
    conf_change[del_v] = 0;
    solution_weight += drop_vertex_score[del_v];
    add_candidate_set.push_back(del_v);
    add_candidate_index[del_v] = add_candidate_set.size()-1;
    is_in_candidates[del_v] = true;
    is_in_solution[del_v] = 0;
    drop_vertex_score[del_v] = 0;

    //cand_neighbor_weight[del_v] = 0;
    long solution_size = solution.size();
    for (auto u : vertex_adjacency_list[del_v]) {
        //conf_change[u.vertex] = 1;
        vertex_num_connect_to_clique[u.vertex].number --;
        vertex_num_connect_to_clique[u.vertex].weight -= u.weight;
        if (is_in_candidates[u.vertex]) {
            add_vertex_score[u.vertex] -= u.weight;
            //cand_neighbor_weight[u.vertex] += u.weight;
            //cand_neighbor_weight[del_v] += u.weight;
        }
        //is_addv_neighbor[u.vertex].is_nei = 1;
        //is_addv_neighbor[u.vertex].weight = u.weight;
        if (is_in_solution[u.vertex]) {
            drop_vertex_score[u.vertex] += u.weight;
        }
    }
    add_vertex_score[del_v] = vertex_num_connect_to_clique[del_v].weight;
    // add add_candidate_set
    for (auto v : vertex_adjacency_list[solution[0]]) {
        if (!is_in_solution[v.vertex] && !is_in_candidates[v.vertex] && solution_size == vertex_num_connect_to_clique[v.vertex].number) {
            add_candidate_set.push_back(v.vertex);
            add_candidate_index[v.vertex] = add_candidate_set.size()-1;
            add_vertex_score[v.vertex] = vertex_num_connect_to_clique[v.vertex].weight; 
            //cand_burden[v.vertex] = add_vertex_score[v.vertex];
            is_in_candidates[v.vertex] = true;
#ifndef greedy_add_construct
            cand_neighbor_weight[v.vertex] = 0;
            for (auto u : vertex_adjacency_list[v.vertex]) {
                if (is_in_candidates[u.vertex]) {
                    cand_neighbor_weight[u.vertex] += u.weight;
                    cand_neighbor_weight[v.vertex] += u.weight;
                }
            }
#endif
        }
    }

    // del swap_candidate_set
    for (vector<swap_element>::size_type i = 0; i < swap_candidate_set.size();) {
        long v2 = swap_candidate_set[i].v2;
        //if (!is_addv_neighbor[v2].is_nei) {
        if (vertex_num_connect_to_clique[v2].number == solution_size) {
            swap_candidate_set[i] = *(swap_candidate_set.end() - 1);
            swap_candidate_set.pop_back();
            tmp_v2_in_swap[v2] = false;
        }
        else {
            // modify swap_score
            long v1 = swap_candidate_set[i].v1;
            long v2 = swap_candidate_set[i].v2;
            swap_candidate_set[i].score = drop_vertex_score[v1] + vertex_num_connect_to_clique[v2].weight;
            i++;
        }
    }
    // add swap_candidate_set
    if (1 != solution_size) {
        long v, v1;
        for (long i=0; i<solution_size; i++) {
            v = solution[i];
            for (auto u : vertex_adjacency_list[v]) {
                long v2 = u.vertex;
                if (solution_size-1 == vertex_num_connect_to_clique[v2].number && !is_in_solution[v2] && !tmp_v2_in_swap[v2]) {
                    v1 = find_v1_from_solution(v2);
                    if (v1) {
                        swap_element se;
                        se.v1 = v1;
                        se.v2 = v2;
                        se.score = vertex_num_connect_to_clique[v2].weight + drop_vertex_score[v1];
                        swap_candidate_set.push_back(se);
                        tmp_v2_in_swap[v2] = true;
                    }
                }
            }
        }
        //sort(swap_candidate_set.begin(), swap_candidate_set.end(), comp);
    }
}

void init_add() {
    long tmp_vertex;
    memset(cand_neighbor_weight, 0, sizeof(long) * (vertex_count+1));
    for (long i=0; i<solution.size(); i++) {
        tmp_vertex = solution[i];
        drop_vertex_score[tmp_vertex] = 0;
        is_in_solution[tmp_vertex] = 0;
        vertex_num_connect_to_clique[tmp_vertex].number = 0;
        vertex_num_connect_to_clique[tmp_vertex].weight = 0;
        for (auto u : vertex_adjacency_list[tmp_vertex]) {
            vertex_num_connect_to_clique[u.vertex].number = 0;
            vertex_num_connect_to_clique[u.vertex].weight = 0;
        }
    }
    for (long i=0; i<add_candidate_set.size(); i++) {
        tmp_vertex = add_candidate_set[i];
        is_in_candidates[tmp_vertex] = false;
        add_vertex_score[tmp_vertex] = 0;
    }
    for (long i=0; i<swap_candidate_set.size(); i++) {
        tmp_v2_in_swap[swap_candidate_set[i].v2] = false; 
    }
    solution.clear();
    add_candidate_set.clear();
    swap_candidate_set.clear();

    solution_weight = 0;
    vector<long>::size_type index, tmp_index;

#ifdef is_simplify
    if (is_new_graph)
    {
        start_vertices = remaining_vertex.vertex; // Vertices from which we choose the start vertex
        untest_pointer = 0;
    }
#endif
    // If all the vertices have been tested, then we should adjust the value of bms
    if(untest_pointer == start_vertices.size())
    {
        untest_pointer = 0;
    }
    
    long best_score = 0, best_start_num = 0, tmp_score, aaa;

    // choose start vertex; Choose randomly k vertices to check
    if (start_vertices.size()-untest_pointer > start_bms_count) {
        for (long i=1; i<=start_bms_count; ++i) {
            tmp_index = untest_pointer+rand()%(start_vertices.size()-untest_pointer);
            if (vertex_upper_bound[start_vertices[tmp_index]]) {
                tmp_score = vertex_upper_bound[start_vertices[tmp_index]];
            }
            else {
                tmp_score = soft_upper_bound(start_vertices[tmp_index]);
                vertex_upper_bound[start_vertices[tmp_index]] = tmp_score;
            }
            if (tmp_score > best_score){
                best_score = tmp_score;
                index = tmp_index;
            }
        }
    }
    else {
        for (tmp_index=untest_pointer; tmp_index<start_vertices.size(); tmp_index++) {
            if (vertex_upper_bound[start_vertices[tmp_index]]) {
                tmp_score = vertex_upper_bound[start_vertices[tmp_index]];
            }
            else {
                tmp_score = soft_upper_bound(start_vertices[tmp_index]);
                vertex_upper_bound[start_vertices[tmp_index]] = tmp_score;
            }
            if (tmp_score > best_score){
                best_score = tmp_score;
                index = tmp_index;
            }
        }
    }
    
    //}
    start_vertex = start_vertices[index];
    
    // Move the chosen start vertex to the untested pointer
    std::swap(start_vertices[index], start_vertices[untest_pointer++]);
    solution.push_back(start_vertex);
    drop_candidate_index[start_vertex] = solution.size()-1;
    is_in_solution[start_vertex] = 1;
    solution_weight = 0;
    total_step ++;
    vertex_step[start_vertex] = total_step;
    
    // Initialize the set of candidate vertices S = N(startv)
    for (auto u : vertex_adjacency_list[start_vertex]) {
        add_candidate_set.push_back(u.vertex);  // start_vertex's neighbor_node u.vertex
        add_candidate_index[u.vertex] = add_candidate_set.size()-1;
        add_vertex_score[u.vertex] = u.weight;
        is_in_candidates[u.vertex] = true;     // put start_vertex's neighbor_node u.vertex into candset
        vertex_num_connect_to_clique[u.vertex].number = 1;
        vertex_num_connect_to_clique[u.vertex].weight = u.weight;
        cand_neighbor_weight[u.vertex] = 0;
        //is_addv_neighbor[u.vertex].is_nei = false;
    }

    for (auto v : add_candidate_set) {
        for (auto n : vertex_adjacency_list[v]){
            if (1 == is_in_candidates[n.vertex]) //
                cand_neighbor_weight[v] += n.weight;
        }
    }
    return;
}

void initGreedyConstruct() {
    init_add();
    long best_score = 0, vertex, tmp_score, modify_vertex, index;
    long best_add_num = 0, tmp_index;

    while (!add_candidate_set.empty()) {
        best_score = 0;
        if (add_candidate_set.size() <= real_bms_count) {
            for (long i=0; i<add_candidate_set.size(); i++) {
                vertex = add_candidate_set[i];
                tmp_score = add_vertex_score[vertex] + cand_neighbor_weight[vertex] / ttt;
                //tmp_score = add_vertex_score[vertex];
                if (tmp_score > best_score) {
                    modify_vertex = vertex;
                    index = i;
                    best_score = tmp_score;
                }
            }
        }
        else {
            for (long i=0; i<real_bms_count; i++) {
                tmp_index = rand()%add_candidate_set.size();
                vertex = add_candidate_set[tmp_index];
                tmp_score = add_vertex_score[vertex] + cand_neighbor_weight[vertex] / ttt;
                //tmp_score = add_vertex_score[vertex];
                if (tmp_score > best_score) {
                    modify_vertex = vertex;
                    index = tmp_index;
                    best_score = tmp_score;
                }
            }
        }
        add(modify_vertex, index, true, 'c');
    }
    return;
}

void local_search() {
    continuous_noimprove = 0;
    continuous_drop = 0;

#ifdef greedy_add_construct
    initGreedyConstruct();
#endif
#ifdef no_greedy_add_construct
    init_add();
#endif
    //if (init_num - last_improved_no >= 10) {
    real_bms_count = real_bms_count + real_bms_count;
    if(real_bms_count > max_bms_count) real_bms_count = min_bms_count;
    if (1 == start_bms_count) {
        start_bms_count = 100;
    }
    else {
        start_bms_count = 1;
    }
    init_num++;
    long best_add_score, tmp_score, best_drop_score, best_swap_score;
    long add_vertex, add_index, vertex, drop_vertex, drop_index, v1, v2;
    long best_add_num, best_swap_num, best_drop_num;
    long tmp_index;
    long step = 0;
    long L = 100;
    unsigned long long cur_start_step = total_step;
    long long circle_step = 0;
    for (step = 0; step < L; step++) {
        total_step ++;
        continuous_noimprove ++;
        if ((double) (clock()) / CLOCKS_PER_SEC > cutoff_time) {
            total_step --;
            return;
        }
        best_add_num = 0; best_swap_num = 0; best_drop_num = 0;
        best_add_score = 0;
        
        if (!add_candidate_set.empty()) {

            long add_candidate_size = add_candidate_set.size();
            if (add_candidate_size < add_bms_count) {
                for (long i=0; i<add_candidate_size; i++) {
                    vertex = add_candidate_set[i];
                    if (conf_change[vertex]) {
#ifndef greedy_add_construct 
                        tmp_score = add_vertex_score[vertex] + cand_neighbor_weight[vertex]/ttt;
#endif
#ifdef greedy_add_construct
                        tmp_score = add_vertex_score[vertex];
#endif
                        if (tmp_score > best_add_score) {
                            best_add_num = 0;
                            best_add_array[best_add_num].vertex_num = vertex;
                            best_add_array[best_add_num].index_num = i;
                            best_add_num ++;
                            best_add_score = tmp_score;
                        }
                        else if (tmp_score == best_add_score) {
#ifdef select_with_step        
                            if (vertex_step[vertex] == vertex_step[best_add_array[0].vertex_num]) {
#endif                            
                                best_add_array[best_add_num].vertex_num = vertex;
                                best_add_array[best_add_num].index_num = i;
                                best_add_num ++;
#ifdef select_with_step                            
                            }
                            else if (vertex_step[vertex] < vertex_step[best_add_array[0].vertex_num]){
                                best_add_num = 0;
                                best_add_array[best_add_num].vertex_num = vertex;
                                best_add_array[best_add_num].index_num = i;
                                best_add_num ++;
                            }
#endif
                        }
                    }
                }
            }
            else {
                long tmp_index = 0;
                for (long i=0; i<add_bms_count; i++) {
                    //for (long i=0; i<start_bms_count; i++) {
                    tmp_index = rand()%add_candidate_set.size();
                    vertex = add_candidate_set[tmp_index];
                    if (conf_change[vertex]) {
#ifndef greedy_add_construct 
                        tmp_score = add_vertex_score[vertex] + cand_neighbor_weight[vertex]/ttt;
#endif
#ifdef greedy_add_construct
                        tmp_score = add_vertex_score[vertex];
#endif
                        if (tmp_score > best_add_score) {
                            best_add_num = 0;
                            best_add_array[best_add_num].vertex_num = vertex;
                            best_add_array[best_add_num].index_num = tmp_index;
                            best_add_num ++;
                            best_add_score = tmp_score;
                        }
                        else if (tmp_score == best_add_score) {
#ifdef select_with_step        
                            if (vertex_step[vertex] == vertex_step[best_add_array[0].vertex_num]) {
#endif                            
                                best_add_array[best_add_num].vertex_num = vertex;
                                best_add_array[best_add_num].index_num = tmp_index;
                                best_add_num ++;
#ifdef select_with_step                            
                            }
                            else if (vertex_step[vertex] < vertex_step[best_add_array[0].vertex_num]){
                                best_add_num = 0;
                                best_add_array[best_add_num].vertex_num = vertex;
                                best_add_array[best_add_num].index_num = tmp_index;
                                best_add_num ++;
                            }
#endif
                        }
                    }
                }
                }
                if (best_add_num > 0) {
                    tmp_index = rand()%best_add_num;
                    add_vertex = best_add_array[tmp_index].vertex_num;
                    add_index = best_add_array[tmp_index].index_num;
#ifndef greedy_add_construct 
                    best_add_score = add_vertex_score[add_vertex]; 
#endif
                }
            }
        best_swap_score = -2147483647;
        if (!swap_candidate_set.empty()) {
            
            long swap_candidate_size = swap_candidate_set.size();
            if (swap_candidate_size <= swap_bms_count) {
                for (long i=0; i<swap_candidate_size; i++) {
                    if (conf_change[swap_candidate_set[i].v2]) {
                        tmp_score = swap_candidate_set[i].score;
                        if (tmp_score > best_swap_score) {
                            best_swap_num = 0;
                            best_swap_array[best_swap_num].vertex_num = swap_candidate_set[i].v1;
                            best_swap_array[best_swap_num].index_num = swap_candidate_set[i].v2;
                            best_swap_num ++;
                            best_swap_score = tmp_score;
                        }
                        else if (tmp_score == best_swap_score) {
#ifdef select_with_step                        
                            if (vertex_step[swap_candidate_set[i].v2] == vertex_step[best_swap_array[0].index_num]) {
#endif
                                best_swap_array[best_swap_num].vertex_num = swap_candidate_set[i].v1;
                                best_swap_array[best_swap_num].index_num = swap_candidate_set[i].v2;
                                best_swap_num ++;
#ifdef select_with_step                
                            }
                            else if (vertex_step[swap_candidate_set[i].v2] < vertex_step[best_swap_array[0].index_num]) {
                                best_swap_num = 0;
                                best_swap_array[best_swap_num].vertex_num = swap_candidate_set[i].v1;
                                best_swap_array[best_swap_num].index_num = swap_candidate_set[i].v2;
                                best_swap_num ++;
                            }
#endif
                        }
                    }
                }
                if (best_swap_num > 0) {
                    tmp_index = rand()%best_swap_num;
                    v1 = best_swap_array[tmp_index].vertex_num;
                    v2 = best_swap_array[tmp_index].index_num;
                }
            }
            else {
                for (long i=0; i<swap_bms_count; i++) {
                    tmp_index = rand()%swap_candidate_size;
                    if (conf_change[swap_candidate_set[tmp_index].v2]) {
                        tmp_score = swap_candidate_set[tmp_index].score;
                        if (tmp_score > best_swap_score) {
                            best_swap_num = 0;
                            best_swap_array[best_swap_num].vertex_num = swap_candidate_set[tmp_index].v1;
                            best_swap_array[best_swap_num].index_num = swap_candidate_set[tmp_index].v2;
                            best_swap_num ++;
                            best_swap_score = tmp_score;
                        }
                        else if (tmp_score == best_swap_score) {
#ifdef select_with_step                        
                            if (vertex_step[swap_candidate_set[tmp_index].v2] == vertex_step[best_swap_array[0].index_num]) {
#endif
                                best_swap_array[best_swap_num].vertex_num = swap_candidate_set[tmp_index].v1;
                                best_swap_array[best_swap_num].index_num = swap_candidate_set[tmp_index].v2;
                                best_swap_num ++;
#ifdef select_with_step                
                            }
                            else if (vertex_step[swap_candidate_set[tmp_index].v2] < vertex_step[best_swap_array[0].index_num]) {
                                best_swap_num = 0;
                                best_swap_array[best_swap_num].vertex_num = swap_candidate_set[tmp_index].v1;
                                best_swap_array[best_swap_num].index_num = swap_candidate_set[tmp_index].v2;
                                best_swap_num ++;
                            }
#endif
                        }
                    }
                }
                if (best_swap_num > 0) {
                    tmp_index = rand()%best_swap_num;
                    v1 = best_swap_array[tmp_index].vertex_num;
                    v2 = best_swap_array[tmp_index].index_num;
                }
            }
        }
        if (best_add_num > 0) {
            continuous_drop = 0;
            if (best_swap_num > 0) {
                if (best_add_score >= best_swap_score) { //add
                    if (vertex_step[add_vertex] > cur_start_step)
                        circle_step ++;
                    add(add_vertex, add_index, true, 'a');
                    add_num++;
                }
                else {  //swap
                    if (vertex_step[v2] > cur_start_step)
                        circle_step ++;
                    
                    drop(v1, drop_candidate_index[v1]);
                    add(v2, add_candidate_index[v2], false, 's');
                    swap_num++;
                }
            }
            else { // add();
                if (vertex_step[add_vertex] > cur_start_step)
                    circle_step ++;
                add(add_vertex, add_index, true, 'a');
                add_num++;
            }
        }
        else if (best_swap_num > 0) {

            if (best_swap_score >= 0) {
                if (vertex_step[v2] > cur_start_step)
                    circle_step ++;
                
                drop(v1, drop_candidate_index[v1]);
                add(v2, add_candidate_index[v2], false, 's');
                
                swap_num++;
                continuous_drop = 0;
            }
            else { //best_swap_score < 0
                if (rand()%MY_RAND_MAX_INT < prob2) {
                    tmp_index = rand()%swap_candidate_set.size();
                    v1 = swap_candidate_set[tmp_index].v1;
                    v2 = swap_candidate_set[tmp_index].v2;
                    best_swap_num = 1;
                    best_swap_score = swap_candidate_set[tmp_index].score;
                }
                best_drop_score = -2147483647;
                if (!solution.empty()) {
                    if (rand()%MY_RAND_MAX_INT < prob3) {
                        drop_index = rand()%solution.size();
                        drop_vertex = solution[drop_index];
                        best_drop_score = drop_vertex_score[drop_vertex];
                    }
                    else {
                        for (long i=0; i<solution.size(); i++) {
                            vertex = solution[i];
                            tmp_score = drop_vertex_score[vertex];
                            if (tmp_score > best_drop_score) {
                                best_drop_num = 0;
                                best_drop_array[best_drop_num].vertex_num = vertex;
                                best_drop_array[best_drop_num].index_num = i;
                                best_drop_num ++;
                                best_drop_score = tmp_score;
                            }
                            else if (tmp_score == best_drop_score) {
#ifdef select_with_step
                                if (vertex_step[vertex] == vertex_step[best_drop_array[0].vertex_num]) {
#endif
                                    best_drop_array[best_drop_num].vertex_num = vertex;
                                    best_drop_array[best_drop_num].index_num = i;
                                    best_drop_num ++;
#ifdef select_with_step
                                }
                                else if (vertex_step[vertex] < vertex_step[best_drop_array[0].vertex_num]) {
                                    best_drop_num = 0;
                                    best_drop_array[best_drop_num].vertex_num = vertex;
                                    best_drop_array[best_drop_num].index_num = i;
                                    best_drop_num ++;
                                    best_drop_score = tmp_score;
                                }
#endif                        
                            }
                        }
                        tmp_index = rand()%best_drop_num;
                        drop_vertex = best_drop_array[tmp_index].vertex_num;
                        drop_index = best_drop_array[tmp_index].index_num;
                    }
                }
                if (best_drop_score <= best_swap_score) {
                    if (vertex_step[v2] > cur_start_step)
                        circle_step ++;
                    
                    drop(v1, drop_candidate_index[v1]);
                    add(v2, add_candidate_index[v2], false, 's');
                    
                    swap_num++;
                    continuous_drop = 0;
                }
                else {
                    continuous_drop ++;
                    
                    drop(drop_vertex, drop_index);
                    drop_num++;
                }
            }
        }
        else if (!solution.empty()){
            continuous_drop ++;

            if (rand()%MY_RAND_MAX_INT < prob3) {
                drop_index = rand()%solution.size();
                drop_vertex = solution[drop_index];
                best_drop_score = drop_vertex_score[drop_vertex];
            }
            else {
                best_drop_score = -2147483647;
                for (long i=0; i<solution.size(); i++) {
                    vertex = solution[i];
                    tmp_score = drop_vertex_score[vertex];
                    if (tmp_score > best_drop_score) {
                        best_drop_num = 0;
                        best_drop_array[best_drop_num].vertex_num = vertex;
                        best_drop_array[best_drop_num].index_num = i;
                        best_drop_num ++;
                        best_drop_score = tmp_score;
                    }
                    else if (tmp_score == best_drop_score) {
#ifdef select_with_step
                        if (vertex_step[vertex] == vertex_step[best_drop_array[0].vertex_num]) {
#endif
                            best_drop_array[best_drop_num].vertex_num = vertex;
                            best_drop_array[best_drop_num].index_num = i;
                            best_drop_num ++;
#ifdef select_with_step
                        }
                        else if (vertex_step[vertex] < vertex_step[best_drop_array[0].vertex_num]) {
                            best_drop_num = 0;
                            best_drop_array[best_drop_num].vertex_num = vertex;
                            best_drop_array[best_drop_num].index_num = i;
                            best_drop_num ++;
                            best_drop_score = tmp_score;
                        }
#endif                        
                    }
                }
                tmp_index = rand()%best_drop_num;
                drop_vertex = best_drop_array[tmp_index].vertex_num;
                drop_index = best_drop_array[tmp_index].index_num;
            }
            drop(drop_vertex, drop_index);
        }
        else {
            total_step --;
            return;
        }
    }
    return;
}
