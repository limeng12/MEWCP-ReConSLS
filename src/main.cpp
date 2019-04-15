/******************
 * Yi Chu: chuyi@ict.ac.cn
 ******************/
#include "local_search.h"
using namespace std;

string input_filename;
int seed;

// Reduction
vector<long> working_vertex;
vector<long> next_working_vertex;
vector<bool> is_pending;
bool* is_removed;
long simplify_threshold = 100;
long size_threshold = 100000;
long no_improved_num = 0;
int* vertex_max_weight;
int* vertex_min_weight;
long* vertex_upper_bound = NULL;

struct input_edge {
    int v1;
    int v2;
    int weight;
};

vector<input_edge> input_edge_cache;
// Generate the edge weight according to some rules
long generate_edge_weight(long v1, long v2) {
    return (v1 + v2) % 200 + 1;
}

#ifdef format_dimacs
// Read the input file and transform it into a structed vector
bool read_input_file(string file_name) {

    cout << file_name << ":" << endl;

    ifstream in_file(file_name);
    if (!in_file.is_open()) {
        cout << "in_file error" << endl;
        return false;
    }

    // Get vertex_count and edge_count
    string line;
    istringstream is;
    string p, tmp;
    do {
        getline(in_file, line);
        is.clear();
        is.str(line);
        is >> p >> tmp >> vertex_count >> edge_count;
    } while (p != "p");

    /*int tmp_a = round(edge_count/vertex_count);
    if (tmp_a > 2)
        degree_weight = 0;
    else if (tmp_a<=2 && tmp_a>1)
        degree_weight = 50;
    else
        degree_weight = 100;
    */    
    input_edge_cache.resize(edge_count);

    int vertex1, vertex2, index = 0;
    while (in_file >> tmp >> vertex1 >> vertex2) {
        
        input_edge_cache[index].v1 = vertex1;
        input_edge_cache[index].v2 = vertex2;
        input_edge_cache[index].weight = generate_edge_weight(vertex1, vertex2);
        index++;
    }
    in_file.close();
    return true;
}
#endif

#ifdef format_konect
// Read the input file and transform it into a structed vector
bool read_input_file(string file_name) {

    cout << file_name << ":" << endl;
    ifstream in_file(file_name);
    if (!in_file.is_open()) {
        cout << "in_file error" << endl;
        return false;
    }
    // Get vertex_count and edge_count
    string line;
    istringstream is;
    string p, tmp;
    do {
        getline(in_file, line);
        is.clear();
        is.str(line);
        is >> p >> edge_count >> vertex_count;
    } while (p != "p");

    cout << edge_count << " " <<vertex_count << endl;
    input_edge_cache.resize(edge_count);

    int vertex1, vertex2, index = 0, weight;
    while (in_file >> vertex1 >> vertex2 >> weight) {
        /*if (index >= input_edge_cache.size()) {
            cout << "The number of edges is wrong!" << endl;
            return false;
        }
        if (vertex1 > vertex_count || vertex2 > vertex_count) {
            cout << "The number of vertices is wrong!" << endl;
            return false;
        }*/
        input_edge_cache[index].v1 = vertex1;
        input_edge_cache[index].v2 = vertex2;
        input_edge_cache[index].weight = weight;
        index++;
    }
    in_file.close();
    return true;
}
#endif

// Build vertex_adjacency_list.
void build(string file_name) {
    // Open the input file
    long list_length = vertex_count + 1;
    vertex_adjacency_list.resize(list_length);
    vertex_edge_weight = (long *)calloc(list_length, sizeof(long));
    vertex_degree = (long *)calloc(list_length, sizeof(long));

//#ifdef is_simplify
    vertex_upper_bound = (long *)calloc(list_length, sizeof(long));
#ifdef is_simplify
    is_pending.resize(list_length);
#endif
    for (auto e : input_edge_cache) {
        // Update the degree of vertices
        vertex_degree[e.v1] ++;
        vertex_degree[e.v2] ++;
    }
   
    for (int v=1; v<=vertex_count; v++) {
        if (vertex_degree[v]) {
            vertex_adjacency_list[v].resize(vertex_degree[v]);
        }
    }
    int* vertex_list = (int *)calloc(list_length, sizeof(int));

//#ifdef is_simplify
    vertex_max_weight = (int *)calloc(list_length, sizeof(int));
    vertex_min_weight = (int *)calloc(list_length, sizeof(int));
//#endif    
    for (auto e : input_edge_cache) {
        weighted_edge tmp_weighted_edge;

        // Construct a list of indent edges and their weights
        tmp_weighted_edge.weight = e.weight;
        
        tmp_weighted_edge.vertex = e.v2;
        vertex_adjacency_list[e.v1][vertex_list[e.v1]++] = tmp_weighted_edge;
        
        tmp_weighted_edge.vertex = e.v1;
        vertex_adjacency_list[e.v2][vertex_list[e.v2]++] = tmp_weighted_edge;
        
        // Sum up the edge weights indent to the vertex
        vertex_edge_weight[e.v1] += tmp_weighted_edge.weight;
        vertex_edge_weight[e.v2] += tmp_weighted_edge.weight;

//#ifdef is_simplify
        if (0 == vertex_min_weight[e.v1] || vertex_min_weight[e.v1] > e.weight) {
            vertex_min_weight[e.v1] = e.weight;
        }
        if (vertex_max_weight[e.v1] < e.weight) {
            vertex_max_weight[e.v1] = e.weight;
        }
        if (0 == vertex_min_weight[e.v2] || vertex_min_weight[e.v2] > e.weight) {
            vertex_min_weight[e.v2] = e.weight;
        }
        if (vertex_max_weight[e.v2] < e.weight) {
            vertex_max_weight[e.v2] = e.weight;
        }
//#endif
    }
    free(vertex_list);
    // Initialize the rest of global variables

    remaining_vertex.init(vertex_count);

    //select vertex;
    best_add_array = (best_element *)malloc(sizeof(best_element)*list_length);
    best_swap_array = (best_element *)malloc(sizeof(best_element)*list_length);
    best_drop_array = (best_element *)malloc(sizeof(best_element)*list_length);

    cand_neighbor_weight = (long *)calloc(list_length, sizeof(long));

    is_in_candidates = (bool*)calloc(list_length, sizeof(bool));
    is_addv_neighbor = (is_add_nei_element*)calloc(list_length, sizeof(is_add_nei_element));
    is_in_neighbor = (is_neighbor_element*)calloc(list_length, sizeof(is_neighbor_element));

    add_vertex_score = (long*)calloc(list_length, sizeof(long));
    drop_vertex_score = (long*)calloc(list_length, sizeof(long));
    is_in_solution = (bool*)calloc(list_length, sizeof(bool));
    vertex_num_connect_to_clique = (vertex_connect_clique*)calloc(list_length, sizeof(vertex_connect_clique));
    tmp_is_v1 = (bool*)calloc(list_length, sizeof(bool));
    tmp_v2_in_swap = (bool*)calloc(list_length, sizeof(bool));
    
    conf_change = (bool *)malloc(list_length * sizeof(bool));
    vertex_step = (unsigned long long *)calloc(list_length, sizeof(unsigned long long));
    add_candidate_index = (long *)calloc(list_length, sizeof(long));
    drop_candidate_index = (long *)calloc(list_length, sizeof(long));

    real_bms_count = min_bms_count;
}

//#ifdef is_simplify
long soft_upper_bound(long vertex) {
    int degree = vertex_degree[vertex];
    if (0 == degree)
        return 0;
    long ans = vertex_edge_weight[vertex];
    if (1 == degree) {
        return ans;
    }
    
    for (auto u : vertex_adjacency_list[vertex]) {
        if (degree >= vertex_degree[u.vertex]) {
            ans += vertex_edge_weight[u.vertex];
        }
        else {
            ans += min_num(degree * vertex_max_weight[u.vertex], vertex_edge_weight[u.vertex] - (vertex_degree[u.vertex] - degree) * vertex_min_weight[u.vertex]);
        }
    }
    return ans/2;
}


#ifdef is_simplify
long simplify_iterative() {
    long ans = 0;
    cout << "simplify_iterative" <<endl;
    //simplify by degree
    long threshold_weight_degree = best_solution_weight;
    working_vertex.clear();
    next_working_vertex.clear();

    for (auto v : remaining_vertex) {
        vertex_num_connect_to_clique[v].number = 0;
        vertex_num_connect_to_clique[v].weight = 0;
        //if(vertex_degree[v] < simplify_threshold) { // && is_in_candidates[v] == 0) {
            working_vertex.push_back(v);
            is_pending[v] = true;
        //}
    }
    long sign_aaa = 0;
    while (!working_vertex.empty()) {
        for (vector<long>::size_type i = 0; i < working_vertex.size(); ++i) {
            auto v = working_vertex[i];
            long aaa;
            if (0 == vertex_upper_bound[v]) {
                aaa = soft_upper_bound(v);
                vertex_upper_bound[v] = aaa;
            }
            else {
                sign_aaa ++;
                aaa = vertex_upper_bound[v];
            }
           
            //delete v
            if (aaa < threshold_weight_degree) {
                // delete edge
                for (auto u : vertex_adjacency_list[v]) {
                    vector<long>::size_type j = 0;
                    for (; j < vertex_adjacency_list[u.vertex].size(); ++j) {
                        if (vertex_adjacency_list[u.vertex][j].vertex == v) {
                            break;
                        }
                    }
                    vertex_adjacency_list[u.vertex][j] = *vertex_adjacency_list[u.vertex].rbegin();
                    vertex_adjacency_list[u.vertex].pop_back();
                    vertex_edge_weight[u.vertex] -= u.weight;
                    vertex_degree[u.vertex] --;
                    edge_count --;
                    vertex_upper_bound[u.vertex] = 0;
                   
                    if (vertex_degree[u.vertex]) {
                        if (vertex_min_weight[u.vertex] == u.weight) {
                            vertex_min_weight[u.vertex] = vertex_max_weight[u.vertex];
                            for (auto t : vertex_adjacency_list[u.vertex]) {
                                vertex_upper_bound[t.vertex] = 0;
                                if (t.weight < vertex_min_weight[u.vertex]) { 
                                    vertex_min_weight[u.vertex] = t.weight;
                                }
                            }
                        }
                        else if (vertex_max_weight[u.vertex] == u.weight) {
                            vertex_max_weight[u.vertex] = vertex_min_weight[u.vertex];
                            for (auto t : vertex_adjacency_list[u.vertex]) {
                                vertex_upper_bound[t.vertex] = 0;
                                if (t.weight > vertex_max_weight[u.vertex]) {
                                    vertex_max_weight[u.vertex] = t.weight;
                                }
                            }
                        }
                    } 
                    
                    //if ((!is_pending[u.vertex]) && vertex_degree[u.vertex] < simplify_threshold && is_in_candidates[u.vertex] == 0) {
                    if ((!is_pending[u.vertex]) && is_in_candidates[u.vertex] == 0) {
                        next_working_vertex.push_back(u.vertex);
                        is_pending[u.vertex] = true;
                    }
                }
                //remove i
                vertex_adjacency_list[v].clear();
                remaining_vertex.remove(v);
                ans ++;
            }
            is_pending[v] = false;
        }
        working_vertex.swap(next_working_vertex);
        //cout << "working_vertex size: " << working_vertex.size() << endl;
        next_working_vertex.clear();
    }
    //cout << "sign aaa:" << sign_aaa << endl;
    return ans;
}

long simplify() {
    long ans = 0;
    cout << "simplify" <<endl;
    //simplify by degree
    long threshold_weight_degree = best_solution_weight;
    working_vertex.clear();
    long working_size = 0;
    working_vertex.resize(remaining_vertex.size());

    for (auto v : remaining_vertex) {
        vertex_num_connect_to_clique[v].number = 0;
        vertex_num_connect_to_clique[v].weight = 0;
        //if (vertex_degree[v] < simplify_threshold) {
            working_vertex[working_size++] = v;
        //}
    }
    long v;
    for (long i=0; i<working_size; i++) {
        v = working_vertex[i];
        //if (!is_removed[v]) {
            long aaa;
            if (0 == vertex_upper_bound[v]) {
                aaa = soft_upper_bound(v);
                vertex_upper_bound[v] = aaa;
            }
            else {
                aaa = vertex_upper_bound[v];
            }

            if (aaa < threshold_weight_degree) {
                if (0 != vertex_degree[v]) {
                    // delete edge
                    for (auto u : vertex_adjacency_list[v]) {
                        vector<long>::size_type j = 0;
                        for (; j < vertex_adjacency_list[u.vertex].size(); ++j) {
                            if (vertex_adjacency_list[u.vertex][j].vertex == v) {
                                break;
                            }
                        }
                        vertex_adjacency_list[u.vertex][j] = *vertex_adjacency_list[u.vertex].rbegin();
                        vertex_adjacency_list[u.vertex].pop_back();
                        vertex_edge_weight[u.vertex] -= u.weight;
                        vertex_degree[u.vertex] --;
                        edge_count--;
                        vertex_upper_bound[u.vertex] = 0;
                        
                        /*if (0 == vertex_degree[u.vertex]) {
                            is_removed[u.vertex] = true;
                            remaining_vertex.remove(u.vertex);
                            ans++;
                        }*/
                        if (0 != vertex_degree[u.vertex]) {
                            if (vertex_min_weight[u.vertex] == u.weight) {
                                vertex_min_weight[u.vertex] = vertex_max_weight[u.vertex];
                                for (auto t : vertex_adjacency_list[u.vertex]) {
                                    vertex_upper_bound[t.vertex] = 0;
                                    if (t.weight < vertex_min_weight[u.vertex]) {
                                        vertex_min_weight[u.vertex] = t.weight;
                                    }
                                }
                            }
                            else if (vertex_max_weight[u.vertex] == u.weight) {
                                vertex_max_weight[u.vertex] = vertex_min_weight[u.vertex];
                                for (auto t : vertex_adjacency_list[u.vertex]) {
                                    vertex_upper_bound[t.vertex] = 0;
                                    if (t.weight > vertex_max_weight[u.vertex]) {
                                        vertex_max_weight[u.vertex] = t.weight;
                                    }
                                }
                            }
                        }
                    }
                }
                //remove i
                vertex_adjacency_list[v].clear();
                //is_removed[v] = true;
                remaining_vertex.remove(v);
                ans ++;
            //}
        }
    }
    return ans;
}
#endif

bool verify_simple(string file_name) {
    cout << "verifying..." << endl;
    ifstream in_file(file_name);
    if (! in_file.is_open()) {
        cout << "in_file error" << endl;
        exit(1);
    }

    //get vertex_count
    string line;
    istringstream is;
    string p, tmp;
    do {
        getline(in_file, line);
        is.clear();
        is.str(line);
        is >> p >> tmp >> vertex_count;
    } while (p != "p");

    // adjacency list
    long weight_p;
#ifdef format_dimacs
    vector<vector<long>> adj_list;
#endif

#ifdef format_konect
    vector<vector<weighted_edge>> adj_list;
#endif

    adj_list.resize(vertex_count + 1);
    long v1, v2;
    cout << vertex_count << endl;

#ifdef format_dimacs
    while (in_file >> tmp >> v1 >> v2) {
        adj_list[v1].push_back(v2);
        adj_list[v2].push_back(v1);
    }
#endif 

#ifdef format_konect
    while (in_file >> v1 >> v2 >> weight_p) {
        weighted_edge tmp_edge;
        tmp_edge.weight = weight_p;
        
        tmp_edge.vertex = v2;
        adj_list[v1].push_back(tmp_edge);
        
        tmp_edge.vertex = v1;
        adj_list[v2].push_back(tmp_edge);
    }
#endif

    in_file.close();

    //check clique
    long sw = 0;
#ifdef format_dimacs
    for (vector<long>::size_type i = 0; i < best_solution.size(); ++i) {
        long v1 = best_solution[i];
        for (vector<long>::size_type j = i + 1; j < best_solution.size(); ++j) {
            long v2 = best_solution[j];
            vector<long>::size_type k = 0;
            for (; k < adj_list[v1].size(); ++k) {
                if (v2 == adj_list[v1][k]) {
                    sw += generate_edge_weight(v1, v2);
                    break;
                }
            }
            if (k >= adj_list[v1].size()) {
                cerr << "wrong anser: " << v1 << " is not neighbor of " << v2 << "!" << endl;
                return false;
            }
        }
    }
#endif

#ifdef format_konect
    for (vector<long>::size_type i = 0; i < best_solution.size(); ++i) {
        long v1 = best_solution[i];
        for (vector<long>::size_type j = i + 1; j < best_solution.size(); ++j) {
            long v2 = best_solution[j];
            vector<long>::size_type k = 0;
            for (; k < adj_list[v1].size(); ++k) {
                if (v2 == adj_list[v1][k].vertex) {
                    sw += adj_list[v1][k].weight;
                    break;
                }
            }
            if (k >= adj_list[v1].size()) {
                cerr << "wrong anser: " << v1 << " is not neighbor of " << v2 << "!" << endl;
                return false;
            }
        }
    }
#endif
    
    cout << "sw: "<< sw<<endl;
    if(sw != best_solution_weight){
        return false;
    }
    cout << "solution verified." << endl;
    return true;
}
bool parse_parameters(int argc, char **argv) {
    int i = 0;
    int temp_para = 0;

    for (i=1; i<argc; i++) {
        if (0 == strcmp(argv[i], "-inst")) {
            i++;
            if (i >= argc)
                return false;
            input_filename = argv[i];
        }
        else if (0 == strcmp(argv[i], "-seed")) {
            i++;
            if (i >= argc)
                return false;
            seed = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-cutoff")) {
            i++;
            if (i >= argc)
                return false;
            cutoff_time = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-min_bms")) {
            i++;
            if (i >= argc)
                return false;
            min_bms_count = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-max_bms")) {
            i++;
            if (i >= argc)
                return false;
            max_bms_count = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-Waim")) {
            i++;
            if (i >= argc)
                return false;
            Waim = atoi(argv[i]);
        }
    }
    return true;
}

int main(int argc, char *argv[]) {

    // The state of final output
    bool exact = false;
    if (argc<7) {
        cout << "./ReConSLS -inst <inst> -seed <seed> -cutoff <cutoff_time>" << endl;
        exit(0);
    }
    // Read the input parameters
    // sscanf(argv[i++],"%ld",&min_bms_count);
    // sscanf(argv[i++],"%ld",&max_bms_count);

    parse_parameters(argc, argv);
    srand(seed);
    clock_t start = clock();
    cout << input_filename << ": " << endl;
    if (!read_input_file(input_filename)) {
        return 0;
    }
    build(input_filename);

    // Main loop
    start_vertices = remaining_vertex.vertex;

#ifdef is_simplify
    better_since_cons = false;
#endif
    //untest_pointer = 0;
    for (tries = 1; tries <= max_tries; tries++) {
        
        //if (tries % 1000000 == 0) srand(seed++);
        
        if ((double) (clock() - start) / CLOCKS_PER_SEC > cutoff_time) break;

        memset(conf_change, 1, sizeof(bool)*(vertex_count+1));
        local_search();

#ifdef is_simplify
        is_new_graph = false;
        if (better_since_cons && init_num - last_improved_no >= 10) {
        //if (better_since_cons) {
            better_since_cons = false;
            long remove_vertex_num = 0;
            if (remaining_vertex.size() > size_threshold) {
            //if (edge_count > size_threshold) {
                remove_vertex_num = simplify();
            }
            else {
                remove_vertex_num = simplify_iterative();
            }
            if (remove_vertex_num > 0) {
                //cout << "remove vertex num:" << remove_vertex_num << " remaining vertex num:" << remaining_vertex.size() << " edge count:" << edge_count << endl;
                
                is_new_graph = true;
                simp_count++;

                if (remaining_vertex.size() == best_solution.size()) {
                    exact = true;
                    break;
                }
            }
        }
#endif

    }

    cout << "uuuu," <<best_solution_weight << ',' << best_solution_time << ',' << best_solution_try;
    if (exact) cout << ",x";
    else cout << ",h";

    cout << ", " << tries << "," << simp_count;
    cout << endl;
    cout << "init num: " << init_num <<endl;

    cout << "add_num: " <<add_num<<"\t" <<"swap_num: " <<swap_num << "\t" <<"drop_num: "<<drop_num <<endl; 
    verify_simple(input_filename);
   
    free(conf_change);
    free(vertex_step);
    free(add_candidate_index);
    free(drop_candidate_index);
    free(vertex_upper_bound);
    //free(is_removed);
    free(best_add_array);
    free(best_swap_array);
    free(best_drop_array);
    free(cand_neighbor_weight);
    free(is_in_candidates);
    free(is_addv_neighbor);
    free(is_in_neighbor);
    free(add_vertex_score);
    free(drop_vertex_score);
    free(is_in_solution);
    free(vertex_num_connect_to_clique);
    free(tmp_is_v1);
    free(tmp_v2_in_swap);
    free(vertex_max_weight);
    free(vertex_min_weight);
    free(vertex_edge_weight);
    free(vertex_degree);
    return 0;
}
