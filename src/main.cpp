#include <iostream>
#include <vector>
#include <fstream>
#include <igraph.h>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <chrono>

#define VERBOSE true

using namespace std;

typedef struct edge_attributes_struct{
    igraph_integer_t edge_idx;
    igraph_real_t label;
    igraph_real_t variant;
}edge_attributes;

void print_variation_graph(const igraph_t *graph){
    igraph_vs_t vs;
    igraph_vit_t vit;

    igraph_vs_all(&vs);
    igraph_vit_create(graph, vs, &vit);
    while (!IGRAPH_VIT_END(vit)) {
        cout << IGRAPH_VIT_GET(vit) << ": " << endl;

        igraph_vector_int_t edges;
        igraph_vector_int_init(&edges, 0);
        igraph_incident(graph, &edges, IGRAPH_VIT_GET(vit), IGRAPH_OUT);
        for(long int i = 0; i < igraph_vector_int_size(&edges); i++){
            igraph_integer_t start, end;
            igraph_edge(graph, igraph_vector_int_get(&edges, i), &start, &end);
            cout << "\t(" << start << ", " << end << ", " 
                << "label: " << igraph_cattribute_EAN(graph, "label", igraph_vector_int_get(&edges, i)) << ", "
                << "variant: " << igraph_cattribute_EAN(graph, "variant" , igraph_vector_int_get(&edges, i)) << ")" << endl; 
            ;
        }
        IGRAPH_VIT_NEXT(vit);
        cout << endl;
    }
    cout<< endl;
}


void print_alignment_graph(const igraph_t *graph){

    cout << "\nPrinting alignment graph:" << endl;
    cout << "pos = " << igraph_cattribute_GAN(graph, "position")
        << ", substring = " << igraph_cattribute_GAS(graph, "substring") << endl;

    igraph_vs_t vs;
    igraph_vit_t vit;

    igraph_vs_all(&vs);
    igraph_vit_create(graph, vs, &vit);
    while (!IGRAPH_VIT_END(vit)) {
        cout << IGRAPH_VIT_GET(vit) << ": " << endl;

        igraph_vector_int_t edges;
        igraph_vector_int_init(&edges, 0);
        igraph_incident(graph, &edges, IGRAPH_VIT_GET(vit), IGRAPH_OUT);
        for(long int i = 0; i < igraph_vector_int_size(&edges); i++){
            igraph_integer_t start, end;
            igraph_edge(graph, igraph_vector_int_get(&edges, i), &start, &end);
            cout << "\t(" << start << ", " << end << ", " 
                << "weight: " << igraph_cattribute_EAN(graph, "weight", igraph_vector_int_get(&edges, i)) << ", "
                << "variant: " << igraph_cattribute_EAN(graph, "variant" , igraph_vector_int_get(&edges, i)) << ")" << endl; 
            ;
        }
        IGRAPH_VIT_NEXT(vit);
        cout << endl;
    }
    cout<< endl;
}


void read_edge_file(string edge_file_name, igraph_t* graph, int* num_variants){
    
    ifstream infile(edge_file_name);

    int num_vertices = -1;

    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);

    int start_vertex, end_vertex;
    char edge_label, variant_number;

    vector<edge_attributes> attributes_vector;
    int i = 0;
    while (infile >> start_vertex >> end_vertex >> edge_label >> variant_number)
    {
        if(start_vertex > num_vertices-1){
            num_vertices = start_vertex + 1;
        }
        if(end_vertex > num_vertices-1){
            num_vertices= end_vertex + 1;
        }

        igraph_real_t label;
        if(edge_label == 'A'){
            label = 0;
        }else if(edge_label == 'T'){
            label = 1;
        }else if(edge_label == 'C'){
            label = 2;
        }else if(edge_label == 'G'){
            label = 3;
        }

        igraph_real_t variant;

        if(variant_number != '-'){
            attributes_vector.push_back({i, label, (igraph_real_t) (variant_number - '0')});
        }else{
            attributes_vector.push_back({i, label, -1});
        }

        //add edge
        igraph_vector_int_push_back(&edges, start_vertex);
        igraph_vector_int_push_back(&edges, end_vertex);
        i++;
    }

    bool directed = true;
    igraph_empty(graph, num_vertices, directed);
    igraph_add_edges(graph, &edges, NULL);

    int max_variant = -1;
    for(edge_attributes a: attributes_vector){
        igraph_cattribute_EAN_set(graph, "label", a.edge_idx, a.label);
        igraph_cattribute_EAN_set(graph, "variant", a.edge_idx, a.variant);
        if(a.variant > max_variant){
            max_variant = a.variant;
        }
    }

    *num_variants = max_variant + 1;

}


void read_pos_substring_file(string pos_string_file_name, vector<int>& positions, vector<string>& substrings){
    ifstream infile(pos_string_file_name);

    int pos;
    string substring;
    while (infile >> pos >> substring){
        positions.push_back(pos);
        substrings.push_back(substring);
    }
    /*
    for(string s : substrings){
        cout << s << endl;
    }

    for(int p : positions){
        cout << p << endl;
    }
    */
}


void get_reachable_vertex(igraph_t* graph, igraph_vs_t vertex_selector, int dist, igraph_vector_int_list_t* reachable_list){
    
    igraph_integer_t order = dist;
    igraph_integer_t mindist = 0;

    igraph_neighborhood(graph, reachable_list, vertex_selector, order, IGRAPH_OUT, mindist);
}


// obtains induced graph in time polylogarithmic in the size of the induced graph
// this is written because the native igraph library requires size proportional to the size of original graph
void get_induced_variation_graph(igraph_t* graph, igraph_vector_int_t* vertex_set, igraph_t* result){
    unordered_map<int, int> umap;
    int n = igraph_vector_int_size(vertex_set);
    for(int i = 0; i < n; i++){
        // unordered map records the existence of a vertex and gives it a new index
        umap[igraph_vector_int_get(vertex_set, i)] = i;
    }

    // tempory vector to incident edges for each vertex
    igraph_vector_int_t incident_edges;
    igraph_vector_int_init(&incident_edges, 0);

    igraph_vector_int_t induced_edges;
    igraph_vector_int_init(&induced_edges, 0);

    vector<edge_attributes> attributes_vector;

    int new_eid = 0;
    for(int i = 0; i < n; i++){
        igraph_incident(graph, &incident_edges, igraph_vector_int_get(vertex_set, i), IGRAPH_OUT);
        
        //cout << igraph_vector_int_get(vertex_set, i) << "~: ";
        
        int m_v = igraph_vector_int_size(&incident_edges);
        for (int j = 0; j < m_v; j++){
            igraph_integer_t start, end;
            int old_eid = igraph_vector_int_get(&incident_edges, j);
            igraph_edge(graph, old_eid, &start, &end);
            
            //cout << "(" << start << ", " << end << ") ";
            
            // check if end is in the vertex_set
            if(umap.find(end) != umap.end()){
                //add edge
                igraph_vector_int_push_back(&induced_edges, umap[start]);
                igraph_vector_int_push_back(&induced_edges, umap[end]);
                attributes_vector.push_back({new_eid, 
                    igraph_cattribute_EAN(graph, "label", old_eid), 
                    igraph_cattribute_EAN(graph, "variant", old_eid)});
                new_eid++;
            }
        }
        //cout << endl;
    }
    bool directed = true;
    igraph_empty(result, n, directed);
    igraph_add_edges(result, &induced_edges, NULL);
    igraph_vector_int_destroy(&induced_edges);
    igraph_vector_int_destroy(&incident_edges);

    for(edge_attributes a: attributes_vector){
        igraph_cattribute_EAN_set(result, "label", a.edge_idx, a.label);
        igraph_cattribute_EAN_set(result, "variant", a.edge_idx, a.variant);
    }

    //print_variation_graph(result);
}


void reverse_weighted_reachable(igraph_t* graph, igraph_integer_t start_vid, igraph_vector_int_t* result, int delta){
    
    // topologically sort the vertices of the reversed graph
    igraph_vector_int_t sorted_vertex;
    igraph_vector_int_init(&sorted_vertex, igraph_vcount(graph));
    igraph_topological_sorting(graph, &sorted_vertex, IGRAPH_IN);
    int n = igraph_vector_int_size(&sorted_vertex);

    //cout << "size = " << n << endl;

    // get start index
    int start_idx = -1;

    //cout << "reverse top order: ";
    for(int i = 0; i < n; i++){
        //cout << i << ":" << igraph_vector_int_get(&sorted_vertex, i) << ", ";
        if(start_vid == igraph_vector_int_get(&sorted_vertex, i)){
            start_idx = i;
            break;
        }
    }
    //cout << endl;

    vector<int> distance(n, INT_MAX-1);
    //cout << "index of start: " << igraph_vector_int_get(&sorted_vertex, start_idx) << endl;
    distance[igraph_vector_int_get(&sorted_vertex, start_idx)] = 0;

    igraph_vector_int_t incident_edges;
    igraph_vector_int_init(&incident_edges, 0);

    for(int i = start_idx; i < n; i++){
        igraph_integer_t vid = igraph_vector_int_get(&sorted_vertex, i);
        // update neightbors
        igraph_incident(graph, &incident_edges, vid, IGRAPH_IN);
        int m = igraph_vector_int_size(&incident_edges);
        //cout << "vid = " << vid  << endl;
        for(int i = 0; i < m; i++){

            igraph_integer_t start, end;
            igraph_edge(graph, igraph_vector_int_get(&incident_edges, i), &start, &end);
            int weight = igraph_cattribute_EAN(graph, "weight", igraph_vector_int_get(&incident_edges, i));
            //cout << "start = " << start << " end = " << end << " weight = " << weight << endl;

            if(distance[start] > distance[end] + weight){
                distance[start] = distance[end] + weight;
            }
        }
    }

    for(int i = 0; i < n; i++){
        //cout << i << ":" << distance[i] << ", ";
        if(distance[i] < delta){
            igraph_vector_int_push_back(result, i);
        }
    }

}


/*
The main idea is that we consider the vertices of the alignment graph as laided out in |P||V| grid.
Each row of length |V| is sorted in the same topological order.
We hash x,y coordinates to track if a vertex as already been encountered and store as value the min dist
result graph should not be initialized
*/
void create_pruned_alignment_graph(igraph_t* graph, string pattern, int delta, igraph_t* result){

    // topologically sort the vertices, this will be used to get the delta_x for every edge considered
    igraph_vector_int_t sorted_vertex;
    igraph_vector_int_init(&sorted_vertex, igraph_vcount(graph));
    igraph_topological_sorting(graph, &sorted_vertex, IGRAPH_OUT);
    int n = igraph_vector_int_size(&sorted_vertex);
    vector<int> inv_sorted_vertex(n);
    for(int i = 0; i < n; i++){
        inv_sorted_vertex[igraph_vector_int_get(&sorted_vertex, i)] = i;
    }

    int sink_id = -1;
    /*
    for(int i = 0; i < n; i++){
        cout << igraph_vector_int_get(&sorted_vertex, i) << ", ";
    }
    cout << endl;
    for(int i = 0; i < n; i++){
        cout << inv_sorted_vertex[i] << ", ";
    }
    cout << endl;
    */

    typedef struct vertex_struct{
        int id;
        int x;
        int y;
        int dist;
    }vertex;

    //  hash to track seen coordinates, distances are stored here, not in queue
    unordered_map<int, vertex> seen_coordinates;

    igraph_vector_int_t alignment_graph_edges; // start, end,   start, end, 
    vector <int> weights; 
    vector <int> variants;

    igraph_vector_int_init(&alignment_graph_edges, 0);

    queue<vertex> q;
    vertex u = {0, inv_sorted_vertex[0], 0, 0};
    q.push(u);

    int new_vid = 1;
    igraph_vector_int_t incident_edges;
    igraph_vector_int_init(&incident_edges, 0);

    while(!q.empty()){

        vertex u = q.front();
        q.pop();
        igraph_integer_t old_uid = igraph_vector_int_get(&sorted_vertex, u.x);
        

        int u_coord = u.y * n + u.x;
        vertex u_ = seen_coordinates[u_coord];
        int u_dist = u_.dist;

        //cout << u.x << ", " << u.y << ": " << u_dist << endl;

        // vertices on same level (insertion)
        if(u_dist + 1 <= delta){

            igraph_incident(graph, &incident_edges, old_uid, IGRAPH_OUT);
            int m = igraph_vector_int_size(&incident_edges);
            //cout << "m = " << m << endl;
            for(int i = 0; i < m; i++){

                igraph_integer_t start, end;
                igraph_edge(graph, igraph_vector_int_get(&incident_edges, i), &start, &end);

                int v_x = inv_sorted_vertex[end];
                int v_y = u.y;
                int v_coordinate = v_y * n + v_x;
                int vid;

                if(seen_coordinates.find(v_coordinate) == seen_coordinates.end()){
                    // never encounter coordinate
                    // add vertex with new vertex id and edge

                    vid = new_vid;
                    seen_coordinates[v_coordinate] = {vid, v_x, v_y, u_dist + 1};
                    new_vid++;
                    
                    
                    // enque vertex
                    q.push({vid, v_x, v_y, -1});

                }else{
                    // coordinate seen before 
                    vertex v = seen_coordinates[v_coordinate];
                    vid = v.id;

                    // update distance if needed
                    if(u_dist + 1 < v.dist){
                        seen_coordinates[v_coordinate] = {vid, v_x, v_y, u_dist + 1};
                    }
                }
                igraph_vector_int_push_back(&alignment_graph_edges, u.id);
                igraph_vector_int_push_back(&alignment_graph_edges, vid);
                weights.push_back(1);
                variants.push_back(igraph_cattribute_EAN(graph, "variant", igraph_vector_int_get(&incident_edges, i)));
            }
        }

        // deletion edges
        if(u.y < pattern.length() && u_dist + 1 <= delta){
            int v_x = u.x;
            int v_y = u.y + 1;
            int v_coordinate = v_y * n + v_x;

            int vid;
            // never seen before
            if(seen_coordinates.find(v_coordinate) == seen_coordinates.end()){

                vid = new_vid;
                seen_coordinates[v_coordinate] = {vid, v_x, v_y, u_dist + 1};
                new_vid++;

                // enque vertex
                q.push({vid, v_x, v_y, -1});

            }else{
                // coordinate seen before 
                vertex v = seen_coordinates[v_coordinate];
                vid = v.id;
                
                // update distance if needed
                if(u_dist + 1 < v.dist){
                    seen_coordinates[v_coordinate] = {vid, v_x, v_y, u_dist + 1};
                }
            }
            igraph_vector_int_push_back(&alignment_graph_edges, u.id);
            igraph_vector_int_push_back(&alignment_graph_edges, vid);
            weights.push_back(1);
            variants.push_back(-2);
        }

        // sub/match edges

        if(u.y < pattern.length()){
            igraph_incident(graph, &incident_edges, old_uid, IGRAPH_OUT);
            int m = igraph_vector_int_size(&incident_edges);

            for(int i = 0; i < m; i++){

                // check if symbol matches
                int sym = (int) igraph_cattribute_EAN(graph, "label", igraph_vector_int_get(&incident_edges, i));

                //cout << "(" << u.x << ", " << u.y << ") sym: " << sym << endl;
                if(    (pattern[u.y] == 'A' && sym == 0) 
                    || (pattern[u.y] == 'T' && sym == 1)
                    || (pattern[u.y] == 'C' && sym == 2)
                    || (pattern[u.y] == 'G' && sym == 3)){

                    igraph_integer_t start, end;
                    igraph_edge(graph, igraph_vector_int_get(&incident_edges, i), &start, &end);

                    int v_x = inv_sorted_vertex[end];
                    int v_y = u.y + 1;
                    int v_coordinate = v_y * n + v_x;
                    int vid;

                    if(seen_coordinates.find(v_coordinate) == seen_coordinates.end()){
                        // never encounter coordinate
                        // add vertex with new vertex id and edge

                        // Adds 0 to the distance from u
                        vid = new_vid;
                        seen_coordinates[v_coordinate] = {vid, v_x, v_y, u_dist};
                        new_vid++;
                        
                        // enque vertex
                        q.push({vid, v_x, v_y, -1});

                    }else{
                        // coordinate seen before 
                        vertex v = seen_coordinates[v_coordinate];
                        vid = v.id;

                        // update distance if needed
                        if(u_dist < v.dist){
                            seen_coordinates[v_coordinate] = {vid, v_x, v_y, u_dist};
                        }
                    }
                    igraph_vector_int_push_back(&alignment_graph_edges, u.id);
                    igraph_vector_int_push_back(&alignment_graph_edges, vid);
                    weights.push_back(0);
                    variants.push_back(igraph_cattribute_EAN(graph, "variant", igraph_vector_int_get(&incident_edges, i)));

                }else if (u_dist + 1 <= delta){

                    igraph_integer_t start, end;
                    igraph_edge(graph, igraph_vector_int_get(&incident_edges, i), &start, &end);

                    int v_x = inv_sorted_vertex[end];
                    int v_y = u.y + 1;
                    int v_coordinate = v_y * n + v_x;
                    int vid;

                    if(seen_coordinates.find(v_coordinate) == seen_coordinates.end()){
                        // never encounter coordinate
                        // add vertex with new vertex id and edge

                        // Adds 1 to the distance from u
                        vid = new_vid;
                        seen_coordinates[v_coordinate] = {vid, v_x, v_y, u_dist + 1};
                        new_vid++;
                        
                        // enque vertex
                        q.push({vid, v_x, v_y, -1});

                    }else{
                        // coordinate seen before 
                        vertex v = seen_coordinates[v_coordinate];
                        vid = v.id;

                        // update distance if needed
                        if(u_dist + 1 < v.dist){
                            seen_coordinates[v_coordinate] = {vid, v_x, v_y, u_dist + 1};
                        }
                    }
                    igraph_vector_int_push_back(&alignment_graph_edges, u.id);
                    igraph_vector_int_push_back(&alignment_graph_edges, vid);
                    weights.push_back(1);
                    variants.push_back(igraph_cattribute_EAN(graph, "variant", igraph_vector_int_get(&incident_edges, i)));
                }
            }
        }


        // edges to sink vertex

        if(u.y == pattern.length()){

            int v_y = u.y+1;
            int v_coordinate = n * v_y;

            int vid;
        
            if(seen_coordinates.find(v_coordinate) == seen_coordinates.end()){
                vid = new_vid;
                seen_coordinates[v_coordinate] = {vid, 0, v_y, u_dist};
                sink_id = vid;
                new_vid++;

                //q.push({vid, 0, n*(pattern.length()+1), -1});

            }else{
                // coordinate seen before 
                vertex v = seen_coordinates[v_coordinate];
                vid = v.id;

                // update distance if needed
                if(u_dist < v.dist){
                    seen_coordinates[v_coordinate] = {vid, 0, v_y, u_dist};
                }
            }

            igraph_vector_int_push_back(&alignment_graph_edges, u.id);
            igraph_vector_int_push_back(&alignment_graph_edges, vid);
            weights.push_back(0);
            variants.push_back(-3);

        }
    }
    // construct temporary forward-pruned only alignment graph
    bool directed = true;
    igraph_t forward_pruned_only_graph;
    igraph_empty(&forward_pruned_only_graph, new_vid, directed);
    igraph_add_edges(&forward_pruned_only_graph, &alignment_graph_edges, NULL);

    for(int i = 0; i < weights.size(); i++){
        igraph_cattribute_EAN_set(&forward_pruned_only_graph, "weight", i, weights.at(i));
        igraph_cattribute_EAN_set(&forward_pruned_only_graph, "variant", i, variants.at(i));
    }

    igraph_vector_int_destroy(&alignment_graph_edges);
    igraph_vector_int_destroy(&sorted_vertex);

    //cout << "============= Forward Prunned Alignment Graph Constructed =============" << endl;
    //print_alignment_graph(&forward_pruned_only_graph);


    // reverse pruning

    // get vertices reachable from the last vertex with distance at most delta
    igraph_vector_int_t reverse_reachable_vertex;
    igraph_vector_int_init(&reverse_reachable_vertex, 0);

    reverse_weighted_reachable(&forward_pruned_only_graph, sink_id, &reverse_reachable_vertex, delta);
    int k = igraph_vector_int_size(&reverse_reachable_vertex);
    
    /*
    cout << "Reverse reachable: ";
    for(int i = 0; i < k; i++){
        cout << igraph_vector_int_get(&reverse_reachable_vertex, i) << ", ";
    }
    cout << endl;
    */

    igraph_vs_t vs;
    igraph_vs_vector(&vs, &reverse_reachable_vertex);

    igraph_induced_subgraph(&forward_pruned_only_graph, result, vs, IGRAPH_SUBGRAPH_AUTO);
    igraph_destroy(&forward_pruned_only_graph);

} 


void create_alignment_graph(igraph_t* variation_graph, int pos, string s, int alpha, int delta, igraph_t* alignment_graph){
    
    igraph_vs_t start_vertex_selector;
    igraph_vs_1(&start_vertex_selector, pos);

    igraph_vector_int_list_t reachable_vertex_list;
    igraph_vector_int_list_init(&reachable_vertex_list, 0);
    
    get_reachable_vertex(variation_graph, start_vertex_selector, alpha + delta, &reachable_vertex_list);
    igraph_vs_destroy(&start_vertex_selector);
    
    // since vertex selector is only for one vertex, we only need the first element
    igraph_vector_int_t* reachable_vertex = igraph_vector_int_list_get_ptr(&reachable_vertex_list, 0);
    
    igraph_t reachable_graph;
    get_induced_variation_graph(variation_graph, reachable_vertex, &reachable_graph);
    
    create_pruned_alignment_graph(&reachable_graph, s, delta, alignment_graph);

    igraph_vector_int_list_destroy(&reachable_vertex_list);
    igraph_destroy(&reachable_graph);
}

void get_variants(igraph_t* g, unordered_set<int>& variants){

    igraph_es_t es;
    igraph_es_all(&es, IGRAPH_EDGEORDER_ID);
    igraph_eit_t eit;
    igraph_eit_create(g, es, &eit);

    while (!IGRAPH_EIT_END(eit)) {
        int variant = igraph_cattribute_EAN(g, "variant" ,IGRAPH_EIT_GET(eit));
        if(variant >= 0){
            variants.insert(variant);
        }

        IGRAPH_EIT_NEXT(eit);
    }
    igraph_es_destroy(&es);
    igraph_eit_destroy(&eit);
}

// returns a vector that assigns each alignment graph to a set using an integer
// also prints the number of sets and size of each corresponding ILP
void analyze_alignment_graph_set(vector<igraph_t>& alignment_graphs, int num_variants){
    

    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);

    int ilp_idx = num_variants;

    for(auto &g: alignment_graphs){
        //print_alignment_graph(&g);

        // obtain set of variants
        unordered_set<int> variants;
        get_variants(&g, variants);

        // add edges for bipartite graph
        for(int v: variants){
            //cout << ilp_idx << ", " << v << endl;
            igraph_vector_int_push_back(&edges, ilp_idx);
            igraph_vector_int_push_back(&edges, v);
        }
        ilp_idx++;
    }
    int n = num_variants + alignment_graphs.size();

    bool directed = false;
    igraph_t bipartite;
    igraph_empty(&bipartite, n, directed);
    igraph_add_edges(&bipartite, &edges, NULL);

    igraph_vector_int_t membership;
    igraph_vector_int_init(&membership, n);
    igraph_connected_components(&bipartite, &membership, NULL, NULL, IGRAPH_STRONG);

    if(VERBOSE){
        cout << "SubILP component assignments" << endl;
        for(int i = 0; i < num_variants; i++){
            cout << "x_" << i << ": " << igraph_vector_int_get(&membership, i) << endl;
        }
        cout << endl;

        for(int i = num_variants; i < n; i++){
            cout << "ILP_" << i << ": " << igraph_vector_int_get(&membership, i) << endl;
        }
        cout << endl;
    }
    // get number of binary variables for each component (size of each sub-ILP)
    unordered_map<int, int> sizes;
    for(int i = 0; i < n; i++){
        int component = igraph_vector_int_get(&membership, i);
        if(i < num_variants){
            if(sizes.find(component) == sizes.end()){
                sizes[component] = 1;
            }else{
                sizes[component] = sizes[component] + 1;
            }
        }else if(i >= num_variants){
            // add number of edges in graph
            igraph_t g = alignment_graphs.at(i - num_variants);
            sizes[component] = sizes[component] + igraph_ecount(&g);
        }
    }

    if(VERBOSE){
        cout << "Component sizes (number of variables added across ILPs in each component)" << endl;
        for (auto size : sizes){
            cout << size.first << ": " << size.second <<  endl;
        }
    }
    // clean up
    igraph_vector_int_destroy(&membership);

}

int main(int argc, char** argv){

    if(argc < 5){
        cout << "Wrong number of arguments provided" << endl;
        return 0;
    }
    string edge_file_name = argv[1];
    string pos_substring_file_name = argv[2];
    int alpha = atoi(argv[3]);
    int delta = atoi(argv[4]);

    igraph_set_attribute_table(&igraph_cattribute_table);

    // construct variation graph
    igraph_t variation_graph;
    int num_variants;

    auto start = chrono::steady_clock::now();
    read_edge_file(edge_file_name, &variation_graph, &num_variants);
    auto stop = chrono::steady_clock::now();
    auto duration = duration_cast<chrono::milliseconds>(stop - start);

    cout << "\n============= Variation Graph Constructed =============\n" << endl;
    cout << "construction time: " << duration.count() << " milliseconds\n" <<endl;
    cout << "number variants: " << num_variants << endl;
    cout << "number vertices: " << igraph_vcount(&variation_graph) << endl;
    cout << "number edges: " << igraph_ecount(&variation_graph) << endl;
    //print_variation_graph(&variation_graph);

    // read locations and strings
    vector<int> positions;
    vector<string> substrings;
    read_pos_substring_file(pos_substring_file_name, positions, substrings);

    int N = positions.size();
    vector<igraph_t> alignment_graphs;

    if(VERBOSE){
        cout << "\nConstructing alignment graphs...\n" << endl;
    }
    for(int i = 0; i < N; i++){
        if(VERBOSE){
            cout << "Constructing alignment graphs for position: " << positions[i] << endl;
        }
        auto start = chrono::steady_clock::now();
        igraph_t alignment_graph;
        create_alignment_graph(&variation_graph, positions[i], substrings[i], alpha, delta, &alignment_graph);
        alignment_graphs.push_back(alignment_graph);
        igraph_cattribute_GAN_set(&alignment_graph, "position", positions[i]);
        igraph_cattribute_GAS_set(&alignment_graph, "substring", substrings[i].c_str());
        auto stop = chrono::steady_clock::now();
        auto duration = duration_cast<chrono::milliseconds>(stop - start);
        cout << "time: " << duration.count() << " milliseconds\n" <<endl;
    }
    cout << "\n============= Alignment Graphs Constructed =============\n" << endl;

    cout << "\nAnalyzing alignment graph set..." << endl;
    analyze_alignment_graph_set(alignment_graphs, num_variants);

    cout << "\n============= Alignment Graphs Analyzed =============\n" << endl;
}
