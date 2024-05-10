#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <igraph.h>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <chrono>
#include <climits>
#include "gurobi_c++.h"

#define VERBOSE false
#define WRITE_ILPS_TO_FILE false
#define PRINT_COMPONENT_SIZES false
#define PRINT_FACTOR 100
#define ILP_TIME_LIMIT 3600 // in seconds

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
        igraph_vector_int_destroy(&edges);
        IGRAPH_VIT_NEXT(vit);
        cout << endl;
    }
    cout<< endl;
}


void print_alignment_graph(igraph_t *graph){

    cout << "\nPrinting alignment graph:" << endl;
    cout << "pos = " << igraph_cattribute_GAN(graph, "position")
        << ", substring = " << igraph_cattribute_GAS(graph, "substring") << endl;

    igraph_vs_t vs;
    igraph_vit_t vit;

    igraph_vs_all(&vs);
    igraph_vit_create(graph, vs, &vit);
    while (!IGRAPH_VIT_END(vit)) {
        cout << IGRAPH_VIT_GET(vit)  
            << " startend: " << igraph_cattribute_VAN(graph, "startend", IGRAPH_VIT_GET(vit))  
            << endl;

        igraph_vector_int_t edges;
        igraph_vector_int_init(&edges, 0);
        igraph_incident(graph, &edges, IGRAPH_VIT_GET(vit), IGRAPH_OUT);
        for(long int i = 0; i < igraph_vector_int_size(&edges); i++){
            igraph_integer_t start, end;
            igraph_edge(graph, igraph_vector_int_get(&edges, i), &start, &end);
            cout << "\t(" << start << ", " << end << ", " 
                << "weight: " << igraph_cattribute_EAN(graph, "weight", igraph_vector_int_get(&edges, i)) << ", "
                << "variant: " << igraph_cattribute_EAN(graph, "variant" , igraph_vector_int_get(&edges, i)) << ", "
                //<< "initsol: " << igraph_cattribute_EAN(graph, "initsol" , igraph_vector_int_get(&edges, i)) 
                << ")" << endl; 
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
    char edge_label;
    string variant_number;

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
        }else if(edge_label == 'N'){
            label = 4;
        }

        igraph_real_t variant;

        if(variant_number != "-"){
            attributes_vector.push_back({i, label, stod(variant_number.c_str())});
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
    igraph_vector_int_destroy(&edges);

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


void write_alignment_graph_to_file(igraph_t *graph, string file_name){
    
    string output = to_string(igraph_ecount(graph)) + "\n";

    igraph_es_t es;
    igraph_eit_t eit;

    igraph_es_all(&es, IGRAPH_EDGEORDER_ID);
    igraph_eit_create(graph, es, &eit);

    int source, sink;

    while (!IGRAPH_EIT_END(eit)) {
        igraph_integer_t start, end;
        igraph_edge(graph, IGRAPH_EIT_GET(eit), &start, &end);
        int weight = igraph_cattribute_EAN(graph, "weight", IGRAPH_EIT_GET(eit));
        int variant = igraph_cattribute_EAN(graph, "variant", IGRAPH_EIT_GET(eit));
        int initsol = igraph_cattribute_EAN(graph, "initsol", IGRAPH_EIT_GET(eit));
        output += to_string(start) + " " + to_string(end) + " " + to_string(weight) + " " + to_string(variant) + " " + to_string(initsol) + "\n";

        if(igraph_cattribute_VAN(graph, "startend", start) == 1){
            source = start;
        }
        if(igraph_cattribute_VAN(graph, "startend", end) == 2){
            sink = end;
        }

        IGRAPH_EIT_NEXT(eit);
    }

    igraph_es_destroy(&es);
    igraph_eit_destroy(&eit);

    output += to_string(source) + " " + to_string(sink);

    ofstream out_file(file_name);
    out_file << output;
    out_file.close();

    
}


void read_alignment_graph_from_file(igraph_t *graph, string file_name){
    
    ifstream in_file(file_name);

    int E;
    in_file >> E;

    if(E == 0){
        bool directed = true;
        igraph_empty(graph, 0, directed);
        return;
    }

    igraph_vector_int_t edges;
    igraph_vector_t weights;
    igraph_vector_t variants;
    igraph_vector_t initsols;

    igraph_vector_int_init(&edges, 2*E);
    igraph_vector_init(&weights, E);
    igraph_vector_init(&variants, E);
    igraph_vector_init(&initsols, E);

    int num_vertices = 0;
    int start, end; 
    float weight, variant, initsol;

    int i = 0;
    int j = 0;
    while (in_file >> start >> end >> weight >> variant >> initsol){
        //cout << start << " " << end << " " << weight << " " << variant << endl;
       
        //edges[i] = start;
        igraph_vector_int_set(&edges, i, start);
        
        //edges[i+1] = end;
        igraph_vector_int_set(&edges, i+1, end);
        
        i += 2;

        //weights[j] = weight;
        igraph_vector_set(&weights, j, weight);

        //variants[j] = variant;
        igraph_vector_set(&variants, j, variant);

        igraph_vector_set(&initsols, j, initsol);

        j++;

        if(start > num_vertices-1){
            num_vertices = start + 1;
        }
        if(end > num_vertices-1){
            num_vertices= end + 1;
        }
    }

    in_file.close();

    bool directed = true;
    igraph_empty(graph, num_vertices, directed);
    igraph_add_edges(graph, &edges, NULL);

    igraph_cattribute_EAN_setv(graph, "weight", &weights);
    igraph_cattribute_EAN_setv(graph, "variant", &variants);
    igraph_cattribute_EAN_setv(graph, "initsol", &initsols);
    igraph_cattribute_VAN_set(graph, "startend", start, 1);
    igraph_cattribute_VAN_set(graph, "startend", end, 2);

    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&variants);
    igraph_vector_destroy(&initsols);

}


void read_pos_substring_file(string pos_string_file_name, vector<int>& positions, vector<string>& substrings){
    ifstream infile(pos_string_file_name);
    string line;

    while (std::getline(infile, line)){

        int pos;
        string substring;
        
        istringstream iss(line); // separates line over spaces
        iss >> pos;

       if (line.size() > 0){
            while (iss >> substring){
                positions.push_back(pos);
                substrings.push_back(substring); 
            }
        }
    }
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

    igraph_vector_int_destroy(&incident_edges);
    igraph_vector_int_destroy(&sorted_vertex);

    for(int i = 0; i < n; i++){
        //cout << i << ":" << distance[i] << ", ";
        if(distance[i] <= delta){
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

    int source_vertex_id = 0;
    int sink_vertex_id = 0;

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
                    || (pattern[u.y] == 'G' && sym == 3)
                    || (pattern[u.y] == 'N' || sym == 4)
                    ){

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
                sink_vertex_id = vid;
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

    igraph_vector_t zeros;
    igraph_vector_init(&zeros, igraph_vcount(&forward_pruned_only_graph));
    igraph_cattribute_VAN_setv(&forward_pruned_only_graph, "startend", &zeros);
    igraph_cattribute_VAN_set(&forward_pruned_only_graph, "startend", source_vertex_id, 1);
    igraph_cattribute_VAN_set(&forward_pruned_only_graph, "startend", sink_vertex_id, 2);


    //cout << "============= Forward Prunned Alignment Graph Constructed =============" << endl;
    //print_alignment_graph(&forward_pruned_only_graph);
    //cout << "Forward pruned alignment graph size: " << igraph_vcount(&forward_pruned_only_graph) << endl;
    // reverse pruning

    // get vertices reachable from the last vertex with distance at most delta
    igraph_vector_int_t reverse_reachable_vertex;
    igraph_vector_int_init(&reverse_reachable_vertex, 0);

    reverse_weighted_reachable(&forward_pruned_only_graph, sink_vertex_id, &reverse_reachable_vertex, delta);

    igraph_vs_t vs;
    igraph_vs_vector(&vs, &reverse_reachable_vertex);

    igraph_induced_subgraph(&forward_pruned_only_graph, result, vs, IGRAPH_SUBGRAPH_AUTO);
    igraph_destroy(&forward_pruned_only_graph);
    igraph_vector_int_destroy(&reverse_reachable_vertex);

} 


void add_initial_solution(igraph_t* alignment_graph){

    // get source and sink
    int start_vertex = 0;
    int end_vertex = 0;
    igraph_vs_t vs;
    igraph_vit_t vit;
    igraph_vs_all(&vs);
    igraph_vit_create(alignment_graph, vs, &vit);
    while (!IGRAPH_VIT_END(vit)) {
        int vid = IGRAPH_VIT_GET(vit);
        int startend = (int)igraph_cattribute_VAN(alignment_graph, "startend", vid);
        if(startend == 1){
            start_vertex = vid;
        }else if (startend == 2){
            end_vertex = vid;
        }
        IGRAPH_VIT_NEXT(vit);
    }
    igraph_vs_destroy(&vs);
    igraph_vit_destroy(&vit);


    // get edge weights
    igraph_es_t es;
    igraph_es_all(&es, IGRAPH_EDGEORDER_ID);
    igraph_vector_t weights;
    igraph_vector_init(&weights, 0);
    igraph_cattribute_EANV(alignment_graph, "weight", es, &weights);
    igraph_es_destroy(&es);

    // get path with minimum weight
    igraph_vector_int_t path_edges;
    igraph_vector_int_init(&path_edges, 0);
    
    //cout << "num vertex in graph: " << igraph_vcount(alignment_graph) << endl;
    //cout << "start:" << start_vertex << " end vertex: " << end_vertex << endl;

    igraph_get_shortest_path_dijkstra(alignment_graph, NULL, &path_edges, start_vertex, end_vertex, &weights, IGRAPH_OUT);
    igraph_vector_destroy(&weights);

    // put it in attribute vector
    igraph_vector_t attribute_vector;
    igraph_vector_init(&attribute_vector, igraph_ecount(alignment_graph));
    int m = igraph_vector_int_size(&path_edges);
    for(int i = 0; i < m; i++){
        int eid = igraph_vector_int_get(&path_edges, i);
        //int weight = igraph_cattribute_EAN(alignment_graph, "weight", eid);
        //cout << "edge id: " << eid << " weight: " << weight << endl;
        igraph_vector_set(&attribute_vector, eid, 1);
    }
    igraph_vector_int_destroy(&path_edges);

    // set edge attrbitue
    igraph_cattribute_EAN_setv(alignment_graph, "initsol", &attribute_vector);
    igraph_vector_destroy(&attribute_vector);

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
    

    //cout << "induced variation graph size: " << igraph_vcount(&reachable_graph) << endl;

    int start_vertex, end_vertex;
    create_pruned_alignment_graph(&reachable_graph, s, delta, alignment_graph);

    //cout << "pruned alignment graph size: " << igraph_vcount(alignment_graph) << endl;

    add_initial_solution(alignment_graph);

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
void analyze_alignment_graph_set(int num_variants,
                                vector<vector<int>*>& alignment_graph_to_variants,
                                vector<int>& global_variable_to_component,
                                vector<int>& subILP_to_component,
                                vector<int>& global_variable_to_local_idx,
                                vector<int>& component_to_global_variable_count,
                                vector<vector<int>>& component_to_subILPs,
                                int* num_components
                                ){
    

    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);

    int ilp_idx = num_variants;

    for(int i = 0; i < alignment_graph_to_variants.size(); i++){
        //print_alignment_graph(&g);

        // obtain set of variants
        //unordered_set<int> variants;
        //get_variants(&g, variants);

        // add edges for bipartite graph
        //cout << i << ": ";
        for(int v: *alignment_graph_to_variants[i]){
            //cout << v << ", ";
            igraph_vector_int_push_back(&edges, ilp_idx);
            igraph_vector_int_push_back(&edges, v);
        }
        //cout << endl;
        ilp_idx++;

    }
    int n = num_variants + alignment_graph_to_variants.size();

    bool directed = false;
    igraph_t bipartite;
    igraph_empty(&bipartite, n, directed);
    igraph_add_edges(&bipartite, &edges, NULL);

    igraph_vector_int_destroy(&edges);

    igraph_vector_int_t membership;
    igraph_vector_int_init(&membership, n);
    igraph_connected_components(&bipartite, &membership, NULL, NULL, IGRAPH_STRONG);


    global_variable_to_component = vector<int>(num_variants);
    subILP_to_component = vector<int>(alignment_graph_to_variants.size());
    global_variable_to_local_idx = vector<int>(num_variants);

    unordered_map<int, int> sizes;

    *num_components = 0;
    for(int i = 0; i < n; i++){

        int component = igraph_vector_int_get(&membership, i);
        if(component + 1 > *num_components){
            *num_components = component + 1;
        }

        if(i < num_variants){
            
            global_variable_to_component[i] = component;

            if(sizes.find(component) == sizes.end()){
                sizes[component] = 1;
            }else{
                sizes[component] = sizes[component] + 1;
            }
            global_variable_to_local_idx[i] = sizes[component] - 1;
            
        }else{
            subILP_to_component[i - num_variants] = component;
        }
    }


    component_to_global_variable_count = vector<int>(*num_components, 0);
    component_to_subILPs = vector<vector<int>>(*num_components);
    for (auto size : sizes){
        component_to_global_variable_count[size.first] = size.second;
        component_to_subILPs[size.first] = vector<int>();
    }

    // clean up
    igraph_vector_int_destroy(&membership);

    
    // print statements for debugging
    if(VERBOSE){
        cout << "Global variable to component" << endl;
        for(int i = 0; i < global_variable_to_component.size(); i++){
            cout << "x_" << i << ": " << global_variable_to_component[i] << ", local:" << global_variable_to_local_idx[i] << endl;
        }

        cout << "\nSubILP to component" << endl;
        for(int i = 0; i < subILP_to_component.size(); i++){
            cout << "ILP_" << i << ": " << subILP_to_component[i] << endl;
        }
    }

    /*
    cout << "\nComponent to global variant count:\n" << endl;
    cout << "[";
    for(int i = 0; i < component_to_global_variable_count.size(); i++){
        //cout << "component_" << i << " has variants: " << component_to_global_variable_count[i] << endl;
        cout << component_to_global_variable_count[i] << ", ";
    }
    cout << "]" << endl;
    */

    unordered_map<int, int> component_to_num_edges;
    for(int i = 0; i < alignment_graph_to_variants.size(); i++){
        
        //igraph_t g = alignment_graphs.at(i);
        int component = subILP_to_component[i];
        component_to_subILPs[component].push_back(i);
        /*
        if(component_to_num_edges.find(component) == component_to_num_edges.end()){
            component_to_num_edges[component] = igraph_ecount(&g);
        }else{
            component_to_num_edges[component] = component_to_num_edges[component] + igraph_ecount(&g);
        }
        */
    }
    /*
    for(int i = 0; i < component_to_global_variable_count.size(); i++){
        cout << "component_" << i << " has total variables: " << component_to_global_variable_count[i] + component_to_num_edges[i] << endl;
    }
    */
}


void construct_ILPs(int component,
                    vector<vector<int>>& component_to_subILPs,
                    vector<int>& global_variable_to_local_idx,
                    vector<int>& component_to_global_variable_count,
                    GRBModel& model,
                    int delta,
                    string graph_directory,
                    vector<igraph_t*>& alignment_graphs,
                    bool use_files
                    ){


    GRBVar* global_var = model.addVars(component_to_global_variable_count[component]);

    GRBLinExpr obj = GRBLinExpr();
    for(int i = 0; i < component_to_global_variable_count[component]; i++){
        obj += global_var[i];

        global_var[i].set(GRB_DoubleAttr_Start, 0);
    }
    model.setObjective(obj, GRB_MAXIMIZE);


    if(PRINT_COMPONENT_SIZES){
        cout << "\nComponent "<< component << " has size (number of alignment graphs): " << component_to_subILPs[component].size() << endl;
    }

    for(int i = 0; i < component_to_subILPs[component].size(); i++){

        auto start = chrono::steady_clock::now();

        int alignment_graph_idx = component_to_subILPs[component][i];

        igraph_t alignment_graph;

        if(use_files == 1){
            string file_name = graph_directory + "g_" + to_string(alignment_graph_idx);
            //cout << file_name << endl;
            read_alignment_graph_from_file(&alignment_graph, file_name);
        }else{
            alignment_graph = *alignment_graphs[alignment_graph_idx];
        }

        if(igraph_ecount(&alignment_graph) == 0){
            cout << "no edges in alignment graph." << endl;
            continue;
        }

        /*
        if(VERBOSE){
            cout << "Alignment_graph_" << i << " ILP construction in component (model): " << component << endl;
        }else if(!VERBOSE && i % 1 == 0){
            cout << "Constructing ILP for alignment_graph: " << i << " out of " << component_to_subILPs[component].size() << " (increments of 1)" << endl;
        }
        */

        GRBVar* local_var = model.addVars(igraph_ecount(&alignment_graph), GRB_BINARY);

        igraph_vs_t vs;
        igraph_vit_t vit;
        igraph_bool_t loops = true;     // there are no loops, this is to get const. query time for degrees

        igraph_vs_all(&vs);
        igraph_vit_create(&alignment_graph, vs, &vit);
        
        
        while (!IGRAPH_VIT_END(vit)) {
            
            int vid = IGRAPH_VIT_GET(vit);
            
            igraph_integer_t in_degree, out_degree;
            igraph_degree_1(&alignment_graph, &in_degree, vid, IGRAPH_IN, loops);
            igraph_degree_1(&alignment_graph, &out_degree, vid, IGRAPH_OUT, loops);

            /*
            cout << "graph: " << i << " vertex: " << vid 
                << " in-degree: " << in_degree 
                << " out-degree: " << out_degree
                << " startend: " << igraph_cattribute_VAN(&alignment_graphs[i], "startend", vid) 
                << endl;
            */

            GRBLinExpr lhs, rhs;
            if(igraph_cattribute_VAN(&alignment_graph, "startend", vid) == 1){
                lhs = GRBLinExpr(1);
                rhs = GRBLinExpr();

                // iterate through outbound edges
                igraph_vector_int_t edges;
                igraph_vector_int_init(&edges, 0);
                igraph_incident(&alignment_graph, &edges, vid, IGRAPH_OUT);
                int m = igraph_vector_int_size(&edges);
                for(int j = 0; j < m; j++){
                    rhs += local_var[igraph_vector_int_get(&edges, j)];
                }
                igraph_vector_int_destroy(&edges);

            }else if(igraph_cattribute_VAN(&alignment_graph, "startend", vid) == 2){
                lhs = GRBLinExpr();
                rhs = GRBLinExpr(1);

                // iterate through inbound edges
                igraph_vector_int_t edges;
                igraph_vector_int_init(&edges, 0);
                igraph_incident(&alignment_graph, &edges, vid, IGRAPH_IN);
                int m = igraph_vector_int_size(&edges);
                for(int j = 0; j < m; j++){
                    lhs += local_var[igraph_vector_int_get(&edges, j)];
                }
                igraph_vector_int_destroy(&edges);                

            }else{
                lhs = GRBLinExpr();
                rhs = GRBLinExpr();

                igraph_vector_int_t edges;
                igraph_vector_int_init(&edges, 0);

                // iterate through in bound edges
                igraph_incident(&alignment_graph, &edges, vid, IGRAPH_IN);
                int m = igraph_vector_int_size(&edges);
                for(int j = 0; j < m; j++){
                    lhs += local_var[igraph_vector_int_get(&edges, j)];
                }

                // iterate through outbound edges
                igraph_incident(&alignment_graph, &edges, vid, IGRAPH_OUT);
                m = igraph_vector_int_size(&edges);
                for(int j = 0; j < m; j++){
                    rhs += local_var[igraph_vector_int_get(&edges, j)];
                }

                igraph_vector_int_destroy(&edges);  

            }

            model.addConstr(lhs, GRB_EQUAL, rhs);
            IGRAPH_VIT_NEXT(vit);
        }

        igraph_vs_destroy(&vs);
        igraph_vit_destroy(&vit);

        // add binding constraints and delta constraint
        igraph_es_t es;
        igraph_eit_t eit;
        igraph_es_all(&es, IGRAPH_EDGEORDER_ID);
        igraph_eit_create(&alignment_graph, es, &eit);

        GRBLinExpr delta_constraint_rhs = GRBLinExpr(delta);
        GRBLinExpr delta_constraint_lhs = GRBLinExpr();

        while(!IGRAPH_EIT_END(eit)){
            int eid = IGRAPH_EIT_GET(eit);
            int variant = igraph_cattribute_EAN(&alignment_graph, "variant", eid);
            if(variant >= 0){
                int idx = global_variable_to_local_idx[variant];

                GRBLinExpr lhs = GRBLinExpr();
                lhs += global_var[idx];
                lhs += local_var[eid];
                GRBLinExpr rhs = GRBLinExpr(1);
                model.addConstr(lhs, GRB_LESS_EQUAL, rhs);
            }

            delta_constraint_lhs += igraph_cattribute_EAN(&alignment_graph, "weight", eid)*local_var[eid];


            // added for feasibility check of initsol, delete afterward
            /*
            GRBLinExpr lhs = GRBLinExpr();
            lhs += local_var[eid];
            GRBLinExpr rhs = GRBLinExpr((int)igraph_cattribute_EAN(&alignment_graph, "initsol", eid));
            model.addConstr(lhs, GRB_EQUAL, rhs);
            */

            local_var[eid].set(GRB_DoubleAttr_Start , (int)igraph_cattribute_EAN(&alignment_graph, "initsol", eid));


            IGRAPH_EIT_NEXT(eit);
        }

        model.addConstr(delta_constraint_lhs, GRB_LESS_EQUAL, delta_constraint_rhs);

        igraph_es_destroy(&es);
        igraph_eit_destroy(&eit);

        /*
        if(VERBOSE){
            auto stop = chrono::steady_clock::now();
            auto duration = duration_cast<chrono::milliseconds>(stop - start);
            cout << "\tTime for alignment_graph_" << alignment_graph_idx << " ILP construction: " << duration.count() << " milliseconds" 
                " component (model): " << component << endl;
        }
        */
        igraph_destroy(&alignment_graph);
    }

    // update and view models for checking
    model.update();
    if(WRITE_ILPS_TO_FILE){
        string file_name = "model_" + to_string(component) + ".lp";
        model.write(file_name);
    }
    
}


int main(int argc, char** argv){

    GRBEnv env = GRBEnv();
    env.set("OutputFlag", "0");
    env.start();

    if(argc < 8){
        cout << "Wrong number of arguments provided" << endl;
        return 0;
    }
    string edge_file_name = argv[1];
    string pos_substring_file_name = argv[2];
    string scratch_directory = argv[3];
    string sol_directory = argv[4];
    int alpha = atoi(argv[5]);
    int delta = atoi(argv[6]);
    int use_files = atoi(argv[7]);

    igraph_set_attribute_table(&igraph_cattribute_table);

    // construct variation graph
    igraph_t variation_graph;
    int num_variants;

    read_edge_file(edge_file_name, &variation_graph, &num_variants);

    cout << "\n============= Variation Graph Constructed =============\n" << endl;
    //cout << "Construction time: " << duration.count() << " milliseconds\n" <<endl;
    cout << "Number of variants: " << num_variants << endl;
    cout << "Number of vertices: " << igraph_vcount(&variation_graph) << endl;
    cout << "Number of edges: " << igraph_ecount(&variation_graph) << endl;

    cout << "\nReading pos_substring file...\n" << endl;
    // read locations and strings
    vector<int> positions;
    vector<string> substrings;

    read_pos_substring_file(pos_substring_file_name, positions, substrings);
    
    int N = positions.size();
    
    cout << "\nConstructing alignment graphs...\n" << endl;
    vector<vector<int>*> alignment_graph_to_variants = vector<vector<int>*>(N);


    vector<igraph_t*> alignment_graphs;
    for(int i = 0; i < N; i++){ 
        alignment_graphs.push_back(new igraph_t); 
    }

    #pragma omp parallel for
    for(int i = 0; i < N; i++){     
        if(VERBOSE){
            cout << "Constructing alignment graphs for position: " << positions[i] << endl;
        }else if (!VERBOSE && i % PRINT_FACTOR == 0){
            cout << "Constructing alignment graph: " << i << " (prints mod " << PRINT_FACTOR << ")" << endl;
        }

        igraph_t alignment_graph;
        create_alignment_graph(&variation_graph, positions[i], substrings[i], alpha, delta, &alignment_graph);

        // obtain set of variants
        alignment_graph_to_variants[i] = new vector<int>;
        unordered_set<int> variants;
        get_variants(&alignment_graph, variants);

        for(int v: variants){
            alignment_graph_to_variants[i]->push_back(v);
        }

        igraph_cattribute_GAN_set(&alignment_graph, "position", positions[i]);
        igraph_cattribute_GAS_set(&alignment_graph, "substring", substrings[i].c_str());
        
        if(use_files == 1){
            string file_name = scratch_directory + "g_" + to_string(i);        
            write_alignment_graph_to_file(&alignment_graph, file_name);
        }else{
            igraph_copy(alignment_graphs[i], &alignment_graph);
        }
        
        
        //sum_of_sizes += igraph_ecount(&alignment_graph);

        igraph_destroy(&alignment_graph);
    }

    //cout << "average alignment graph size: " << (float)sum_of_sizes/(float)N << endl;
    
    
    positions.clear();
    positions.shrink_to_fit();
    substrings.clear();
    substrings.shrink_to_fit();

    igraph_destroy(&variation_graph);
    cout << "\n============= Alignment Graphs Constructed =============\n" << endl;

    cout << "Analyzing alignment graph set..." << endl;

    vector<int> global_variable_to_component;
    vector<int> subILP_to_component;
    vector<int> global_variable_to_local_idx;
    vector<int> component_to_variant_count;
    vector<vector<int>> component_to_subILPs;
    int num_components;
    analyze_alignment_graph_set(num_variants,
                                alignment_graph_to_variants, 
                                global_variable_to_component, 
                                subILP_to_component, 
                                global_variable_to_local_idx, 
                                component_to_variant_count,
                                component_to_subILPs,
                                &num_components);

    alignment_graph_to_variants.clear();
    alignment_graph_to_variants.shrink_to_fit();

    cout << "\n============= Alignment Graphs Analyzed =============\n" << endl;
    cout << "number of components: " << num_components << endl;

    // we need to invert the map from separate ILPs
    map<pair<int, int>, int> component_index_pair_to_variant;
    for(int i = 0; i < num_variants; i++){
        
        auto p = make_pair(global_variable_to_component[i], global_variable_to_local_idx[i]);

        if(component_index_pair_to_variant.find(p) == component_index_pair_to_variant.end()){
            component_index_pair_to_variant[p] = i;
        }
    }

    
    cout << "\nConstructing ILP models from alignment graphs...\n" << endl;
    
    vector<int> solution = vector<int>(num_variants, 1);

    bool some_ILP_timed_out = false;

    #pragma omp parallel for
    for(int component = 0; component < num_components; component++){
        
        if(component % PRINT_FACTOR == 0){
            cout << "Solving ILP_" << component << " out of " <<  num_components << " (prints mod " << PRINT_FACTOR << ")" <<endl;
        }
        
        //cout << "Solving ILP_" << component << " out of " <<  num_components << endl;


        //cout << "component: " << component << endl; 
        if(component_to_subILPs[component].size() > 0){

            GRBModel model = GRBModel(env);

            construct_ILPs(component,
                            component_to_subILPs,
                            global_variable_to_local_idx, 
                            component_to_variant_count,
                            model,
                            delta,
                            scratch_directory,
                            alignment_graphs,
                            use_files
                            );
            //cout << "\n============= ILP model "<< component << " constructed =============\n" << endl;
            
            
            model.getEnv().set(GRB_DoubleParam_TimeLimit, ILP_TIME_LIMIT);
            auto start = chrono::steady_clock::now();
            model.optimize();
            auto stop = chrono::steady_clock::now();
            auto duration = duration_cast<chrono::seconds>(stop - start);
            
            if(duration.count() > .95 * ILP_TIME_LIMIT){
                cout << "ILP was timed out" << endl;
                some_ILP_timed_out = true;
            }
            //if(VERBOSE){
            //    cout << "time: " << duration.count() << " milliseconds\n" <<endl;
            //}
            
            //cout << "\n============= ILP models solved =============\n" << endl;

            GRBVar* vars = NULL;
            vars = model.getVars();
            //cout << "ILP_" << component << " local assignments: " << endl;
            //cout << "variable count: " << component_to_variant_count[component] << endl;
            for (int j = 0; j < component_to_variant_count[component]; j++) {
                //cout << vars[j].get(GRB_StringAttr_VarName) << " = " << vars[j].get(GRB_DoubleAttr_X) << endl;
                solution[component_index_pair_to_variant[make_pair(component, j)]] = vars[j].get(GRB_DoubleAttr_X);
            }
        }

    }
        
    int sum = 0;
    cout << "\nFinal assignment: " << endl;


    string sol_file_name = sol_directory + "ILP_sol_" + to_string(alpha) + "_" + to_string(delta) + "_" + to_string(N);
    ofstream out_file(sol_file_name);
    for(int i = 0; i < solution.size(); i++){
        cout << "x_" << i << " = " << solution[i] << endl;
        out_file << solution[i] << endl;
        sum += solution[i];
    }
    out_file.close();
    cout << "\nNumber of variants deleted: " << sum << " out of " << num_variants <<endl;
    if(some_ILP_timed_out){
        cout << "Some ILP was timed-out, solution may be sub-optimal" << endl;
    }else{
        cout << "No ILPs were timed-out, solution is optimal" << endl;
    }
    

}
