#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include </usr/local/include/igraph/igraph.h>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <chrono>
#include <climits>

using namespace std;

typedef struct edge_attributes_struct{
    igraph_integer_t edge_idx;
    igraph_real_t label;
}edge_attributes;

void print_variation_graph(const igraph_t* graph){
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
                << "label: " << igraph_cattribute_EAN(graph, "label", igraph_vector_int_get(&edges, i))
                << ")" << endl;   

        }
        igraph_vector_int_destroy(&edges);
        IGRAPH_VIT_NEXT(vit);
        cout << endl;
    }
    cout<< endl;
}


void print_alignment_graph(igraph_t *graph){

    cout << "\nPrinting alignment graph:" << endl;
    //cout << "pos = " << igraph_cattribute_GAN(graph, "position")
    //    << ", substring = " << igraph_cattribute_GAS(graph, "substring") << endl;

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
                //<< "variant: " << igraph_cattribute_EAN(graph, "variant" , igraph_vector_int_get(&edges, i)) << ", "
                //<< "initsol: " << igraph_cattribute_EAN(graph, "initsol" , igraph_vector_int_get(&edges, i)) 
                << ")" << endl; 
        }
        IGRAPH_VIT_NEXT(vit);
        cout << endl;
    }
    cout<< endl;
}


void read_edge_file(string edge_file_name, igraph_t* graph){
   
    ifstream infile(edge_file_name);

    int num_vertices = -1;

    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);

    int start_vertex, end_vertex;
    char edge_label;
    string variant_number;

    vector<edge_attributes> attributes_vector;
    int i = 0;
    while (infile >> start_vertex >> end_vertex >> edge_label)
    {
        cout << start_vertex << " " << end_vertex << " " << edge_label << endl;
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

        attributes_vector.push_back({i, label});

        //add edge
        igraph_vector_int_push_back(&edges, start_vertex);
        igraph_vector_int_push_back(&edges, end_vertex);
        i++;
    }

    //cout << num_vertices << endl;
    bool directed = true;
    igraph_empty(graph, num_vertices, directed);
    igraph_add_edges(graph, &edges, NULL);
    igraph_vector_int_destroy(&edges);

    for(edge_attributes a: attributes_vector){
        cout << a.edge_idx << " " << a.label << endl;
        igraph_cattribute_EAN_set(graph, "label", a.edge_idx, a.label);
    }

}


void read_read_file(string read_file_name, vector<string>& reads){
    ifstream infile(read_file_name);
    string line;

    while (std::getline(infile, line)){

        string substring;
        
        istringstream iss(line); // separates line over spaces
        //iss >> pos;

       if (line.size() > 0){
            while (iss >> substring){
                //positions.push_back(pos);
                reads.push_back(substring); 
            }
        }
    }
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

    igraph_vector_int_init(&alignment_graph_edges, 0);

    int source_vertex_id = 0;
    int sink_vertex_id = 0;

    int new_vid = 1;
    igraph_vector_int_t incident_edges;
    igraph_vector_int_init(&incident_edges, 0);


    queue<vertex> q;
    vertex s = {0, -1, 0, 0};


    // edges from source
    for(int v_x = 0; v_x < n; v_x++){

        int v_y = 0;
        int v_coordinate = v_y * n + v_x;

        int vid = new_vid;
        seen_coordinates[v_coordinate] = {vid, v_x, v_y, 0};
        new_vid++;

        igraph_vector_int_push_back(&alignment_graph_edges, 0);
        igraph_vector_int_push_back(&alignment_graph_edges, vid);
        weights.push_back(0);

        q.push({vid, v_x, v_y, 0});
    }


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

        }
    }
    // construct temporary forward-pruned only alignment graph
    bool directed = true;
    igraph_t forward_pruned_only_graph;
    igraph_empty(&forward_pruned_only_graph, new_vid, directed);
    igraph_add_edges(&forward_pruned_only_graph, &alignment_graph_edges, NULL);

    for(int i = 0; i < weights.size(); i++){
        igraph_cattribute_EAN_set(&forward_pruned_only_graph, "weight", i, weights.at(i));
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

int main(int argc, char** argv){
    cout << "hi" << endl;

    /*
    GRBEnv env = GRBEnv();
    env.set("OutputFlag", "0");
    env.start();
    */

    if(argc < 3){
        cout << "Wrong number of arguments provided" << endl;
        return 0;
    }
    string graph_file_name = argv[1];
    string read_file_name = argv[2];

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_t variation_graph;

    read_edge_file(graph_file_name, &variation_graph);

    print_variation_graph(&variation_graph);

    vector<string> reads;
    read_read_file(read_file_name, reads);

    /*
    for(int  i = 0; i < reads.size(); i++){
        cout << reads[i] << endl;
    }
    */

    int N = reads.size();

    vector<igraph_t*> alignment_graphs;
    for(int i = 0; i < N; i++){ 
        alignment_graphs.push_back(new igraph_t); 
    }

    int delta = 0;

    for(int i = 0; i < N; i++){  
        igraph_t alignment_graph;
        create_pruned_alignment_graph(&variation_graph, reads[i], delta, &alignment_graph);
        print_alignment_graph(&alignment_graph);
        igraph_copy(alignment_graphs[i], &alignment_graph);
        igraph_destroy(&alignment_graph);
    }

    return 0;
}

    //alignment_graph_to_variants.clear();
    //alignment_graph_to_variants.shrink_to_fit();

    //cout << "\n============= Alignment Graphs Analyzed =============\n" << endl;
    //cout << "number of components: " << num_components << endl;

    /*we need to invert the map from separate ILPs
    map<pair<int, int>, int> component_index_pair_to_variant;
    for(int i = 0; i < num_variants; i++){
        
        auto p = make_pair(global_variable_to_component[i], global_variable_to_local_idx[i]);

        if(component_index_pair_to_variant.find(p) == component_index_pair_to_variant.end()){
            component_index_pair_to_variant[p] = i;
        }
    }*/
