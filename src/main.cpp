#include <iostream>
#include <vector>
#include <fstream>
#include <igraph.h>
#include <unordered_map>
#include <queue>

using namespace std;

typedef struct edge_attributes_struct{
    igraph_integer_t edge_idx;
    igraph_real_t label;
    igraph_real_t variant;
}edge_attributes;

void print_variation_graph(igraph_t *graph){
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


void print_alignment_graph(igraph_t *graph){
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


void read_edge_file(string edge_file_name, igraph_t* graph){
    
    ifstream infile(edge_file_name);

    int num_vertices;
    infile >> num_vertices;
    bool directed = true;
    igraph_empty(graph, num_vertices, directed);

    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);

    int start_vertex, end_vertex;
    char edge_label, variant_number;

    vector<edge_attributes> attributes_vector;
    int i = 0;
    while (infile >> start_vertex >> end_vertex >> edge_label >> variant_number)
    {
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

    igraph_add_edges(graph, &edges, NULL);

    for(edge_attributes a: attributes_vector){
        igraph_cattribute_EAN_set(graph, "label", a.edge_idx, a.label);
        igraph_cattribute_EAN_set(graph, "variant", a.edge_idx, a.variant);
    }

    print_variation_graph(graph);
    cout << "============= Graph Constructed =============" << endl;
    
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



/*
The main idea is that we consider the vertices of the alignment graph as laided out in |P||V| grid.
Each row of length |V| is sorted in the same topological order.
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

    //  hash to track seen coordinates
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
        

        int u_coord = u.y * n + inv_sorted_vertex[u.id];
        vertex u_ = seen_coordinates[u_coord];
        int u_dist = u_.dist;

        // vertices on same level (insertion)
        if(u_dist + 1 <= delta){

            igraph_incident(graph, &incident_edges, old_uid, IGRAPH_OUT);
            int m = igraph_vector_int_size(&incident_edges);

            for(int i = 0; i < m; i++){

                igraph_integer_t start, end;
                igraph_edge(graph, igraph_vector_int_get(&incident_edges, i), &start, &end);

                int x = inv_sorted_vertex[end];
                int y = u.y;
                int v_coordinate = y*n + x;
                int vid;

                if(seen_coordinates.find(v_coordinate) == seen_coordinates.end()){
                    // never encounter coordinate within reachable distance
                    // add vertex with new vertex id and edge

                    vid = new_vid;
                    seen_coordinates[v_coordinate] = {vid, x, y, u_dist + 1};
                    new_vid++;
                    
                    
                    // enque vertex
                    q.push({vid, x, y, -1});

                }else{
                    // coordinate seen before 
                    vertex v = seen_coordinates[v_coordinate];
                    vid = v.id;

                    if(u_dist + 1 < v.dist){
                        seen_coordinates[v_coordinate] = {vid, x, y, u_dist + 1};
                    }
                }
                igraph_vector_int_push_back(&alignment_graph_edges, u.id);
                igraph_vector_int_push_back(&alignment_graph_edges, vid);
                weights.push_back(1);
                variants.push_back(igraph_cattribute_EAN(graph, "variant", igraph_vector_int_get(&incident_edges, i)));
            }
        }

        // deletion edges

        // sub/match edges
    }

    bool directed = true;
    igraph_empty(result, new_vid, directed);
    igraph_add_edges(result, &alignment_graph_edges, NULL);

    for(int i = 0; i < weights.size(); i++){
        igraph_cattribute_EAN_set(result, "weight", i, weights.at(i));
        igraph_cattribute_EAN_set(result, "variant", i, variants.at(i));
    }

    igraph_vector_int_destroy(&alignment_graph_edges);
    igraph_vector_int_destroy(&sorted_vertex);

    print_alignment_graph(result);
} 


int main(){

    igraph_set_attribute_table(&igraph_cattribute_table);

    // construct variation graph
    string edge_file_name = "edge_file.txt";
    igraph_t variation_graph;
    read_edge_file(edge_file_name, &variation_graph);

    // read locations and strings
    string pos_substring_file_name = "loc_substring.txt";
    vector<int> positions;
    vector<string> substrings;
    read_pos_substring_file(pos_substring_file_name, positions, substrings);

    // how to get reachable subgraph for each start vertex
    int dist = 5;
    int start_vertex = 2;
    igraph_vs_t vertex_selector;
    igraph_vector_int_list_t reachable_list;

    igraph_vs_1(&vertex_selector, start_vertex);
    igraph_vector_int_list_init(&reachable_list, 0);
    
    get_reachable_vertex(&variation_graph, vertex_selector, dist, &reachable_list);
    igraph_vs_destroy(&vertex_selector);

    // assuming only one vertex is specified by vertex selector, we get first vector
    igraph_vector_int_t* reachable_vertex = igraph_vector_int_list_get_ptr(&reachable_list, 0);

    igraph_t reachable_graph;
    get_induced_variation_graph(&variation_graph, reachable_vertex, &reachable_graph);

    igraph_t alignment_graph;
    string P= "TCA";
    int delta = 1;
    create_pruned_alignment_graph(&reachable_graph, P, delta, &alignment_graph);

    /*    
    for(long int i = 0; i < igraph_vector_int_size(reachable); i++){
        cout << igraph_vector_int_get(reachable, i) << ", ";
    }
    cout << endl;
    */
}