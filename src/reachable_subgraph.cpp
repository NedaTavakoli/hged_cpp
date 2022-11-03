//https://stackoverflow.com/questions/59397776/iterating-over-boost-graph-vertexes-based-on-maximal-depth
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <iostream>
#include <map>

using boost::make_iterator_range;
using Name   = char;
using Graph  = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, Name>;
using Edge   = Graph::edge_descriptor;
using Vertex = Graph::vertex_descriptor;

using Depths = std::map<Vertex, size_t>;

struct Recorder : public boost::base_visitor<Recorder> {
    using event_filter = boost::on_tree_edge;

    Depths& _ref;
    Recorder(Depths& r):_ref(r){}

    void operator()(Edge e, const Graph& g) const {
        auto depth = 1 + _ref[source(e, g)];

        if (auto [it, isnew] = _ref.emplace(target(e, g), depth); !isnew) {
            it->second = std::max(it->second, depth);
        }
    }
};

int main() {
    enum : Vertex { A, B, C, D, E, F, G, H, I, COUNT };
    Graph g(COUNT);

    // give names
    for (auto v : make_iterator_range(vertices(g)))
        g[v] = 'A' + v;

    // add edges
    for (auto [s,t] : {
            std::pair(D,G), {D, G}, {C, G}, {C, F}, {B, F}, {B, E}, {A, E},
            {G, H}, {F, I}, {E, H},
            {H, I} })
    {
        add_edge(s, t, g);
    }

    boost::queue<Vertex> queue;
    std::vector<boost::default_color_type> colors(num_vertices(g));
    auto color_map = boost::make_iterator_property_map(begin(colors), get(boost::vertex_index, g));

    Depths depths;
    Recorder r{depths};

    auto all = vertices(g);
    boost::breadth_first_search(g, all.first, all.second, queue, boost::make_bfs_visitor(r), color_map);

    for (auto v : make_iterator_range(vertices(g))) {
        boost::breadth_first_search(g, v, queue, boost::make_bfs_visitor(r), color_map);
    }

    std::map<size_t, std::set<Vertex> > by_depth;
    for (auto [v,d] : depths)
        by_depth[d].insert(v);

    for (auto& [d,vs] : by_depth) {
        std::cout << "depth:" << d;
        for (auto v : vs) std::cout << " " << g[v];
        std::cout << "\n";
    }
}
