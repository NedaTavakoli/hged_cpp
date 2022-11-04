#include <boost/graph/filtered_graph.hpp>

template <typename Graph, typename V = typename Graph::vertex_descriptor>
using Filtered = boost::filtered_graph<Graph, boost::keep_all, std::function<bool(V)>>;

template <typename Graph, typename Nodes, typename V = typename Graph::vertex_descriptor>
Filtered<Graph> induced_subgraph(Graph& g, Nodes nodes) {
    auto n = std::make_shared<Nodes>(std::move(nodes));
    return {g, boost::keep_all{}, [n](V v) { return n->count(v); }};
}

#include <boost/graph/breadth_first_search.hpp>
template <typename Graph, typename V = typename Graph::vertex_descriptor>
Filtered<Graph> reachable_subgraph(Graph& g, V source, unsigned depth) {
    std::vector<unsigned> dist(num_vertices(g));

    auto vis = make_bfs_visitor(record_distances(dist.data(), boost::on_tree_edge{}));
    breadth_first_search(g, source, visitor(vis));

    std::set<V> N;
    typename Graph::vertex_iterator b, e;
    for (boost::tie(b, e) = vertices(g); b != e; ++b) {
        auto d = dist[*b];
        if (d && d < depth)
            N.insert(*b);
    }

    return induced_subgraph(g, N);
}

// DEMONSTRATION
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/random.hpp>
#include <fmt/ranges.h>
#include <random>

int main() {
    boost::adjacency_list<> g;
    std::mt19937 prng(77); // manually selected for example
    generate_random_graph(g, 100, 150, prng);

    print_graph(g);
    fmt::print("g has vertices [0..{}]\n", num_vertices(g) - 1);

    auto sub = reachable_subgraph(g, 5, 3);
    fmt::print("nodes reachable < 3 from 5: {}\n", boost::make_iterator_range(vertices(sub)));

    print_graph(sub);
}
