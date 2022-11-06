#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/filtered_graph.hpp>
template <typename Graph, typename V = typename Graph::vertex_descriptor>
auto reachable_subgraph(Graph& g, V source, unsigned depth) {
    std::vector<unsigned>                  dist(num_vertices(g));
    std::vector<boost::default_color_type> colors(num_vertices(g));

    auto stop_when = [&dist, depth](auto v, auto const&) { return dist.at(v) >= depth; };

    depth_first_visit(g, source,
                      make_dfs_visitor(record_distances(dist.data(), boost::on_tree_edge())),
                      colors.data(), stop_when);

    return boost::filtered_graph( //
        g, boost::keep_all{}, std::function{[depth, dist = std::move(dist)](V v) {
            auto d = dist.at(v);
            return d > 0 && d < depth;
        }});
}

// DEMONSTRATION
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/random.hpp>
#include <fmt/ranges.h>
#include <random>

int main() {
    boost::adjacency_list g;
    std::mt19937 prng(77); // manually selected for example
    generate_random_graph(g, 100, 150, prng);

    // print_graph(g);
    fmt::print("g has vertices [0..{}]\n", num_vertices(g) - 1);

    auto sub = reachable_subgraph(g, 5, 3);
    //fmt::print("nodes reachable < 3 from 5: {}\n", boost::make_iterator_range(vertices(sub)));

    print_graph(sub);
}