#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>


// Edge weight. is a labled string
typedef boost::property<boost::edge_weight_t, std::string> EdgeWeightProperty;

typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::no_property, EdgeWeightProperty> G;
typedef boost::property<boost::edge_weight_t, int> EdgeWeightProperty;

using V = G::vertex_descriptor;
using E = G::edge_descriptor;

int main() {
    G g;


    V v1 = add_vertex(g);
    V v2 = add_vertex(g);
    E e1 = boost::add_edge(v1, v2, EdgeWeightProperty(2), g);


        // The property map associated with the weights.
        boost::property_map < G,
                      boost::edge_weight_t >::type EdgeWeightMap = get(boost::edge_weight, g);


    print_graph(g);
}