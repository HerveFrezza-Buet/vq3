#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>


// This example shows the use of the dijkstra algorithm.


// Graph definition
//
///////////////

// Vertex values have to be instrumented with shortest path
// computation structure.

//                                              ## Node properties :
using layer_0 = prototype;                      // prototypes are 2D points (this is the "user defined" value).
using layer_1 = vq3::shortest_path<layer_0>;    // This holds informations built by dijkstra.
using vertex  = layer_1;

// Here, each edge hosts its cost value.  It is the optional value
// (not mandatory) so that it is not recomputed once it has been
// calculated first.

//                                              ## Node properties :
using edge   = vq3::decorator::tagged<void>;    // We need tags on the edges for CHL


using graph  = vq3::graph<vertex, edge>; 

// Distance
//
////////////////


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v, p);}

