#include <opencv2/opencv.hpp>
#include <vq3.hpp>
#include <vq3demo.hpp>

/*
  git means "Graph Induced Topology". It constists in using an
  auxiliary graph to compute distances between elements and perform
  linear combination of elements. The auxiliary graph gives the
  topological structure of the space the data samples are living in.
*/


// Let us consider 2D points.

using sample = demo2d::Point;

// Let us define the auxiliary graph.



//                                                         ## Node properties :
using vlayer_0 = demo2d::Point;                            // The graph nodes are points.
using vlayer_1 = vq3::decorator::path::shortest<vlayer_0>; // This holds informations built by dijkstra.
using vertex   = vlayer_1;
  


using graph = vq3::graph<vertex, void>;


// Let us define some functions in order to set up dijkstra
// computation.

// This is the cost of an edge, constant here.
double edge_cost(const graph::ref_edge& ref_e) {return 1.;}

// This is the distance between a node value and a sample.
double d2(const vertex& v, const sample& p) {return demo2d::d2(v.vq3_value, p);}

// This is the linear interpolation between samples, when no graph is used.
sample interpolate(const sample& a, const sample& b, double lambda) {return (1-lambda)*a + lambda*b;}

// This is the  way a shortest path is computed in our graph.
void shortest_path(graph& g, graph::ref_vertex start, graph::ref_vertex dest) {
  vq3::path::a_star<false, false>(g, start, dest, edge_cost,
				  [start](const graph::ref_vertex& ref_v){ // This is the heuristic
				    return demo2d::d((*start)().vq3_value, (*ref_v)().vq3_value); 
				  });
}

// We now are ready to start actual computation.




int main(int argc, char* argv[]) {
  
  graph g;

  // Let us build a very simple auxiliary graph.
  auto A = g += demo2d::Point(-.5, -.5);
  auto B = g += demo2d::Point( .5, -.5);
  auto C = g += demo2d::Point( .5,  .5);
  auto D = g += demo2d::Point(-.5,  .5);

  g.connect(A, B);
  g.connect(B, C);
  g.connect(C, D);

  // Let us consider two samples.
  auto x = demo2d::Point(-.5,  .1);
  auto y = demo2d::Point(-.5, -.1);

  // They are usual Eucidian points, distance and linear
  // interpolations between them are linear, ignoring the auxiliary
  // graph.


  // The idea is to set up "mirror" points, which are the same points,
  // but living in a world where the distances are graph-related (or
  // graph-induced).
  
  auto traits = vq3::topology::gi::traits<sample>(g, d2, interpolate, shortest_path); // This gathers all the required definitions.
  auto X      = vq3::topology::gi::value(traits, x);                                  // X is x, but topology is no more Eucldian, it is graph-induced.
  auto Y      = vq3::topology::gi::value(traits, y);                                  // The same for Y and y.

  // x == X(), y == Y()   The () operator allows to retrieve the value.

  // Let us set up a display.
  
  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  auto image = cv::Mat(480, 640, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = demo2d::opencv::direct_orthonormal_frame(image.size(), .8*image.size().height, true);
  int slider = 500;
  cv::createTrackbar("walk ratio", "image", &slider, 1000, nullptr);

  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
                                                                           [](const vertex& v) {return                      true;},  // always draw
                                                                           [](const vertex& v) {return               v.vq3_value;},  // position
                                                                           [](const vertex& v) {return                         6;},  // radius
                                                                           [](const vertex& v) {return cv::Scalar(200,  80,  80);},  // color
                                                                           [](const vertex& v) {return                        -1;}); // thickness

  auto draw_edge = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
  								     [](const vertex&, const vertex&) {return         true;},  // always draw
  								     [](const vertex& v) {return               v.vq3_value;},  // position
  								     []()                {return cv::Scalar(255, 200, 200);},  // color
  								     []()                {return                         1;}); // thickness
                    

  // This is the loop.
  
  std::cout << std::endl
            << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl
            << "Press ESC to quit." << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl;
  
  int keycode = 0;
  while(keycode != 27) {
    
    auto T = X;                    // T is another value with graph-induced topology.
    double lambda = slider*.001;   // lambda (in [0, 1]) is an interpolation coefficient.
    T += lambda*(Y - T);           // T <- (1 - lambda) * T + lambda * Y

    // T is updated as in usual Kohonen rule. When lambda varies from
    // 0 to 1 (use the slider), T goes from X (initial T value) to Y,
    // gradually. As you can see on the figure, this gradual walk
    // follows the graph, rather than following a straight line.
    
    image = cv::Scalar(255, 255, 255);
    cv::circle(image, frame(T()), 11, cv::Scalar(0, 200, 200), -1); 
    
    g.foreach_edge(draw_edge);
    g.foreach_vertex(draw_vertex);
    cv::circle(image, frame(X()), 3, cv::Scalar(0, 0, 200), -1);
    cv::circle(image, frame(Y()), 3, cv::Scalar(0, 0, 200), -1);
      
    cv::imshow("image", image);
    keycode = cv::waitKey(50) & 0xFF;
  }

  return 0;
}

