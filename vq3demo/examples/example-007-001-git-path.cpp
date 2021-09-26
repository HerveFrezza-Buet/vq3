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

//                                                        ## Edge properties :
using elayer_0 = vq3::decorator::optional_cost<void>;     // Edge cost.
using edge     = elayer_0;

using graph = vq3::graph<vertex, edge>;


// Let us define some functions in order to set up dijkstra
// computation.

// This is the cost of an edge, constant here.
double edge_cost(const graph::ref_edge& ref_e) {
  auto& opt_cost = (*ref_e)().vq3_cost; 
  if(!opt_cost) { // if cost not calculated yet.
    auto extr_pair = ref_e->extremities();
    const auto& pt1 = (*(extr_pair.first))().vq3_value;
    const auto& pt2 = (*(extr_pair.second))().vq3_value;
    opt_cost = demo2d::d(pt1, pt2);
  }
  return *opt_cost;
}

// This is the distance between a node value and a sample.
double d(const sample& p1, const sample& p2) {return demo2d::d(p1, p2);}

// This is the distance between a node value and a sample.
double D(const vertex& v, const sample& p) {return demo2d::d(v.vq3_value, p);}

// This is the linear interpolation between samples, when no graph is used.
sample interpolate(const sample& a, const sample& b, double lambda) {return (1-lambda)*a + lambda*b;}

// This is the  way a shortest path is computed in our graph.
void shortest_path(graph& g, graph::ref_vertex start, graph::ref_vertex dest) {
  vq3::path::a_star<false, false>(g, start, dest, edge_cost,
				  [start](const graph::ref_vertex& ref_v){ // This is the heuristic
				    return demo2d::d((*start)().vq3_value, (*ref_v)().vq3_value); 
				  });
}

// We now are ready to start actual computation....

// ... but first, let us add cv stuff to be able to place our 2
// samples.


void on_mouse( int event, int x, int y, int, void* user_data);
struct callback_data {
  demo2d::Point start, end;

private:
  demo2d::opencv::Frame& frame;
  bool start_changed;
  bool end_changed;

  friend void on_mouse( int event, int x, int y, int, void* user_data);
  
public:
  callback_data(demo2d::opencv::Frame& frame)
    : start({-.5,  .1}),
      end({-.5, -.1}),
      frame(frame),
      start_changed(false), end_changed(false) {}
  callback_data()                                = delete;
  callback_data(const callback_data&)            = delete;
  callback_data& operator=(const callback_data&) = delete;
  void set_start(const demo2d::Point& pt) {start_changed = true; start = pt;}
  void set_end(const demo2d::Point& pt)   {end_changed   = true; end   = pt;}
  bool start_has_changed()                {auto res = start_changed; start_changed = false; return res;}
  bool end_has_changed()                  {auto res =   end_changed; end_changed   = false; return res;}
};


// Ok, now we start.

int main(int argc, char* argv[]) {
  
  graph g;

  // Let us build a very simple auxiliary graph.
  auto V1 = g += demo2d::Point(-.5, -.5);
  auto V2 = g += demo2d::Point( .5, -.5);
  auto V3 = g += demo2d::Point( .5,  .5);
  auto V4 = g += demo2d::Point(-.5,  .5);

  g.connect(V1, V2);
  g.connect(V2, V3);
  g.connect(V3, V4);


  // Let us set up a display.
  
  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  auto image = cv::Mat(480, 640, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = demo2d::opencv::direct_orthonormal_frame(image.size(), .8*image.size().height, true);
  int slider = 500;
  cv::createTrackbar("walk ratio", "image", &slider, 1000, nullptr);
  callback_data udata(frame);
  cv::setMouseCallback("image", on_mouse, reinterpret_cast<void*>(&udata));

  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
                                                                           [](const vertex& v) {return                      true;},  // always draw
                                                                           [](const vertex& v) {return               v.vq3_value;},  // position
                                                                           [](const vertex& v) {return                         6;},  // radius
                                                                           [](const vertex& v) {return cv::Scalar(200,  80,  80);},  // color
                                                                           [](const vertex& v) {return                        -1;}); // thickness

  auto draw_edge = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
  								     [](const vertex&, const vertex&, const edge&) {return         true;},  // always draw
  								     [](const vertex& v) {return               v.vq3_value;},  // position
  								     [](const edge&)     {return cv::Scalar(255, 200, 200);},  // color
  								     [](const edge&)     {return                         1;}); // thickness

  
  // Let us consider two samples, udata.start and udata.end
  // (they can be changed with the mouse).

  // They are usual Euclidian points, distance and linear
  // interpolations between them are linear, ignoring the auxiliary
  // graph.

  // The idea is to set up "mirror" points, which are the same points,
  // but living in a world where the distances are graph-related (or
  // graph-induced).
  
  auto traits = vq3::topology::gi::traits<sample>(g, d, D, interpolate, shortest_path);   // This gathers all the required definitions.
  auto X      = vq3::topology::gi::value(traits, udata.start);                            // X is start, but topology is no more Eucldian, it is graph-induced.
  auto Y      = vq3::topology::gi::value(traits, udata.end);                              // The same for Y and end.

  // start == X(), end == Y()   The () operator allows to retrieve the value.

  // This is the loop.
  
  std::cout << std::endl
            << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl
	    << "left  click sets the starting point." << std::endl
	    << "right click sets the ending point." << std::endl
            << "Press ESC to quit." << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl;
  
  int keycode = 0;
  while(keycode != 27) {

    if(udata.start_has_changed()) X = udata.start; // This re-initialize X/Y from another value...
    if(udata.end_has_changed())   Y = udata.end;   // ... this triggers computation of closest node.
    
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
    cv::circle(image, frame(X()), 3, cv::Scalar(0, 200, 0), -1);
    cv::circle(image, frame(Y()), 3, cv::Scalar(0, 0, 200), -1);
      
    cv::imshow("image", image);
    keycode = cv::waitKey(50) & 0xFF;
  }

  return 0;
}


void on_mouse( int event, int x, int y, int, void* udata) {
  if(event != cv::EVENT_LBUTTONDOWN && event != cv::EVENT_RBUTTONDOWN)
    return;
  auto& data = *(reinterpret_cast<callback_data*>(udata));
  if(event == cv::EVENT_LBUTTONDOWN)
    data.set_start(data.frame(cv::Point(x, y)));
  else
    data.set_end(data.frame(cv::Point(x, y)));
}
