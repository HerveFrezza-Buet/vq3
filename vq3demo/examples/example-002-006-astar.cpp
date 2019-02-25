#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>
#include <algorithm>
#include <vector>
#include <iterator>


// This example shows the use of the A* algorithm.


// Graph definition
//
///////////////

// Vertex values have to be instrumented with shortest path
// computation structure.
//                                                         ## Node properties :
using vlayer_0 = vq3::demo2d::Point;                       // prototypes are 2D points (this is the "user defined" value).
using vlayer_1 = vq3::decorator::path::shortest<vlayer_0>; // This holds informations built by dijkstra.
using vlayer_2 = vq3::decorator::tagged<vlayer_1>;         // We will tag extermities of the shortest path for display.
using vertex   = vlayer_2;

// Here, each edge hosts its cost value.  It is the optional value
// (not mandatory) so that it is not recomputed once it has been
// calculated first.
//                                                        ## Node properties :
using elayer_0 = vq3::decorator::efficiency<void>;        // We consider only efficient edges for paths.
using elayer_1 = vq3::decorator::optional_cost<elayer_0>; // Edge cost.
using elayer_2 = vq3::decorator::tagged<elayer_1>;        // We will tag edges belonging to a shortest path for display.
using edge     = elayer_2;


using graph   = vq3::graph<vertex, edge>;


// Distance
//
////////////////


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v.vq3_value, p);}

// We compute the edge cost as being its length. As we use an optional
// cost, we can check if the computation is necessary.
double edge_cost(const typename graph::ref_edge& ref_e) {
  auto& opt_cost = (*ref_e)().vq3_cost; 
  if(!opt_cost) { // if cost not calculated yet.
    auto extr_pair = ref_e->extremities();
    const auto& pt1 = (*(extr_pair.first))().vq3_value;
    const auto& pt2 = (*(extr_pair.second))().vq3_value;
    opt_cost = vq3::demo2d::d(pt1, pt2);
  }
  return *opt_cost;
}


// Callback
//
//////////////

struct callback_data {
  graph& g;
  std::mt19937& rd;
  vq3::demo2d::opencv::Frame& frame;
  graph::ref_vertex start, dest;
  callback_data(graph& g, std::mt19937& rd, vq3::demo2d::opencv::Frame& frame) : g(g), rd(rd), frame(frame) {}
  callback_data()                                = delete;
  callback_data(const callback_data&)            = delete;
  callback_data& operator=(const callback_data&) = delete;
};

void on_mouse( int event, int x, int y, int, void* user_data) {
  // A mouse click reinitializes the graph.
  if(event != cv::EVENT_LBUTTONDOWN )
    return;
  
  auto& data = *(reinterpret_cast<callback_data*>(user_data));
  
  vq3::utils::clear_all_tags(data.g, false);
  data.start = data.dest;
  data.dest  = vq3::utils::closest(data.g, data.frame(cv::Point(x,y)), d2);
  if(data.start == data.dest) data.start = nullptr;

  if(data.start) (*(data.start))().vq3_tag = true;
  if(data.dest)  (*(data.dest))().vq3_tag  = true;

  // false true : we do not consider vertex efficiency, but only edge efficiency.
  vq3::path::a_star<false,true>(data.g, data.start, data.dest, edge_cost,
				[start = data.start](const graph::ref_vertex& ref_v){ // This is the heuristic
				  if(start)
				    return vq3::demo2d::d((*start)().vq3_value, (*ref_v)().vq3_value); // direct distance is lower than real cost.
				  else
				    return 0.0;
				});

  // This draws the path by tagging the edges belonging to it.
  auto end = vq3::path::end(data.g);
  if(data.start == nullptr) {
    // Dijkstra has computed the full tree.
    data.g.foreach_vertex([end](const graph::ref_vertex ref_v) {
	for(auto it = vq3::path::begin(ref_v); it != end; ++it) {
	  // *it is a ref_vertex
	  if(auto ref_e = it.get_edge(); ref_e)
	    (*ref_e)().vq3_tag = true;
	}
      });
  }
  else {
    for(auto it = vq3::path::begin(data.start); it != end; ++it) {
      // *it is a ref_vertex
      if(auto ref_e = it.get_edge(); ref_e)
	(*ref_e)().vq3_tag = true;
    }
  }
  
}


// Main
//
//////////////

#define NB_VERTICES_PER_M2      200
#define NB_SAMPLES_PER_M2    100000
#define NB_UNEFFICIENT_EDGES    100

int main(int argc, char* argv[]) {

  std::random_device rd;  
  std::mt19937 random_device(rd());

  // Let us build the graph.
  //
  //////////////

  double intensity = 1.;
  double side      = 1.;
  auto   density   = vq3::demo2d::sample::rectangle(side, side, intensity);

  graph g;

  // We add vertices
  for(auto& location : vq3::demo2d::sample::sample_set(random_device, density, NB_VERTICES_PER_M2)) g += location;

  // We add edges
  for(auto& sample : vq3::demo2d::sample::sample_set(random_device, density, NB_SAMPLES_PER_M2)) {
    auto closest = vq3::utils::two_closest(g, sample, d2);
    if(g.get_edge(closest.first, closest.second) == nullptr) 
      g.connect(closest.first, closest.second);
  }

  // We handle edge efficiencies. Only NB_UNEFFICIENT_EDGES are unefficient.
  vq3::utils::clear_edge_efficiencies(g, true);
  std::vector<graph::ref_edge> edges;
  vq3::utils::collect_edges(g, std::back_inserter(edges));
  std::shuffle(edges.begin(), edges.end(), random_device);
  auto end = edges.begin() + NB_UNEFFICIENT_EDGES;
  for(auto it = edges.begin(); it != end; ++it)
    (*(*it))().vq3_efficient = false;

  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  auto image = cv::Mat(480, 640, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .9*image.size().height, true);
  callback_data user_data(g, random_device, frame);
  cv::setMouseCallback("image", on_mouse, reinterpret_cast<void*>(&user_data));
  
  auto draw_edge = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
  								     [](const vertex&, const vertex&, const edge&) {return   true;},  // always draw
  								     [](const vertex& v) {return                      v.vq3_value;},  // position
  								     [](const edge& e)   {
								       if(e.vq3_tag)            return cv::Scalar(000, 000, 000);
								       else if(e.vq3_efficient) return cv::Scalar(255, 100, 100);
								       else                     return cv::Scalar(200, 200, 200);},   // color
  								     [](const edge& e)   {
								       if(e.vq3_tag)            return 3;
								       else if(e.vq3_efficient) return 2;
								       else                     return 2;});                          // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return                  true;},  // always draw
									   [](const vertex& v) {return           v.vq3_value;},  // position
									   [](const vertex& v) {
									     if(v.vq3_tag) return 6;
									     else return          4;},                           // radius
									   [](const vertex& v) {
									     if(v.vq3_tag) return cv::Scalar(000, 000, 200);
									     else          return cv::Scalar(200,  80,  80);},   // color
									   [](const vertex& v) {return                   -1;});  // thickness
  
    
  std::cout << std::endl
            << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl
            << "Click on vertices to compute a path. If you click twice on the same" << std::endl
	    << "vertex, you will get a full Dijkstra computation." << std::endl
	    << "Press ESC to quit." << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl;
  int keycode = 0;
  while(keycode != 27) {    
    image = cv::Scalar(255, 255, 255);
    g.foreach_edge(draw_edge); 
    g.foreach_vertex(draw_vertex);
    cv::imshow("image", image);
    keycode = cv::waitKey(100) & 0xFF;
  }

  return 0;
}
