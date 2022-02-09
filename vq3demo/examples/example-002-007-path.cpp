#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>
#include <optional>
#include <algorithm>
#include <vector>
#include <iterator>
#include <stack>


// This example shows the use of path walking.


// Graph definition
//
///////////////

// Vertex values have to be instrumented with shortest path
// computation structure.
//                                                               ## Node properties :
using vlayer_0 = demo2d::Point;                                  // prototypes are 2D points (this is the "user defined" value).
using vlayer_1 = vq3::decorator::path::shortest<vlayer_0>;       // This holds informations built by dijkstra.
using vlayer_2 = vq3::decorator::path::length<vlayer_1, double>; // This holds accumulation along paths travels (needed if non-default behavior is implemented).
using vlayer_3 = vq3::decorator::tagged<vlayer_2>;               // We will tag extermities of the shortest path for display.
using vlayer_4 = vq3::decorator::efficiency<vlayer_3>;           // We consider only efficient edges for paths.
using vertex   = vlayer_4;

// Here, each edge hosts its cost value.  It is the optional 
// value here, so that it is not recomputed once it has been
// calculated first.
//                                                      ## Node properties :
using cost_type = std::optional<double>;                // The edge cost is a double, that may not be computed yet.  
using elayer_0 = vq3::decorator::cost<void, cost_type>; // Edge cost.
using elayer_1 = vq3::decorator::tagged<elayer_0>;      // We will tag edges belonging to a shortest path for display.
using edge     = elayer_1;

using graph   = vq3::graph<vertex, edge>;


// Distance
//
////////////////


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const demo2d::Point& p) {return demo2d::d2(v.vq3_value, p);}

// We compute the edge cost as being its length. As we use an optional
// cost, we can check if the computation is necessary.
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

// This is how we interpolate between 2 vertices.
auto vertex_interpolation(const graph::ref_vertex& start, const graph::ref_vertex& dest, double lambda) {
  if(start == nullptr) return (*dest )().vq3_value;
  if(dest  == nullptr) return (*start)().vq3_value;
  auto& start_value = (*start)().vq3_value;
  auto& dest_value  = (*dest )().vq3_value;
  return (1 - lambda) * start_value + lambda * dest_value;
}

// We will travel along some shortest paths. Default behavior is that
// the length of edges, consider in the travel interpolation, is the
// cost of the edges ised in shortest path implementation. In this
// case, non need to define anything more, and the
// vq3::decorator::path::length layer of our node could be
// removed. Let us implement a non default behavoir here. Let us
// consider that the length of each edge is 1. So traveling half of
// the edge (0.5) would lead to a position such as there is the same
// number of edges for reaching both path extremities, not the same
// path length. Two custom functions are needed.

void compute_cumulated_costs(graph::ref_vertex start, graph::ref_vertex dest) {
  // We can iterate from start to dest, but accumulation is to be made
  // from dest to start.. this is why we stack the vertices.
  std::stack<graph::ref_vertex> visited;
  auto it  = vq3::path::begin(start);
  auto end = vq3::path::end<graph::ref_vertex>();
  while(it != end) visited.push(*it++);
  double length = 0;
  while(!visited.empty()) {
    (*(visited.top()))().vq3_length = length++;
    visited.pop();
  }
}
	
double cumulated_cost_of(graph::ref_vertex ref_v) {
  // We only ave to return the lengths we have accumulated at that
  // vertex. Here, we have used the vq3::decorator::path::length
  // decoration to hold that value.
  return (*ref_v)().vq3_length;
}
				  

// Callback
//
//////////////

struct callback_data {
  graph& g;
  demo2d::opencv::Frame& frame;
  graph::ref_vertex start, dest;
  callback_data(graph& g, demo2d::opencv::Frame& frame) : g(g), frame(frame) {}
  callback_data()                                = delete;
  callback_data(const callback_data&)            = delete;
  callback_data& operator=(const callback_data&) = delete;
  					
};

void on_mouse( int event, int x, int y, int, void* user_data) {
  // A mouse click reinitializes the graph.
  if(event != cv::EVENT_LBUTTONDOWN )
    return;
  
  auto& data   = *(reinterpret_cast<callback_data*>(user_data));
  auto closest = vq3::utils::closest(data.g, data.frame(cv::Point(x,y)), d2);
    
  if(closest != data.start && closest != data.dest) {
    auto& efficient = (*closest)().vq3_efficient;
    efficient = !efficient;
  }
    
  // true false: we do not consider edge efficiency, but only vertex efficiency.
  vq3::path::dijkstra<true, false>(data.g, data.start, data.dest, edge_cost);
}


// Main
//
//////////////

#define NB_VERTICES_PER_M2      200
#define NB_SAMPLES_PER_M2    100000

int main(int argc, char* argv[]) {

  std::random_device rd;  
  std::mt19937 random_device(rd());

  // Let us build the graph.
  //
  //////////////

  double intensity = 1.;
  double radius    = 0.5;
  auto density     = demo2d::sample::disk(radius, intensity);

  graph g;
  
  auto sampler = demo2d::sample::base_sampler::random(random_device, NB_VERTICES_PER_M2);

  // We add vertices
  for(auto& location : demo2d::sample::sample_set(random_device, sampler, density)) g += location;

  // We add edges
  sampler = NB_SAMPLES_PER_M2;
  for(auto& sample : demo2d::sample::sample_set(random_device, sampler, density)) {
    auto closest = vq3::utils::two_closest(g, sample, d2);
    if(g.get_edge(closest.first, closest.second) == nullptr) 
      g.connect(closest.first, closest.second);
  }

  // We handle vertex efficiencies. 
  vq3::utils::clear_vertex_efficiencies(g, true);
  
  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  auto image = cv::Mat(800, 800, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = demo2d::opencv::direct_orthonormal_frame(image.size(), .9*image.size().height, true);  
  callback_data user_data(g, frame);

  user_data.start = vq3::utils::closest(g, demo2d::Point(-100, 0), d2);
  user_data.dest  = vq3::utils::closest(g, demo2d::Point( 100, 0), d2);

  // This is for path drawing.
  vq3::utils::clear_vertex_tags(g, false);
  (*(user_data.start))().vq3_tag = true;
  (*(user_data.dest ))().vq3_tag = true;
  
  // true false: we do not consider edge efficiency, but only vertex efficiency.
  vq3::path::dijkstra<true, false>(g, user_data.start, user_data.dest, edge_cost);

  cv::setMouseCallback("image", on_mouse, reinterpret_cast<void*>(&user_data));

  int slider = 500;
  cv::createTrackbar("walk ratio", "image", &slider, 1000, nullptr);
  
  
  auto draw_edge = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
  								     [](const vertex&, const vertex&, const edge&) {return   true;},  // always draw
  								     [](const vertex& v) {return                      v.vq3_value;},  // position
  								     [](const edge& e)   {
								       if(e.vq3_tag)            return cv::Scalar(000, 000, 000);
								       else                     return cv::Scalar(200,  80,  80);},   // color
  								     [](const edge& e)   {
								       if(e.vq3_tag)            return 3;
								       else                     return 1;});                          // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return                  true;},           // always draw
									   [](const vertex& v) {return           v.vq3_value;},           // position
									   [](const vertex& v) {
									     if(v.vq3_tag) return 6;
									     else return          4;},                                    // radius
									   [](const vertex& v) {
									     if(v.vq3_tag)            return cv::Scalar(000, 000, 200);
									     else if(v.vq3_efficient) return cv::Scalar(200,  80,  80);
									     else                     return cv::Scalar(255, 200, 200);}, // color
									   [](const vertex& v) {return -1;});                             // thickness
  
    
  std::cout << std::endl
            << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl
            << "Click on vertices to toggle their efficiency." << std::endl
	    << "Press ESC to quit." << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl;
  int keycode = 0;
  while(keycode != 27) {    

    // Let us recompute the path drawing info
    
    vq3::utils::clear_edge_tags(g, false);
    auto end = vq3::path::end(g);
    for(auto it = vq3::path::begin(user_data.start); it != end; ++it)
      if(auto ref_e = it.get_edge(); ref_e)
	(*ref_e)().vq3_tag = true;

    // Let us draw
    
    image = cv::Scalar(255, 255, 255);
    g.foreach_edge(draw_edge); 
    g.foreach_vertex(draw_vertex);

    // Let us compute and draw the walking position on the path.
    // Let us use default functions for dealing with accumulations.
    auto position = vq3::path::travel(user_data.start, user_data.dest,
				      slider*.001,
				      vertex_interpolation,
				      vq3::path::travel_defaults::compute_cumulated_costs<graph::ref_vertex>,
				      vq3::path::travel_defaults::cumulated_cost_of<graph::ref_vertex>);
    if(position) // position is an optional
      cv::circle(image, frame(*position), 7, cv::Scalar(000, 200, 200), -1); // Filled circle.

    // Now, we can do the same, but using custom measurement of path
    // length for interpolation. Here, we consider each edge haveing a
    // length of 1.
    position = vq3::path::travel(user_data.start, user_data.dest,
				 slider*.001,
				 vertex_interpolation,
				 compute_cumulated_costs, // We use our functions rather...
				 cumulated_cost_of);      // ... than the default ones.
    if(position) 
      cv::circle(image, frame(*position), 7, cv::Scalar(000, 200, 200),  3); // Outlined circle.
    
    cv::imshow("image", image);
    keycode = cv::waitKey(50) & 0xFF;
  }

  return 0;
}
