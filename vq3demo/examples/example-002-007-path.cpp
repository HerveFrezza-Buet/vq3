#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>
#include <algorithm>
#include <vector>
#include <iterator>


// This example shows the use of path walking.


// Graph definition
//
///////////////

// Vertex values have to be instrumented with shortest path
// computation structure.
//                                                         ## Node properties :
using vlayer_0 = demo2d::Point;                            // prototypes are 2D points (this is the "user defined" value).
using vlayer_1 = vq3::decorator::path::shortest<vlayer_0>; // This holds informations built by dijkstra.
using vlayer_2 = vq3::decorator::tagged<vlayer_1>;         // We will tag extermities of the shortest path for display.
using vlayer_3 = vq3::decorator::efficiency<vlayer_2>;     // We consider only efficient edges for paths.
using vertex   = vlayer_3;

// Here, each edge hosts its cost value.  It is the optional value
// (not mandatory) so that it is not recomputed once it has been
// calculated first.
//                                                        ## Node properties :
using elayer_0 = vq3::decorator::optional_cost<void>;     // Edge cost.
using elayer_1 = vq3::decorator::tagged<elayer_0>;        // We will tag edges belonging to a shortest path for display.
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
  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  auto image = cv::Mat(480, 640, CV_8UC3, cv::Scalar(255,255,255));
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
    auto position = vq3::path::travel(user_data.start, user_data.dest,
				      slider*.001,
				      [](const graph::ref_vertex& start,
					 const graph::ref_vertex& dest,
					 double lambda) {
					// This is interpolation.
					if(start == nullptr) return (*dest )().vq3_value;
					if(dest  == nullptr) return (*start)().vq3_value;
					auto& start_value = (*start)().vq3_value;
					auto& dest_value  = (*dest )().vq3_value;
					return (1 - lambda) * start_value + lambda * dest_value;
				      });
    if(position) // position is an optional
      cv::circle(image, frame(*position), 5, cv::Scalar(000, 200, 200), -1);
    
    
    cv::imshow("image", image);
    keycode = cv::waitKey(50) & 0xFF;
  }

  return 0;
}
