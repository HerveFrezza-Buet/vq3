#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <vector>
#include <iterator>
#include <algorithm>
#include <random>


/*

  Here, we build up and run a batch Kohonen self-organizing map from
  scratch, using general vq3 features. Next examples show the use of processors 
  in order to run this in parallel.

*/



#define NB_SAMPLES_PER_M2   10000
#define GRID_WIDTH             20
#define GRID_HEIGHT            20

#define SOM_H_RADIUS            5.1
#define SOM_MAX_DIST            (unsigned int)(SOM_H_RADIUS)

using sample    = demo2d::Point;
using prototype = demo2d::Point;

//                                                           ## Node properties :
using layer_0 = prototype;                                   // prototypes are 2D points (this is the "user defined" value).
using layer_1 = vq3::decorator::tagged<layer_0>;             // we add a tag for topology computation.
using layer_2 = vq3::decorator::sum<layer_1, demo2d::Point>; // we add the ability to hande a sum of points.
using layer_3 = vq3::decorator::grid_pos<layer_2>;           // we add the ability to register a grid position.
using vertex  = layer_3;

using graph  = vq3::graph<vertex, void>;

using topology_key_type = int;

// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const demo2d::Point& p) {return demo2d::d2(v.vq3_value, p);}

struct callback_data {
  graph& g;
  std::mt19937& rd;

  callback_data(graph& g, std::mt19937& rd) : g(g), rd(rd) {}
  callback_data()                                = delete;
  callback_data(const callback_data&)            = delete;
  callback_data& operator=(const callback_data&) = delete;
};

void on_mouse( int event, int x, int y, int, void* user_data) {
  // A mouse click reinitializes the graph.
  if(event != cv::EVENT_LBUTTONDOWN )
    return;
  auto& data = *(reinterpret_cast<callback_data*>(user_data));
  data.g.foreach_vertex([&data](const graph::ref_vertex& v_ref) {(*v_ref)().vq3_value = demo2d::uniform(data.rd, {-.5, -.5}, {.5, .5});});
}

int main(int argc, char* argv[]) {
  std::random_device rd;  
  std::mt19937 random_device(rd());

  graph g;
  
  vq3::utils::make_grid(g, GRID_WIDTH, GRID_HEIGHT,
			[&random_device](unsigned int w, unsigned int h) {
			  graph::vertex_value_type value(demo2d::uniform(random_device, {-.5, -.5}, {.5, .5}));
			  value.vq3_gridpos = {w, h};
			  return value;
			});

  auto topology = vq3::topology::table<topology_key_type>(g);
  topology.declare_distance(0,
			    [](unsigned int edge_distance) {return std::max(0., 1 - edge_distance/double(SOM_H_RADIUS));},
			    SOM_MAX_DIST,
			    1e-3); // We consider node and edge-based neihborhoods.
  topology.update_full();


  double intensity = 1. ;
  double radius    =  .5;
  double hole      =  .3;
  auto density     = demo2d::sample::disk(radius, intensity) - demo2d::sample::disk(hole, intensity);

  
  auto sampler = demo2d::sample::base_sampler::random(random_device, NB_SAMPLES_PER_M2);
  auto S_      = demo2d::sample::sample_set(random_device, sampler, density); // This re-toss points at each begin...end iteration.
  
  std::vector<demo2d::Point> S; // Let us use a single sample of S_.
  auto out = std::back_inserter(S);
  std::copy(S_.begin(), S_.end(), out);
  
  
  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  auto image       = cv::Mat(480, 640, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = demo2d::opencv::direct_orthonormal_frame(image.size(), .9*image.size().height, true);
  callback_data user_data(g, random_device);
  cv::setMouseCallback("image", on_mouse, reinterpret_cast<void*>(&user_data));
  auto dd          = demo2d::opencv::dot_drawer<demo2d::Point>(image, frame,
							       [](const demo2d::Point& pt) {return                      true;},
							       [](const demo2d::Point& pt) {return                        pt;},
							       [](const demo2d::Point& pt) {return                         1;},
							       [](const demo2d::Point& pt) {return cv::Scalar(200, 200, 200);},
							       [](const demo2d::Point& pt) {return                        -1;});
  auto draw_edge   = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
								       [](const vertex& v1, const vertex& v2) {return true;}, // always draw
								       [](const vertex& v)     {return         v.vq3_value;}, // position
								       []()                    {return cv::Scalar(0, 0, 0);}, // color
								       []()                    {return                  1;}); // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return                                                   true;},
									   [](const vertex& v) {return                                            v.vq3_value;},  // position
									   [](const vertex& v) {return                                                      3;},  // radius
									   [](const vertex& v) {return cv::Scalar(0,
														  255*v.vq3_gridpos.first/(GRID_WIDTH - 1),
														  255*v.vq3_gridpos.second/(GRID_HEIGHT - 1));},  // color
									   [](const vertex& v) {return                                                     -1;}); // thickness
  

  
  image = cv::Scalar(255, 255, 255);
  std::copy(S.begin(), S.end(), dd);
  g.foreach_edge(draw_edge); 
  g.foreach_vertex(draw_vertex);
  cv::imshow("image", image);
  cv::waitKey(1000);

  
  
  std::cout << std::endl
	    << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << "click in the image to restart, press ESC to quit." << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl;
  int keycode = 0;
  while(keycode != 27) {

    // This is the SOM algorithm (with minibatch) : 2 steps

    // Step #1 : submit all samples.
    for(auto& sample : S)
      for(auto& topo_info : topology.neighborhood(vq3::utils::closest(g, sample, d2), 0))
	(*(topology(topo_info.index)))().vq3_sum.increment(topo_info.value, sample);
    
    // Step #2 : update prototypes.
    g.foreach_vertex([](const graph::ref_vertex& ref_v) {
	auto& v = (*ref_v)();
	if(v.vq3_sum.nb != 0)
	  v.vq3_value = v.vq3_sum.average<double>();
	v.vq3_sum.clear();
      });
    
    image = cv::Scalar(255, 255, 255);
    std::copy(S.begin(), S.end(), dd);
    g.foreach_edge(draw_edge); 
    g.foreach_vertex(draw_vertex);
    cv::imshow("image", image);
    keycode = cv::waitKey(1) & 0xFF;
  }
  
  return 0;
}
