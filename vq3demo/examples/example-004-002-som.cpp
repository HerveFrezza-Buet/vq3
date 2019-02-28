#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <iterator>
#include <chrono>

/*

  Here, we build up and run a Kohonen self-organizing map from vq3
  dedicated tools. We will use multi-threading in order to speed-up
  the computation. Moreover, just for fun, we will not use a grid but
  a crown-shaped mesh.

*/


#define NB_CHL_SAMPLES        10000
#define SOM_NB_PROTOTYPES       400
#define SOM_H_RADIUS            5.1
#define SOM_MAX_DIST            (unsigned int)(SOM_H_RADIUS)

using sample    = vq3::demo2d::Point;
using prototype = vq3::demo2d::Point;

//                                                      ## Node properties :
using layer_0 = prototype;                              // prototypes are 2D points (this is the "user defined" value).
using layer_1 = vq3::decorator::tagged<layer_0>;        // we add a tag for topology computation.
using layer_2 = vq3::demo::decorator::colored<layer_1>; // we add a color to the nodes.
using vertex  = layer_2;

using graph  = vq3::graph<vertex, void>;

using topology_key_type = int;

// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v.vq3_value, p);}

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
  data.g.foreach_vertex([&data](const graph::ref_vertex& v_ref) {(*v_ref)().vq3_value = vq3::demo2d::uniform(data.rd, {-.5, -.5}, {.5, .5});});
}

// At each epoch, the sample are presented and some data is collected,
// for each vertex. From each data epoch, the new vertex value can be
// computed. vq3 offers classes for epoch data, they are intended to
// be stacked in order to define an agglomerate epoch data class.

using epoch_data_0 = vq3::epoch::data::none<sample, vertex, prototype>; // This is the root of the stack.
using epoch_data_1 = vq3::epoch::data::wtm<epoch_data_0>;               // This gathers computation for batch winner-take-most.
using epoch_data   = epoch_data_1;

int main(int argc, char* argv[]) {
  if(argc != 3) {
    std::cout << "Usage : " << argv[0] << " nb_threads nb_samples" << std::endl;
    return 0;
  }

  unsigned int nb_threads = std::atoi(argv[1]);
  unsigned int nb_samples = std::atoi(argv[2]);
  
  std::random_device rd;  
  std::mt19937 random_device(rd());

  // Let us build up a dataset.
  //
  ///////////////////
  
  double intensity   = 1. ;
  double radius      =  .5;
  double hole        =  .3;
  auto density       = vq3::demo2d::sample::disk(radius, intensity) - vq3::demo2d::sample::disk(hole, intensity);

  std::vector<vq3::demo2d::Point> data;
  data.reserve(nb_samples);
  auto out = std::back_inserter(data);
  for(unsigned int i = 0; i < nb_samples; ++i)
    *(out++) = vq3::demo2d::sample::get_one_sample(random_device, density);

  // Let us build up a grid-like graph, from a Competitive Hebbian Learning (CHL) process.
  //
  ///////////////////
  
  double intensity_         =   1;
  double side_              = 255;
  vq3::demo2d::Point center = {.5*side_, .5*side_};
  auto density_             = vq3::demo2d::sample::rectangle(side_, side_, intensity_) + center;

  graph g;

  // we define the nodes
  for(unsigned int i = 0; i < SOM_NB_PROTOTYPES; ++i) {
    auto pos = vq3::demo2d::sample::get_one_sample(random_device, density_);
    auto ref_v = g += pos;
    (*ref_v)().vq3_color = cv::Scalar(0, pos.x, pos.y);
  }

  // we define the edges
  for(unsigned int i = 0; i < NB_CHL_SAMPLES; ++i) {
    auto closest = vq3::utils::two_closest(g, vq3::demo2d::sample::get_one_sample(random_device, density_), d2);
    if(g.get_edge(closest.first, closest.second) == nullptr) 
      g.connect(closest.first, closest.second);
  }

  // we initialize the prototypes to a random value.
  g.foreach_vertex([&random_device](const graph::ref_vertex& v_ref) {(*v_ref)().vq3_value = vq3::demo2d::uniform(random_device, {-.5, -.5}, {.5, .5});});


  // Image settings
  //
  ///////////////////

  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  auto image       = cv::Mat(480, 640, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .9*image.size().height, true);
  callback_data user_data(g, random_device);
  cv::setMouseCallback("image", on_mouse, reinterpret_cast<void*>(&user_data));
  auto draw_edge   = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
								       [](const vertex& v1, const vertex& v2) {return true;}, // always draw
								       [](const vertex& v)     {return         v.vq3_value;}, // position
								       []()                    {return cv::Scalar(0, 0, 0);}, // color
								       []()                    {return                  1;}); // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return        true;},  // always draw
									   [](const vertex& v) {return v.vq3_value;},  // position
									   [](const vertex& v) {return           3;},  // radius
									   [](const vertex& v) {return v.vq3_color;},  // color
									   [](const vertex& v) {return          -1;}); // thickness
  
    

  // The is the loop....
  //
  ///////////////////

  image = cv::Scalar(255, 255, 255);
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

  // First, we need a structure for handling SOM-like computation. The
  // template argument is the type of input samples.

  auto topology = vq3::topology::table<int>(g);
  topology.declare_distance(0,
			    [](unsigned int edge_distance) {return std::max(0., 1 - edge_distance/double(SOM_H_RADIUS));},
			    SOM_MAX_DIST,
			    1e-3); // We consider node and edge-based neihborhoods.
  topology.update_full();
  // This is the winner-take-most parallel processor.
  auto winner_take_most = vq3::epoch::wtm::processor(topology);
			    

  
  while(keycode != 27) {

    auto t_start = std::chrono::high_resolution_clock::now();
    // Learning : the returned value of thus function is ignored here. See next examples.
    winner_take_most.process<epoch_data>(nb_threads, 0,
					 data.begin(), data.end(),
					 [](const vq3::demo2d::Point& p) -> const vq3::demo2d::Point& {return p;},         // how to get the sample from a data content (*it).
					 [](vertex& vertex_value) -> vq3::demo2d::Point& {return vertex_value.vq3_value;}, // how to get a **reference** to the prototype from the vertex value.
					 d2);
    auto t_end = std::chrono::high_resolution_clock::now();
    
    int duration =  std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
    std::cout << "Step duration (" << nb_threads << " threads, " << nb_samples << " samples) : " << std::setw(6) << duration << " ms.    \r"<< std::flush;
    
    image = cv::Scalar(255, 255, 255);
    g.foreach_edge(draw_edge); 
    g.foreach_vertex(draw_vertex);
    cv::imshow("image", image);
    keycode = cv::waitKey(1) & 0xFF;
  }

  std::cout << std::endl;

  return 0;
}
