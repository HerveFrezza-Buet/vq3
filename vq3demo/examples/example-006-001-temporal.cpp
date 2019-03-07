#include <vq3demo.hpp>
#include <random>
#include <string>

#define GRID_WIDTH             15
#define GRID_HEIGHT             8

#define SOM_H_RADIUS            5.1
#define SOM_MAX_DIST            (unsigned int)(SOM_H_RADIUS)

//                                                                               ## Node properties :
using layer_0 = vq3::demo2d::Point;                                              // prototypes are 2D points (this is the "user defined" value).
using layer_1 = vq3::decorator::tagged<layer_0>;                                 // we add a tag for topology computation.
using layer_2 = vq3::decorator::smoother<layer_1, vq3::demo2d::Point, 1, 21, 2>; // smoothing of prototypes.
using vertex  = layer_2;

using graph  = vq3::graph<vertex, void>;

using neighbour_key_type = std::string;

// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v.vq3_value, p);}

using epoch_data_0 = vq3::epoch::data::none<vq3::demo2d::Point, vertex, vq3::demo2d::Point>; // This is the root of the stack.
using epoch_data_1 = vq3::epoch::data::wtm<epoch_data_0>;                                    // This gathers computation for batch winner-take-most.
using epoch_data   = epoch_data_1;

#define SPEED_TO_METER .5

int main(int argc, char* argv[]) {
  
  unsigned int nb_threads = std::thread::hardware_concurrency();
  if(nb_threads == 0) {
    nb_threads = 1;
    std::cout << "The harware multi-threading capabilities cannot be queried. I use a single thread." << std::endl;
  }
  else
    std::cout << "I use " << nb_threads << " thread(s)." << std::endl;
  std::cout << std::endl;

  std::random_device rd;  
  std::mt19937 random_device(rd());

  
  // Distribution and display
  //
  ///////////////////
  
  int N_slider =  5000;

  vq3::demo2d::opencv::HueSelector selector;
  auto video_data = vq3::demo2d::opencv::sample::video_data(0, selector.build_pixel_test());

  auto density = vq3::demo2d::opencv::sample::webcam(video_data);

  auto input_size  = video_data.image.size();
  video_data.frame = vq3::demo2d::opencv::direct_orthonormal_frame(input_size, .5*input_size.width, true);

  
  cv::namedWindow("video", CV_WINDOW_AUTOSIZE);
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  cv::createTrackbar("nb/m^2",           "image", &N_slider, 10000, nullptr);
  selector.build_sliders("video");

  
  auto image = cv::Mat(600, 800, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .4*image.size().width, true);
  
  auto dd    = vq3::demo2d::opencv::dot_drawer<vq3::demo2d::Point>(image, frame,
								   [](const vq3::demo2d::Point& pt) {return                      true;},
								   [](const vq3::demo2d::Point& pt) {return                        pt;},
								   [](const vq3::demo2d::Point& pt) {return                         1;},
								   [](const vq3::demo2d::Point& pt) {return cv::Scalar(255, 120, 120);},
								   [](const vq3::demo2d::Point& pt) {return                        -1;});

  auto draw_edge   = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
								       [](const vertex& v1, const vertex& v2) {return   true;}, // always drawing
  								       [](const vertex& v) {return               v.vq3_value;}, // position
  								       []()                {return cv::Scalar(127, 127, 127);}, // color
  								       []()                {return                        1;}); // thickness
  
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return                      true;},  // always drawing
  									   [](const vertex& v) {return               v.vq3_value;},  // position
  									   [](const vertex& v) {return                         3;},  // radius
  									   [](const vertex& v) {return cv::Scalar(127, 127, 127);},  // color
  									   [](const vertex& v) {return                        -1;}); // thickness

  auto smooth_edge   = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
									 [](const vertex& v1, const vertex& v2) {
									   return v1.vq3_smoother.get<0>() && v2.vq3_smoother.get<0>();}, // draw anly if smoothed vertices are available.
									 [](const vertex& v) {return   v.vq3_smoother.get<0>().value();}, // position
									 []()                {return         cv::Scalar(  0,   0, 200);}, // color
									 []()                {return                                1;}); // thickness
  
  auto smooth_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									     [](const vertex& v) {return (bool)(v.vq3_smoother.get<0>());},  // draw only is the smoothed vertex is available.
									     [](const vertex& v) {return v.vq3_smoother.get<0>().value();},  // position
									     [](const vertex& v) {return                               3;},  // radius
									     [](const vertex& v) {return       cv::Scalar(  0,   0, 200);},  // color
									     [](const vertex& v) {return                              -1;}); // thickness
  
  auto smooth_speed = vq3::demo2d::opencv::segment_at_vertex_drawer<graph::ref_vertex>(image, frame,
										       [](const vertex& v) {return (bool)(v.vq3_smoother.get<0>());},  // draw only is the smoothed vertex is available.
										       [](const vertex& v) {return v.vq3_smoother.get<0>().value();},  // position
										       [](const vertex& v) {return v.vq3_smoother.get<1>().value()
													    *                      -SPEED_TO_METER;},  // speed
										       [](const vertex& v) {return       cv::Scalar(  0, 200,   0);},  // color
										       [](const vertex& v) {return                               3;}); // thickness

  // Data for computation
  //
  ///////////////////

  graph g;
  
  vq3::utils::make_grid(g, GRID_WIDTH, GRID_HEIGHT,
			[&random_device](unsigned int w, unsigned int h) {
			  graph::vertex_value_type value(vq3::demo2d::uniform(random_device, {-.5, -.5}, {.5, .5}));
			  return value;
			});
  

  // We need to register the input samples in a vector since we want
  // to both use and display them.
  std::vector<vq3::demo2d::Point> S;
  
  // This keeps up to date information about the graph topology.
  auto topology = vq3::topology::table<neighbour_key_type>(g);
  topology.declare_distance("som",   [](unsigned int edge_distance) {return std::max(0., 1 - edge_distance/double(SOM_H_RADIUS));}, SOM_MAX_DIST, 1e-3);

  // This is the SOM algorithm.
  auto som = vq3::algo::som::processor(topology);

  // We compute all the neighborhoods once, beforehand.
  topology.update_full();

  
  // This is the loop
  //
  //////////////////

  std::cout << std::endl
	    << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << "click in the image to restart, press ESC to quit." << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl;

  vq3::temporal::dt_averager frame_delay(.05);
  
  int keycode = 0;
  while(keycode != 27) {
    ++video_data; // get next frame.
    
    // Get the samples
    
    auto S_ = vq3::demo2d::sample::sample_set(random_device, density, N_slider);
    S.clear();
    std::copy(S_.begin(), S_.end(), std::back_inserter(S));

    // SOM update
    
    som.process<epoch_data>(nb_threads, "som",
			    S.begin(), S.end(),
			    [](const vq3::demo2d::Point& p) -> const vq3::demo2d::Point& {return p;},
			    [](vertex& vertex_value) -> vq3::demo2d::Point& {return vertex_value.vq3_value;},
			    d2);
    
    // Temporal update
    
    frame_delay.tick();
    double delay = frame_delay().value_or(1.);
    
    g.foreach_vertex([delay](graph::ref_vertex ref_v) {
	auto& value = (*ref_v)();
	value.vq3_smoother += value.vq3_value;
	value.vq3_smoother.set_timestep(delay);
      });
    
    // Display
    
    image = cv::Scalar(255, 255, 255);
    
    std::copy(S.begin(), S.end(), dd);
    g.foreach_edge(draw_edge); 
    g.foreach_vertex(draw_vertex);
    g.foreach_edge(smooth_edge); 
    g.foreach_vertex(smooth_speed);
    g.foreach_vertex(smooth_vertex);
    
    cv::imshow("image", image);
    selector.build_image(video_data.image);
    cv::imshow("video", selector.image);
    keycode = cv::waitKey(1) & 0xFF;
  }
  
  return 0;
}
