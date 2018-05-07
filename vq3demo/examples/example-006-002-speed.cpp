#include <vq3demo.hpp>
#include <random>
#include <algorithm>

#define ANGLE_PERIOD 500
#define D_ANGLE (360./ANGLE_PERIOD)
#define SPEED_TO_METER .5

// Graph definition
//
///////////////

using sample    = vq3::demo2d::Point;
using prototype = vq3::demo2d::Point;

struct Param {
  double alpha() {return .1;}
  unsigned int min_updates() {return 3;}
};

//                                                                                 ## Node properties :
using vlayer_0 = prototype;                                                        // prototypes are 2D points (this is the "user defined" value).
using vlayer_1 = vq3::decorator::tagged<vlayer_0>;                                 // we add a tag for topology computation.
using vlayer_2 = vq3::decorator::online::mean_std<vlayer_1, double, Param>;        // we add distortion statistics over time. 
using vlayer_3 = vq3::decorator::smoother<vlayer_2, vq3::demo2d::Point, 1, 21, 2>; // we smooth of prototypes.
using vertex   = vlayer_3;

//                                                        ## Edge properties :
using elayer_0 = vq3::decorator::tagged<void>;            // we add a tag for CHL computation.
using edge     = elayer_0;

using graph  = vq3::graph<vertex, edge>;
  

// Epoch data for SOM-like pass
//
////////////////

using epoch_wtm = vq3::epoch::data::wtm<vq3::epoch::data::none<sample>>;


// Distance
//
////////////////


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double dist(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v.vq3_value, p);}


// Here is the main.

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
  
  int N_slider =   500;
  int T_slider =   400;
  int S_slider =   170;
  int H_slider =   250;
  int W_slider =     1;
  int K_slider =     2;
  
  int Z_slider =  2000;
  
  
  // Input distribution
  //
  ///////////////////

#define BAR_SIDE 2.
  
  auto bar_theta     = vq3::demo::dyn::linear<vq3::demo::dyn::bound::wrap>(0, 0, 360, D_ANGLE);
  auto bar_pos_x     = vq3::demo::dyn::sin(-2, 2,  0,     D_ANGLE);
  auto bar_size_x    = vq3::demo::dyn::cos(BAR_SIDE, 5,  0, 1.25*D_ANGLE);
  auto bar_pos       = vq3::demo2d::dyn::point(bar_pos_x, 0.);
  auto bar_size      = vq3::demo2d::dyn::point(bar_size_x, BAR_SIDE);
  double bar_density =   1;
  double bar_width   =   1;
  double bar_height  =   1;
  auto bar = (vq3::demo2d::sample::rectangle(bar_width, bar_height, bar_density) * bar_size()) % bar_theta() + bar_pos();
  
  
  auto density = bar;
  

  
  
  // Display
  //
  ///////////////////

  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  cv::createTrackbar("nb/m^2",              "image", &N_slider,   2000, nullptr);
  cv::createTrackbar("T",                   "image", &T_slider,   1000, nullptr);
  cv::createTrackbar("100*sigma_coef",      "image", &S_slider,    300, nullptr);
  cv::createTrackbar("100*h_radius",        "image", &H_slider,   1000, nullptr);
  cv::createTrackbar("nb wide SOM steps",   "image", &W_slider,     20, nullptr);
  cv::createTrackbar("nb narrow SOM steps", "image", &K_slider,     20, nullptr);
  
  auto image = cv::Mat(600, 800, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .1*image.size().width, true);


  cv::namedWindow("speed", CV_WINDOW_AUTOSIZE);
  cv::createTrackbar("zoom", "speed", &Z_slider, 3000, nullptr);
  auto speed_image = cv::Mat(500, 500, CV_8UC3, cv::Scalar(255,255,255));

  
  auto dd = vq3::demo2d::opencv::dot_drawer<vq3::demo2d::Point>(image, frame,
								[](const vq3::demo2d::Point& pt) {return                      true;},
								[](const vq3::demo2d::Point& pt) {return                        pt;},
								[](const vq3::demo2d::Point& pt) {return                         1;},
								[](const vq3::demo2d::Point& pt) {return cv::Scalar(255, 120, 120);},
								[](const vq3::demo2d::Point& pt) {return                        -1;});
  
  auto smooth_edge = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
								       [](const vertex& v1, const vertex& v2, const edge& e) {
									 return v1.vq3_smoother.get<0>() && v2.vq3_smoother.get<0>();},   // draw anly if smoothed vertices are available.
								       [](const vertex& v) {return   v.vq3_smoother.get<0>().value();},   // position
								       [](const edge& e)     {return         cv::Scalar(  0,   0, 200);}, // color
								       [](const edge& e)     {return                                1;}); // thickness
  
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


  // We need to register the input samples in a vector since we want
  // to both use and display them.
  std::vector<vq3::demo2d::Point> S;

  graph g;
  
  auto vertices  = vq3::utils::vertices(g);
  auto gngt      = vq3::algo::gngt::processor<prototype, sample>(g, vertices);
  auto wtm       = vq3::epoch::wtm::processor(g, vertices);
  auto evolution = vq3::algo::gngt::by_default::evolution(random_device);
  
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

  unsigned int step = 0;
  unsigned int mode = 0;
  
  int keycode = 0;
  while(keycode != 27) {

    // Ubdate the distribution
    switch(mode) {
    case 0: // Rotation
      ++bar_theta; 
      break;
    case 1: // stretch
      ++bar_size_x;
      ++bar_size;
      break;
    case 2: // translation
      ++bar_pos_x;
      ++bar_pos;
      break;
    default:
      break;
    }
    
    // Get the samples and plot them.
    
    auto S_ = vq3::demo2d::sample::sample_set(random_device, density, N_slider);
    S.clear();
    std::copy(S_.begin(), S_.end(), std::back_inserter(S));

    // Step

    // The input has evolved. We update the graph so that it fits the
    // new distribution, keeping the previously computed topology.

    vertices.update_topology(g);
    
    // First, wide kernel fit, for a quicker update.
    double h = std::max(H_slider,1)*.01;
    wtm.update_topology([h](unsigned int edge_distance) {return std::max(0., 1.-edge_distance/h);}, (unsigned int)h, 0);
    for(int wta_step = 0; wta_step < W_slider; ++wta_step)
      wtm.update_prototypes<epoch_wtm>(nb_threads,
    				       S.begin(), S.end(),
    				       [](const sample& s) {return s;},
    				       [](vertex& v) -> prototype& {return v.vq3_value;},
    				       dist);

    // Second, narrow kernel fit, allowing the graph to spread over the samples.
    wtm.update_topology([](unsigned int edge_distance) {return edge_distance == 0 ? 1.0 : 0.1;}, 1, 0);
    for(int wta_step = 0; wta_step < K_slider; ++wta_step)
      wtm.update_prototypes<epoch_wtm>(nb_threads,
    				       S.begin(), S.end(),
    				       [](const sample& s) {return s;},
    				       [](vertex& v) -> prototype& {return v.vq3_value;},
    				       dist);
    

    // Now, we can compute GNG-T topology evolution.
    
    double e = T_slider/1000.0;
    double expo_min = -5;
    double expo_max = -1;
    
    evolution.density    = N_slider;
    evolution.T          = std::pow(10, expo_min*(1-e) + expo_max*e);
    evolution.sigma_coef = S_slider*.01;
    
    gngt.epoch(nb_threads,
	       S.begin(), S.end(),
	       [](const sample& s) {return s;},
	       [](vertex& v) -> prototype& {return v.vq3_value;},
	       [](const prototype& p) {return p + vq3::demo2d::Point(-1e-5,1e-5);},
	       dist,
	       evolution);
    
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
    std::copy(S.begin(),  S.end(), dd);
    g.foreach_edge(smooth_edge); 
    g.foreach_vertex(smooth_speed);
    g.foreach_vertex(smooth_vertex);
    

    speed_image = cv::Scalar(255, 255, 255);
    double max_speed = Z_slider*.001;
    auto speed_frame = vq3::demo2d::opencv::direct_orthonormal_frame(speed_image.size(), .5*speed_image.size().width/std::max(max_speed, 1e-3), true);
    cv::line(speed_image,
	     speed_frame(vq3::demo2d::Point(0, -max_speed)),
	     speed_frame(vq3::demo2d::Point(0,  max_speed)),
	     cv::Scalar(0,0,0), 1);
    cv::line(speed_image,
	     speed_frame(vq3::demo2d::Point(-max_speed, 0)),
	     speed_frame(vq3::demo2d::Point( max_speed, 0)),
	     cv::Scalar(0,0,0), 1);
    g.foreach_vertex([&speed_image, &speed_frame](graph::ref_vertex ref_v) {
	auto& vertex = (*ref_v)();
	if(vertex.vq3_smoother.get<1>()) // If speed is available
	  cv::circle(speed_image, speed_frame(vertex.vq3_smoother.get<1>().value()), 3, cv::Scalar(0,0,0), -1);
      });
    
    cv::imshow("image", image);
    cv::imshow("speed", speed_image);
    keycode = cv::waitKey(1) & 0xFF;

    // Input mode
    ++step;
    if(step >= ANGLE_PERIOD) {
      step = 0;
      ++mode;
      if(mode == 3)
	mode = 0;
    }
  }
  
  return 0;
}
