#include <vq3demo.hpp>
#include <random>

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
using vlayer_1 = vq3::decorator::efficiency<vlayer_0>;                             // for connected components
using vlayer_2 = vq3::decorator::tagged<vlayer_1>;                                 // we add a tag for topology computation and connected components.
using vlayer_3 = vq3::decorator::online::mean_std<vlayer_2, double, Param>;        // we add distortion statistics over time. 
using vlayer_4 = vq3::decorator::smoother<vlayer_3, vq3::demo2d::Point, 1, 21, 2>; // we smooth of prototypes.
using vlayer_5 = vq3::demo::decorator::colored<vlayer_4>;
using vertex   = vlayer_5;

//                                                        ## Edge properties :
using elayer_0 = vq3::decorator::tagged<void>;            // we add a tag for CHL computation.
using elayer_1 = vq3::demo::decorator::colored<elayer_0>;
using edge     = elayer_1;

using graph  = vq3::graph<vertex, edge>;
  

// Distance
//
////////////////


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double dist(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v.vq3_value, p);}




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
  
  auto color_of_label = vq3::demo2d::opencv::colormap::random(random_device);

  vq3::demo2d::opencv::HueSelector selector;
  auto video_data = vq3::demo2d::opencv::sample::video_data(0, selector.build_pixel_test());
  
  int N_slider =   500;
  int T_slider =   500;
  int S_slider =   170;
  int K_slider =     5;
  
  
  // Input distribution
  //
  ///////////////////

  auto density = vq3::demo2d::opencv::sample::webcam(video_data);

  auto input_size  = video_data.image.size();
  video_data.frame = vq3::demo2d::opencv::direct_orthonormal_frame(input_size, .5*input_size.width, true);
  
  
  // Display
  //
  ///////////////////

  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  cv::createTrackbar("nb/m^2",         "image", &N_slider,   2000, nullptr);
  cv::createTrackbar("T",              "image", &T_slider,   1000, nullptr);
  cv::createTrackbar("100*sigma_coef", "image", &S_slider,    300, nullptr);
  cv::createTrackbar("nb SOM steps",   "image", &K_slider,     20, nullptr);
  
  cv::namedWindow("video", CV_WINDOW_AUTOSIZE);
  selector.build_sliders("video");
  
  auto image = cv::Mat(600, 800, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .4*image.size().width, true);
  
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
  auto evolution = vq3::algo::gngt::by_default::evolution(random_device);
  
  gngt.nb_threads       = nb_threads;
  gngt.neighbour_weight = .1;
  gngt.distance         = dist;
  gngt.prototype        = [](vertex& v) -> prototype& {return v.vq3_value;};
  
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
    
    // Get the samples and plot them.
    
    auto S_ = vq3::demo2d::sample::sample_set(random_device, density, N_slider);
    S.clear();
    std::copy(S_.begin(), S_.end(), std::back_inserter(S));

    // Step

    double e = T_slider/1000.0;
    double expo_min = -5;
    double expo_max = -1;
    
    evolution.density    = N_slider;
    evolution.T          = std::pow(10, expo_min*(1-e) + expo_max*e);
    evolution.sigma_coef = S_slider*.01;
    
    gngt.epoch(K_slider, 1,
	       S.begin(), S.end(),
	       [](const sample& s) {return s;},
	       [](const prototype& p) {return p + vq3::demo2d::Point(-1e-5,1e-5);},
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
    
    cv::imshow("image", image);
    selector.build_image(video_data.image);
    cv::imshow("video", selector.image);
    keycode = cv::waitKey(1) & 0xFF;
  }
  
  return 0;
}
