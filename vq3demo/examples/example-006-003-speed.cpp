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
using vlayer_5 = vq3::decorator::labelled<vlayer_4>;                               // for vertex labelling
using vertex   = vlayer_5;

//                                                        ## Edge properties :
using elayer_0 = vq3::decorator::tagged<void>;            // we add a tag for CHL computation.
using elayer_1 = vq3::decorator::efficiency<elayer_0>;    // for connected components
using elayer_2 = vq3::decorator::labelled<elayer_1>;      // for edge labelling     
using edge     = elayer_2;

using graph  = vq3::graph<vertex, edge>;
  

// Epoch data for SOM-like pass
//
////////////////

using epoch_wtm = vq3::epoch::data::wtm<vq3::epoch::data::none<sample, vertex, prototype>>;
  

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
  
  auto color_of_label = vq3::demo2d::opencv::colormap::random(random_device);

  vq3::demo2d::opencv::HueSelector selector;
  selector.H_slider =   0;
  selector.T_slider =  20;
  selector.S_slider =  80;
  selector.V_slider = 200;
  auto video_data = vq3::demo2d::opencv::sample::video_data(0, selector.build_pixel_test());
  
  int N_slider =  5000;
  int T_slider =   400;
  int S_slider =   170;
  int E_slider =   100;
  
  int Z_slider =  2000;
  
  
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
  cv::createTrackbar("nb/m^2",               "image", &N_slider, 10000, nullptr);
  cv::createTrackbar("1000*min_edge_length", "image", &E_slider,   500, nullptr);
  cv::createTrackbar("T",                    "image", &T_slider,  1000, nullptr);
  cv::createTrackbar("100*sigma_coef",       "image", &S_slider,   300, nullptr);
  
  cv::namedWindow("video", CV_WINDOW_AUTOSIZE);
  selector.build_sliders("video");
  
  auto image = cv::Mat(600, 800, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .4*image.size().width, true);
  
  cv::namedWindow("speed", CV_WINDOW_AUTOSIZE);
  cv::createTrackbar("zoom", "speed", &Z_slider, 3000, nullptr);
  auto speed_image = cv::Mat(500, 500, CV_8UC3, cv::Scalar(255,255,255));
  
  auto dd = vq3::demo2d::opencv::dot_drawer<vq3::demo2d::Point>(image, frame,
								[](const vq3::demo2d::Point& pt) {return                      true;},
								[](const vq3::demo2d::Point& pt) {return                        pt;},
								[](const vq3::demo2d::Point& pt) {return                         1;},
								[](const vq3::demo2d::Point& pt) {return cv::Scalar(230, 230, 230);},
								[](const vq3::demo2d::Point& pt) {return                        -1;});
  
  auto smooth_edge = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
								       [](const vertex& v1, const vertex& v2, const edge& e) {
									 return v1.vq3_smoother.get<0>() && v2.vq3_smoother.get<0>();},   // draw anly if smoothed vertices are available.
								       [](const vertex& v) {return   v.vq3_smoother.get<0>().value();},   // position
								       [&color_of_label](const edge& e)   {
									 if(e.vq3_efficient)
									   return color_of_label(e.vq3_label);
									 else
									   return cv::Scalar(200, 200, 200);},                            // color
								       [](const edge& e)     {return                                1;}); // thickness
  
  auto smooth_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									     [](const vertex& v) {return (bool)(v.vq3_smoother.get<0>());},  // draw only is the smoothed vertex is available.
									     [](const vertex& v) {return v.vq3_smoother.get<0>().value();},  // position
									     [](const vertex& v) {return                               3;},  // radius
									     [&color_of_label](const vertex& v) {
									     if(v.vq3_efficient)
									       return color_of_label(v.vq3_label);
									     else
									       return cv::Scalar(200, 200, 200);},                           // color
									     [](const vertex& v) {return                              -1;}); // thickness
  
  auto smooth_speed = vq3::demo2d::opencv::segment_at_vertex_drawer<graph::ref_vertex>(image, frame,
										       [](const vertex& v) {return (bool)(v.vq3_smoother.get<0>());},  // draw only is the smoothed vertex is available.
										       [](const vertex& v) {return v.vq3_smoother.get<0>().value();},  // position
										       [](const vertex& v) {return v.vq3_smoother.get<1>().value()
													    *                      -SPEED_TO_METER;},  // speed
										       [](const vertex& v) {return       cv::Scalar(  0, 0, 0);},      // color
										       [](const vertex& v) {return                               1;}); // thickness
  
  
  // Data for computation
  //
  ///////////////////


  // We need to register the input samples in a vector since we want
  // to both use and display them.
  std::vector<vq3::demo2d::Point> S;

  
  graph g;
  auto topology  = vq3::topology::table(g);
  auto gngt      = vq3::algo::gngt::processor<prototype, sample>(topology);
  auto wtm       = vq3::epoch::wtm::processor(topology);
  auto evolution = vq3::algo::gngt::by_default::evolution();
  
  // This is the loop
  //
  //////////////////

  std::cout << std::endl
	    << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << "press ESC to quit." << std::endl
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

    topology([](unsigned int edge_distance) {return edge_distance == 0 ? 1.0 : 0.1;}, 1, 0);
    for(int wta_step = 0; wta_step < 2; ++wta_step)
      wtm.process<epoch_wtm>(nb_threads,
			     S.begin(), S.end(),
			     [](const sample& s) {return s;},
			     [](vertex& v) -> prototype& {return v.vq3_value;},
			     dist);
    
    double e = T_slider/1000.0;
    double expo_min = -7;
    double expo_max = -3;
    
    evolution.density    = N_slider;
    evolution.T          = std::pow(10, expo_min*(1-e) + expo_max*e);
    evolution.sigma_coef = S_slider*.01;
    
    gngt.process(nb_threads,
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

    // Let us label the connected components

    vq3::utils::clear_vertex_efficiencies(g, true); // All vertices are considered for connected components.
    g.foreach_edge([thresh = E_slider*E_slider*1e-6](graph::ref_edge ref_e) {
	  auto extr = ref_e->extremities();
	  if(vq3::invalid_extremities(extr))
	    ref_e->kill();
	  else {
	    auto& A =  (*(extr.first) )().vq3_value;
	    auto& B =  (*(extr.second))().vq3_value;
	    // Long edge are not considered for connected components.
	    (*ref_e)().vq3_efficient = vq3::demo2d::d2(A,B) < thresh;
	  }
      });
    auto components = vq3::connected_components::make(g);
    vq3::labelling::conservative(components.begin(), components.end());
    vq3::labelling::edges_from_vertices(g);
    
    
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
    g.foreach_vertex([&speed_image, &speed_frame, &color_of_label](graph::ref_vertex ref_v) {
	auto& vertex = (*ref_v)();
	if(vertex.vq3_smoother.get<1>()) // If speed is available
	  cv::circle(speed_image, speed_frame(vertex.vq3_smoother.get<1>().value()), 3, color_of_label(vertex.vq3_label), -1);
      });
								
    
    
    cv::imshow("image", image);
    selector.build_image(video_data.image);
    cv::imshow("video", selector.image);
    cv::imshow("speed", speed_image);
    keycode = cv::waitKey(1) & 0xFF;
  }
  
  return 0;
}
