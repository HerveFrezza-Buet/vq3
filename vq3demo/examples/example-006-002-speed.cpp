#include <vq3demo.hpp>
#include <random>
#include <algorithm>
#include <string>

#define ANGLE_PERIOD 500
#define D_ANGLE (360./ANGLE_PERIOD)
#define SPEED_TO_METER .5

#define EVOLUTION_MARGIN_ABOVE        .20
#define EVOLUTION_MARGIN_BELOW        .30
#define EVOLUTION_TOPOLOGICAL_RATIO   .15

#define GNGT_ALPHA                    .05
#define GNGT_NB_SAMPLES_PER_PROTOTYPE  10
#define GNGT_NB_WTA_1                   5
#define GNGT_NB_WTA_2                   2
#define GNGT_NB_WTA_3                   0


#define SOM_H_RADIUS                  5.1
#define SOM_MAX_DIST                  (unsigned int)(SOM_H_RADIUS)
#define NARROW_SOM_COEF               .02
#define AVERAGE_RADIUS                5

#define FIXED_FRAME_DELAY             .03

#define N_SLIDER_INIT  300
#define T_SLIDER_INIT  500
#define Z_SLIDER_INIT 2000


// Graph definition
//
///////////////

using sample    = vq3::demo2d::Point;
using prototype = vq3::demo2d::Point;

//                                                                                 ## Node properties :
using vlayer_0 = prototype;                                                        // prototypes are 2D points (this is the "user defined" value).
using vlayer_1 = vq3::decorator::tagged<vlayer_0>;                                 // we add a tag for topology computation.
using vlayer_2 = vq3::decorator::smoother<vlayer_1, vq3::demo2d::Point, 1, 21, 2>; // we smooth of prototypes.
using vertex   = vlayer_2;

//                                                        ## Edge properties :
using elayer_0 = vq3::decorator::tagged<void>;            // we add a tag for CHL computation.
using edge     = elayer_0;

using graph    = vq3::graph<vertex, edge>;

using neighbour_key_type = std::string;


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
  
  int N_slider = N_SLIDER_INIT;
  int T_slider = T_SLIDER_INIT;
  int Z_slider = Z_SLIDER_INIT;
  
  
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
  cv::createTrackbar("nb/m^2",              "image", &N_slider, 1000, nullptr);
  cv::createTrackbar("T",                   "image", &T_slider, 1000, nullptr);
  
  auto image = cv::Mat(600, 800, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .1*image.size().width, true);


  cv::namedWindow("speed", CV_WINDOW_AUTOSIZE);
  cv::createTrackbar("zoom", "speed", &Z_slider, 5000, nullptr);
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
  
  auto topology = vq3::topology::table<neighbour_key_type>(g);
  topology.declare_distance("wide som",   [](unsigned int edge_distance) {return std::max(0., 1 - edge_distance/double(SOM_H_RADIUS));},   SOM_MAX_DIST, 1e-3);
  topology.declare_distance("narrow som", [](unsigned int edge_distance) {return edge_distance == 0 ? 1 : NARROW_SOM_COEF ;           },              1,  0.0);
  topology.declare_distance("avg",        [](unsigned int edge_distance) {return 1;                                                   }, AVERAGE_RADIUS,  0.0);
  
  auto gngt = vq3::algo::gngt::processor<sample>(topology);
  
  auto evolution = vq3::algo::gngt::by_default::evolution();

  gngt.alpha              = GNGT_ALPHA;
  gngt.samples_per_vertex = GNGT_NB_SAMPLES_PER_PROTOTYPE;
  gngt.nb_wta_1           = GNGT_NB_WTA_1;
  gngt.nb_wta_2           = GNGT_NB_WTA_2;
  gngt.nb_wta_3           = GNGT_NB_WTA_3;
  
  evolution.margin_above  = EVOLUTION_MARGIN_ABOVE;
  evolution.margin_below  = EVOLUTION_MARGIN_BELOW;
  evolution.topo_ratio    = EVOLUTION_TOPOLOGICAL_RATIO;
  
  
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

  unsigned int step = 0;
  unsigned int mode = 0;
  
  auto sampler = vq3::demo2d::sample::base_sampler::random(random_device, N_slider);
  
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
    
    sampler = N_slider;
    auto S_ = vq3::demo2d::sample::sample_set(random_device, sampler, density);
    S.clear();
    std::copy(S_.begin(), S_.end(), std::back_inserter(S));

    // Step


    double e = T_slider/1000.0;
    double expo_min = -2;
    double expo_max =  0;
    
    evolution.density    = N_slider;
    evolution.T          = std::pow(10, expo_min*(1-e) + expo_max*e);

    // We compute the topology evolution of the graph...
    gngt.process(nb_threads,
		 S.begin(), S.end(),                                                    // The sample set. Shuffle if the dataser is not sampled randomly.
		 [](const sample& s) {return s;},                                       // get sample from *iter (identity here).
		 [](vertex& v) -> prototype& {return v.vq3_value;},                     // get a prototype reference from the vertex value.
		 [](const prototype& p) {return p + vq3::demo2d::Point(-1e-5,1e-5);},   // get a point close to a prototype.
		 dist,                                                                  // The squared distance, faster, used for bmu-related stuff.
		 "wide som", "narrow som", "avg",                                       // Neighborhood keys.
		 evolution,
		 true);
    
    // Temporal update
    g.foreach_vertex([](graph::ref_vertex ref_v) {
	auto& value = (*ref_v)();
	value.vq3_smoother += value.vq3_value;
	value.vq3_smoother.set_timestep(FIXED_FRAME_DELAY);
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
