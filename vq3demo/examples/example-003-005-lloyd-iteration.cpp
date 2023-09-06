
#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <set>
#include <cmath>
#include <limits>
#include <vector>

// This example illustrates step by step lloyd iteration.

// Graph definition
//
///////////////
using sample    = demo2d::Point;
using prototype = demo2d::Point;

using vlayer_0  = prototype;
using vlayer_1  = vq3::demo::decorator::colored<vlayer_0>; 
using vertex    = vlayer_1;
using graph     = vq3::graph<vertex, void>;


using epoch_data_0 = vq3::epoch::data::none<sample, vertex, prototype>;  // This is the root of the stack.
using epoch_data_1 = vq3::epoch::data::wta<epoch_data_0>;                // This gathers computation for batch winner-take-all.
using epoch_data_2 = vq3::epoch::data::delta<epoch_data_1>;              // This records variations of the prototype (for further convergence test).
using epoch_data_3 = vq3::epoch::data::bmu<epoch_data_2,
					   vq3::epoch::data::bmu_dist_accum<prototype,
									    sample>>;    // This enables the computation of distortion. This is not mandatory for finding the prototype positions.         
using epoch_data   = epoch_data_3;

using topology_key_type = int;


// Distance
//
////////////////


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double dist2(const vertex& v, const demo2d::Point& p) {return demo2d::d2(v.vq3_value, p);}


// Main
//
////////////////

#define NB_SAMPLES_PER_M2 5e2
#define NB_VERTICES 12

int main(int argc, char* argv[]) {
  std::random_device rd;  
  std::mt19937 random_device(rd());

  graph g;
  

  std::vector<cv::Scalar> color_of_sample;
  auto color_it = color_of_sample.begin();
  
  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  auto image       = cv::Mat(360, 1024, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = demo2d::opencv::direct_orthonormal_frame(image.size(), .325*image.size().width, true);
  auto dd          = demo2d::opencv::dot_drawer<demo2d::Point>(image, frame,
							       [](const demo2d::Point& pt) {return                      true;},
							       [](const demo2d::Point& pt) {return                        pt;},
							       [](const demo2d::Point& pt) {return                         3;},
							       [&color_it](const demo2d::Point& pt) {return    *(color_it++);},
							       [](const demo2d::Point& pt) {return                        -1;});
  
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return                true;},  // always draw
									   [](const vertex& v) {return         v.vq3_value;},  // position
									   [](const vertex& v) {return                   7;},  // radius
									   [](const vertex& v) {return         v.vq3_color;},  // color
									   [](const vertex& v) {return                  -1;}); // thickness
  
  auto draw_vertex_contour = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
										   [](const vertex& v) {return                true;},  // always draw
										   [](const vertex& v) {return         v.vq3_value;},  // position
										   [](const vertex& v) {return                   9;},  // radius
										   [](const vertex& v) {return  cv::Scalar(0, 0, 0);},  // color
										   [](const vertex& v) {return                  -1;}); // thickness



  // Distribution settings
  //
  ///////////////////

  
  double thickness = .2;
  double sep       =  1;
  double i         =  1;

  double w         = 1;
  double h         = 1;
  demo2d::Point p1 = {-sep, 0};
  auto rect        = demo2d::sample::rectangle(w, h, i) + p1;

  double R         = .5;
  double r         = R - thickness;
  demo2d::Point p2 = {sep, 0};
  auto crown       = (demo2d::sample::disk(R, i) - demo2d::sample::disk(r, i)) + p2;

  demo2d::Point A  = {-sep + w/2, -thickness/2};
  demo2d::Point B  = { sep - r,    thickness/2};
  auto bar         = demo2d::sample::rectangle(A, B, i);
  
  auto density = rect || bar || crown;
  
  std::vector<demo2d::Point> S;
  auto out = std::back_inserter(S);

  auto sampler = demo2d::sample::base_sampler::random(random_device, NB_SAMPLES_PER_M2);
  auto SS = demo2d::sample::sample_set(random_device, sampler, density);
  for(auto ss_it = SS.begin(); ss_it != SS.end();) {
    *out++ = *ss_it++;
    color_of_sample.push_back(cv::Scalar(200, 200, 200));
  }

  std::array<cv::Scalar, NB_VERTICES> colors {
    cv::Scalar(150,   0,   0),
    cv::Scalar(150, 150,   0),
    cv::Scalar(  0, 150,   0),
    cv::Scalar(  0, 150, 150),
    cv::Scalar(  0,   0, 150),
    cv::Scalar(150,   0, 150),
    cv::Scalar(255, 200, 200),
    cv::Scalar(255, 255, 200),
    cv::Scalar(200, 255, 200),
    cv::Scalar(200, 255, 255),
    cv::Scalar(200, 200, 255),
    cv::Scalar(255, 200, 255)};
  for(auto& color: colors) {
    auto ref_v = g += demo2d::sample::get_one_sample(random_device, density);
    (*ref_v)().vq3_color = color;
  }
  
  auto topology = vq3::topology::table<topology_key_type>(g);
  topology.update(); // We update the topology without considering edge-based neighborhood.
  auto wta = vq3::epoch::wta::processor(topology);
  
  // This is the loop
  //
  //////////////////
  auto stop = false;
  
  image = cv::Scalar(255, 255, 255);
  color_it = color_of_sample.begin();
  std::copy(S.begin(), S.end(), dd);
  g.foreach_vertex(draw_vertex_contour);
  g.foreach_vertex(draw_vertex);
  cv::imshow("image", image);
  if((cv::waitKey(0) & 0xFF) == 27) stop = true;
  
  while(!stop) {
    
    color_it = color_of_sample.begin();
    for(const auto& pt : S) {
      auto ref_v = vq3::utils::closest(g, pt, dist2);
      *(color_it++) = (*ref_v)().vq3_color;
    }
    
    image = cv::Scalar(255, 255, 255);
    color_it = color_of_sample.begin();
    std::copy(S.begin(), S.end(), dd);
    g.foreach_vertex(draw_vertex_contour);
    g.foreach_vertex(draw_vertex);
    cv::imshow("image", image);
    if((cv::waitKey(0) & 0xFF) == 27) stop = true;
    
    auto epoch_result = wta.process<epoch_data>(1,
						S.begin(), S.end(), 
						[](const demo2d::Point& s) {return s;},           // Gets the sample from *it.
						[](vertex& v) -> prototype& {return v.vq3_value;},// Gets the prototype ***reference*** from the vertex value.
						dist2);           
    image = cv::Scalar(255, 255, 255);
    color_it = color_of_sample.begin();
    std::copy(S.begin(), S.end(), dd);
    g.foreach_vertex(draw_vertex_contour);
    g.foreach_vertex(draw_vertex);
    cv::imshow("image", image);
    if((cv::waitKey(0) & 0xFF) == 27) stop = true;   
    
  }
  
  return 0;
}
