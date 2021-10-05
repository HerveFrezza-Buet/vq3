#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <iostream>
#include <iomanip>
#include <chrono>

// This example shows how to perform multi-threaded competitive
// Hebbian learning... with the addition ouf counting features.



// Graph definition
//
///////////////

using vertex = demo2d::Point;

using elayer_0 = vq3::decorator::tagged<void>;        // We need tags on the edges.
using elayer_1 = vq3::decorator::chl_count<elayer_0>; // We add an atomic counter on the edges.
using edge     = elayer_1;

using graph  = vq3::graph<vertex, edge>;


// Distance
//
////////////////

// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const demo2d::Point& p) {return demo2d::d2(v, p);}

auto color_of_count(std::size_t counts, std::size_t max) {
  unsigned char level = (unsigned char)(255*(1-counts/double(max)+.5));
  return cv::Scalar(255, level, level);
}

// Main
//
////////////////

#define NB_VERTICES_PER_M2     50
#define NB_SAMPLES_PER_M2   10000


int main(int argc, char* argv[]) {
  if(argc != 2) {
    std::cout << "Usage : " << argv[0] << " nb_threads" << std::endl;
    return 0;
  }

  unsigned int nb_threads = std::atoi(argv[1]);
  
  std::random_device rd;  
  std::mt19937 random_device(rd());


  // Image settings
  //
  ///////////////////

  
  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  auto image       = cv::Mat(800, 800, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = demo2d::opencv::direct_orthonormal_frame(image.size(), .48*image.size().width, true);
  auto dd          = demo2d::opencv::dot_drawer<demo2d::Point>(image, frame,
									 [](const demo2d::Point& pt) {return                      true;},
									 [](const demo2d::Point& pt) {return                        pt;},
									 [](const demo2d::Point& pt) {return                         1;},
									 [](const demo2d::Point& pt) {return cv::Scalar(200, 200, 200);},
									 [](const demo2d::Point& pt) {return                        -1;});
  std::size_t max_count = 0;
  auto draw_edge   = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
								       [](const vertex& v1, const vertex& v2, const edge& e) {return                                       true;},  // always draw
								       [](const vertex& v)                                   {return                                          v;},  // position
								       [&max_count](const edge&   e)                         {return color_of_count(e.vq3_chl_count, max_count);},  // color
								       [](const edge&   e)                                   {return                                          3;}); // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return                true;},  // always draw
									   [](const vertex& v) {return                   v;},  // position
									   [](const vertex& v) {return                   3;},  // radius
									   [](const vertex& v) {return cv::Scalar(0, 0, 0);},  // color
									   [](const vertex& v) {return                  -1;}); // thickness

  
  // Initializations
  //
  //////////////////

  double side = 2;
  double intensity = 1.0;
  auto background = demo2d::sample::rectangle(side, side, intensity);

  graph g;
  auto v_sampler = demo2d::sample::base_sampler::random(random_device, NB_VERTICES_PER_M2);
  for(auto v : demo2d::sample::sample_set(random_device, v_sampler, background)) g += v;

  
  auto samples = demo2d::sample::custom(demo2d::sample::BBox(-1, -1, 2, 2), [](const demo2d::Point& p) {return .5*(p.x + 1.0);});
  std::vector<demo2d::Point> S;
  auto s_sampler = demo2d::sample::base_sampler::random(random_device, NB_SAMPLES_PER_M2);
  auto S_ = demo2d::sample::sample_set(random_device, s_sampler, samples);
  auto out = std::back_inserter(S);
  std::copy(S_.begin(), S_.end(), out);

  // Computation
  //
  //////////////

  auto chl = vq3::epoch::chl::processor_with_counts(g);
  chl.process(nb_threads,
	      S.begin(), S.end(),
	      [](const demo2d::Point& s) {return s;}, // Gets the sample from *it.
	      d2,                                     // d2(prototype, sample).
	      edge());                                // New edge initialization value.

  g.foreach_edge([&max_count](auto ref_e) {
		   auto extr = ref_e->extremities();
		   if(vq3::invalid_extremities(extr))
		     ref_e->kill();
		   else if(std::size_t c = (*ref_e)().vq3_chl_count; c > max_count)
		     max_count = c;
		 });
    
  
  // Display
  //
  //////////
    
  
  image = cv::Scalar(255, 255, 255);
  std::copy(S.begin(), S.end(), dd);
  g.foreach_edge(draw_edge); 
  g.foreach_vertex(draw_vertex);
  cv::imshow("image", image);
  cv::waitKey();
    
  return 0;
}
  
