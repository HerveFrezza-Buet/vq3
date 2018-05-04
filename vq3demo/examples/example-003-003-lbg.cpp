#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <iostream>
#include <iomanip>
#include <chrono>


// This example shows how to perform a k-means, i.e the
// Linde-Buzo-Gray (LBG) algorithm.



// Graph definition
//
///////////////

using vertex = vq3::demo2d::Point;
using graph  = vq3::graph<vertex, void>;

using sample    = vq3::demo2d::Point;
using prototype = vq3::demo2d::Point;


// Distance
//
////////////////


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double dist2(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v, p);}

// Statistics
//
////////////////

// We need an epoch data stack in order to compute statistics.
using epoch_data_0 = vq3::epoch::data::none<sample>;          // This is the root of the stack, the sample type has to be provided.
using epoch_data_1 = vq3::epoch::data::bmu<epoch_data_0>;     // We collect distortion for each best-matching vertex.
using epoch_data   = epoch_data_1;


// Main
//
////////////////

#define NB_SAMPLES_PER_M2   50000
#define MAX_DIST2           1e-10  // If a prototypes changes less that this squared distance, it is considered as constant.


int main(int argc, char* argv[]) {
  if(argc != 3) {
    std::cout << "Usage : " << argv[0] << " nb_threads k" << std::endl
	      << "try k=500" << std::endl;
    return 0;
  }

  unsigned int nb_threads = std::atoi(argv[1]);
  unsigned int K          = std::atoi(argv[2]);
  
  std::random_device rd;  
  std::mt19937 random_device(rd());

  graph g;


  // Image settings
  //
  ///////////////////

  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  auto image       = cv::Mat(768, 1024, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .2*image.size().width, true);
  auto dd          = vq3::demo2d::opencv::dot_drawer<vq3::demo2d::Point>(image, frame,
									 [](const vq3::demo2d::Point& pt) {return                      true;},
									 [](const vq3::demo2d::Point& pt) {return                        pt;},
									 [](const vq3::demo2d::Point& pt) {return                         1;},
									 [](const vq3::demo2d::Point& pt) {return cv::Scalar(200, 200, 200);},
									 [](const vq3::demo2d::Point& pt) {return                        -1;});
  auto draw_edge   = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
								       [](const vertex& v1, const vertex& v2) {return true;},  // always draw
								       [](const vertex& v)   {return                     v;},  // position
								       []()                  {return cv::Scalar(255, 0, 0);},  // color
								       []()                  {return                     3;}); // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return                true;},  // always draw
									   [](const vertex& v) {return                   v;},  // position
									   [](const vertex& v) {return                   3;},  // radius
									   [](const vertex& v) {return cv::Scalar(0, 0, 0);},  // color
									   [](const vertex& v) {return                  -1;}); // thickness
   
  // Distribution settings
  //
  ///////////////////

  
  double thickness = .2;
  double sep       =  1;
  double i         =  1;

  double w              = 1;
  double h              = 1;
  vq3::demo2d::Point p1 = {-sep, 0};
  auto rect             = vq3::demo2d::sample::rectangle(w, h, i) + p1;

  double R              = .5;
  double r              = R - thickness;
  vq3::demo2d::Point p2 = {sep, 0};
  auto crown            = (vq3::demo2d::sample::disk(R, i) - vq3::demo2d::sample::disk(r, i)) + p2;

  vq3::demo2d::Point A = {-sep + w/2, -thickness/2};
  vq3::demo2d::Point B = { sep - r,    thickness/2};
  auto bar             = vq3::demo2d::sample::rectangle(A, B, i);

  auto density = rect || bar || crown;
  
  // Some initializations
  //
  ///////////////////
  
 
  // We need to register the input samples in a vector since we want
  // to both use and display them.
  std::vector<vq3::demo2d::Point> S;
  
  auto S_ = vq3::demo2d::sample::sample_set(random_device, density, NB_SAMPLES_PER_M2);
  std::copy(S_.begin(), S_.end(), std::back_inserter(S));

  
  // This is the algorithm
  //
  //////////////////
  
  
  auto t_start = std::chrono::high_resolution_clock::now();
  vq3::algo::lbg<prototype>(random_device,
			    nb_threads, g, K,
			    S.begin(), S.end(),
			    [](const vq3::demo2d::Point& s) {return s;},                    // Gets the sample from *it.
			    [](vertex& v) -> vertex& {return v;},                           // Gets the prototype ***reference*** from the vertex value.
			    dist2,                                                          // dist2(prototype, sample).
			    [&random_device](const vq3::demo2d::Point& proto) {             // Makes a prototype nearly similar to proto.
			      return vq3::demo2d::alter(random_device, proto, MAX_DIST2);},  
			    [](const vertex& previous, const vertex& current) {             // check function : if check(previous, current) is false for each vertex, the convergence is achieved.
			      return dist2(previous, current) > MAX_DIST2;
			    },
			    true);                                                          // Verbosity             
  auto t_end = std::chrono::high_resolution_clock::now();

    
  int duration =  std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
  std::cout << std::endl
	    << "convergence duration (" << nb_threads << " threads, " << S.size() << " samples, k = " << K << ") : " << std::setw(6) << duration << " ms." << std::endl
	    << std::endl
	    << "press any key to exit." << std::endl;

  // Let is display the result
  //
  //////////////////
  
  image = cv::Scalar(255, 255, 255);
  std::copy(S.begin(), S.end(), dd);
  g.foreach_edge(draw_edge); 
  g.foreach_vertex(draw_vertex);
  cv::imshow("image", image);
  cv::waitKey(0);

  // Now, let us plot the histogram of local distortions. This is why
  // we have built up an epoch data stack.<graph::ref_vertex>
  auto vertices = vq3::utils::vertices(g);
  auto wta = vq3::epoch::wta::processor(g, vertices);
  vertices.update_topology(g);
  auto epoch_result = wta.update_prototypes<epoch_data>(nb_threads,
							S.begin(), S.end(), 
							[](const vq3::demo2d::Point& s) {return s;}, // Gets the sample from *it.
							[](vertex& v) -> vertex& {return v;},        // Gets the prototype ***reference*** from the vertex value.
							dist2);                                      // dist2(prototype, sample).
  
  auto disto_histo = vq3::demo::gnuplot::histogram("distortion distribution", "histo");
  disto_histo.set_xlabel("distortion values");
  disto_histo.set_term_pdf();
  disto_histo.set_sci(.75);
  disto_histo.set_bins(90,20); // 20 bins covers the smallest interval containing 90% of the values.
  auto out_histo = disto_histo.output_iterator();
  for(auto& data : epoch_result) 
    if(data.vq3_bmu_accum.nb > 0)
      *(out_histo++) = data.vq3_bmu_accum.value;
  disto_histo.make();


  return 0;
}
