#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <iostream>
#include <iomanip>
#include <chrono>


// This example shows how to perform winner-take-all. We will
// implement a k-means (not LBG, see next example for that). We will
// use multi-threading.



// Graph definition
//
///////////////

using vertex = vq3::demo2d::Point;
using graph  = vq3::graph<vertex, void>; 

// This defines the data used during the epochs in order to collect
// and gather information from threads.
using epoch_data_0 = vq3::epoch::data::none<vq3::demo2d::Point>;         // This is the root of the stack, the sample type has to be provided.
using epoch_data_1 = vq3::epoch::data::wta_check<epoch_data_0, vertex>;  // This gathers computation for batch winner-take-all and provides arguments for convergence test.
using epoch_data   = epoch_data_1;


// Distance
//
////////////////


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double dist2(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v, p);}


// Main
//
////////////////

#define NB_SAMPLES_PER_M2   50000
#define NB_VERTICES           200  // The k of k-means...
#define MAX_DIST2           1e-10  // If a prototypes changes less that this squared distance, it is considered as constant.


int main(int argc, char* argv[]) {
  if(argc != 3) {
    std::cout << "Usage : " << argv[0] << " <uniform|unbalanced> nb_threads" << std::endl;
    return 0;
  }

  bool uniform = std::string(argv[1]) == std::string("uniform");
  unsigned int nb_threads = std::atoi(argv[2]);
  
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

  double w3             = .1;
  double h3             = .1;
  auto source           = vq3::demo2d::sample::rectangle(w3, h3, i) + p1;
  

  auto density = rect || bar || crown;
  
  // Some initializations
  //
  ///////////////////
  
  // Let us generate the graph as random unconnected vertices, taken from the distribution.
  
  for(unsigned int i=0; i < NB_VERTICES; ++i)
    if(uniform)
      g += vq3::demo2d::sample::get_one_sample(random_device, density);
    else
      g += vq3::demo2d::sample::get_one_sample(random_device, source);
 
  // We need to register the input samples in a vector since we want
  // to both use and display them.
  std::vector<vq3::demo2d::Point> S;
  
  auto S_ = vq3::demo2d::sample::sample_set(random_device, density, NB_SAMPLES_PER_M2);
  std::copy(S_.begin(), S_.end(), std::back_inserter(S));

  bool stop = false;
  unsigned int step = 0;

  // This is an structure that stores vertex-related computation. It
  // can be shared by several processors, this is why it is allocated
  // first and then passed to each processor. Here, we have only one
  // processor, so this external allocation may sound artificial.
  auto vertices = vq3::utils::vertices(g);
  vertices.update_topology(g);
  auto wta = vq3::epoch::wta::processor(g, vertices);
  
  // This is the loop
  //
  //////////////////
  
  
  while(!stop) {
    ++step;

    // Vertex update
    
    auto t_start = std::chrono::high_resolution_clock::now();
    auto epoch_result = wta.update_prototypes<epoch_data>(nb_threads,
							  S.begin(), S.end(), 
							  [](const vq3::demo2d::Point& s) {return s;}, // Gets the sample from *it.
							  [](vertex& v) -> vertex& {return v;},        // Gets the prototype ***reference*** from the vertex value.
							  dist2);                                      // dist2(prototype, sample).
    auto t_end = std::chrono::high_resolution_clock::now();

    
    int duration =  std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
    std::cout << "Step duration (" << nb_threads << " threads) : " << std::setw(6) << duration << " ms.    \r" << std::flush;

    // Let us check the convergence for each vertex.
    stop = true;
    for(auto& d : epoch_result)
      if(dist2(d.wq3_wta_previous_prototype, d.wq3_wta_current_prototype) > MAX_DIST2) {
	stop = false;
	break;
      }
				  
    
    // Display
    
    image = cv::Scalar(255, 255, 255);
    std::copy(S.begin(), S.end(), dd);
    g.foreach_edge(draw_edge); 
    g.foreach_vertex(draw_vertex);
    cv::imshow("image", image);
    cv::waitKey(1);
    
  }

  

  std::cout << std::endl
	    << std::endl
	    << std::endl
	    << "##########" << std::endl
	    << std::endl
	    << "Convergence reached after " << step << " steps. Press any key to end." << std::endl
	    << std::endl
	    << "##########" << std::endl
	    << std::endl;
  cv::imshow("image", image);
  cv::waitKey(0);

  return 0;
}
