#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <set>
#include <cmath>


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
using epoch_data_2 = vq3::epoch::data::bmu<epoch_data_1>;                // This enables the computation of distortion. This is not mandatory for finding the prototype positions.         
using epoch_data   = epoch_data_2;


// Distance
//
////////////////


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double dist2(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v, p);}


// Main
//
////////////////

#define NB_SAMPLES          20000
#define NB_VERTICES           200  // The k of k-means...
#define MAX_DIST2           1e-10  // If a prototypes changes less that this squared distance, it is considered as constant.


int main(int argc, char* argv[]) {
  if(argc < 3) {
    std::cout << "Usage : " << argv[0] << " <uniform|unbalanced> nb_threads [i1 | [i2 | ...] ]" << std::endl
	      << "    i1 i2 ... : Successive steps where a snaphot is taken. if i1=-1, a snapshot is taken at each step." << std::endl
	      << "                If snapshots are asked, the last step is recordered systematically." << std::endl;
    return 0;
  }

  bool uniform = std::string(argv[1]) == std::string("uniform");
  unsigned int nb_threads = std::atoi(argv[2]);

  std::set<unsigned int> snapshots;
  bool snap_all = false;

  if(argc > 3) {
    int i = std::atoi(argv[3]);
    if(i < 0) snap_all = true;
    else 
      for(int arg = 3; arg < argc; ++arg)  
	snapshots.insert((unsigned int)(std::atoi(argv[arg])));
  }
  
  
  std::mt19937 random_device(0);

  graph g;


  // Image settings
  //
  ///////////////////

  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  auto image       = cv::Mat(360, 1024, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .325*image.size().width, true);
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
 
  // We need to register the input samples in a vector since we want
  // to both use and display them.
  std::vector<vq3::demo2d::Point> S;
  auto out = std::back_inserter(S);
  for(unsigned int i = 0; i < NB_SAMPLES; ++i)
    *out++ = vq3::demo2d::sample::get_one_sample(random_device, density);
  
  // Let us generate the graph as random unconnected vertices, taken from the distribution.
  
  for(unsigned int i=0; i < NB_VERTICES; ++i)
    if(uniform)
      g += vq3::demo2d::sample::get_one_sample(random_device, density);
    else
      g += vq3::demo2d::sample::get_one_sample(random_device, source);

  bool stop = false;
  unsigned int step = 0;
  double disto;
  double sqrt_disto;

  // This is an structure that stores vertex-related computation. It
  // can be shared by several processors, this is why it is allocated
  // first and then passed to each processor. Here, we have only one
  // processor, so this external allocation may sound artificial.
  auto vertices = vq3::utils::vertices(g);
  vertices.update_topology(g);
  auto wta = vq3::epoch::wta::processor(g, vertices);

  std::string prefix = "wta";
  if(uniform) prefix += "-uniform";
  else        prefix += "-unbalanced";
  auto filename = vq3::demo::videoframe_name(prefix, "png");


  std::string data_filename;
  std::ofstream data(prefix+".data");
  if(!data) {
    std::cerr << "cannot write \"" << data_filename << "\"." << std::endl;
    return 1;
  }
  
  
  // This is the loop
  //
  //////////////////
  
  image = cv::Scalar(255, 255, 255);
  std::copy(S.begin(), S.end(), dd);
  g.foreach_edge(draw_edge); 
  g.foreach_vertex(draw_vertex);
    
  if(snap_all)
    cv::imwrite(filename(), image);
  else if(snapshots.find(0) != snapshots.end())
    cv::imwrite(filename(0), image);
  cv::imshow("image", image);
  cv::waitKey(10);
  
  while(!stop) {

    // Vertex update
    ++step;
    
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
    
    // Let us compute and save the distortion
    auto distortion_accum = vq3::epoch::distortion(epoch_result.begin(), epoch_result.end());
    disto      = distortion_accum.average<>();
    sqrt_disto = std::sqrt(disto);
    data << step << ' ' << disto << ' ' << sqrt_disto << std::endl;

    // Display
    
    image = cv::Scalar(255, 255, 255);
    std::copy(S.begin(), S.end(), dd);
    g.foreach_edge(draw_edge); 
    g.foreach_vertex(draw_vertex);

    if(snap_all)
      cv::imwrite(filename(), image);
    else if(stop)
      cv::imwrite(filename("final"), image);
    else if(snapshots.find(step) != snapshots.end())
      cv::imwrite(filename(step), image);
    
    cv::imshow("image", image);
    cv::waitKey(10);
    
  }

  

  cv::imshow("image", image);

  if(snap_all || snapshots.size() > 0)
    cv::waitKey(2000);
  else {
    std::cout << std::endl
	      << std::endl
	      << std::endl
	      << "##########" << std::endl
	      << std::endl
	      << "Convergence reached after " << step << " steps. Press any key to end." << std::endl
	      << std::endl
	      << "##########" << std::endl
	      << std::endl;
    cv::waitKey(0);
  }

  

  
  // This file can be included in your latex code.
  std::string latex_filename;
  std::ofstream latex(prefix+".tex");
  if(!latex) {
    std::cerr << "cannot write \"" << latex_filename << "\"." << std::endl;
    return 1;
  }

  std::string varname = "\\Wta";
  if(uniform) varname += "Uniform";
  else        varname += "Unbalanced";
      
  latex << "% generated by " << argv[0] << std::endl
	<< "\\newcommand{" << varname << "FinalStep}[0]{" << step << '}' << std::endl
	<< "\\newcommand{" << varname << "FinalDisto}[0]{" << std::fixed << std::setprecision(3) << disto*1000 << "\\times 10^{-3}}" << std::endl
	<< "\\newcommand{" << varname << "FinalSqrtDisto}[0]{" << std::fixed << std::setprecision(3) << sqrt_disto*100 << "\\times 10^{-2}}" << std::endl
	<< "\\newcommand{" << varname << "NbSamples}[0]{" << NB_SAMPLES << '}' << std::endl
	<< "\\newcommand{" << varname << "NbPrototypes}[0]{" << NB_VERTICES << '}' << std::endl;
  

  return 0;
}
