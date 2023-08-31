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
using sample    = demo2d::Point;
using prototype = demo2d::Point;
using vertex    = prototype;
using graph     = vq3::graph<vertex, void>; 

// This defines the data used during the epochs in order to collect
// and gather information from threads.
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
double dist2(const vertex& v, const demo2d::Point& p) {return demo2d::d2(v, p);}


// Main
//
////////////////
#define NB_SAMPLES_COEF 1e4
#define NB_SAMPLES_PER_M2 NB_SAMPLES_COEF*1.6 // for hexagrid
unsigned int NB_SAMPLES; 
unsigned int NB_VERTICES;   // The k of k-means...
#define MAX_DIST2   1e-10  // If a prototypes changes less that this squared distance, it is considered as constant.

#define NB_BINS 500
#define MAX_HISTO 5.

int main(int argc, char* argv[]) {
  if(argc < 5) {
    std::cout << "Usage : " << argv[0] << " <uniform|unbalanced> <rectangle|multidim|multidensity> <uniform|hexagonal|resample> nb_threads [i1 | [i2 | ...] ]" << std::endl
	      << "    i1 i2 ... : Successive steps where a snaphot is taken. if i1=-1, a snapshot is taken at each step." << std::endl
	      << "                If snapshots are asked, the last step is recordered systematically." << std::endl
	      << "    uniform|unbalanced : the initialization of prototypes." << std::endl
	      << "    rectangle, .... : The shape of the distribution." << std::endl
	      << "    uniform|hexagrid|resample : " << std::endl
	      << "      - uniform : samples are taken uniformly." << std::endl
	      << "      - hexagonal: samples are taken in an hewagonal mesh." << std::endl
	      << "      - resample: uniform re-sample of the dataset at each Lloyd iteration." << std::endl
	      << std::endl;
      
    return 0;
  }

  bool uniform = std::string(argv[1]) == std::string("uniform");
  std::string distrib {argv[2]};
  bool resample = std::string(argv[3]) == std::string("resample");
  bool hexagrid = std::string(argv[3]) == std::string("hexagonal");
  unsigned int nb_threads = std::atoi(argv[4]);

  std::set<unsigned int> snapshots;
  bool snap_all = false;

  if(argc > 5) {
    int i = std::atoi(argv[5]);
    if(i < 0) snap_all = true;
    else 
      for(int arg = 5; arg < argc; ++arg)  
	snapshots.insert((unsigned int)(std::atoi(argv[arg])));
  }
  
  
  std::random_device rd;  
  std::mt19937 random_device(rd());

  graph g;


  // Image settings
  //
  ///////////////////

  
  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  auto image       = cv::Mat(720, 1024, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = demo2d::opencv::direct_orthonormal_frame(image.size(), .325*image.size().width, true);
  auto dd          = demo2d::opencv::dot_drawer<demo2d::Point>(image, frame,
									 [](const demo2d::Point& pt) {return                      true;},
									 [](const demo2d::Point& pt) {return                        pt;},
									 [](const demo2d::Point& pt) {return                         1;},
									 [](const demo2d::Point& pt) {return cv::Scalar(200, 200, 200);},
									 [](const demo2d::Point& pt) {return                        -1;});
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

  demo2d::Point offset {0., 0.5};
  
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

  double w3        = .1;
  double h3        = .1;
  auto source      = demo2d::sample::rectangle(w3, h3, i) + p1 + offset;

  auto multidim = (rect || bar || crown) + offset;

  auto bb = multidim->bbox();
  auto min_pt = bb.bottom_left();
  auto max_pt = bb.top_right();
  auto rectangle = demo2d::sample::rectangle(min_pt, max_pt, i);

  auto f = [](auto pt) {
    if (pt.x < 0) return -pt.x;
    if (pt.x < .8) return 1.0;
    return .2;
  };

  auto multidensity = demo2d::sample::custom(bb, f);
  
  
  demo2d::sample::density density;
  
  
  if(distrib == "multidim") {
    density     = multidim;
    NB_SAMPLES  = 2*NB_SAMPLES_COEF;
    NB_VERTICES = 200;
  }
  else if(distrib == "rectangle") {
    density     = rectangle;
    NB_SAMPLES  = 5*NB_SAMPLES_COEF;
    NB_VERTICES = 300;
  }
  else if(distrib == "multidensity") {
    density     = multidensity;
    NB_SAMPLES  = 5*NB_SAMPLES_COEF;
    NB_VERTICES = 300;
  } 
  else {
    std::cout << distrib << " is not a valid distribution name" << std::endl;
    ::exit(0);
  }
    
  
  // Some initializations
  //
  ///////////////////

  
  vq3::demo2d::opencv::histogram histo {{-1.5, -1.}, {1.5, -.2}};
  histo.frame_margin = .1;
  histo.title        = "local distortions";
  
  // We need to register the input samples in a vector since we want
  // to both use and display them.
  std::vector<demo2d::Point> S;
  auto out = std::back_inserter(S);

  if(hexagrid) {
    auto sampler = demo2d::sample::base_sampler::triangles(random_device, NB_SAMPLES_PER_M2);
    auto SS = demo2d::sample::sample_set(random_device, sampler, density);
    std::copy(SS.begin(), SS.end(), out);
  }
  else
    for(unsigned int i = 0; i < NB_SAMPLES; ++i)
      *out++ = demo2d::sample::get_one_sample(random_device, density);

  std::cout << std::endl
	    << "Using " << S.size() << " samples." << std::endl
	    << std::endl;
  
  // Let us generate the graph as random unconnected vertices, taken from the distribution.
  
  for(unsigned int i=0; i < NB_VERTICES; ++i)
    if(uniform)
      g += demo2d::sample::get_one_sample(random_device, density);
    else
      g += demo2d::sample::get_one_sample(random_device, source);

  bool stop = false;
  unsigned int step = 0;
  double disto;
  double sqrt_disto;

  // This is an structure that stores vertex-related computation. It
  // can be shared by several processors, this is why it is allocated
  // first and then passed to each processor. 
  auto topology = vq3::topology::table<topology_key_type>(g);
  topology.update(); // We update the topology without considering edge-based neighborhood.
  auto wta = vq3::epoch::wta::processor(topology);

  std::string prefix = "wta";
  if(uniform) prefix += "-uniform";
  else        prefix += "-unbalanced";
  auto filename = demo::videoframe_name(prefix, "png");


  std::string data_filename;
  std::ofstream data(prefix+".data");
  if(!data) {
    std::cerr << "cannot write \"" << data_filename << "\"." << std::endl;
    return 1;
  }

  if(resample)
    std::cout << std::endl
	      << "We are in resample mode, so press ESC" << std::endl
	      << "to interrupt when you consider the convergence" << std::endl
	      << "is achieved." << std::endl
	      << std::endl;
  
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

    if(resample) {
      S.clear();
      for(unsigned int i = 0; i < NB_SAMPLES; ++i)
	*out++ = demo2d::sample::get_one_sample(random_device, density);
    }
    
    auto t_start = std::chrono::high_resolution_clock::now();
    auto epoch_result = wta.process<epoch_data>(nb_threads,
						S.begin(), S.end(), 
						[](const demo2d::Point& s) {return s;}, // Gets the sample from *it.
						[](vertex& v) -> vertex& {return v;},   // Gets the prototype ***reference*** from the vertex value.
						dist2);                                 // dist2(prototype, sample).
    auto t_end = std::chrono::high_resolution_clock::now();

    auto hout = histo.output_iterator();
    for(auto& data : epoch_result)
      if(data.vq3_bmu_accum.nb > 0) 
	*(hout++) = data.vq3_bmu_accum.value;
    histo.set_bins(0., MAX_HISTO, NB_BINS);
    histo.make();

    
    int duration =  std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
    std::cout << "Step duration (" << nb_threads << " threads) : " << std::setw(6) << duration << " ms.    \r" << std::flush;

    // Let us check the convergence for each vertex.
    stop = true;
    for(auto& d : epoch_result)
      if(dist2(d.vq3_previous_prototype, d.vq3_current_prototype) > MAX_DIST2) {
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
    histo.draw(image, frame);
    histo.vline(image, frame, 2.5, {0., 0., 0.}, 1);

    if(snap_all)
      cv::imwrite(filename(), image);
    else if(stop)
      cv::imwrite(filename("final"), image);
    else if(snapshots.find(step) != snapshots.end())
      cv::imwrite(filename(step), image);
    
    cv::imshow("image", image);
    if((cv::waitKey(1) & 0xFF) == 27) stop = true;
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
