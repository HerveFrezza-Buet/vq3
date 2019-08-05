#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <iterator>
#include <cstdlib>
#include <thread>

/*

  Here, we build up and run a Kohonen seolf-organizing map. We use vq3
  tools for handling measures and statistics.

*/

#define NB_SAMPLES_PER_M2   20000
#define GRID_WIDTH             30
#define GRID_HEIGHT            30
#define NB_STEPS              250

#define SOM_H_RADIUS            5.1
#define SOM_MAX_DIST            (unsigned int)(SOM_H_RADIUS)

using sample    = vq3::demo2d::Point;
using prototype = vq3::demo2d::Point;

//                                                                ## Node properties :
using layer_0 = prototype;                                        // prototypes are 2D points (this is the "user defined" value).
using layer_1 = vq3::decorator::tagged<layer_0>;                  // we add a tag for topology computation.
using layer_2 = vq3::demo::decorator::colored<layer_1>;           // we add a color to the nodes.
using vertex  = layer_2;

using graph  = vq3::graph<vertex, void>;

// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v.vq3_value, p);}

using epoch_data_0 = vq3::epoch::data::none<sample, vertex, prototype>; // This is the root of the stack.
using epoch_data_1 = vq3::epoch::data::wtm<epoch_data_0>;               // This gathers computation for batch winner-take-most.
using epoch_data_2 = vq3::epoch::data::bmu<epoch_data_1,
					   vq3::epoch::data::bmu_sqrt_dist_accum<prototype,
										 sample>>;  // This gathers computation for the BMU (sum of d1(sample, proto))
using epoch_data   = epoch_data_2;

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

  // Let us build the graph.
  //
  //////////////

  graph g;
  
  vq3::utils::make_grid(g, GRID_WIDTH, GRID_HEIGHT,
			[&random_device](unsigned int w, unsigned int h) {
			  graph::vertex_value_type value(vq3::demo2d::uniform(random_device, {-.5, -.5}, {.5, .5}));
			  value.vq3_color = cv::Scalar(0, 255*w/double(GRID_WIDTH-1), 255*h/double(GRID_HEIGHT-1));
			  return value;
			});

  // Let us build the dataset.
  //
  //////////////

  
  double intensity = 1. ;
  double radius    =  .5;
  double hole      =  .3;
  auto density     = vq3::demo2d::sample::disk(radius, intensity) - vq3::demo2d::sample::disk(hole, intensity);

  auto sampler = vq3::demo2d::sample::base_sampler::random(random_device, NB_SAMPLES_PER_M2);
  auto S_      = vq3::demo2d::sample::sample_set(random_device, sampler, density); // This re-toss points at each begin...end iteration.
  
  std::vector<vq3::demo2d::Point> S; // Let us use a single sample of S_.
  auto out = std::back_inserter(S);
  std::copy(S_.begin(), S_.end(), out);

  

  // Image settings
  //
  ///////////////////

  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  auto image       = cv::Mat(480, 640, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .9*image.size().height, true);
  auto dd          = vq3::demo2d::opencv::dot_drawer<vq3::demo2d::Point>(image, frame,
									 [](const vq3::demo2d::Point& pt) {return                      true;},
									 [](const vq3::demo2d::Point& pt) {return                        pt;},
									 [](const vq3::demo2d::Point& pt) {return                         1;},
									 [](const vq3::demo2d::Point& pt) {return cv::Scalar(200, 200, 200);},
									 [](const vq3::demo2d::Point& pt) {return                        -1;});
  auto draw_edge   = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
								       [](const vertex& v1, const vertex& v2) {return true;}, // always draw
								       [](const vertex& v)     {return         v.vq3_value;}, // position
								       []()                    {return cv::Scalar(0, 0, 0);}, // color
								       []()                    {return                  1;}); // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return        true;},  // always draw
									   [](const vertex& v) {return v.vq3_value;},  // position
									   [](const vertex& v) {return           3;},  // radius
									   [](const vertex& v) {return v.vq3_color;},  // color
									   [](const vertex& v) {return          -1;}); // thickness

  // Gnuplot charts
  //
  ///////////////////

  auto disto_chart = vq3::demo::gnuplot::curve("distortion evolution", "disto");
  disto_chart.set_term_pdf();
  disto_chart.set_xlabel("steps");
  disto_chart.set_ylabel("distortion");
  disto_chart.set_yrange(0, 0.25);
  auto out_disto = disto_chart.output_iterator();

  auto disto_histo = vq3::demo::gnuplot::histogram("distortion distribution", "histo");
  disto_histo.set_xlabel("distortion values");
  disto_histo.set_term_pdf();
  disto_histo.set_sci(.5);
  disto_histo.set_bins(0, 0.004, 20);
  disto_histo.set_yrange(0, 150);
  
  
  // The is the loop....
  //
  ///////////////////

  image = cv::Scalar(255, 255, 255);
  std::copy(S.begin(), S.end(), dd);
  g.foreach_edge(draw_edge); 
  g.foreach_vertex(draw_vertex);
  cv::imshow("image", image);
  cv::waitKey(1000);

  // First, we need a structure for handling SOM-like computation. The
  // template argument is the type of input samples. For the sake of
  // function notation homogeneity, let us use the alias algo::som
  // indtead of epoch::wtm.
  auto topology = vq3::topology::table<int>(g);
  topology.declare_distance(0,
			    [](unsigned int edge_distance) {return std::max(0., 1 - edge_distance/double(SOM_H_RADIUS));},
			    SOM_MAX_DIST,
			    1e-3); // We consider node and edge-based neihborhoods.
  topology.update_full();
  auto som = vq3::algo::som::processor(topology);


  std::cout << std::endl;

  std::vector<epoch_data> epoch_results;
  for(unsigned int step = 0; step < NB_STEPS; ++step) {
    // Learning : the returned value of this function is ignored here. See next examples.
    epoch_results = som.process<epoch_data>(nb_threads, 0,
					    S.begin(), S.end(),
					    [](const vq3::demo2d::Point& p) -> const vq3::demo2d::Point& {return p;},
					    [](vertex& vertex_value) -> vq3::demo2d::Point& {return vertex_value.vq3_value;},
					    d2);
    
    // since our epoch data stack has a bmu layer, we can get the
    // accumulator handling the distortion for all nodes. The distance
    // we use is the squared Euclidian distance (d2). The distortion
    // is taken from the epoch data of all vertices.
    auto distortion_accum = vq3::epoch::distortion(epoch_results.begin(), epoch_results.end());
    // Let us plot the square root of the average d2.
    *(out_disto++) = {double(step), std::sqrt(distortion_accum.average<>())};

    std::cout << "step " << std::setw(3) << step+1 << '/' << NB_STEPS << "   \r" << std::flush;
    
    image = cv::Scalar(255, 255, 255);
    std::copy(S.begin(), S.end(), dd);
    g.foreach_edge(draw_edge);
    g.foreach_vertex(draw_vertex);
    cv::imshow("image", image);
    cv::waitKey(1);
  }
  std::cout << std::endl << std::endl;
  // Let us generate the gnuplot files.
  
  disto_chart.make();

  auto out_histo = disto_histo.output_iterator();
  for(auto& data : epoch_results) 
    if(data.vq3_bmu_accum.nb > 0)
      *(out_histo++) = data.vq3_bmu_accum.value;
  disto_histo.make();
  
  

  return 0;
}
