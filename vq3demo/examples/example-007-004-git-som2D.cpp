#include <vector>
#include <iterator>
#include <algorithm>

#include <opencv2/opencv.hpp>
#include <vq3.hpp>
#include <vq3demo.hpp>

/*

  This examples show that processors cannot really be exploited with
  GIT values. Indeed, each distance computation involves a A*
  computation on the support graph. On that graph, vertices are tagged
  during A*, so many A* instances cannot run in parallel since they
  cannot benefit from a dedicated tagging process (tags are shared).

  This example shows how to perform with such constraints.

*/

// This contains the definition for the auxiliary graph.

namespace aux { 
  
  using sample    = demo2d::Point;
  using prototype = demo2d::Point;
  
  //                                                         ## Node properties :
  using vlayer_0 = prototype;                                // prototypes are 2D points (this is the "user defined" value).
  using vlayer_1 = vq3::decorator::path::shortest<vlayer_0>; // This holds informations built by dijkstra.
  using vlayer_2 = vq3::demo::decorator::colored<vlayer_1>;  // We add a color to the vertices.
  using vertex   = vlayer_2;
  
  //                                                        ## Node properties :
  using elayer_0 = vq3::decorator::optional_cost<void>;     // Edge cost.
  using elayer_1 = vq3::decorator::tagged<elayer_0>;        // We need a tag for CHL.
  using edge     = elayer_1;
  
  
  using graph = vq3::graph<vertex, edge>;

  // We compute the edge cost as being its length. As we use an optional
  // cost, we can check if the computation is necessary.
  double edge_cost(const graph::ref_edge& ref_e) {
    auto& opt_cost = (*ref_e)().vq3_cost; 
    if(!opt_cost) { // if cost not calculated yet.
      auto extr_pair = ref_e->extremities();
      const auto& pt1 = (*(extr_pair.first))().vq3_value;
      const auto& pt2 = (*(extr_pair.second))().vq3_value;
      opt_cost = demo2d::d(pt1, pt2);
    }
    return *opt_cost;
  }

  // This is the distance between a node value and a sample.
  double d2(const sample& p1, const sample& p2) {return demo2d::d2(p1, p2);}

  // This is the distance between a node value and a sample.
  double D2(const vertex& v, const sample& p) {return demo2d::d2(v.vq3_value, p);}

  // The linear interpolation
  sample interpolate(const sample& a, const sample& b, double lambda) {return (1-lambda)*a + lambda*b;}

  // The shortest path function
  void shortest_path(graph& g, graph::ref_vertex start, graph::ref_vertex dest) {
    vq3::path::a_star<false, false>(g, start, dest, edge_cost,
    				    [start](const graph::ref_vertex& ref_v){ // This is the heuristic
    				      if(start) return demo2d::d((*start)().vq3_value, (*ref_v)().vq3_value); 
    				      else      return 0.0;
    				    });
  }

  // This is the traits type for building up values related to the
  // support graph. The function used inside decltype(...) is only
  // usefull (and very convenient) for this "traits" type definition.
  using traits = decltype(vq3::topology::gi::traits_val<sample, graph>(d2, D2, interpolate, shortest_path));
}



// This namespace defines the SOM graph

namespace som { 
  using sample    = vq3::topology::gi::Value<aux::traits>;          // The samples' topology is dependent on the support graph.
  using prototype = sample;                                         // prototypes and samples are the same here.
  
  //                                                                ## Node properties :
  using layer_0 = prototype;                                        // 
  using layer_1 = vq3::decorator::tagged<layer_0>;                  // we add a tag for topology computation.
  using vertex  = layer_1;
  
  using graph  = vq3::graph<vertex, void>;
  using topology_key_type = int;
}

#define SLIDER_INIT                1000
#define NB_SAMPLES_PER_M2_SUPPORT 20000
#define GRID_WIDTH                   30
#define GRID_HEIGHT                   5
#define ALPHA                        .1
#define CHUNK_SIZE                   20


int main(int argc, char* argv[]) {
  
  std::random_device rd;  
  std::mt19937 random_device(rd());
  
  unsigned int nb_threads = std::thread::hardware_concurrency(); // We use parallelisme for CHL here.

  /////
  //
  // Setting up the input density
  //
  /////
  
  double intensity = 1. ;
  double radius    =  .5;
  double hole      =  .3;
  auto density     = demo2d::sample::disk(radius, intensity) - demo2d::sample::disk(hole, intensity);

  
  /////
  //
  // Setting up the support graph
  //
  /////
  
  aux::graph g_aux; // this is the graph

  
  // Setting vertices of the support graph
  {
    auto sampler_random = demo2d::sample::base_sampler::random(random_device, NB_SAMPLES_PER_M2_SUPPORT);
    auto S              = demo2d::sample::sample_set(random_device, sampler_random, density);
    for(auto pt : S)    g_aux += pt;
  }

  
  // Setting edges of the support graph with CHL.
  // (a Delaunay triangulation would have been better but it is not available in vq3)
  {
    std::vector<demo2d::Point> S;
    auto sampler_triangles = demo2d::sample::base_sampler::triangles(random_device, 20*NB_SAMPLES_PER_M2_SUPPORT);
    auto S_                = demo2d::sample::sample_set(random_device, sampler_triangles, density);
    auto out               = std::back_inserter(S);
    std::copy(S_.begin(), S_.end(), out);
    
    auto chl = vq3::epoch::chl::processor(g_aux);
    chl.process(nb_threads,
                S.begin(), S.end(),
                [](const demo2d::Point& s) {return s;},      // Gets the sample from *it.
                aux::D2,                                     // D2(prototype, sample).
                aux::edge());                                // New edge initialization value.

    // We remove some edges afterwards.
    g_aux.foreach_edge([](auto& ref_e) {
	auto extr_pair = ref_e->extremities();           
	if(vq3::invalid_extremities(extr_pair)) {
	  ref_e->kill();
	  return;
	}
	auto& A = (*(extr_pair.first))().vq3_value;
	auto& B = (*(extr_pair.second))().vq3_value;
	if(A.y > 0 && B.y > 0 && (A.x * B.x < 0))
	  ref_e->kill();
      });
  }



  /////
  //
  // Setting up the SOM (2D)
  //
  /////

  som::graph g_som;

  auto traits = vq3::topology::gi::traits<aux::sample>(g_aux, aux::d2, aux::D2, aux::interpolate, aux::shortest_path);
  
  vq3::utils::make_grid(g_som, GRID_WIDTH, GRID_HEIGHT,
			[&random_device, &traits, bbox = density->bbox()](unsigned int w, unsigned int h) {
			  return vq3::topology::gi::value(traits, bbox.uniform(random_device));
			});

  
  // This is the topology for WTM computation.
  auto topology = vq3::topology::table<som::topology_key_type>(g_som);
  
  // We will need a distance for selecting the closest prototype. It
  // is easily available from the traits instance.
  auto som_d2 = vq3::topology::gi::distance<som::vertex>(traits,
							 [](const som::vertex& vertex) -> const som::prototype& {return vertex.vq3_value;});

  
  /////
  //
  // Setting up the display
  //
  /////
  
  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  auto image = cv::Mat(800, 800, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = demo2d::opencv::direct_orthonormal_frame(image.size(), .9*image.size().width, true);

  int old_slider = -1;
  int slider = SLIDER_INIT;
  cv::createTrackbar("100 * kernel radius", "image", &slider, 3000, nullptr);


  auto draw_edge_aux = vq3::demo2d::opencv::edge_drawer<aux::graph::ref_edge>(image, frame,
  									      [](const aux::vertex& v1, const aux::vertex& v2, const aux::edge& e) {return true;}, // always draw
  									      [](const aux::vertex& v)  {return                                     v.vq3_value;}, // position
  									      [](const aux::edge&   e)  {return                       cv::Scalar(210, 210, 210);}, // color
  									      [](const aux::edge&   e)  {return                                              1;}); // thickness

  auto draw_vertex_som = vq3::demo2d::opencv::vertex_drawer<som::graph::ref_vertex>(image, frame,
  										    [](const som::vertex& v) {return                      true;},  // always draw
  										    [](const som::vertex& v) {return             v.vq3_value();},  // position
  										    [](const som::vertex& v) {return                         5;},  // radius
  										    [](const som::vertex& v) {return cv::Scalar(  0,   0, 180);},  // color
  										    [](const som::vertex& v) {return                        -1;}); // thickness

  auto draw_edge_som = vq3::demo2d::opencv::edge_drawer<som::graph::ref_edge>(image, frame,
  									      [](const som::vertex& v1, const som::vertex& v2) {return true;},  // always draw
  									      [](const som::vertex& v)  {return               v.vq3_value();},  // position
  									      []()                      {return   cv::Scalar(  0,   0,   0);},  // color
  									      []()                      {return                           3;}); // thickness
 

  /////
  //
  // Running loop
  //
  /////
  
  std::cout << std::endl
            << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl
            << "Press ESC to quit." << std::endl
            << "Press 'return' to reset the map." << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl;
  
  int keycode = 0;
  while(keycode != 27) {

    if(keycode == 10)
      g_som.foreach_vertex([&random_device, bbox=density->bbox()](auto& ref_v){(*ref_v)().vq3_value = bbox.uniform(random_device);});
    
    if(slider != old_slider) {
      // We change the WTM kernel.
      topology.declare_distance(0, // This the key of the single neighborhood declared in this example.
				[radius = slider * .01](unsigned int edge_distance) {return std::max(0., 1 - edge_distance / radius);},
				slider * .01,
				1e-3);
      topology.update_full();
      old_slider = slider;
    }

    for(unsigned int i = 0; i < CHUNK_SIZE; ++i) {
      auto sample_point = demo2d::sample::get_one_sample(random_device, density);
      vq3::online::wtm::learn(topology, 0,
			      [](som::vertex& vertex) -> som::prototype& {return vertex.vq3_value;},
			      som_d2, vq3::topology::gi::value(traits, sample_point), ALPHA); // Our sample is a GIT value.
    }
    
    
    image = cv::Scalar(255, 255, 255);
    
    g_aux.foreach_edge(draw_edge_aux);

    g_som.foreach_edge(draw_edge_som);
    g_som.foreach_vertex(draw_vertex_som);

    
    cv::imshow("image", image);
    keycode = cv::waitKey(1) & 0xFF;
  }

  
  return 0;
}
