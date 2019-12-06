#include <vector>
#include <iterator>
#include <algorithm>

#include <opencv2/opencv.hpp>
#include <vq3.hpp>
#include <vq3demo.hpp>


/* 
   1D SOMs are known to approximate the travelling salesman problem
   (TSP). The effect is much more relavant when the AD som is fed with
   values whose topology reflects the "visitable" places of the input
   space. The visitable places are covered by a graph, which gives
   this space its topology. This graph is the auxiliary graph for GIT
   values. The whole SOM computes in that topological space. 
*/




namespace aux { // This contains the definition for the auxiliary graph.
  
  using sample    = vq3::demo2d::Point;
  using prototype = vq3::demo2d::Point;
  
  //                                                         ## Node properties :
  using vlayer_0 = prototype;                                // prototypes are 2D points (this is the "user defined" value).
  using vlayer_1 = vq3::decorator::path::shortest<vlayer_0>; // This holds informations built by dijkstra.
  using vertex   = vlayer_1;
  
  //                                                        ## Edge properties :
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
      opt_cost = vq3::demo2d::d(pt1, pt2);
    }
    return *opt_cost;
  }

  // The distance between a node and a sample.
  double d2(const vertex& v, const sample& p) {return vq3::demo2d::d2(v.vq3_value, p);}

  // The linear interpolation
  sample interpolate(const sample& a, const sample& b, double lambda) {return (1-lambda)*a + lambda*b;}

  // The shortest path function
  void shortest_path(graph& g, graph::ref_vertex start, graph::ref_vertex dest) {
    vq3::path::a_star<false, false>(g, start, dest, edge_cost,
    				    [start](const graph::ref_vertex& ref_v){ // This is the heuristic
    				      if(start) return vq3::demo2d::d((*start)().vq3_value, (*ref_v)().vq3_value); 
    				      else      return 0.0;
    				    });
  }

  // This is the traits type for building up values related to the
  // support graph. The function used inside decltype(...) is only
  // usefull (and very convenient) for this "traits" type definition.
  using traits = decltype(vq3::topology::gi::traits_val<sample, graph>(d2, interpolate, shortest_path));
}

namespace som { // This namespace defines the SOM graph
  using sample    = vq3::topology::gi::Value<aux::traits>;          // The samples' topology is dependent on the support graph.
  using prototype = sample;                                         // prototypes and samples are the same here.
  
  //                                                                ## Node properties :
  using layer_0 = prototype;                                        // 
  using layer_1 = vq3::decorator::tagged<layer_0>;                  // we add a tag for topology computation.
  using vertex  = layer_1;
  
  using graph  = vq3::graph<vertex, void>;
  using topology_key_type = int;
}


void rebuild_support_graph(aux::graph& g_aux,
			   std::mt19937& random_device,
			   vq3::demo2d::sample::density& density,
			   unsigned int nb_threads);

// Simulation parameters.

#define SLIDER_INIT               1000
#define NB_SAMPLES_PER_M2_SUPPORT 5000
#define NB_SOM_VERTICES            100 
#define ALPHA                       .1

int main(int argc, char* argv[]) {
  
  std::random_device rd;  
  std::mt19937 random_device(rd());
  
  unsigned int nb_threads = std::thread::hardware_concurrency(); // We use parallelisme for CHL here.

  
  
  /////
  //
  // Setting up the support graph
  //
  /////
  
  aux::graph g_aux; // this is the graph

  // We sill build it up from a density with a fancy shape.

  // wire shapes thickness
  double thickness = .1;
  double fig_pos   = .60;
  double up_bar    = .70;
    
  // Left square
  double i1             = 1;
  double w1             = 1;
  double h1             = 1;
  vq3::demo2d::Point o1 = {-fig_pos, 0};
  auto shape1           = vq3::demo2d::sample::rectangle(w1, h1, i1) + o1;
    
  // Right crown
  double i2             =  1;
  double r2             = .5;
  double h2             = r2 - thickness;
  vq3::demo2d::Point o2 = { fig_pos, 0};
  auto shape2 = (vq3::demo2d::sample::disk(r2, i2) - vq3::demo2d::sample::disk(h2, i2)) + o2;
  
  // Bars
  
  double i3               =  1;
  vq3::demo2d::Point min3 = {-fig_pos - thickness*.5, .5};
  vq3::demo2d::Point max3 = {-fig_pos + thickness*.5,  up_bar + thickness*.5};
  auto shape3             = vq3::demo2d::sample::rectangle(min3, max3, i3);
  
  double i4               =  1;
  vq3::demo2d::Point min4 = {fig_pos - thickness*.5, .5};
  vq3::demo2d::Point max4 = {fig_pos + thickness*.5,  up_bar + thickness*.5};
  auto shape4             = vq3::demo2d::sample::rectangle(min4, max4, i4);
  
  double i5               =  1;
  vq3::demo2d::Point min5 = {-fig_pos, up_bar - thickness*.5};
  vq3::demo2d::Point max5 = { fig_pos, up_bar + thickness*.5};
  auto shape5             = vq3::demo2d::sample::rectangle(min5, max5, i5);

  // Bridge

  
  double i6             =  1;
  double w6             =  .4;
  double h6             =  thickness;
  auto shape6           = vq3::demo2d::sample::rectangle(w6, h6, i6);
  
  
  // All
  auto density = shape1 || shape2 || shape3 || shape4 || shape5 || shape6;

  // See the code below.
  rebuild_support_graph(g_aux, random_device, density, nb_threads);  

  /////
  //
  // Setting up the SOM (1D)
  //
  /////

  som::graph g_som;

  auto traits = vq3::topology::gi::traits<aux::sample>(g_aux, aux::d2, aux::interpolate, aux::shortest_path);

  // This tosses a random value in the distribution bounding box.
    
  // We build up a ring of random graph-induced values.
  auto prev  = g_som += vq3::topology::gi::value(traits, density->bbox().uniform(random_device));
  auto curr  = prev;
  auto first = prev;
  for(unsigned int i = 1; i < NB_SOM_VERTICES; ++i, prev = curr) {
    curr = g_som += vq3::topology::gi::value(traits, density->bbox().uniform(random_device));
    g_som.connect(prev, curr);
  }
  g_som.connect(prev, first);

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
  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  auto image = cv::Mat(750, 1500, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .3*image.size().width, true);

  int old_slider = -1;
  int slider = SLIDER_INIT;
  cv::createTrackbar("100 * kernel radius", "image", &slider, 5000, nullptr);

  auto draw_edge_aux = vq3::demo2d::opencv::edge_drawer<aux::graph::ref_edge>(image, frame,
  									      [](const aux::vertex& v1, const aux::vertex& v2, const aux::edge& e) {return true;}, // always draw
  									      [](const aux::vertex& v)  {return                                     v.vq3_value;}, // position
  									      [](const aux::edge&   e)  {return                       cv::Scalar(255, 210, 210);}, // color
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
            << "Press space to toggle the bridge." << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl;
  
  int keycode = 0;
  while(keycode != 27) {

    if(slider != old_slider) {
      // We change the WTM kernel.
      topology.declare_distance(0, // This the key of the single neighborhood declared in this example.
				[radius = slider * .01](unsigned int edge_distance) {return std::max(0., 1 - edge_distance / radius);},
				slider * .01,
				1e-3);
      topology.update_full();
      old_slider = slider;
    }

    auto sample_point = vq3::demo2d::sample::get_one_sample(random_device, density);
    vq3::online::wtm::learn(topology, 0,
			    [](som::vertex& vertex) -> som::prototype& {return vertex.vq3_value;},
			    som_d2, vq3::topology::gi::value(traits, sample_point), ALPHA); // Our sample is a GIT value.

    
    
    image = cv::Scalar(255, 255, 255);
    
    g_aux.foreach_edge(draw_edge_aux);

    g_som.foreach_edge(draw_edge_som);
    g_som.foreach_vertex(draw_vertex_som);

    
    cv::imshow("image", image);
    keycode = cv::waitKey(1) & 0xFF;
    if(keycode == 32) {
      i6 = 1 - i6;
      rebuild_support_graph(g_aux, random_device, density, nb_threads);
      // The auxiliary graph has changed. We need to reset all the GIT
      // values that we handle. Indeed, they are the g_som vertices.
      g_som.foreach_vertex([](auto& ref_v){(*ref_v)().vq3_value.reset();});
    }
  }

  
  return 0;
}


void rebuild_support_graph(aux::graph& g_aux,
			   std::mt19937& random_device,
			   vq3::demo2d::sample::density& density,
			   unsigned int nb_threads) {
  
  g_aux.foreach_vertex([](auto ref_v){ref_v->kill();});
  
  // Setting vertices of the support graph
  {
    auto sampler_triangles = vq3::demo2d::sample::base_sampler::triangles(random_device, NB_SAMPLES_PER_M2_SUPPORT);
    auto S                 = vq3::demo2d::sample::sample_set(random_device, sampler_triangles, density);
    for(auto pt : S)       g_aux += pt;
  }
  
  // Setting edges of the support graph with CHL.
  {
    std::vector<vq3::demo2d::Point> S;
    auto sampler_triangles = vq3::demo2d::sample::base_sampler::triangles(random_device, 5*NB_SAMPLES_PER_M2_SUPPORT);
    auto S_                = vq3::demo2d::sample::sample_set(random_device, sampler_triangles, density);
    auto out               = std::back_inserter(S);
    std::copy(S_.begin(), S_.end(), out);
    
    auto chl = vq3::epoch::chl::processor(g_aux);
    chl.process(nb_threads,
                S.begin(), S.end(),
                [](const vq3::demo2d::Point& s) {return s;},      // Gets the sample from *it.
                aux::d2,                                          // d2(prototype, sample).
                aux::edge());                                     // New edge initialization value.
  }
}
