#include <vector>
#include <iterator>
#include <algorithm>
#include <optional>

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
  
  //                                                         ## Vertex properties :
  using vlayer_0 = prototype;                                // prototypes are 2D points (this is the "user defined" value).
  using vlayer_1 = vq3::decorator::path::shortest<vlayer_0>; // This holds informations built by dijkstra.
  using vlayer_2 = vq3::demo::decorator::colored<vlayer_1>;  // We add a color to the vertices.
  using vertex   = vlayer_2;
  
  //                                                                  ## Edge properties :
  using elayer_0 = vq3::decorator::cost<void, std::optional<double>>; // Edge cost.
  using elayer_1 = vq3::decorator::tagged<elayer_0>;                  // We need a tag for CHL.
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
  double d(const sample& p1, const sample& p2) {return demo2d::d(p1, p2);}

  // This is the distance between a node value and a sample.
  double D(const vertex& v, const sample& p) {return demo2d::d(v.vq3_value, p);}

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
  using traits = decltype(vq3::topology::gi::traits_val<sample, graph>(d, D, interpolate, shortest_path,
								       vq3::path::travel_defaults::compute_cumulated_costs<graph::ref_vertex>,
								       vq3::path::travel_defaults::cumulated_cost_of<graph::ref_vertex>));
}


// This namespace defines the k-means graph

namespace kmeans {
  using sample    = vq3::topology::gi::Value<aux::traits>;   // The samples' topology is dependent on the support graph.
  using prototype = sample;                                  // prototypes and samples are the same here.
  
  //                                                         ## Node properties :
  using vlayer_0 = prototype;                                // prototypes are 2D points sensitive to the auxiliary graph topology.
  using vlayer_1 = vq3::demo::decorator::colored<vlayer_0>;  // We add a color to the vertices.
  using vertex   = vlayer_1;
  using graph    = vq3::graph<vertex, void>;
}

// Let measure the distortion thanks to epochs computation.

using epoch_data_0 = vq3::epoch::data::none<kmeans::sample, kmeans::vertex, kmeans::prototype>; // This is the root of the stack.
using epoch_data_1 = vq3::epoch::data::bmu<epoch_data_0,
					   vq3::epoch::data::bmu_sqrt_dist_accum<kmeans::prototype,kmeans::sample>>;
using epoch_data   = epoch_data_1;


#define NB_SAMPLES_PER_M2_SUPPORT   500
#define K                            10
#define CHUNK_SIZE                  100
#define ALPHA_SLIDER_INIT           100

// Execution mode
enum class Mode : char {Cont = 'c', Step = 's'};

int main(int argc, char* argv[]) {
  
  std::random_device rd;  
  std::mt19937 random_device(rd());
  
  unsigned int nb_threads = std::thread::hardware_concurrency(); // We use parallelisme for CHL only here.

  /////
  //
  // Setting up the support graph
  //
  /////
  
  aux::graph g_aux; // this is the graph
  
  // We sill build it up from a density with a fancy shape.
  double intensity     =  1;
  double thickness     = .2;
  double radius        = .5;
  double hole_radius   = radius - thickness;
  double break_width   = .1;
  double break_height  = 3*thickness;
  double bar_width     = .75*thickness;
  double bar_height    = .5;
  double bbar_width    = .9*bar_height;
  
  auto crown        = demo2d::sample::disk(radius, intensity) - demo2d::sample::disk(hole_radius, intensity);
  auto middle_break = demo2d::sample::rectangle(break_width, break_height, intensity);
  auto bar          = demo2d::sample::rectangle(bar_width, bar_height, intensity);
  auto bbar         = demo2d::sample::rectangle(bbar_width, bar_width, intensity);

  demo2d::Point up      = {0., hole_radius + .5*thickness};
  demo2d::Point barpos  = {bar_width, bar_height*.5};
  demo2d::Point bbarpos = {0, bar_height*.9};
  
  // All
  auto density = (((crown + up) || (crown - up)) - middle_break) || (bar + barpos) || (bar - barpos)|| (bbar + bbarpos) || (bbar - bbarpos);
  
  // Setting vertices of the support graph
  {
    auto sampler_triangles = demo2d::sample::base_sampler::triangles(random_device, NB_SAMPLES_PER_M2_SUPPORT);
    auto S                 = demo2d::sample::sample_set(random_device, sampler_triangles, density);
    for(auto pt : S) {
      auto ref_v = g_aux += pt;
      (*ref_v)().vq3_color = cv::Scalar(180, 180, 180);
    }
      
  }
  
  // Setting edges of the support graph with CHL.
  {
    std::vector<demo2d::Point> S;
    auto sampler_triangles = demo2d::sample::base_sampler::triangles(random_device, 5*NB_SAMPLES_PER_M2_SUPPORT);
    auto S_                = demo2d::sample::sample_set(random_device, sampler_triangles, density);
    auto out               = std::back_inserter(S);
    std::copy(S_.begin(), S_.end(), out);
    
    auto chl = vq3::epoch::chl::processor(g_aux);
    chl.process(nb_threads,
                S.begin(), S.end(),
                [](const demo2d::Point& s) {return s;},      // Gets the sample from *it.
                aux::D,                                     // D(prototype, sample).
                aux::edge());                                // New edge initialization value.
    
  }
  
  /////
  //
  // Setting up the k_means
  //
  /////
  
  kmeans::graph g_kmeans;
  auto traits = vq3::topology::gi::traits<aux::sample>(g_aux, aux::d, aux::D, aux::interpolate, aux::shortest_path,
						       vq3::path::travel_defaults::compute_cumulated_costs<aux::graph::ref_vertex>,
						       vq3::path::travel_defaults::cumulated_cost_of<aux::graph::ref_vertex>);

  {
    auto colormap = demo2d::opencv::colormap::random(random_device, .1, 50);
    for(unsigned int i = 0; i < K; ++i) {
      auto ref_v = g_kmeans += vq3::topology::gi::value(traits, demo2d::sample::get_one_sample(random_device, density));
      (*ref_v)().vq3_color = colormap(i);
    }
  }
  
  // We will need a distance for selecting the closest prototype. It
  // is easily available from the traits instance.
  auto kmeans_d = vq3::topology::gi::distance<kmeans::vertex>(traits,
							      [](const kmeans::vertex& vertex) -> const kmeans::prototype& {return vertex.vq3_value;});

  ////
  //
  // Set up the distortion computation
  //
  ////

  auto topology   = vq3::topology::table<int>(g_kmeans);
  topology.update(); // We do not need a full update since we have no edges. Moreover, a full update would have required a "tag" decoration of the vertices.
  auto distortion = vq3::epoch::wta::processor(topology); 

  // This is the dataset used to measure the distortion. We use the quxiliary graph vertices positions.
  std::vector<kmeans::sample> S;
  auto out = std::back_inserter(S);
  g_aux.foreach_vertex([&out, &traits](const auto& ref_v){*(out++) = vq3::topology::gi::value(traits, (*ref_v)().vq3_value);});

  /////
  //
  // Setting up the display
  //
  /////

  int slider_alpha = ALPHA_SLIDER_INIT;
  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  cv::createTrackbar("1000 * alpha", "image", &slider_alpha, 500, nullptr);
  auto image = cv::Mat(600, 350, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = demo2d::opencv::direct_orthonormal_frame(image.size(), .8*image.size().width, true);

  auto draw_vertex_aux = vq3::demo2d::opencv::vertex_drawer<aux::graph::ref_vertex>(image, frame,
  										    [](const aux::vertex& v) {return        true;},  // always draw
  										    [](const aux::vertex& v) {return v.vq3_value;},  // position
  										    [](const aux::vertex& v) {return           3;},  // radius
  										    [](const aux::vertex& v) {return v.vq3_color;},  // color	
  										    [](const aux::vertex& v) {return          -1;}); // thickness

  auto draw_edge_aux = vq3::demo2d::opencv::edge_drawer<aux::graph::ref_edge>(image, frame,
  									      [](const aux::vertex& v1, const aux::vertex& v2, const aux::edge& e) {return true;}, // always draw
  									      [](const aux::vertex& v)  {return                                     v.vq3_value;}, // position
  									      [](const aux::edge&   e)  {return                       cv::Scalar(210, 210, 210);}, // color
  									      [](const aux::edge&   e)  {return                                              1;}); // thickness

  auto draw_vertex_kmeans_bg = vq3::demo2d::opencv::vertex_drawer<kmeans::graph::ref_vertex>(image, frame,
  											     [](const kmeans::vertex& v) {return          true;},  // always draw
  											     [](const kmeans::vertex& v) {return v.vq3_value();},  // position
  											     [](const kmeans::vertex& v) {return             5;},  // radius
  											     [](const kmeans::vertex& v) {return   v.vq3_color;},  // color
  											     [](const kmeans::vertex& v) {return            -1;}); // thickness
  
  auto draw_vertex_kmeans_fg = vq3::demo2d::opencv::vertex_drawer<kmeans::graph::ref_vertex>(image, frame,
  											     [](const kmeans::vertex& v) {return                true;},  // always draw
  											     [](const kmeans::vertex& v) {return       v.vq3_value();},  // position
  											     [](const kmeans::vertex& v) {return                   5;},  // radius
  											     [](const kmeans::vertex& v) {return cv::Scalar(0, 0, 0);},  // color
  											     [](const kmeans::vertex& v) {return                   2;}); // thickness
  
  /////
  //
  // Running loop
  //
  /////

  std::cout << std::endl
	    << std::endl
	    << "ESC           - quit." << std::endl
	    << "return key    - restart." << std::endl
	    << "space         - one step." << std::endl
	    << "c             - running mode." << std::endl
	    << std::endl;

  Mode mode = Mode::Step;
  bool compute = true;
  int wait_ms;
  int keycode = 0;
  while(keycode != 27 /* ESC */) {

    
    compute = mode != Mode::Step;
    
    if(keycode == 32 || keycode == 10) { // space or return key pressed.
      mode = Mode::Step;
      compute = true;
    }
    else if(keycode == 99) {   // 'c' key pressed.
      mode = Mode::Cont;
      compute = true;
    }
    
    wait_ms = 500;
    if(compute) {
      wait_ms = 1;

      if(keycode == 10) { // return key was pressed, we reset the prototypes
	auto colormap = demo2d::opencv::colormap::random(random_device, .1, 50);
	unsigned int i=0;
	g_kmeans.foreach_vertex([&random_device, &density, &colormap, &i](kmeans::graph::ref_vertex ref_v){
	    (*ref_v)().vq3_value = demo2d::sample::get_one_sample(random_device, density); // GIT values can be initialized from a regular value.
	    (*ref_v)().vq3_color = colormap(i++);
	  });
      }
      else // We compute the usual k-mean update.
	for(unsigned int i = 0; i < CHUNK_SIZE; ++i) {
	  auto sample_point = demo2d::sample::get_one_sample(random_device, density);
	  vq3::online::wta::learn(g_kmeans,
	  			  [](kmeans::vertex& vertex_value) -> kmeans::prototype& {return vertex_value.vq3_value;},
	  			  kmeans_d, vq3::topology::gi::value(traits, sample_point),
				  slider_alpha*.001); // Our sample is a GIT value.
	}

      // Let us colorize the auxiliary graph according to the color of the closest vertex in the kmeans graph.
      g_aux.foreach_vertex([&g_kmeans, &kmeans_d, &traits](auto& ref_v) {
	  auto closest = vq3::utils::closest(g_kmeans, vq3::topology::gi::value(traits, (*ref_v)().vq3_value), kmeans_d);
	  (*ref_v)().vq3_color = (*closest)().vq3_color;
	});

      // Let us measure the distortion. Indeed, we collect for each
      // k-mean vertex the average distance with the closest
      // sample. Then, we average this for all vertices.
      
      auto epoch_results = distortion.process<epoch_data>(1, // A single thread is mandatory here since GIT is not thread safe.
      							  S.begin(), S.end(),
      							  [](const kmeans::sample& s) -> const kmeans::sample& {return s;},                        // Retrieves the sample from *it
      							  [](kmeans::vertex& vertex_value) -> kmeans::prototype& {return vertex_value.vq3_value;}, // Retrieves the prototype from the vertex value.
      							  kmeans_d);
      double average = 0;
      for(auto& data : epoch_results)
      	average += data.vq3_bmu_accum.average();
      average /= epoch_results.size();
      std::cout << "average per vertex distortion : " << average << std::endl;
      
    }
    
    image = cv::Scalar(255, 255, 255);
    
    g_aux.foreach_edge(draw_edge_aux);
    g_aux.foreach_vertex(draw_vertex_aux);
    g_kmeans.foreach_vertex(draw_vertex_kmeans_bg);
    g_kmeans.foreach_vertex(draw_vertex_kmeans_fg);

    
    cv::imshow("image", image);
    keycode = cv::waitKey(wait_ms) & 0xFF;
  }

  
  return 0;
}
