
#include <tuple>
#include <thread>
#include <random>
#include <vector>
#include <iterator>

#include <demo2d.hpp>
#include <vq3.hpp>
#include <vq3demo.hpp>

// Input distribution settings
#define NB_SAMPLES_PER_M2          7000
#define NOISE_RADIUS                1.2
#define NOISE_DENSITY               .33
#define BIG_RADIUS         NOISE_RADIUS
#define MEDIUM_RADIUS                .8
#define SMALL_RADIUS                 .4
#define ARC_THICKNESS                .1
#define ARC_SEP_THICKNESS            .1
#define BAR_OFFSET .02

// K-means settings
#define K 6
#define COLORMAP				\
  { cv::Scalar(255, 100, 100),			\
    cv::Scalar(100, 255, 100),			\
    cv::Scalar(100, 100, 255),			\
    cv::Scalar(100, 255, 255),			\
    cv::Scalar(255, 255, 100),			\
    cv::Scalar(255, 100, 255)}	
#define EDGE_THICKNESS .05                                         // The radius of the tube around the edge for counting samples.
#define NB_AROUND_PER_M2  (2*(EDGE_THICKNESS)*(NB_SAMPLES_PER_M2)) // The number of samples in a 1-length edge placed in a 1-density region.
#define CHL_FRACTION .1                                            // The fraction of samples used for teh auxiliary graph.

#define IMG_SIDE 800
#define ZOOM     .4


namespace aux {

  // This hosts cost information
  struct Cost {
  private:
    mutable std::optional<double> cost;
    
  public:
    double edge_length = 0;
    double relative_density = 0;

    Cost()                       = default;
    Cost(const Cost&)            = default;
    Cost& operator=(const Cost&) = default;

    // This computes the edge density from the count of samples in the surrounding tube.
    void operator=(std::size_t count) {
      relative_density = std::min(1.0, count/ (edge_length * NB_AROUND_PER_M2));
    }

#define D0  0
#define C0 1000
#define D1 .4
#define C1  900
#define D2 .6
#define C2    2
#define D3   1
#define C3    1
    // This is a piecewise linear function, transforming density into a dijkstra cost.
    // Low density implies very hight cost.
    double operator()() const {
      if(!cost) { // We set the cost if it has never been computed so far.
	double coef;
	if(relative_density < D1)
	  coef = C0 +  (C1 - C0) * relative_density/D1;
	else if(relative_density < D2)
	  coef = C1 + (C2 - C1) * (relative_density - D1)/(D2 - D1);
	else
	  coef = C2 + (C3 - C2) * (relative_density - D2)/(1 - D2);
     
	cost = edge_length*coef;
      }
      return *cost;
    }
  };
  
  using sample    = demo2d::Point;
  using prototype = sample;
  
  using vlayer_0 = prototype;        
  using vlayer_1 = vq3::decorator::path::shortest<vlayer_0>;       // We host Dijkstra information.
  using vlayer_2 = vq3::decorator::path::length<vlayer_1, double>; // We will host the path length, that will be distinct
  //                                                                  from the cumulated cost computed by Dijkstra.
  using vertex   = vlayer_2;
  
  using elayer_0 = vq3::decorator::cost<void, Cost>;   // Our specific cost is hosted here.
  using elayer_1 = vq3::decorator::tagged<elayer_0>;   // We need tags for CHL.
  using elayer_2 = vq3::decorator::counter<elayer_1>;  // We need a counter for counting the samples close to the edge.
  using edge     = elayer_2;
  
  using graph = vq3::graph<vertex, edge>;

  // These functions are required for auxiliary graph definitions and traits.
  double edge_cost(const graph::ref_edge& ref_e) {return (*ref_e)().vq3_cost();}
  double d(const sample& p1, const sample& p2)   {return demo2d::d(p1,          p2);}
  double D(const vertex& v, const sample& p)     {return demo2d::d(v.vq3_value, p );}
  sample interpolate(const sample& a, const sample& b, double lambda) {return (1-lambda)*a + lambda*b;}
  
  void shortest_path(graph& g, graph::ref_vertex start, graph::ref_vertex dest) {
    vq3::path::a_star<false, false>(g, start, dest, edge_cost,
				    [start](const graph::ref_vertex& ref_v){ // This is the A* heuristics.
				      if(start) return demo2d::d((*start)().vq3_value, (*ref_v)().vq3_value); 
				      else      return 0.0;
				    });
  }

  // This is for learning, once a path is found. We consider all edge costs as the edge length.
  void compute_cumulated_costs(graph::ref_vertex start, graph::ref_vertex dest) {
    // Cumulated costs are from destination to start, but the
    // available path iteration goes from start to destination. We
    // need a stack.
    std::stack<std::pair<graph::ref_vertex, double>> visited;
    auto end = vq3::path::begin(dest); ++end; // end is not reached, it is the iterator after dest.
    for(auto it  = vq3::path::begin(start); it != end; ++it) {
      auto ref_v = *it;
      if(ref_v == dest)
	visited.push({ref_v, 0});
      else 
	visited.push({ref_v, (*(it.get_edge()))().vq3_cost.edge_length});
    }

    // Now we accumate the lengths.
    double length = 0;
    while(!visited.empty()) {
      auto [ref_v, l] = visited.top();
      length += l;
      (*ref_v)().vq3_length = length;
      visited.pop();
    }
  }

  // This function says hox to get the cumulated cost for travelling paths when learning is performed.
  double cumulated_cost_of(graph::ref_vertex ref_v) {return (*ref_v)().vq3_length;}

  // Now, the traits can be defined.
  using traits = decltype(vq3::topology::gi::traits_val<sample, graph>(d, D, interpolate, shortest_path,
								       compute_cumulated_costs, cumulated_cost_of));

  // We measure the edge density from the counter.
  void compute_edge_relative_density(const graph::ref_edge& ref_e) {
    auto& edge_data = (*ref_e)();
    edge_data.vq3_cost = edge_data.vq3_counter; // Performs relative-density cost computation.
  }
  
}

// This is the k-means graph (no egdes there actually...)
namespace kmeans {
  using sample    = vq3::topology::gi::Value<aux::traits>;
  using prototype = sample;
  
  using vlayer_0  = prototype;                                
  using vlayer_1  = vq3::demo::decorator::colored<vlayer_0>;  // We will color the prototypes.
  using vertex    = vlayer_1;
  
  using graph     = vq3::graph<vertex, void>;
}


int main(int argc, char* argv[]) {

  unsigned int nb_threads = std::thread::hardware_concurrency();
  std::mt19937 random_device(0);

  // Learning rates and nb steps.
  std::vector<std::pair<double, std::size_t>> alphas = {{.2, 500}, {.1, 500}, {.05, 500}, {.02, 500}};

  // The dataset.
  std::vector<demo2d::Point> S;
  {
    double side              = 2*NOISE_RADIUS;
    double intensity         = 1;
    double noise_density     = NOISE_DENSITY;
    double bar_thickness     = 2* ARC_SEP_THICKNESS;
    double big_ext           = BIG_RADIUS;
    double big_int           = BIG_RADIUS - 2*ARC_THICKNESS;
    double big_bar_length    = 2*(BIG_RADIUS + BAR_OFFSET);
    double medium_ext        = MEDIUM_RADIUS;
    double medium_int        = MEDIUM_RADIUS - 2*ARC_THICKNESS;
    double medium_bar_length = 2*(MEDIUM_RADIUS + BAR_OFFSET);  
    double small_ext         = SMALL_RADIUS;
    double small_bar_length  = 2*(SMALL_RADIUS + BAR_OFFSET);
    double angle             = 90;
    auto density = demo2d::sample::rectangle(side, side, noise_density);
    auto crown   = demo2d::sample::disk(big_ext, intensity) - demo2d::sample::disk(big_int, intensity);
    auto pieces  = crown - demo2d::sample::rectangle(big_bar_length, bar_thickness, intensity);
    density      = density || pieces;
    crown        = demo2d::sample::disk(medium_ext, intensity) - demo2d::sample::disk(medium_int, intensity);
    pieces       = crown - demo2d::sample::rectangle(medium_bar_length, bar_thickness, intensity);
    density      = density || (pieces % angle);
    crown        = demo2d::sample::disk(small_ext, intensity);
    pieces       = crown - demo2d::sample::rectangle(small_bar_length, bar_thickness, intensity);
    density      = density || crown;
    auto sampler = demo2d::sample::base_sampler::random(random_device, NB_SAMPLES_PER_M2);
    auto S_      = demo2d::sample::sample_set(random_device, sampler, density);
    std::copy(S_.begin(), S_.end(), std::back_inserter(S));
  }

  // Setting up the auxiliary graph

  aux::graph g_aux;
  
  std::cout << "Starting auxiliary graph..." << std::endl;
  auto aux_end = S.begin() + (std::size_t)(S.size() * CHL_FRACTION);
  for(auto it = S.begin(); it != aux_end; ++it) g_aux += *it;
  auto chl = vq3::epoch::chl::processor(g_aux);
  chl.process(nb_threads, S.begin(), S.end(), [](const demo2d::Point& s) {return s;}, aux::D, aux::edge());  
  g_aux.foreach_edge([](auto ref_e){ // We compute edge lengths
		       auto [ref_a, ref_b] = ref_e->extremities();
		       (*ref_e)().vq3_cost.edge_length = aux::d((*ref_a)().vq3_value, (*ref_b)().vq3_value);
		     });
  
  auto e_counter = vq3::epoch::count::edges::processor(g_aux);
  vq3::utils::clear_edge_counters(g_aux, 0);
  // We ccount the number of examples in a cylinder around each edge.
  e_counter.process(nb_threads,
		    S.begin(), S.end(),
		    [](auto& content) {return content;}, 
		    [](auto& ref_e, const demo2d::Point& sample) { 
		      auto [ref_a, ref_b] = ref_e->extremities(); 
		      return demo2d::cylinder_d2(sample, {(*ref_a)().vq3_value, (*ref_b)().vq3_value}) < EDGE_THICKNESS * EDGE_THICKNESS;
		    });
  // We compute the edge densities.
  g_aux.foreach_edge(aux::compute_edge_relative_density);
  std::cout << "... Done." << std::endl;


  // Setting up te k-means graph
  
  kmeans::graph g_kmeans;
  auto traits = vq3::topology::gi::traits<aux::sample>(g_aux, aux::d, aux::D, aux::interpolate, aux::shortest_path, aux::compute_cumulated_costs, aux::cumulated_cost_of);
  auto kmeans_d = vq3::topology::gi::distance<kmeans::vertex>(traits, [](const kmeans::vertex& vertex) -> const kmeans::prototype& {return vertex.vq3_value;});
  std::array<cv::Scalar, K> colormap = COLORMAP;
  auto end = S.begin() + K;
  unsigned int i = 0;
  for(auto it = S.begin(); it != end; ++it) {
    auto ref_v = g_kmeans += vq3::topology::gi::value(traits, *it);
    (*ref_v)().vq3_color = colormap[i++];
  }

  // Learning
  
  std::cout << "Starting k-means..." << std::endl;
  auto it = S.end();
  for(auto [alpha, nb] : alphas) {
    for(std::size_t i = 0; i < nb; ++i) {
      if(it == S.end()) {
	std::shuffle(S.begin(), S.end(), random_device);
	it = S.begin();
      }
      vq3::online::wta::learn(g_kmeans,
			      [](kmeans::vertex& vertex_value) -> kmeans::prototype& {return vertex_value.vq3_value;},
			      kmeans_d, vq3::topology::gi::value(traits, *(it++)),
			      alpha); 
    }
  }
  std::cout << "... Done." << std::endl;

  // Display
  auto image = cv::Mat(IMG_SIDE, IMG_SIDE, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = demo2d::opencv::direct_orthonormal_frame(cv::Size(IMG_SIDE, IMG_SIDE), ZOOM*IMG_SIDE, true);

  cv::Scalar sample_color {200, 200, 200};
  auto draw_samples = demo2d::opencv::dot_drawer<demo2d::Point>(image, frame,
								[             ](const demo2d::Point& pt) {return                      true;},
								[             ](const demo2d::Point& pt) {return                        pt;},
								[             ](const demo2d::Point& pt) {return                         1;},
								[&sample_color](const demo2d::Point& pt) {return              sample_color;},
								[             ](const demo2d::Point& pt) {return                        -1;});
  

  image = cv::Scalar(255, 255, 255);
  std::copy(S.begin(), S.end(), draw_samples);
  cv::imwrite("noisy-git-inputs.png", image);
  std::cout << "File \"noisy-git-inputs.png\" generated." << std::endl;
  
  auto draw_edge_aux = vq3::demo2d::opencv::edge_drawer<aux::graph::ref_edge>(image, frame,
									      [](const aux::vertex& v1, const aux::vertex& v2, const aux::edge& e) {return true;}, 
									      [](const aux::vertex& v)  {return v.vq3_value;}, 
									      [](const aux::edge&   e)  {
										double component = 255*(.95*(1 - e.vq3_cost.relative_density));
										return cv::Scalar(component, component, component);
									      }, 
									      [](const aux::edge&   e)  {return 2;});
  cv::imwrite("noisy-git-auxiliary.png", image);
  std::cout << "File \"noisy-git-auxiliary.png\" generated." << std::endl;
  g_aux.foreach_edge(draw_edge_aux);

  image = cv::Scalar(255, 255, 255);
  std::vector<demo2d::Point> SS;
  {
    double rect_width     = NOISE_RADIUS*2.1;
    double rect_height    = NOISE_RADIUS*2.1;
    double rect_intensity = 1.0;
    auto density          = demo2d::sample::rectangle(rect_width, rect_height, rect_intensity);
    auto sampler          = demo2d::sample::base_sampler::grid(random_device, NB_SAMPLES_PER_M2);
    auto S_               = demo2d::sample::sample_set(random_device, sampler, density);
    std::copy(S_.begin(), S_.end(), std::back_inserter(SS));
  }
  auto draw_vertex_white = vq3::demo2d::opencv::vertex_drawer<kmeans::graph::ref_vertex>(image, frame,
											 [](const kmeans::vertex& v) {return                      true;}, 
											 [](const kmeans::vertex& v) {return             v.vq3_value();}, 
											 [](const kmeans::vertex& v) {return                        15;}, 
											 [](const kmeans::vertex& v) {return cv::Scalar(255, 255, 255);}, 
											 [](const kmeans::vertex& v) {return                        -1;});
  auto draw_vertex_color = vq3::demo2d::opencv::vertex_drawer<kmeans::graph::ref_vertex>(image, frame,
											 [](const kmeans::vertex& v) {return          true;}, 
											 [](const kmeans::vertex& v) {return v.vq3_value();}, 
											 [](const kmeans::vertex& v) {return            10;}, 
											 [](const kmeans::vertex& v) {return   v.vq3_color;}, 
											 [](const kmeans::vertex& v) {return            -1;});
  auto draw_vertex_black = vq3::demo2d::opencv::vertex_drawer<kmeans::graph::ref_vertex>(image, frame,
											 [](const kmeans::vertex& v) {return                true;}, 
											 [](const kmeans::vertex& v) {return       v.vq3_value();},  
											 [](const kmeans::vertex& v) {return                  10;},
											 [](const kmeans::vertex& v) {return cv::Scalar(0, 0, 0);}, 
											 [](const kmeans::vertex& v) {return                   2;});
  auto draw_voronoi_sample = demo2d::opencv::dot_drawer<demo2d::Point>(image, frame,
								       [](const demo2d::Point& pt) {return                      true;},
								       [](const demo2d::Point& pt) {return                        pt;},
								       [](const demo2d::Point& pt) {return                         2;},
								       [&g_kmeans, &traits, &kmeans_d](const demo2d::Point& pt) {
									 if(auto closest = vq3::utils::closest(g_kmeans, vq3::topology::gi::value(traits, pt), kmeans_d); closest)
									   return (*closest)().vq3_color;
									 else
									   return cv::Scalar(0, 0, 0);
								       },
								       [](const demo2d::Point& pt) {return                        3;});
  std::cout << "Drawing voronoi cells..." << std::endl;
  image = cv::Scalar(255, 255, 255);
  std::copy(SS.begin(), SS.end(), draw_voronoi_sample);
  sample_color = cv::Scalar(0, 0, 0);
  std::copy(S.begin(), S.end(), draw_samples);
  g_kmeans.foreach_vertex(draw_vertex_white);
  g_kmeans.foreach_vertex(draw_vertex_color);
  g_kmeans.foreach_vertex(draw_vertex_black);
  std::cout << "... Done." << std::endl;
  cv::imwrite("noisy-git-kmeans.png", image);
  std::cout << "File \"noisy-git-kmeans.png\" generated." << std::endl;
  return 0;
}
