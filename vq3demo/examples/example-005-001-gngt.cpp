#include <vector>
#include <iterator>
#include <tuple>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <string>

#include <opencv2/opencv.hpp>

#include <vq3.hpp>
#include <vq3demo.hpp>

// Histograms configuration.

#define NB_BINS 20

#define SOM_H_RADIUS 3.1
#define SOM_MAX_DIST (unsigned int)(SOM_H_RADIUS)

#define NARROW_SOM_COEF .02


// Vertex colors

#define COLOR_BELOW    cv::Scalar(000, 000, 255) // Red, too crowdy nodes here.
#define COLOR_ABOVE    cv::Scalar(255, 000, 000) // Blue, nodes are too sparse here.
#define COLOR_INBOUNDS cv::Scalar(000, 136, 246) // Orange


// GUI sliders

#define INIT_SLIDER_N               5000
#define INIT_SLIDER_T                150
#define INIT_SLIDER_DENSITY           25
#define INIT_SLIDER_MARGIN_ABOVE      35
#define INIT_SLIDER_MARGIN_BELOW      20
#define INIT_SLIDER_AVERAGE_RADIUS     8
#define INIT_SLIDER_EVOLUTION_RATIO   30
#define INIT_SLIDER_NB_WTA_AFTER       1
#define INIT_SLIDER_ALPHA             50
#define INIT_SLIDER_SAMPLE_PER_VERTEX 10


// Execution mode

enum class Mode : char {Cont = 'c', Step = 's'};


// GUI Callback data

struct callback_data {
  double high_density  =   1;
  double noise_density = .05;

  double i1, i2, i3, i4; // object densities.
  bool b1 = true;    // visibility of object 1
  bool b2 = true;    // visibility of object 2
  bool b3 = true;    // visibility of object 3
  bool b4 = false;   // visibility of object 4

  vq3::demo2d::sample::density d1;
  vq3::demo2d::sample::density d2;
  vq3::demo2d::sample::density d3;
  vq3::demo2d::sample::density d4;

  const vq3::demo2d::opencv::Frame& frame;
  int& d_slider;

  callback_data(const vq3::demo2d::opencv::Frame& frame, int& d_slider)
    : frame(frame), d_slider(d_slider) {
    update();
  }

  void toggle1() {b1 = !b1;}
  void toggle2() {b2 = !b2;}
  void toggle3() {b3 = !b3;}
  void toggle4() {b4 = !b4;}

  void update() {
    if(b1) i1 = high_density;   else i1 = 0;
    if(b2) i2 = high_density;   else i2 = 0;
    if(b3) i3 = .01*d_slider;   else i3 = 0;
    if(b4) i4 = noise_density;  else i4 = 0;
  }
};

void on_mouse( int event, int x, int y, int, void* user_data) {
  auto& data = *(reinterpret_cast<callback_data*>(user_data));
    
  // A mouse click reinitializes the graph.
  if(event == cv::EVENT_LBUTTONDOWN) {
    auto click = data.frame(cv::Point(x,y));
    if(data.d1->bbox().contains(click))      data.toggle1();
    else if(data.d2->bbox().contains(click)) data.toggle2();
    else if(data.d3->bbox().contains(click)) data.toggle3();
    else                                     data.toggle4();
  }
}


// Graph definition

using sample    = vq3::demo2d::Point;
using prototype = vq3::demo2d::Point;

//                                                        ## Node properties :
using vlayer_0 = prototype;                               // prototypes are 2D points (this is the "user defined" value).
using vlayer_1 = vq3::decorator::tagged<vlayer_0>;        // we add a tag for topology computation.
using vlayer_2 = vq3::demo::decorator::colored<vlayer_1>; // we add a color for nice drawings.
using vertex   = vlayer_2;

//                                             ## Edge properties :
using elayer_0 = vq3::decorator::tagged<void>; // we add a tag for CHL computation.
using edge     = elayer_0;

using graph  = vq3::graph<vertex, edge>;

using neighbour_key_type = std::string;


// Distance

// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double dist(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v.vq3_value, p);}


// Network evolution rule.

// A default one is provided by vq3. Let us set up an evolution rule
// for the sake of illustration. It corresponds to the default one
// with extra addings for display.


struct Evolution {

  // These attributes are set from sliders.
  
  double T;
  double density;
  double margin_above = .01*INIT_SLIDER_MARGIN_ABOVE;
  double margin_below = .01*INIT_SLIDER_MARGIN_BELOW;
  double topo_ratio   = .01*INIT_SLIDER_EVOLUTION_RATIO;
  
  
  double NT           = 0;    // It is computed as density*T.
  bool unfreezed      = true; // Tells wether evolution is actually allowed (the user can feeze it).

  // Collect the (vertex_idx, bmu_error) whose BMU error is too high and too_low.
  
  std::vector<std::pair<std::size_t, double>> above;
  std::vector<std::pair<std::size_t, double>> below;

  // This shows the distribution of vertex errors.
  
  vq3::demo2d::opencv::histogram histo;

	  	  
  Evolution()
    : histo({2.4, -0.8}, {4.4, -0.1}) {
    histo.frame_margin = .1;
    histo.title        = "Averaged distortions";
  }

  // GNG-T calls this method when it considers to perform a
  // modification of the number of vertices.
  template<typename TABLE, typename BMU_RESULT, typename CLONE_PROTOTYPE>
  bool operator()(TABLE&                 topology,
		  const BMU_RESULT&      neighboring_bmu_epoch_result,
		  const CLONE_PROTOTYPE& clone_prototype) {
    double topology_changed = false;

    above.clear();
    auto above_out = std::back_inserter(above);

    below.clear();
    auto below_out = std::back_inserter(below);
    
    NT = density*T;

    double above_bound = NT*(1+margin_above);
    double below_bound = NT*(1-margin_below);

    // We prepare the histogram
    histo.NT           = NT;                         // display a red bar for the value NT
    histo.value_bounds = {0, 2*NT};                  // value horizontal axis has the range [0, 2*NT]
    histo.range        = {below_bound, above_bound}; // display a horizontal range corresponding to inbound values.
    auto hout          = histo.output_iterator();    // provide an output iterator for feeding the histogram.

    // Let us consider all the errors.
    std::size_t vertex_idx = 0;
    for(auto& res : neighboring_bmu_epoch_result) {
      if(res.vq3_bmu_accum.nb != 0) { // We consider only the vertices which have been a bmu at least once.
	double error  = res.vq3_bmu_accum.average();
	auto& color = (*(topology(vertex_idx)))().vq3_color; // This is a reference to the vertex color.
	*(hout++) = error;
	if(error > above_bound) {
	  *(above_out++) = {vertex_idx, error};
	  color = COLOR_ABOVE;
	}
	else if(error < below_bound) {
	  *(below_out++) = {vertex_idx, error};
	  color = COLOR_BELOW;
	}
	else
	  color = COLOR_INBOUNDS;
      }
      else {
	// We remove from the graph the vertices that have never been selected as a BMU.
	if(unfreezed) {
	  topology(vertex_idx)->kill();
	  topology_changed = true;
	}
      }
      ++vertex_idx;
    }

    // We sort out of bounds vertices (small error first for below, big error first for above).
    std::sort(above.begin(), above.end(),
	      [](const std::pair<std::size_t, double>& p1, const std::pair<std::size_t, double>& p2) {
		return p1.second > p2.second;});
    std::sort(below.begin(), below.end(),
	      [](const std::pair<std::size_t, double>& p1, const std::pair<std::size_t, double>& p2) {
		return p1.second < p2.second;});

    if(unfreezed) {

      // We clone a topo_ratio fraction of the above vertices.
      
      auto above_end = above.begin();
      if(above_end != above.end()) {// if not empty
	std::advance(above_end, std::max((std::size_t)(above.size()*topo_ratio), (size_t)1));
	for(auto it = above.begin(); it != above_end; ++it) {
	  auto ref_v = topology.g += clone_prototype((*(topology(it->first)))().vq3_value);
	  (*ref_v)().vq3_color = COLOR_INBOUNDS;
	}
      }

      // we delete a topo_ration fraction of the below vertices.
      auto below_end = below.begin();
      if(below_end != below.end()) {// if not empty
	std::advance(below_end, std::max((std::size_t)(below.size()*topo_ratio), (size_t)1));
	for(auto it = below.begin(); it != below_end; ++it) 
	  topology(it->first)->kill();
      }
      
      topology_changed = (above.begin() != above_end) || (below.begin() != below_end);
    }

    // We tell GNG-T if any change in the vertices (add/remove) has occurred.
    return topology_changed;
  }
};

inline Evolution make_evolution() {
  return Evolution();
}


// Main


int main(int argc, char* argv[]) {
  
  unsigned int nb_threads = std::thread::hardware_concurrency();
  if(nb_threads == 0) {
    nb_threads = 1;
    std::cout << "The harware multi-threading capabilities cannot be queried. I use a single thread." << std::endl;
  }
  else
    std::cout << std::endl
	      << "I use " << nb_threads << " thread(s)." << std::endl;


  
  std::random_device rd;  
  std::mt19937 random_device(rd());

  graph g;
  
  int slider_N                 = INIT_SLIDER_N;
  int slider_T                 = INIT_SLIDER_T;
  int slider_density           = INIT_SLIDER_DENSITY;
  int slider_margin_above      = INIT_SLIDER_MARGIN_ABOVE;
  int slider_margin_below      = INIT_SLIDER_MARGIN_BELOW;
  int slider_average_radius    = INIT_SLIDER_AVERAGE_RADIUS;
  int slider_evolution_ratio   = INIT_SLIDER_EVOLUTION_RATIO;
  int slider_nb_wta_after      = INIT_SLIDER_NB_WTA_AFTER;
  int slider_alpha             = INIT_SLIDER_ALPHA;
  int slider_sample_per_vertex = INIT_SLIDER_SAMPLE_PER_VERTEX;
  
  int old_slider_average_radius = slider_average_radius;

  // Let us used an histogram for non-averaged error distribution.
  vq3::demo2d::opencv::histogram histo({2.4, 0.2}, {4.4, 0.9});
  histo.frame_margin = .1;
  histo.title = "Raw distortions";
  
  // Image settings
  
  cv::namedWindow("algorithm", CV_WINDOW_AUTOSIZE);
  cv::namedWindow("params", CV_WINDOW_AUTOSIZE);
  
  cv::createTrackbar("nb/m^2",                 "params", &slider_N,               50000, nullptr);
  cv::createTrackbar("T",                      "params", &slider_T,                1000, nullptr);
  cv::createTrackbar("above_margin*100",       "params", &slider_margin_above,      100, nullptr);
  cv::createTrackbar("below_margin*100",       "params", &slider_margin_below,      100, nullptr);
  cv::createTrackbar("right square density",   "params", &slider_density,           100, nullptr);
  cv::createTrackbar("statial average radius", "params", &slider_average_radius,     10, nullptr);
  cv::createTrackbar("growth/shrink ratio",    "params", &slider_evolution_ratio,   100, nullptr);
  cv::createTrackbar("nb_wta_after",           "params", &slider_nb_wta_after,       20, nullptr);
  cv::createTrackbar("alpha",                  "params", &slider_alpha,             200, nullptr);
  cv::createTrackbar("online samples/vertex",  "params", &slider_sample_per_vertex, 100, nullptr);
  
  auto image       = cv::Mat(600, 1500, CV_8UC3, cv::Scalar(255,255,255));
  auto params      = cv::Mat(1, 600, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = vq3::demo2d::opencv::direct_orthonormal_frame(220, 220, {512, 288});
  callback_data cb(frame, slider_density);
  cv::setMouseCallback("algorithm", on_mouse, reinterpret_cast<void*>(&cb));
  
  auto dd          = vq3::demo2d::opencv::dot_drawer<vq3::demo2d::Point>(image, frame,
									 [](const vq3::demo2d::Point& pt) {return                      true;},
									 [](const vq3::demo2d::Point& pt) {return                        pt;},
									 [](const vq3::demo2d::Point& pt) {return                         1;},
									 [](const vq3::demo2d::Point& pt) {return cv::Scalar(200, 200, 200);},
									 [](const vq3::demo2d::Point& pt) {return                        -1;});
  auto draw_edge   = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
								       [](const vertex& v1, const vertex& v2, const edge& e) {return true;}, // always draw
								       [](const vertex& v)              {return               v.vq3_value;}, // position
								       [](const edge&   e)              {return cv::Scalar(150, 150, 150);}, // color
								       [](const edge&   e)              {return                        1;}); // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return        true;},  // always draw
									   [](const vertex& v) {return v.vq3_value;},  // position
									   [](const vertex& v) {return           3;},  // radius
									   [](const vertex& v) {return v.vq3_color;},  // color
									   [](const vertex& v) {return          -1;}); // thickness
  
  // Distribution settings

  // wire shapes thickness
  double thickness = .1;
  double rect_dist = 1.5;
    
  // Left square
  double w1             = 1;
  double h1             = 1;
  vq3::demo2d::Point o1 = {-rect_dist, 0};
  auto shape1           = vq3::demo2d::sample::rectangle(w1, h1, cb.i1) + o1;
  cb.d1                 = shape1;
    
  // Middle crown
  double i2   =  1;
  double r2   = .5;
  double h2   = r2 - thickness;
  auto shape2 = vq3::demo2d::sample::disk(r2, cb.i2) - vq3::demo2d::sample::disk(h2, i2);

  // Right square
  double w3             = 1;
  double h3             = 1;
  vq3::demo2d::Point o3 = {rect_dist, 0};
  auto shape3           = vq3::demo2d::sample::rectangle(w3, h3, cb.i3) + o3;
  cb.d3                 = shape3;

  // Bar
  vq3::demo2d::Point min4 = {-rect_dist +  .5, -thickness*.5};
  vq3::demo2d::Point max4 = {             -h2,  thickness*.5};
  auto shape4             = vq3::demo2d::sample::rectangle(min4, max4, cb.i2);
  cb.d2                   = shape2 || shape4;

  // Noise
  double w5             = 4.5;
  double h5             = 2;
  vq3::demo2d::Point o5 = {0, 0};
  auto shape5           = vq3::demo2d::sample::rectangle(w5, h5, cb.i4) + o5;
  cb.d4                 = shape5;

  // All
  auto density = shape1 || shape2 || shape3 || shape4 || shape5;

  

  // Some initializations
 
  std::cout << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << "click (left button) on the shape locations to toggle them." << std::endl
	    << "click (right button) on a prototype to display its statistics." << std::endl
	    << std::endl
	    << "press <space> to compute step by step." << std::endl
	    << "press c to exit step by step mode (cont)." << std::endl
	    << "press ESC to quit." << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl;
  int keycode = 0;
 
  // We need to register the input samples in a vector since we want
  // to both use and display them. Moreover, for gngt, the sample set
  // is iterated several times.
  std::vector<vq3::demo2d::Point> S;

  // This keeps up to date information about the graph topology.
  auto topology = vq3::topology::table<neighbour_key_type>(g);
  topology.declare_distance("wide som",   [](unsigned int edge_distance) {return std::max(0., 1 - edge_distance/double(SOM_H_RADIUS));},          SOM_MAX_DIST, 1e-3);
  topology.declare_distance("narrow som", [](unsigned int edge_distance) {return edge_distance == 0 ? 1 : NARROW_SOM_COEF ;           },                     1,  0.0);
  topology.declare_distance("avg",        [](unsigned int edge_distance) {return 1;                                                   }, slider_average_radius,  0.0);

  // This processes the topology evolution (number of vertices and edges)
  auto gngt     = vq3::algo::gngt::processor<sample>(topology);
  
  // This is how the default evolution would have been obtained.
  //   auto evolution = vq3::algo::gngt::by_default::evolution();
  // But we use here the one that we have handcrafted.
  auto evolution = make_evolution();

  
  // This is the loop

  Mode mode = Mode::Step;
  bool compute = true;
  int wait_ms;

  while(keycode != 27) {       // no ESC key pressed.
    if(slider_average_radius != old_slider_average_radius) {
      // we need to redeclare avg neighborhoods.
      topology.declare_distance("avg", [](unsigned int edge_distance) {return 1;}, slider_average_radius, 0.0);
      old_slider_average_radius = slider_average_radius;
    }

    
    compute = mode != Mode::Step;
    
    if(keycode == 32) {        // space key pressed.
      mode = Mode::Step;
      compute = true;
    }
    else if(keycode == 99) {   // 'c' key pressed.
      mode = Mode::Cont;
      compute = true;
    }
    else if(keycode == 102)    // 'f' key pressed.
      evolution.unfreezed = !evolution.unfreezed;
    

    cb.update();
    wait_ms = 500;
    if(compute) {
      wait_ms = 1;
    
      // Get the samples
      
      auto S_ = vq3::demo2d::sample::sample_set(random_device, density, slider_N);
      S.clear();
      auto out = std::back_inserter(S);
      std::copy(S_.begin(), S_.end(), out);

      // Step

      // T do not evolve proportionally to the slider.
      double e = slider_T/1000.0;
      double expo_min = -5;
      double expo_max = -1;
      evolution.T             = std::pow(10, expo_min*(1-e) + expo_max*e);
      

      evolution.density       =     slider_N;
      evolution.margin_above  = .01*slider_margin_above;
      evolution.margin_below  = .01*slider_margin_below;
      evolution.topo_ratio    = .01*slider_evolution_ratio;

      gngt.nb_wta_after       = slider_nb_wta_after;
      gngt.alpha              = .001*slider_alpha;
      gngt.samples_per_vertex = slider_sample_per_vertex;

      // We compute the topology evolution of the graph...
      gngt.process(nb_threads,
		   S.begin(), S.end(),                                                             // The sample set. Shuffle if the dataser is not sampled randomly.
		   [](const sample& s) {return s;},                                                // get sample from *iter (identity here).
		   [](vertex& v) -> prototype& {return v.vq3_value;},                              // get a prototype reference from the vertex value.
		   [](const prototype& p) {return p + vq3::demo2d::Point(-1e-5,1e-5);},            // get a point close to a prototype.
		   dist,                      
		   "wide som", "narrow som", "avg",                                                // Neighborhood keys.
		   evolution);

      // Display
    
      image = cv::Scalar(255, 255, 255);
      
      // Display the points and the graph
      std::copy(S.begin(), S.end(), dd);
      g.foreach_edge(draw_edge); 
      g.foreach_vertex(draw_vertex);

      // Display the frozen status.
      if(!evolution.unfreezed)
	cv::putText(image, "frozen", {5, 30}, cv::FONT_HERSHEY_PLAIN, 2., cv::Scalar(0,0,0), 3);

      // We display the 2 histograms.
      try {

	// This is the histogram of the errors used in the evolution algorithm.
	evolution.histo.set_bins(NB_BINS);
	evolution.histo.make();
	evolution.histo.draw(image, frame);

	// This is the istogram of actual vertex errors, with settings
	// similar to the previous one.
	histo.NT           = evolution.histo.NT;
	histo.value_bounds = evolution.histo.value_bounds;
	histo.range        = evolution.histo.range;
	auto hout = histo.output_iterator();
	for(auto& data : gngt.bmu_results)
	  if(data.vq3_bmu_accum.nb > 0)
	    *(hout++) = data.vq3_bmu_accum.value;
	histo.set_bins(NB_BINS);
	histo.make();
	histo.draw(image, frame);
      }
      catch(std::runtime_error& e) {}
    }
    
    
    cv::imshow("algorithm", image );
    cv::imshow("params",    params);
    keycode = cv::waitKey(wait_ms) & 0xFF;
  }
    

  return 0;
}

