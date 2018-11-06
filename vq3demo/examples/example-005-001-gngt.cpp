#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <cmath>


#define NT_NB_BINS       20


enum class Mode : char {Cont = 'c', Step = 's'};

struct Param {
  double alpha() {return .1;}
  unsigned int min_updates() {return 3;}
};

// Graph definition
//
///////////////

using sample    = vq3::demo2d::Point;
using prototype = vq3::demo2d::Point;

//                                                                          ## Node properties :
using vlayer_0 = prototype;                                                 // prototypes are 2D points (this is the "user defined" value).
using vlayer_1 = vq3::decorator::tagged<vlayer_0>;                          // we add a tag for topology computation.
using vlayer_2 = vq3::decorator::online::mean_std<vlayer_1, double, Param>; // we add distortion statistics over time. 
using vlayer_3 = vq3::demo::decorator::colored<vlayer_2>;                   // we add a color for nice drawings.
using vertex   = vlayer_3;

//                                                        ## Edge properties :
using elayer_0 = vq3::decorator::tagged<void>;            // we add a tag for CHL computation.
using edge     = elayer_0;

using graph  = vq3::graph<vertex, edge>;


// Epoch data for SOM-like pass
//
////////////////

using epoch_wtm = vq3::epoch::data::wtm<vq3::epoch::data::none<sample, vertex, prototype>>;

// Distance
//
////////////////


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double dist(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v.vq3_value, p);}


// Network evolution rule
//
////////////////

// There is a default evolution rule (see comments in the
// main). Nevertheless, let us customize the evolution, for the sake
// of illustration of what it is. Our evolution is similar to the
// default one, but it adds colors to the vertices.

struct Evolution {
  double T;
  double density;
  double sigma_coef;
  vq3::demo2d::opencv::colormap::jet colormap;
  vq3::demo2d::opencv::histogram histo; // This collects information for display.
  vq3::stats::MeanStd mean_std;
  double radius = 0;
  double NT = 0;
	  
  Evolution()
    : histo({0, -1.3}, {2.4, -0.6}) {
    histo.frame_margin = .1;
  }
  
  Evolution(const Evolution&) = default;


  template<typename TABLE, typename BMU_RESULT, typename CLONE_PROTOTYPE>
  void operator()(TABLE& table,
		  const BMU_RESULT& bmu_epoch_result,
		  const CLONE_PROTOTYPE& clone_prototype) {

    std::vector<typename TABLE::graph_type::ref_vertex> above;
    std::vector<typename TABLE::graph_type::ref_vertex> below;
    auto out_above = std::back_inserter(above);
    auto out_below = std::back_inserter(below);
    
    NT = density*T;

    histo.range.reset(); // we do not display the range except if we are told to afterwards.
    histo.NT  = NT;
    auto hout = histo.output_iterator();
    
    mean_std.clear();
    
    for(auto& res : bmu_epoch_result)
      if(res.vq3_bmu_accum.nb != 0)
	mean_std  = res.vq3_bmu_accum.value;
    double spatial_dmean = std::sqrt(mean_std.variance())*sigma_coef;
    
    colormap = {NT - spatial_dmean, NT + spatial_dmean};
    
    unsigned int idx = 0;
    for(auto& res : bmu_epoch_result) {
      auto& ref_v = table(idx++);
      if(res.vq3_bmu_accum.nb == 0)
	ref_v->kill(); // We kill a vertex which has never won the competition.
      else {
	*(hout++) = res.vq3_bmu_accum.value;
	auto& vertex = (*ref_v)();
	if(auto& omstd = vertex.vq3_online_mean_std; omstd) {
	  auto [m, std] = omstd(); // We get the vertex distortion statistics.
	  auto dm = sigma_coef*std;
	  
	  vertex.vq3_color = colormap(m);
	  radius = dm + spatial_dmean;
	  
	  if(NT < m - radius) // There are not enough nodes (cold color).
	    *(out_above++) = ref_v;
	  else if(m + radius < NT)  // There are too many nodes (hot color).
	    *(out_below++) = ref_v;
	}
	else 
	  vertex.vq3_color = cv::Scalar(0,0,0);
      }
    }

    if(above.size() > 0)
      for(auto it = above.begin(); it != above.end(); ++it)
	table.g += clone_prototype((*(*it))().vq3_value);
    
    if(below.size() > 0)
      for(auto it = below.begin(); it != below.end(); ++it)
	(*it)->kill();
  }
};

inline Evolution make_evolution() {
  return Evolution();
}


// GUI Callback data
//
////////////////

struct callback_data {
  double high_density  =   1;
  double low_density   = .25;
  double noise_density = .05;

  double i1, i2, i3, i4; // object densities.
  bool b1 = true;    // visibility of object 1
  bool b2 = true;    // visibility of object 2
  bool b3 = true;    // visibility of object 3
  bool b4 = true;    // visibility of object 3

  vq3::demo2d::sample::density d1;
  vq3::demo2d::sample::density d2;
  vq3::demo2d::sample::density d3;
  vq3::demo2d::sample::density d4;

  const vq3::demo2d::opencv::Frame& frame;
  graph& g;
  graph::ref_vertex ref_v;

  callback_data(const vq3::demo2d::opencv::Frame& frame, graph& g)
    : frame(frame), g(g) {
    update();
  }

  void toggle1() {b1 = !b1;}
  void toggle2() {b2 = !b2;}
  void toggle3() {b3 = !b3;}
  void toggle4() {b4 = !b4;}

  void update() {
    if(b1) i1 = high_density;   else i1 = 0;
    if(b2) i2 = high_density;   else i2 = 0;
    if(b3) i3 = low_density;    else i3 = 0;
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
    
    data.update();
  }
  
  if(event == cv::EVENT_RBUTTONDOWN)
    data.ref_v = vq3::utils::closest(data.g, data.frame(cv::Point(x,y)), dist);
}


// Main
//
////////////////

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
  int N_slider =  5000;
  int T_slider =   150;
  int S_slider =   150;
  

  // Image settings
  //
  ///////////////////

  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  
  cv::createTrackbar("nb/m^2",         "image", &N_slider, 50000, nullptr);
  cv::createTrackbar("T",              "image", &T_slider,  1000, nullptr);
  cv::createTrackbar("100*sigma_coef", "image", &S_slider,   300, nullptr);
  auto image       = cv::Mat(576, 1024, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .2*image.size().width, true);
  callback_data cb(frame, g);
  cv::setMouseCallback("image", on_mouse, reinterpret_cast<void*>(&cb));
  
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
  //
  ///////////////////

  // Thread shapes thickness
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
  //
  ///////////////////
  

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

  auto topology = vq3::topology::table(g);

  // This processes the topology evolution (number of vertices and edges)
  auto gngt     = vq3::algo::gngt::processor<prototype, sample>(topology);

  // This updates the prototypes (i.e. the vertex positions).
  auto wtm      = vq3::epoch::wtm::processor(topology);
  
  // This is how the default evolution would have been obtained.
  // auto evolution = vq3::algo::gngt::by_default::evolution(random_device);

  auto evolution = make_evolution();
  
  // This is the loop
  //
  //////////////////

  Mode mode = Mode::Cont;
  bool compute = true;
  int wait_ms;

  while(keycode != 27) {
    compute = mode != Mode::Step;
    
    if(keycode == 32) {
      mode = Mode::Step;
      compute = true;
    }
    else if(keycode == 99) {
      mode = Mode::Cont;
      compute = true;
    }

    wait_ms = 500;
    if(compute) {
      wait_ms = 1;
    
      // Get the samples
      
      auto S_ = vq3::demo2d::sample::sample_set(random_device, density, N_slider);
      S.clear();
      auto out = std::back_inserter(S);
      std::copy(S_.begin(), S_.end(), out);

      // Step

      double e = T_slider/1000.0;
      double expo_min = -5;
      double expo_max = -1;
      

      evolution.density    = N_slider;
      evolution.T          = std::pow(10, expo_min*(1-e) + expo_max*e);
      evolution.sigma_coef = S_slider*.01;

      // We compute the topology evolution of the graph...
      gngt.process(nb_threads,
		   S.begin(), S.end(),
		   [](const sample& s) {return s;},
		   [](vertex& v) -> prototype& {return v.vq3_value;},
		   [](const prototype& p) {return p + vq3::demo2d::Point(-1e-5,1e-5);},
		   dist,
		   evolution);

      // ... and update the nodes with a SOM-like pass. Here,
      // classically (for GNG), only the immediate neighbours are
      // updated.
      topology([](unsigned int edge_distance) {return edge_distance == 0 ? 1.0 : 0.1;}, 1, 0);
      wtm.process<epoch_wtm>(nb_threads,
			     S.begin(), S.end(),
			     [](const sample& s) {return s;},
			     [](vertex& v) -> prototype& {return v.vq3_value;},
			     dist);
							    

      // Display
    
      image = cv::Scalar(255, 255, 255);

      // Display selected node information.
      if(cb.ref_v) {// If there is a selected vertex to display
	if(cb.ref_v->is_killed()) // If the vertex has died...
	  cb.ref_v = nullptr;     // ... we release the pointer to it.
	else 
	  cv::circle(image, frame((*(cb.ref_v))().vq3_value), 7, cv::Scalar(0,0,0), -1);
      }
      
      // Display the points and the graph
      std::copy(S.begin(), S.end(), dd);
      g.foreach_edge(draw_edge); 
      g.foreach_vertex(draw_vertex);
    

      try {
	evolution.histo.set_bins(NT_NB_BINS);
	evolution.histo.make();
	evolution.histo.draw(image, frame);
	
	auto [mm, ss, aa] = evolution.histo.get_normal_fit();
	evolution.histo.normal_var(image, frame, mm, ss, aa, 50, cv::Scalar(255,0,0), 1, false);
	evolution.histo.interval(image, frame,
				 evolution.NT - evolution.radius,
				 evolution.NT + evolution.radius,
				 6, cv::Scalar(0,0,0), 3);

	if(cb.ref_v) {
	  auto& vertex = (*(cb.ref_v))();
	  if(vertex.vq3_online_mean_std) {
	    auto [m, s] =  vertex.vq3_online_mean_std();
	    evolution.histo.normal_var(image, frame, m, s, aa, 50, cv::Scalar(0,0,0),   5, false);
	    evolution.histo.normal_var(image, frame, m, s, aa, 50, vertex.vq3_color,    3, false);
	  }
	}
      }
      catch(std::runtime_error& e) {}
    }
    
    
    cv::imshow("image", image);
    keycode = cv::waitKey(wait_ms) & 0xFF;
  }
    

  return 0;
}
