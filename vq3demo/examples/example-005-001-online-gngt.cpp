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


#define INIT_SLIDER_N               5000
#define INIT_SLIDER_DENSITY           25
#define INIT_SLIDER_T                 50
#define INIT_SLIDER_ALPHA             10

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

  demo2d::sample::density d1;
  demo2d::sample::density d2;
  demo2d::sample::density d3;
  demo2d::sample::density d4;

  const demo2d::opencv::Frame& frame;
  int& d_slider;

  callback_data(const demo2d::opencv::Frame& frame, int& d_slider)
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

using sample    = demo2d::Point;
using prototype = demo2d::Point;

//                                                        ## Node properties :
using vlayer_0 = prototype;                               // prototypes are 2D points (this is the "user defined" value).
using vlayer_1 = vq3::decorator::tagged<vlayer_0>;        // we add a tag for topology computation.
using vlayer_2 = vq3::demo::decorator::colored<vlayer_1>; // we add a color for nice drawings.
using vertex   = vlayer_2;

//                                                        ## Edge properties :
using elayer_0 = vq3::decorator::tagged<void>;            // we add a tag for CHL computation.
using edge     = elayer_0;

using graph  = vq3::graph<vertex, edge>;

using neighbour_key_type = int;

// Distance

// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
// d2 is faster, but d is more stable in our 2D case.
double dist(const vertex& v, const sample& p) {return demo2d::d(v.vq3_value, p);}


// Network evolution rule.
struct Evolution {
  // here, evolution consists of adding or removing prototypes so as
  // to have a determined number of prototypes.
  unsigned int T = 0;
  
  // GNG-T calls this method when it considers to perform a
  // modification of the number of vertices.
  template<typename TABLE, typename BMU_RESULT, typename CLONE_PROTOTYPE, typename FCT_ERROR_OF_ACCUM>
  bool operator()(TABLE&                    topology,
		  const BMU_RESULT&         neighboring_bmu_epoch_result,
		  const CLONE_PROTOTYPE&    clone_prototype,
		  const FCT_ERROR_OF_ACCUM& error_of_accum) {

    // If some node do not win, we remove them.
    bool removed = false;
    std::size_t vertex_idx = 0;
    for(auto& res : neighboring_bmu_epoch_result) {
      if(res.vq3_bmu_accum.nb == 0) {
	topology(vertex_idx)->kill();
	removed = true;
      }
      ++vertex_idx;
    }
    if(removed)
      return true;
    
    if(topology.size() < T) {
      auto max_iter = std::max_element(neighboring_bmu_epoch_result.begin(), neighboring_bmu_epoch_result.end(),
					 [&error_of_accum](auto& content1, auto& content2){
					 return error_of_accum(content1.vq3_bmu_accum) < error_of_accum(content2.vq3_bmu_accum);});
      auto ref_vertex = topology(std::distance(neighboring_bmu_epoch_result.begin(), max_iter));
      topology.g += clone_prototype((*ref_vertex)().vq3_value);
      return true;
    }

    if(topology.size() > T) {
      auto min_iter = std::max_element(neighboring_bmu_epoch_result.begin(), neighboring_bmu_epoch_result.end(),
					 [&error_of_accum](auto& content1, auto& content2){
					 return error_of_accum(content1.vq3_bmu_accum) < error_of_accum(content2.vq3_bmu_accum);});
      auto ref_vertex = topology(std::distance(neighboring_bmu_epoch_result.begin(), min_iter));
      ref_vertex->kill();
      return true;
    }
    
    return false;
  }
};

inline Evolution make_evolution() {return Evolution();}

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
  int slider_density           = INIT_SLIDER_DENSITY;
  int slider_T                 = INIT_SLIDER_T;
  int slider_alpha             = INIT_SLIDER_ALPHA;


  // Image settings
  
  cv::namedWindow("algorithm", CV_WINDOW_AUTOSIZE);
  cv::namedWindow("params", CV_WINDOW_AUTOSIZE);
  
  cv::createTrackbar("nb/m^2",                 "params", &slider_N,               50000, nullptr);
  cv::createTrackbar("right square density",   "params", &slider_density,           100, nullptr);
  cv::createTrackbar("T",                      "params", &slider_T,                 200, nullptr);
  cv::createTrackbar("100*alpha",              "params", &slider_alpha,             200, nullptr);


  auto image       = cv::Mat(600, 1500, CV_8UC3, cv::Scalar(255,255,255));
  auto params      = cv::Mat(1, 600, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = demo2d::opencv::direct_orthonormal_frame(image.size(), .5*image.size().height, true);
  callback_data cb(frame, slider_density);
  cv::setMouseCallback("algorithm", on_mouse, reinterpret_cast<void*>(&cb));
  
  
  auto dd          = demo2d::opencv::dot_drawer<demo2d::Point>(image, frame,
							       [](const demo2d::Point& pt) {return                      true;},
							       [](const demo2d::Point& pt) {return                        pt;},
							       [](const demo2d::Point& pt) {return                         1;},
							       [](const demo2d::Point& pt) {return cv::Scalar(200, 200, 200);},
							       [](const demo2d::Point& pt) {return                        -1;});
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
  double w1        = 1;
  double h1        = 1;
  demo2d::Point o1 = {-rect_dist, 0};
  auto shape1      = demo2d::sample::rectangle(w1, h1, cb.i1) + o1;
  cb.d1            = shape1;
    
  // Middle crown
  double i2   =  1;
  double r2   = .5;
  double h2   = r2 - thickness;
  auto shape2 = demo2d::sample::disk(r2, cb.i2) - demo2d::sample::disk(h2, i2);

  // Right square
  double w3        = 1;
  double h3        = 1;
  demo2d::Point o3 = {rect_dist, 0};
  auto shape3      = demo2d::sample::rectangle(w3, h3, cb.i3) + o3;
  cb.d3            = shape3;

  // Bar
  demo2d::Point min4 = {-rect_dist +  .5, -thickness*.5};
  demo2d::Point max4 = {             -h2,  thickness*.5};
  auto shape4        = demo2d::sample::rectangle(min4, max4, cb.i2);
  cb.d2              = shape2 || shape4;

  // Noise
  double w5        = 4.5;
  double h5        = 2;
  demo2d::Point o5 = {0, 0};
  auto shape5      = demo2d::sample::rectangle(w5, h5, cb.i4) + o5;
  cb.d4            = shape5;

  // All
  auto density = shape1 || shape2 || shape3 || shape4 || shape5;

  

  // Some initializations
 
  std::cout << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << "click (left button) on the shape locations to toggle them." << std::endl
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
  std::vector<demo2d::Point> S;
  
  // This keeps up to date information about the graph topology.
  auto topology = vq3::topology::table<neighbour_key_type>(g);
  topology.declare_distance(0,
			    [](unsigned int edge_distance) {return edge_distance == 0 ? 1 : .2;},
			    1,  0.0);

  // This processes the topology evolution (number of vertices and edges)
  auto gngt = vq3::algo::online::gngt::processor<sample>(topology);
  
  // This is how the default evolution would have been obtained.
  //   auto evolution = vq3::algo::gngt::by_default::evolution();
  // But we use here the one that we have handcrafted.
  auto evolution = make_evolution();

  
  auto sampler = demo2d::sample::base_sampler::random(random_device, slider_N);
  
  Mode mode = Mode::Cont;
  bool compute = true;
  int wait_ms;

  while(keycode != 27) {       // no ESC key pressed.
    compute = mode != Mode::Step;
    
    if(keycode == 32) {        // space key pressed.
      mode = Mode::Step;
      compute = true;
    }
    else if(keycode == 99) {   // 'c' key pressed.
      mode = Mode::Cont;
      compute = true;
    }
    

    cb.update();
    wait_ms = 500;
    if(compute) {
      wait_ms = 1;
      
      // Get the samples
      
      sampler = slider_N;
      auto S_ = demo2d::sample::sample_set(random_device, sampler, density);
      S.clear();
      auto out = std::back_inserter(S);
      std::copy(S_.begin(), S_.end(), out);

      // GNG-T online
      evolution.T = slider_T;
      gngt.alpha  = .001*slider_alpha;
      // We compute the topology evolution of the graph...
      
      gngt.process(nb_threads,
		   S.begin(), S.end(),                                              // The sample set. Shuffle if the dataset is not sampled randomly.
		   [](const sample& s) {return s;},                                 // get sample from *iter (identity here).
		   [](vertex& v) -> prototype& {return v.vq3_value;},               // get a prototype reference from the vertex value.
		   [](const prototype& p) {return p + demo2d::Point(-1e-5,1e-5);},  // get a point close to a prototype.
		   dist,                                                            // The distance used for bmu-related stuff.
		   0,                                                               // Neighborhood key.
		   evolution);  
    
      // Display
    
      image = cv::Scalar(255, 255, 255);
      
      // Display the points and the graph
      std::copy(S.begin(), S.end(), dd);
      g.foreach_edge(draw_edge); 
      g.foreach_vertex(draw_vertex);
    }
    
    cv::imshow("algorithm", image );
    cv::imshow("params",    params);
    keycode = cv::waitKey(wait_ms) & 0xFF;
  }
  
  return 0;
}

