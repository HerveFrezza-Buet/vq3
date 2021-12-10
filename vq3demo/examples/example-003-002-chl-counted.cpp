#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <optional>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <limits>

// This example shows how to perform multi-threaded competitive
// Hebbian learning... with the addition ouf counting features.



// Graph definition
//
///////////////

using vertex = demo2d::Point;

using elayer_0 = vq3::decorator::tagged<void>;      // We need tags on the edges.
using elayer_1 = vq3::decorator::counter<elayer_0>; // We add an atomic counter on the edges.
using edge     = elayer_1;

using graph  = vq3::graph<vertex, edge>;


// Distance
//
////////////////

// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const demo2d::Point& p) {return demo2d::d2(v, p);}

auto color_of_count(std::size_t counts, std::size_t max) {
  unsigned char level = (unsigned char)(255*(1-counts/double(max))+.5);
  return cv::Scalar(255, level, level);
}

struct callback_data {
private:
  
  mutable bool new_pointer;
  const demo2d::opencv::Frame& frame;

public:
  
  std::optional<demo2d::Point> pointer;
  callback_data(const demo2d::opencv::Frame& frame) : new_pointer(false), frame(frame), pointer() {}

  void unclick()           {if(pointer) {pointer = std::nullopt; new_pointer = true;}            }
  void click(int x, int y) {pointer  = frame(cv::Point(x,y));    new_pointer = true;             }
  operator bool() const    {auto res = new_pointer;              new_pointer = false; return res;}
};

void on_mouse( int event, int x, int y, int, void* user_data) {
  auto& data = *(reinterpret_cast<callback_data*>(user_data));
  if(event == cv::EVENT_LBUTTONDOWN) data.click(x, y);
  if(event == cv::EVENT_RBUTTONDOWN) data.unclick();
}

// Main
//
////////////////

#define NB_VERTICES_PER_M2     50
#define NB_SAMPLES_PER_M2   50000

graph::ref_edge select_edge(graph& g, const demo2d::Point location);
void select_samples(graph& g, const std::vector<demo2d::Point>& S,
		    graph::ref_edge selected_edge,
		    std::vector<demo2d::Point>& selected_samples);

int main(int argc, char* argv[]) {
  if(argc != 4) {
    std::cout << "Usage : " << argv[0] << " nb_threads <random | triangles | grid> <uniform | gradient>" << std::endl;
    return 0;
  }

  unsigned int nb_threads = std::atoi(argv[1]);
  std::string graph_kind  = argv[2];
  std::string dist_kind   = argv[3];
  
  std::random_device rd;  
  std::mt19937 random_device(rd());
  
  graph::ref_edge selected_edge = nullptr;

  // Image settings
  //
  ///////////////////

  
  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  auto image       = cv::Mat(800, 800, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = demo2d::opencv::direct_orthonormal_frame(image.size(), .48*image.size().width, true);
  callback_data cb {frame};
  cv::setMouseCallback("image", on_mouse, reinterpret_cast<void*>(&cb));
  auto dd          = demo2d::opencv::dot_drawer<demo2d::Point>(image, frame,
									 [](const demo2d::Point& pt) {return                      true;},
									 [](const demo2d::Point& pt) {return                        pt;},
									 [](const demo2d::Point& pt) {return                         1;},
									 [](const demo2d::Point& pt) {return cv::Scalar(200, 200, 200);},
									 [](const demo2d::Point& pt) {return                        -1;});
  
  auto dd_selected = demo2d::opencv::dot_drawer<demo2d::Point>(image, frame,
							       [](const demo2d::Point& pt) {return                  true;},
							       [](const demo2d::Point& pt) {return                    pt;},
							       [](const demo2d::Point& pt) {return                     2;},
							       [](const demo2d::Point& pt) {return cv::Scalar(0, 0, 200);},
							       [](const demo2d::Point& pt) {return                    -1;});
  std::size_t max_count  = 0;
  auto draw_edge         = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
									     [](const vertex& v1, const vertex& v2, const edge& e) {return                                       true;},  // always draw
									     [](const vertex& v)                                   {return                                          v;},  // position
									     [&max_count](const edge&   e)                         {return color_of_count(e.vq3_counter, max_count);},  // color
									     [](const edge&   e)                                   {return                                          3;}); // thickness
  auto draw_edge_contour = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
									     [](const vertex& v1, const vertex& v2, const edge& e) {return                 true;},  // always draw
									     [](const vertex& v)                                   {return                    v;},  // position
									     [&max_count](const edge&   e)                         {return cv::Scalar(0, 0 , 0);},  // color
									     [](const edge&   e)                                   {return                    5;}); // thickness
  auto draw_edge_selected_contour = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
										      [](const vertex& v1, const vertex& v2, const edge& e) {return                 true;},  // always draw
										      [](const vertex& v)                                   {return                    v;},  // position
										      [&max_count](const edge&   e)                         {return cv::Scalar(0, 0 , 0);},  // color
										      [](const edge&   e)                                   {return                    9;}); // thickness
  
  // Initializations
  //
  //////////////////

  double side = 2;
  double intensity = 1.0;
  auto background = demo2d::sample::rectangle(side, side, intensity);

  graph g;
  if(graph_kind == "random") {
    auto v_sampler = demo2d::sample::base_sampler::random(random_device, NB_VERTICES_PER_M2);
    for(auto v : demo2d::sample::sample_set(random_device, v_sampler, background)) g += v;
  }
  else if(graph_kind == "triangles") {
    auto v_sampler = demo2d::sample::base_sampler::triangles(random_device, NB_VERTICES_PER_M2);
    for(auto v : demo2d::sample::sample_set(random_device, v_sampler, background)) g += v;
  }
  else if(graph_kind == "grid") {
    auto v_sampler = demo2d::sample::base_sampler::grid(random_device, NB_VERTICES_PER_M2);
    for(auto v : demo2d::sample::sample_set(random_device, v_sampler, background)) g += v;
  }
  else {
    std::cout << "Graph kind must be in {random, triangles, grid}, \"" << graph_kind << "\" provided." << std::endl;
    ::exit(0);
  }

  auto samples = background;
  if(dist_kind == "gradient")
    samples = demo2d::sample::custom(demo2d::sample::BBox(-1, -1, 2, 2), [](const demo2d::Point& p) {return .5*(p.x + 1.0);});
  std::vector<demo2d::Point> S;
  auto s_sampler = demo2d::sample::base_sampler::random(random_device, NB_SAMPLES_PER_M2);
  auto S_ = demo2d::sample::sample_set(random_device, s_sampler, samples);
  auto out = std::back_inserter(S);
  std::copy(S_.begin(), S_.end(), out);

  // Computation
  //
  //////////////

  auto chl = vq3::epoch::chl::processor_with_counts(g);
  chl.process(nb_threads,
	      S.begin(), S.end(),
	      [](const demo2d::Point& s) {return s;}, // Gets the sample from *it.
	      d2,                                     // d2(prototype, sample).
	      edge());                                // New edge initialization value.

  g.foreach_edge([&max_count](auto ref_e) {
		   auto extr = ref_e->extremities();
		   if(vq3::invalid_extremities(extr))
		     ref_e->kill();
		   else if(std::size_t c = (*ref_e)().vq3_counter; c > max_count)
		     max_count = c;
		 });
    
  
  // Display
  //
  //////////
  
  std::cout << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << "left-click near the edges." << std::endl
	    << "right-click to cancel." << std::endl
	    << "press ESC to quit." << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << std::endl;
    
  
  int keycode = 0;
  std::vector<demo2d::Point> selected_samples;
  while(keycode != 27) {     
    image = cv::Scalar(255, 255, 255);

    if(cb) {
      if(cb.pointer) {
	selected_edge = select_edge(g, *(cb.pointer));
	select_samples(g, S, selected_edge, selected_samples);
      }
      else {
	selected_edge = nullptr;
	selected_samples.clear();
      }
    }
    
    std::copy(S.begin(), S.end(), dd);
    std::copy(selected_samples.begin(), selected_samples.end(), dd_selected);
    g.foreach_edge(draw_edge_contour);
    if(selected_edge) draw_edge_selected_contour(selected_edge);
    g.foreach_edge(draw_edge); 
    cv::imshow("image", image);
    keycode = cv::waitKey(100) & 0xFF;
  }
    
  return 0;
}
  
graph::ref_edge select_edge(graph& g, const demo2d::Point location) {
  graph::ref_edge res;
  double min_dist = std::numeric_limits<double>::max();

  g.foreach_edge([&res, &min_dist, &location](auto ref_e) {
		   auto [A, B] = ref_e->extremities();
		   if(auto d = demo2d::d2(location, {(*A)(), (*B)()}); d < min_dist) {
		     min_dist = d;
		     res = ref_e;
		   }
		 });
  return res;
}

void select_samples(graph& g,
		    const std::vector<demo2d::Point>& S,
		    graph::ref_edge selected_edge,
		    std::vector<demo2d::Point>& selected_samples) {
  selected_samples.clear();
  for(auto& M : S)
    if(auto [ref_A, ref_B] = vq3::utils::two_closest(g, M, d2); g.get_edge(ref_A, ref_B) == selected_edge)
      selected_samples.push_back(M);
}
