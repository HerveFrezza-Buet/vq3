#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <vector>
#include <iterator>
#include <algorithm>
#include <random>

/*

  Here, we build up and run a Kohonen self-organizing map from
  scratch, using general vq3 features.

  In this introductive example, we use vq3 for implementing the online
  version of SOMs. Keep in mind that the library is rather oriented to
  batch mode, as next examples show.

*/



#define NB_SAMPLES_PER_M2   10000
#define GRID_WIDTH             20
#define GRID_HEIGHT            20

#define SOM_H_RADIUS            5.1
#define SOM_MAX_DIST            (unsigned int)(SOM_H_RADIUS)

#define ALPHA                  .05

#define NB_SAMPLES_PER_FRAME   10


using sample    = demo2d::Point;
using prototype = demo2d::Point;

//                                                                ## Node properties :
using layer_0 = prototype;                                        // prototypes are 2D points (this is the "user defined" value).
using layer_1 = vq3::decorator::tagged<layer_0>;                  // we add a tag for topology computation.
using layer_2 = vq3::decorator::grid_pos<layer_1>;                // we add the ability to register a grid position.
using vertex  = layer_2;

using graph  = vq3::graph<vertex, void>;

using topology_key_type = int;


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const demo2d::Point& p) {return demo2d::d2(v.vq3_value, p);}


struct callback_data {
  graph& g;
  std::mt19937& rd;

  callback_data(graph& g, std::mt19937& rd) : g(g), rd(rd) {}
  callback_data()                                = delete;
  callback_data(const callback_data&)            = delete;
  callback_data& operator=(const callback_data&) = delete;
};

void on_mouse( int event, int x, int y, int, void* user_data) {
  // A mouse click reinitializes the graph.
  if(event != cv::EVENT_LBUTTONDOWN )
    return;
  auto& data = *(reinterpret_cast<callback_data*>(user_data));
  data.g.foreach_vertex([&data](const graph::ref_vertex& v_ref) {(*v_ref)().vq3_value = demo2d::uniform(data.rd, {-.5, -.5}, {.5, .5});});
}

int main(int argc, char* argv[]) {
  std::random_device rd;  
  std::mt19937 random_device(rd());

  graph g;
  
  vq3::utils::make_grid(g, GRID_WIDTH, GRID_HEIGHT,
			[&random_device](unsigned int w, unsigned int h) {
			  graph::vertex_value_type value(demo2d::uniform(random_device, {-.5, -.5}, {.5, .5}));
			  value.vq3_gridpos = {w, h};
			  return value;
			});

  auto topology = vq3::topology::table<topology_key_type>(g);
  topology.declare_distance(0, // This the key of the single neighborhood declared in this example.
			    [](unsigned int edge_distance) {return std::max(0., 1 - edge_distance/double(SOM_H_RADIUS));},
			    SOM_MAX_DIST,
			    1e-3); // We consider node and edge-based neihborhoods.
  topology.update_full(); // Update the neighborhoods.... neighborhood #0 here.

  bool wtm_mode = true; // false means wta.
  

  double intensity = 1. ;
  double radius    =  .5;
  double hole      =  .3;
  auto density     = demo2d::sample::disk(radius, intensity) - demo2d::sample::disk(hole, intensity);

  
  auto sampler = demo2d::sample::base_sampler::random(random_device, NB_SAMPLES_PER_M2);
  auto S_      = demo2d::sample::sample_set(random_device, sampler, density); // This re-toss points at each begin...end iteration.
  std::vector<demo2d::Point> S; // Let us use a single sample of S_.
  auto out = std::back_inserter(S);
  std::copy(S_.begin(), S_.end(), out);

  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  auto image       = cv::Mat(480, 640, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = demo2d::opencv::direct_orthonormal_frame(image.size(), .9*image.size().height, true);

  callback_data user_data(g, random_device);
  cv::setMouseCallback("image", on_mouse, reinterpret_cast<void*>(&user_data));
  auto dd          = demo2d::opencv::dot_drawer<demo2d::Point>(image, frame,
							       [](const demo2d::Point& pt) {return                      true;},
							       [](const demo2d::Point& pt) {return                        pt;},
							       [](const demo2d::Point& pt) {return                         1;},
							       [](const demo2d::Point& pt) {return cv::Scalar(200, 200, 200);},
							       [](const demo2d::Point& pt) {return                        -1;});
  auto draw_edge   = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
								       [](const vertex& v1, const vertex& v2) {return true;}, // always draw
								       [](const vertex& v)     {return         v.vq3_value;}, // position
								       []()                    {return cv::Scalar(0, 0, 0);}, // color
								       []()                    {return                  1;}); // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return                                                   true;},
									   [](const vertex& v) {return                                            v.vq3_value;},  // position
									   [](const vertex& v) {return                                                      3;},  // radius
									   [](const vertex& v) {return cv::Scalar(0,
														  255*v.vq3_gridpos.first/(GRID_WIDTH - 1),
														  255*v.vq3_gridpos.second/(GRID_HEIGHT - 1));},  // color
									   [](const vertex& v) {return                                                     -1;}); // thickness


  
  std::cout << std::endl
	    << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << "click in the image to restart, press ESC to quit. 'space' toggles wta/wtm mode." << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl;
  int keycode = 0;
  auto sample_it = S.begin();
  graph::ref_vertex bmu;
  
  while(keycode != 27) {

    if(wtm_mode)
      for(unsigned int i=0; i < NB_SAMPLES_PER_FRAME; ++i) {
	bmu = vq3::online::wtm::learn(topology, 0, // 0 is our neighbourhood key.
				      [](vertex& vertex_value) -> demo2d::Point& {return vertex_value.vq3_value;},
				      d2,
				      *(sample_it++),
				      ALPHA);
	if(sample_it == S.end()) sample_it=S.begin();
      }
    else
      for(unsigned int i=0; i < NB_SAMPLES_PER_FRAME; ++i) {
	bmu = vq3::online::wta::learn(g, 
				      [](vertex& vertex_value) -> demo2d::Point& {return vertex_value.vq3_value;},
				      d2,
				      *(sample_it++),
				      ALPHA);
	if(sample_it == S.end()) sample_it=S.begin();
      }

    
    image = cv::Scalar(255, 255, 255);

    cv::circle(image, frame((*bmu)().vq3_value), 10, cv::Scalar(0,0,0), -1);
    
    std::copy(S.begin(), S.end(), dd);
    if(wtm_mode)
      g.foreach_edge(draw_edge); 
    g.foreach_vertex(draw_vertex);
    cv::imshow("image", image);
    keycode = cv::waitKey(1) & 0xFF;
    if(keycode == 32)
      wtm_mode = !wtm_mode;

    if(++sample_it == S.end())
      sample_it = S.begin();
  }
  
  
  
  return 0;
}
