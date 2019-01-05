#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>
#include <iostream>
#include <cstdlib>
#include <string>

#define NB_VERTICES_PER_M2    500
#define NB_SAMPLES_PER_M2   10000

#define GRID_WIDTH        28
#define GRID_HEIGHT       21
#define GRID_SQUARE_SIZE .04

// This is what we want for vertex value.
struct ScalarAt {
  vq3::demo2d::Point pos;
  double value;

  ScalarAt() : pos(), value(0) {}
  ScalarAt(const vq3::demo2d::Point& pos) : pos(pos), value(0) {}
  ScalarAt(const ScalarAt&)            = default;
  ScalarAt& operator=(const ScalarAt&) = default;
};

// The topological manipulations on the graph requires te vertices to
// be tagged. This is why our actual vertex values need to be
// decorated with a tag. vq3 offer templates for incremental
// decorations of the value type you initially provide.

// decoration...
using user_vertex   = ScalarAt;
using tagged_vertex = vq3::decorator::tagged<user_vertex>;
// ... so finally
using vertex        = tagged_vertex; 

/* Usage :
 
   graph::ref_vertex ref_v;
   vertex&           vertex_full_value = (*ref_v)(); 
   ScalarAt&         user_content      = vertex_full_value.vq3_value;
*/

using graph = vq3::graph<vertex, void>;

// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v.vq3_value.pos, p);}

struct callback_data {
  vq3::demo2d::opencv::Frame& frame;
  vq3::topology::Table<graph>& topology;

  callback_data(vq3::demo2d::opencv::Frame& frame,
		vq3::topology::Table<graph>& topology)
    : frame(frame),
      topology(topology) {}
  callback_data()                                = delete;
  callback_data(const callback_data&)            = delete;
  callback_data& operator=(const callback_data&) = delete;
  
};

void on_mouse( int event, int x, int y, int, void* user_data);

int main(int argc, char* argv[]) {
  std::random_device rd;  
  std::mt19937 random_device(rd());

  // Let us define a graph with two connected component. The first one
  // is obtained by using the Competitive Hebbian Learning algorithm
  // on a coronal shaped distribution. The second one is a grid.

  graph g; 

  // Component #1
  
  double intensity    = 1. ;
  double crown_radius =  .55;
  double hole_radius  =  .3;
  vq3::demo2d::Point crown_center{-.6, 0};

  auto density = (vq3::demo2d::sample::disk(crown_radius, intensity) - vq3::demo2d::sample::disk(hole_radius, intensity)) + crown_center;
  
  // First, we add vertices
  for(auto& prototype : vq3::demo2d::sample::sample_set(random_device, density, NB_VERTICES_PER_M2)) g += ScalarAt(prototype);

  // Then, we sample points, and connect the two closest prototypes (if not connected yet).
  for(auto& sample : vq3::demo2d::sample::sample_set(random_device, density, NB_SAMPLES_PER_M2)) {
    auto closest = vq3::utils::two_closest(g, sample, d2);
    if(g.get_edge(closest.first, closest.second) == nullptr) 
      g.connect(closest.first, closest.second);
  }

  // Component #2

  vq3::utils::make_grid(g, GRID_WIDTH, GRID_HEIGHT,
			[](unsigned int w, unsigned int h) {return ScalarAt(vq3::demo2d::Point(w *GRID_SQUARE_SIZE, h * GRID_SQUARE_SIZE));});
  
  // If edges had values, we would have used
  /*
    vq3::utils::make_grid(g, GRID_WIDTH, GRID_HEIGHT,
    [](unsigned int w, unsigned int h) {return ScalarAt(vq3::demo2d::Point(w *GRID_SQUARE_SIZE, h * GRID_SQUARE_SIZE));},
    [](unsigned int w, unsigned int h, unsigned int ww, unsigned int hh) {return some_edge_value;});
  */

  auto topology = vq3::topology::table(g);
  topology(); // We only inform the topology table about vertices, ignoring edge-based neighborhoods.
  
  // Let us draw the graph
  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  auto image        = cv::Mat(480, 640, CV_8UC3, cv::Scalar(255,255,255));
  auto frame        = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .4*image.size().width, true);
  callback_data user_data(frame, topology);
  cv::setMouseCallback("image", on_mouse, reinterpret_cast<void*>(&user_data));
  
  auto draw_edge = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
  								     [](const vertex& v1, const vertex& v2) {return   true;},  // always draw
  								     [](const vertex& v) {return           v.vq3_value.pos;},  // position
  								     []()                {return cv::Scalar(200, 200, 200);},  // color
  								     []()                {return                         1;}); // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return                                      true;},  // always draw
									   [](const vertex& v) {return                           v.vq3_value.pos;},  // position
									   [](const vertex& v) {return                                         4;},  // radius
									   [](const vertex& v) {return cv::Scalar(0, 255 * v.vq3_value.value, 0);},  // color
									   [](const vertex& v) {return                                        -1;}); // thickness
  




  
  std::cout << std::endl
	    << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << "click in the image, press ESC to quit." << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl;
  int keycode = 0;
  while(keycode != 27) {
    g.foreach_edge(draw_edge); 
    g.foreach_vertex(draw_vertex);
    cv::imshow     ("image", image);
    keycode = cv::waitKey(100) & 0xFF;
  }

  return 0;

}

void on_mouse( int event, int x, int y, int, void* user_data) {
  if(event != cv::EVENT_LBUTTONDOWN )
    return;
  auto& data = *(reinterpret_cast<callback_data*>(user_data));

  // Get the coordinates of the clicked point.
  auto click_pos = data.frame(cv::Point(x,y));

  // Let us clear the graph values.
  data.topology.g.foreach_vertex([](const graph::ref_vertex& ref_v) {(*ref_v)().vq3_value.value = 0;});

  // Now, let us find the vertex which is the closest to click_pos.
  auto ref_v = vq3::utils::closest(data.topology.g, click_pos, d2);

  // Now, let us compute its neighborhood. This is why we have
  // tag-decorated our node values. 
  vq3::utils::clear_vertex_tags(data.topology.g, false);
  auto n = data.topology.neighborhood(ref_v,
				      [](unsigned int edge_distance) {return std::pow(.9, edge_distance);},
				      0,
				      0.0);

  // Now, from value/index pairs in the neighborhood, let us change the nodes.
  for(auto& info : n) {
    auto& ref_v = data.topology(info.index);
    (*(ref_v))().vq3_value.value = info.value;
  }
  

}
