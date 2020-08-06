#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>

#define NB_VERTICES_PER_M2    500
#define NB_SAMPLES_PER_M2   10000

#define GRID_WIDTH        40
#define GRID_HEIGHT       30

// This is what we want for vertex value.
struct ScalarAt {
  demo2d::Point pos;
  double value;

  ScalarAt() : pos(), value(0) {}
  ScalarAt(const demo2d::Point& pos) : pos(pos), value(0) {}
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


using graph = vq3::graph<vertex, void>;

/* We can register global neighborhood computations in a single
   topology table.  This is why each computation is referred by a
   key. Here, let us use strings as identifiers for distances. */
using key_type = std::string;


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const demo2d::Point& p) {return demo2d::d2(v.vq3_value.pos, p);}

struct callback_data {
  demo2d::opencv::Frame&            frame;
  vq3::topology::Table<graph, key_type>& topology;

  callback_data(demo2d::opencv::Frame& frame,
		vq3::topology::Table<graph, key_type>& topology)
    : frame(frame),
      topology(topology) {}
  callback_data()                                = delete;
  callback_data(const callback_data&)            = delete;
  callback_data& operator=(const callback_data&) = delete;
};

void on_mouse(int event, int x, int y, int, void* user_data);

int main(int argc, char* argv[]) {
  std::random_device rd;  
  std::mt19937 random_device(rd());

  // Let us define a graph with two connected component. The first one
  // is obtained by using the Competitive Hebbian Learning algorithm
  // on a coronal shaped distribution. The second one is a grid.

  graph g; 

  vq3::utils::make_grid(g, GRID_WIDTH, GRID_HEIGHT,
			[](unsigned int w, unsigned int h) {return ScalarAt(demo2d::Point(w, h));},
			false, true); // Do not loop on with, loop on height.
  

  auto topology = vq3::topology::table<key_type>(g);

  // let us register different neighbourhood distances. Keys are strings here (see key_type definition).
  topology.declare_distance("linear truncated",   [](unsigned int edge_distance) {return         1-edge_distance*.05;}, 10, 0.00);
  topology.declare_distance("exponential",        [](unsigned int edge_distance) {return std::pow(.9, edge_distance);}, std::numeric_limits<unsigned int>::max(), 0.01);
  topology.declare_distance("piecewise constant", [](unsigned int edge_distance) {return       .25*(edge_distance/5);}, 20, 0.00);

  // As the graph topology won't change, we can compute all
  // neighborhoods for all distances at once, and then use it withaout
  // any further computation.
  topology.update_full(); 
  
  // Let us draw the graph
  
  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  auto image        = cv::Mat(480, 640, CV_8UC3, cv::Scalar(50, 50, 50));
  double unit_size  = 640/(GRID_WIDTH+1.0);
  auto frame        = demo2d::opencv::direct_orthogonal_frame(unit_size, unit_size,
							      demo2d::Point(unit_size, 480 - unit_size));
  callback_data user_data(frame, topology);
  cv::setMouseCallback("image", on_mouse, reinterpret_cast<void*>(&user_data));
  
  auto draw_edge = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
  								     [](const vertex& v1, const vertex& v2) {return   true;},  // always draw
  								     [](const vertex& v) {return           v.vq3_value.pos;},  // position
  								     []()                {return       cv::Scalar(0, 0, 0);},  // color
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
	    << "click in the image (try left, middle and right button), press ESC to quit." << std::endl
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

void on_mouse(int event, int x, int y, int, void* user_data) {
  std::string key;

  switch(event) {
  case cv::EVENT_LBUTTONDOWN : key = "linear truncated";   break;
  case cv::EVENT_MBUTTONDOWN : key = "exponential";        break;
  case cv::EVENT_RBUTTONDOWN : key = "piecewise constant"; break;
  default: return;
  }

  std::cout << "Using distance \"" << key << "\"." << std::endl;
  
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
  auto n = data.topology.neighborhood(ref_v, key); // This is an acces to the pre-computed neighborhood information.

  // Now, from value/index pairs in the neighborhood, let us change the nodes.
  for(auto& info : n) {
    auto& ref_v = data.topology(info.index);
    (*(ref_v))().vq3_value.value = info.value;
  }
}
