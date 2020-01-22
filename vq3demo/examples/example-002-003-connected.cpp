#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>
#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>

// Let us define the graph-related types. We need decorators to
// perform connected component computation.

using user_vertex    = demo2d::Point;
using vertex_layer_1 = vq3::decorator::efficiency<user_vertex>;       // for connected components
using vertex_layer_2 = vq3::decorator::tagged<vertex_layer_1>;        // for connected components
using vertex_layer_3 = vq3::decorator::labelled<vertex_layer_2>;      // for vertex labelling
using vertex_layer_4 = vq3::demo::decorator::colored<vertex_layer_3>;
using vertex         = vertex_layer_4;

using edge_layer_1   = vq3::decorator::efficiency<void>;              // for connected components
using edge_layer_2   = vq3::decorator::labelled<edge_layer_1>;        // for edge labelling
using edge_layer_3   = vq3::demo::decorator::colored<edge_layer_2>;
using edge           = edge_layer_3;

using graph = vq3::graph<vertex, edge>;


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const demo2d::Point& p) {return demo2d::d2(v.vq3_value, p);}


auto on_circle(double degree) {
  double rad = degree*(3.141592653589793/180.0);
  return demo2d::Point(std::cos(rad), std::sin(rad));
}

#define NB_POINTS 75

void labelling(graph& g, bool conservative) {
  auto components = vq3::connected_components::make(g);
  if(conservative)
    vq3::labelling::conservative(components.begin(), components.end());
  else
    vq3::labelling::basic(components.begin(), components.end());
  vq3::labelling::edges_from_vertices(g);
}

int main(int argc, char* argv[]) {
  if(argc != 2) {
    std::cout << std::endl
	      << "Usage : " << argv[0] << " [basic | conservative]" << std::endl
	      << std::endl;
    return 0;
  }

  std::string mode = argv[1];
  bool conservative;
  if(mode == "conservative")
    conservative = true;
  else if(mode == "basic")
    conservative = false;
  else {
    std::cout << "Argument \"" << mode << "\" is not allowed. Use \"basic\" or \"conservative\"." << std::endl;
    return 1;
  }

  
  std::random_device rd;  
  std::mt19937 random_device(rd());

  auto color_of_label = demo2d::opencv::colormap::random(random_device);

  // Let us define a crown-shaped graph.

  graph g; 

  double coef = 360.0/NB_POINTS;
  auto p0   = g += on_circle(0);
  auto prev = p0;
  for(unsigned int i=1; i < NB_POINTS; ++i) {
    auto current = g += on_circle(i*coef);
    g.connect(prev, current);
    prev = current;
  }
  g.connect(prev,p0);

  // We set everything efficient.
  vq3::utils::clear_all_efficiencies(g, true);

  // Let us draw the graph
  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  auto image       = cv::Mat(480, 640, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = demo2d::opencv::direct_orthonormal_frame(image.size(), .3*image.size().width, true);
  
  auto draw_edge   = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
  								       [](const vertex& v1, const vertex& v2, const edge& e) {return true;}, // always draw
  								       [](const vertex& v) {return v.vq3_value;},                            // position
  								       [&color_of_label](const edge& e)   {
									 if(e.vq3_efficient)
									   return color_of_label(e.vq3_label);
									 else
									   return cv::Scalar(200, 200, 200);},                               // color
  								       [](const edge& e)   {
									 if(e.vq3_efficient)
									   return 3;
									 else
									   return 1;});                                                      // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
  									   [](const vertex& v) {return true;},        // always draw
  									   [](const vertex& v) {return v.vq3_value;}, // position
  									   [](const vertex& v) {return 5;},           // radius
  									   [&color_of_label](const vertex& v) {
									     if(v.vq3_efficient)
									       return color_of_label(v.vq3_label);
									     else
									       return cv::Scalar(200, 200, 200);},    // color
  									   [](const vertex& v) {return -1;});         // thickness
  


  std::cout << std::endl
	    << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << "press ESC to quit," << std::endl
	    << "press <space> successively to perform steps." << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl;
  

  int keycode = 0;

  // Let us fill the list of all efficient edges.
  std::list<graph::ref_edge> efficient_edges;
  vq3::utils::collect_edges(g, std::back_inserter(efficient_edges));

  // Let us label the graph at first.
  labelling(g, conservative);
  
  while(keycode != 27) {
    if(keycode == 32) {
      if(efficient_edges.size() > 0) {
	// Let is pick an efficient edge randomly.
	auto it = efficient_edges.begin();
	std::advance(it, std::uniform_int_distribution<unsigned int>(0, efficient_edges.size()-1)(random_device));
	(*(*it))().vq3_efficient = false;
	efficient_edges.erase(it);

	labelling(g, conservative);
      }
    }
    
    image = cv::Scalar(255,255,255);
    g.foreach_edge(draw_edge); 
    g.foreach_vertex(draw_vertex);
    cv::imshow("image", image);
    keycode = cv::waitKey(100) & 0xFF;

  }


  

  return 0;

}

