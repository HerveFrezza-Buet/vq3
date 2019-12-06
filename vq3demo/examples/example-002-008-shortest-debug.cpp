#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <utility>
#include <iostream>
#include <fstream>

#define NB_VERTICES                 200
#define NB_SAMPLES_PER_M2_SUPPORT 10000

// Graph definition
//
///////////////

enum class Status : char {Source = 's', Destination = 'd', Computed = 'c'};

// Vertex values have to be instrumented with shortest path
// computation structure.
//                                                         ## Node properties :
using vlayer_0 = vq3::demo2d::Point;                       // prototypes are 2D points (this is the "user defined" value).
using vlayer_1 = vq3::decorator::path::shortest<vlayer_0>; // This holds informations built by dijkstra.
using vlayer_2 = vq3::decorator::custom<vlayer_1, Status>; // We will tag some nodes for display.
using vertex   = vlayer_2;

//                                                         ## Edge properties :
using elayer_0 = double;                                   // Edge cost.
using elayer_1 = vq3::decorator::tagged<elayer_0>;         // We need a tag for CHL.
using edge     = elayer_1;

using graph   = vq3::graph<vertex, edge>;


// Distance
//
////////////////


// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with points.
double d2(const vertex& v, const vq3::demo2d::Point& p) {return vq3::demo2d::d2(v.vq3_value, p);}




// Main
//
//////////////

#define NB_VERTICES_PER_M2      200
#define NB_SAMPLES_PER_M2    100000

int main(int argc, char* argv[]) {

  std::random_device rd;  
  std::mt19937 random_device(rd());

  // Let us build the graph.
  //
  //////////////

  double intensity = 1.;
  double side      = 1.;
  double hole_side = .5;
  
  auto density = vq3::demo2d::sample::rectangle(side, side, intensity)
    -  vq3::demo2d::sample::rectangle(hole_side, hole_side, intensity);
  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  auto image = cv::Mat(480, 1280, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .9*image.size().height, true);

  auto draw_edge = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
  								     [](const vertex&, const vertex&, const edge&) {return   true;},  // always draw
  								     [](const vertex& v) {return                      v.vq3_value;},  // position
  								     [](const edge& e)   {return        cv::Scalar(200, 200, 200);},  // color
  								     [](const edge& e)   {return                                1;}); // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return                  true;},  // always draw
									   [](const vertex& v) {return           v.vq3_value;},  // position
									   [](const vertex& v) {return                     3;},  // radius
									   [](const vertex& v) {return cv::Scalar(200,  80,  80);},   // color
									   [](const vertex& v) {return                   -1;});  // thickness


  graph g1,g2;

  // We need to specify the serialization functions at graph creation.

  g1.vertex_to_stream = [](std::ostream& os, const vertex& v) {os << v.vq3_value << std::endl;};
  g1.edge_to_stream   = [](std::ostream& os, const edge&   e) {os << e.vq3_value << std::endl;};
  
  g2.vertex_from_stream = [](std::istream& is, vertex& v) {char c; is >> v.vq3_value; is.get(c);};
  g2.edge_from_stream   = [](std::istream& is, edge&   e) {char c; is >> e.vq3_value; is.get(c);;};
  

  for(unsigned int i = 0; i < NB_VERTICES; ++i)
    g1 += vq3::demo2d::sample::get_one_sample(random_device, density);

  
  // Setting edges of the support graph with CHL.
  {
    std::vector<vq3::demo2d::Point> S;
    auto sampler = vq3::demo2d::sample::base_sampler::triangles(random_device, NB_SAMPLES_PER_M2_SUPPORT);
    auto S_      = vq3::demo2d::sample::sample_set(random_device, sampler, density);
    auto out     = std::back_inserter(S);
    std::copy(S_.begin(), S_.end(), out);
    
    auto chl = vq3::epoch::chl::processor(g1);
    chl.process(std::thread::hardware_concurrency(),
                S.begin(), S.end(),
                [](const vq3::demo2d::Point& s) {return s;}, // Gets the point from *it.
                d2,                                          // d2(vertex, sample).
		0);                                          // New edge initialization value.

    // We set the edge cost to the Euclidian distance between vertices.
    g1.foreach_edge([](auto ref_e) {
	auto [v1, v2] = ref_e->extremities();
	(*ref_e)().vq3_value = vq3::demo2d::d((*v1)().vq3_value, (*v2)().vq3_value);
      });
  }

  // we create g1 as a copy of g1.... we use the serialization into a
  // file for that, which is a bit tricky.
  {
    std::ofstream file("tmp.gph");
    file << g1;
  }
  
  // {
  //   std::ifstream file("tmp.gph");
  //   file >> g2;
  // }

  // Then, we move all the vertices of the two graphs.
  auto shift = vq3::demo2d::Point(.75, 0);
  g1.foreach_vertex([s = -shift](auto ref_v){(*ref_v)().vq3_value += s;});
  g2.foreach_vertex([s =  shift](auto ref_v){(*ref_v)().vq3_value += s;});

  
  
    
  std::cout << std::endl
            << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl
	    << "Press any key to display the next step." << std::endl
	    << "Press r to restart." << std::endl
	    << "Press ESC to quit." << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl;
  int keycode = 0;
  while(keycode != 27) {    
    image = cv::Scalar(255, 255, 255);
    g1.foreach_edge(draw_edge); 
    g1.foreach_vertex(draw_vertex);
    g2.foreach_edge(draw_edge); 
    g2.foreach_vertex(draw_vertex);
    cv::imshow("image", image);
    keycode = cv::waitKey(0) & 0xFF;
  }
  return 0;
}
