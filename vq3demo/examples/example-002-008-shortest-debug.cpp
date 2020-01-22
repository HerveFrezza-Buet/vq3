#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <iterator>
#include <utility>
#include <iostream>
#include <fstream>

#define NB_VERTICES                5000
#define NB_SAMPLES_PER_M2_SUPPORT 50000

// Graph definition
//
///////////////

enum class Status : char {Source = 's', Destination = 'd', Computed = 'c', None = 'n'};

auto color_of(Status s) {
  switch(s) {
  case Status::Source :
    return cv::Scalar( 80,  80, 200);
  case Status::Destination :
    return cv::Scalar( 80, 200,  80);
  case Status::Computed :
    return cv::Scalar(200,  80,  80);
  case Status::None :
    return cv::Scalar(0, 0, 0);
  }
  return cv::Scalar(0, 0, 0);
}

auto radius_of(Status s) {
  switch(s) {
  case Status::Source :
    return 11;
  case Status::Destination :
    return 11;
  case Status::Computed :
    return 3;
  case Status::None :
    return 0;
  }
  return 0;
}

// Vertex values have to be instrumented with shortest path
// computation structure.
//                                                         ## Node properties :
using vlayer_0 = demo2d::Point;                            // prototypes are 2D points (this is the "user defined" value).
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
double d2(const vertex& v, const demo2d::Point& p) {return demo2d::d2(v.vq3_value, p);}




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
  
  auto density = demo2d::sample::rectangle(side, side, intensity)
    -  demo2d::sample::rectangle(hole_side, hole_side, intensity);
  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  auto image = cv::Mat(480, 1280, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = demo2d::opencv::direct_orthonormal_frame(image.size(), .9*image.size().height, true);

  auto draw_edge = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
  								     [](const vertex&, const vertex&, const edge&) {return   true;},  // always draw
  								     [](const vertex& v) {return                      v.vq3_value;},  // position
  								     [](const edge& e)   {return        cv::Scalar(200, 200, 200);},  // color
  								     [](const edge& e)   {return                                1;}); // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return                       true;},  // always draw
									   [](const vertex& v) {return                v.vq3_value;},  // position
									   [](const vertex& v) {return    radius_of(v.vq3_custom);},  // radius
									   [](const vertex& v) {return     color_of(v.vq3_custom);},  // color
									   [](const vertex& v) {return                         -1;}); // thickness


  graph g1,g2;

  // We need to specify the serialization functions at graph creation.

  g1.vertex_to_stream = [](std::ostream& os, const vertex& v) {os << v.vq3_value << std::endl;};
  g1.edge_to_stream   = [](std::ostream& os, const edge&   e) {os << e.vq3_value << std::endl;};
  
  g2.vertex_from_stream = [](std::istream& is, vertex& v) {char c; is >> v.vq3_value; is.get(c);};
  g2.edge_from_stream   = [](std::istream& is, edge&   e) {char c; is >> e.vq3_value; is.get(c);;};
  

  for(unsigned int i = 0; i < NB_VERTICES; ++i)
    g1 += demo2d::sample::get_one_sample(random_device, density);

  
  // Setting edges of the support graph with CHL.
  {
    std::vector<demo2d::Point> S;
    auto sampler = demo2d::sample::base_sampler::triangles(random_device, NB_SAMPLES_PER_M2_SUPPORT);
    auto S_      = demo2d::sample::sample_set(random_device, sampler, density);
    auto out     = std::back_inserter(S);
    std::copy(S_.begin(), S_.end(), out);
    
    auto chl = vq3::epoch::chl::processor(g1);
    chl.process(std::thread::hardware_concurrency(),
                S.begin(), S.end(),
                [](const demo2d::Point& s) {return s;}, // Gets the point from *it.
                d2,                                          // d2(vertex, sample).
		0);                                          // New edge initialization value.

    // We set the edge cost to the Euclidian distance between vertices.
    g1.foreach_edge([](auto ref_e) {
	auto [v1, v2] = ref_e->extremities();
	(*ref_e)().vq3_value = demo2d::d((*v1)().vq3_value, (*v2)().vq3_value);
      });
  }

  // we create g1 as a copy of g1.... we use the serialization into a
  // file for that, which is a bit tricky.
  {
    std::ofstream file("tmp.gph");
    file << g1;
  }
  
  {
    std::ifstream file("tmp.gph");
    file >> g2;
  }

  // Then, we move all the vertices of the two graphs.
  auto shift = demo2d::Point(.75, 0);
  g1.foreach_vertex([s = -shift](auto ref_v){(*ref_v)().vq3_value += s;});
  g2.foreach_vertex([s =  shift](auto ref_v){(*ref_v)().vq3_value += s;});
  
  
  // These vectors store the vertices in the order of their
  // computation by the two shortest path algoritms.
  std::vector<graph::ref_vertex> dijkstra_vertices;
  std::vector<graph::ref_vertex> a_star_vertices;

  // Let us apply dijkstra and a* to the graphs, with similar nodes.
  auto source_locus             = demo2d::Point(0,   1);
  auto destination_locus        = demo2d::Point(-.3, -.3);
  auto source_1                 = vq3::utils::closest(g1, source_locus - shift, d2);
  auto source_2                 = vq3::utils::closest(g2, source_locus + shift, d2);
  auto destination_1            = vq3::utils::closest(g1, destination_locus - shift, d2);
  auto destination_2            = vq3::utils::closest(g2, destination_locus + shift, d2);
  (*source_1)().vq3_custom      = Status::Source;
  (*source_2)().vq3_custom      = Status::Source;
  (*destination_1)().vq3_custom = Status::Destination;
  (*destination_2)().vq3_custom = Status::Destination;

  // We collect the vertices in their order of computation by the shortest path algorithms.
  vq3::path::dijkstra<false, false>(g1, source_1, destination_1,
				    [](auto& ref_e){return (*ref_e)().vq3_value;},
				    std::back_inserter(dijkstra_vertices));
  
  vq3::path::a_star<false, false>(g2, source_2, destination_2,
				  [](auto& ref_e){return (*ref_e)().vq3_value;},
				  [start = source_2](const auto& ref_v){ // This is the heuristic
				   if(start)
				     return demo2d::d((*start)().vq3_value, (*ref_v)().vq3_value); // direct distance is lower than real cost.
				   else
				     return 0.0;
				  },
				  std::back_inserter(a_star_vertices));
 
  
  // The iterators used for display.
  auto current_dijkstra = dijkstra_vertices.begin();
  auto current_a_star    = a_star_vertices.begin();
  g1.foreach_vertex([](auto ref_v){(*ref_v)().vq3_custom = Status::None;});
  g2.foreach_vertex([](auto ref_v){(*ref_v)().vq3_custom = Status::None;});

  
  std::cout << std::endl
            << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl
	    << "Press ESC to quit." << std::endl
            << std::endl
            << "##################" << std::endl
            << std::endl;
  int keycode = 0;

  while(keycode != 27) {    

    if(current_dijkstra == dijkstra_vertices.end() && current_a_star == a_star_vertices.end()) {
      // We initialize the statue used for display.
      g1.foreach_vertex([](auto ref_v){(*ref_v)().vq3_custom = Status::None;});
      g2.foreach_vertex([](auto ref_v){(*ref_v)().vq3_custom = Status::None;});
      current_dijkstra = dijkstra_vertices.begin();
      current_a_star = a_star_vertices.begin();
    }    
    if(current_dijkstra !=  dijkstra_vertices.end())
      (*(*(current_dijkstra++)))().vq3_custom = Status::Computed;
    if(current_a_star !=  a_star_vertices.end())
      (*(*(current_a_star++)))().vq3_custom = Status::Computed;
    
    (*source_1)().vq3_custom = Status::Source;
    (*source_2)().vq3_custom = Status::Source;
    (*destination_1)().vq3_custom = Status::Destination;
    (*destination_2)().vq3_custom = Status::Destination;
    
    image = cv::Scalar(255, 255, 255);
    g1.foreach_edge(draw_edge); 
    g1.foreach_vertex(draw_vertex);
    g2.foreach_edge(draw_edge); 
    g2.foreach_vertex(draw_vertex);
    cv::imshow("image", image);
    keycode = cv::waitKey(1) & 0xFF;
  }
  return 0;
}
