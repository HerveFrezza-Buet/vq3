#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <vector>
#include <iterator>
#include <random>
#include <fstream>

using vertex = vq3::demo2d::Point;
using edge   = double;
using graph  = vq3::graph<vertex, edge>;

/* 
   Warning :

   The vertex or edge type is required to permitt a serialization by <<
   and >> operators. You may need to implement these, when values are
   decorated for example.

*/

#define NB_VERTICES  20
#define NB_EDGES     50
#define OFFSET       .02
int main(int argc, char* argv[]) {

  std::random_device rd;  
  std::mt19937 random_device(rd());
  
  vq3::demo2d::sample::BBox bbox(1, 1);
  std::vector<graph::ref_vertex> vertices;
  graph g1, g2;

  auto out = std::back_inserter(vertices);
  for(unsigned int i=0; i < NB_VERTICES; ++i)
    *(out++) = g1 += bbox.uniform(rd); // We add vertices made of random points.

  auto random_value = std::uniform_real_distribution<double>(0., 1.);
  for(unsigned int i=0; i < NB_EDGES; ++i) {
    auto [ref_v1, ref_v2] = vq3::utils::two_closest(g1, bbox.uniform(rd), vq3::demo2d::d2);
    g1.connect(ref_v1, ref_v2, random_value(random_device));
  }

  {
    std::ofstream file("example.gph");
    file << g1;
  }

  {
    std::ifstream file("example.gph");
    file >> g2;
  }

  g2.foreach_vertex([](auto ref_v){(*ref_v)() += vq3::demo2d::Point(OFFSET, -OFFSET);});
  
  
  auto image = cv::Mat(600, 600, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .9*image.size().width, true);
  
  auto draw_vertex_1 = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									     [](const vertex& v) {return                  true;},  // always draw
									     [](const vertex& v) {return                     v;},  // position
									     [](const vertex& v) {return                     7;},  // radius
									     [](const vertex& v) {return cv::Scalar(0, 0, 200);},  // color
									     [](const vertex& v) {return                    -1;}); // thickness
  
  auto draw_edge_1 = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
								       [](const vertex&, const vertex&, const edge&) {return       true;},  // always draw
								       [](const vertex& v) {return                                    v;},  // position
								       [](const edge  e)   {int x = 255*e; return cv::Scalar(x, x, 255);},  // color
								       [](const edge  e)   {return                                    5;}); // thickness
  
  auto draw_vertex_2 = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									     [](const vertex& v) {return                  true;},  // always draw
									     [](const vertex& v) {return                     v;},  // position
									     [](const vertex& v) {return                     7;},  // radius
									     [](const vertex& v) {return cv::Scalar(200, 0, 0);},  // color
									     [](const vertex& v) {return                    -1;}); // thickness
  
  auto draw_edge_2 = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
								       [](const vertex&, const vertex&, const edge&) {return       true;},  // always draw
								       [](const vertex& v) {return                                    v;},  // position
								       [](const edge  e)   {int x = 255*e; return cv::Scalar(255, x, x);},  // color
								       [](const edge  e)   {return                                    5;}); // thickness
  

  g1.foreach_edge(draw_edge_1);
  g1.foreach_vertex(draw_vertex_1);
  g2.foreach_edge(draw_edge_2);
  g2.foreach_vertex(draw_vertex_2);
  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  cv::imshow     ("image", image);
  cv::waitKey(0);
  
  return 0;
}
