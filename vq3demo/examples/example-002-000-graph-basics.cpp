#include <utility>
#include <string>

#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>


using vertex = std::pair<char, demo2d::Point>;
using edge   = double;
using graph             = vq3::graph<vertex, edge>;

int main(int argc, char* argv[]) {
  graph g;

  auto A = g += std::make_pair('A', demo2d::Point(-.5, -.5));
  auto B = g += std::make_pair('B', demo2d::Point( .5, -.5));
  auto C = g += std::make_pair('C', demo2d::Point( .5,  .5));

  auto length = [](auto a, auto b) {return demo2d::d((*a)().second, (*b)().second);};
  
  g.connect(A, B, length(A, B));
  g.connect(A, C, length(A, C));
  g.connect(B, C, length(B, C));

  
  auto image = cv::Mat(600, 600, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = demo2d::opencv::direct_orthonormal_frame(image.size(), .6*image.size().width, true);
  
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return                  true;},  // always draw
									   [](const vertex& v) {return              v.second;},  // position
									   [](const vertex& v) {return                     5;},  // radius
									   [](const vertex& v) {return cv::Scalar(200, 0, 0);},  // color
									   [](const vertex& v) {return                    -1;}); // thickness
  
  auto print_vertex = vq3::demo2d::opencv::vertex_printer<graph::ref_vertex>(image, frame,
									     [](const vertex& v) {return                    true;},  // always draw
									     [](const vertex& v) {return std::string(1, v.first);},  // text
									     [](const vertex& v) {return                v.second;},  // position
									     [](const vertex& v) {return   cv::Scalar(0, 0, 255);},  // color
									     [](const vertex& v) {return                       2;},  // thickness
									     [](const vertex& v) {return std::make_pair(10, -10);},  // text offset
									     [](const vertex& v) {return                      1.;}); // font scale
  
  auto draw_edge = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
  								     [](const vertex&, const vertex&, const edge&) {return   true;},  // always draw
  								     [](const vertex& v) {return                         v.second;},  // position
  								     [](const edge& e)   {return              cv::Scalar(0, 0, 0);},  // color
  								     [](const edge& e)   {return                                1;}); // thickness
  
  auto print_edge = vq3::demo2d::opencv::edge_printer<graph::ref_edge>(image, frame,
								       [](const vertex&, const vertex&, const edge&) {return   true;},  // always draw
								       [](const edge& e)   {return                std::to_string(e);},  // text
								       [](const vertex& v) {return                         v.second;},  // position
								       [](const edge& e)   {return        cv::Scalar(255, 180, 100);},  // color
								       [](const edge& e)   {return                                1;},  // thickness
								       [](const edge& e)   {return          std::make_pair(10, -10);},  // text offset
								       [](const edge& e)   {return                               .5;}); // font scale
  
  
  // This is the drawing....
  g.foreach_edge(draw_edge);
  g.foreach_vertex(draw_vertex);
  g.foreach_edge(print_edge);
  g.foreach_vertex(print_vertex);
  
  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  cv::imshow ("image", image);
  cv::waitKey(0);


  return 0;
}
  
