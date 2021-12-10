#include <iostream>
#include <iomanip>
#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>



// This example shows how to perform multi-threaded counts on vertices
// or edges.  The purpose is to match each sample againts all the
// vertices (resp. edges) of the graph with a "matters" boolean
// function. At the end, the vertices (resp. edges) have counted the
// number of positive match.

// The parallelization is made on the samples.



// Graph definition
//
///////////////


using vertex = vq3::decorator::counter<demo2d::Point>;
using edge   = vq3::decorator::counter<void>;
using graph  = vq3::graph<vertex, edge>;

#define NB_SAMPLES_PER_M2   1000
#define RADIUS                .5

void print_counts(graph::ref_vertex ref_O,
		  graph::ref_vertex ref_A,
		  graph::ref_vertex ref_B,
		  graph::ref_vertex ref_C,
		  graph::ref_vertex ref_D) {
  std::cout << std::setw(5) << 'O'
	    << std::setw(5) << 'A'
	    << std::setw(5) << 'B'
	    << std::setw(5) << 'C'
	    << std::setw(5) << 'D'
	    << std::endl;
  std::cout << std::setw(5) << (*ref_O)().vq3_counter
	    << std::setw(5) << (*ref_A)().vq3_counter
	    << std::setw(5) << (*ref_B)().vq3_counter
	    << std::setw(5) << (*ref_C)().vq3_counter
	    << std::setw(5) << (*ref_D)().vq3_counter
	    << std::endl;
  std::cout << std::endl;
}

int main(int argc, char* argv[]) {
  if(argc != 2) {
    std::cout << "Usage : " << argv[0] << " nb_threads." << std::endl;
    return 0;
  }

  unsigned int nb_threads = std::atoi(argv[1]);
  
  std::random_device rd;  
  std::mt19937 random_device(rd());
  
  cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
  auto image       = cv::Mat(800, 800, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = demo2d::opencv::direct_orthonormal_frame(image.size(), .48*image.size().width, true);

  double side = 2;
  double intensity = 1.0;
  auto distrib = demo2d::sample::rectangle(side, side, intensity);
  auto sampler = demo2d::sample::base_sampler::random(random_device, NB_SAMPLES_PER_M2);
  std::vector<demo2d::Point> S;

  
  graph g;
  auto O = (g += demo2d::Point(  0,   0));
  auto A = (g += demo2d::Point(-.8,  .8));
  auto B = (g += demo2d::Point( .8,  .8));
  auto C = (g += demo2d::Point( .8, -.8));
  auto D = (g += demo2d::Point(-.8, -.8));
  
  auto dd = demo2d::opencv::dot_drawer<demo2d::Point>(image, frame,
						      [](const demo2d::Point& pt) {return                      true;},  // always draw
						      [](const demo2d::Point& pt) {return                        pt;},  // position
						      [](const demo2d::Point& pt) {return                         1;},  // radius
						      [](const demo2d::Point& pt) {return cv::Scalar(200, 200, 200);},  // color
						      [](const demo2d::Point& pt) {return                        -1;}); // thickness
  
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return                true;},  // always draw
									   [](const vertex& v) {return         v.vq3_value;},  // position
									   [](const vertex& v) {return                   5;},  // radius
									   [](const vertex& v) {return cv::Scalar(0, 0, 0);},  // color
									   [](const vertex& v) {return                  -1;}); // thickness
  
  auto draw_area = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									 [](const vertex& v)                {return                true;},  // always draw
									 [](const vertex& v)                {return         v.vq3_value;},  // position
									 [r=frame(RADIUS)](const vertex& v) {return                   r;},  // radius
									 [](const vertex& v)                {return cv::Scalar(0, 0, 0);},  // color
									 [](const vertex& v)                {return                   1;}); // thickness

  // Counting processors
  auto v_counter = vq3::epoch::count::vertices::processor(g);
  
  // Display
  //
  //////////
  
  std::cout << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << "press ESC to quit." << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << std::endl;

  int keycode = 0;
  
  while(keycode != 27) {     
    image = cv::Scalar(255, 255, 255);
    S.clear();
    auto S_ = demo2d::sample::sample_set(random_device, sampler, distrib);
    auto out = std::back_inserter(S);
    std::copy(S_.begin(), S_.end(), out);

    vq3::utils::clear_all_counters(g, 0);
    v_counter.process(nb_threads,
		      S.begin(), S.end(),
		      [](auto& content) {return content;}, // sample_of(*it) --> sample
		      [](auto& ref_v, const demo2d::Point& sample) { // do the sample matter for that vertex ?
			return demo2d::d2((*ref_v)().vq3_value, sample) < RADIUS * RADIUS;
		      });

    print_counts(O, A, B, C, D);
    
    std::copy(S.begin(), S.end(), dd);
    g.foreach_vertex(draw_vertex);
    g.foreach_vertex(draw_area);
    cv::imshow("image", image);
    keycode = cv::waitKey(100) & 0xFF;
  }
    
  return 0;
}
  
