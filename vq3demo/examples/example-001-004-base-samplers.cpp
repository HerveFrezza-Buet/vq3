#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <random>
#include <algorithm>

#include <iostream>

#define NB_SAMPLES_PER_M2 1000

int main(int argc, char* argv[]) {
  std::random_device rd;  
  std::mt19937 random_device(rd());


  auto theta = vq3::demo::dyn::linear<vq3::demo::dyn::bound::wrap>(0, 0, 360, .5);

#define SIZE 1
#define THICKNESS .2

  double fill_1   =  1;
  double fill_2   = .25;
  double radius_1 = SIZE * .7;
  double radius_2 = radius_1 - 1.5*THICKNESS;

  double rec_width  = SIZE*.5;
  double rec_height = SIZE;

  double bar_width  = SIZE+THICKNESS;
  double bar_height = THICKNESS;

  vq3::demo2d::Point o_left  = {-1.25*SIZE,    0};
  vq3::demo2d::Point o_right = {-0.75*SIZE,    0};
  vq3::demo2d::Point h       = {.5 + radius_1, 0};

  
  auto left  = vq3::demo2d::sample::rectangle(rec_width, rec_height, fill_2) + o_left;
  auto right = vq3::demo2d::sample::rectangle(rec_width, rec_height, fill_1) + o_right;
  auto bar   = vq3::demo2d::sample::rectangle(bar_width, bar_height, fill_1);
  auto disk  = vq3::demo2d::sample::disk(radius_1, fill_1);
  auto hole  = vq3::demo2d::sample::disk(radius_2, fill_1);

  auto density = (left || right || bar || ((disk - hole) + h)) % theta();
  
  auto image = cv::Mat(1000, 1000, CV_8UC3, cv::Scalar(255,255,255));
  
  auto frame = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), image.size().width/4.0, true);
  auto dd    = vq3::demo2d::opencv::dot_drawer<vq3::demo2d::Point>(image, frame,
								   [](const vq3::demo2d::Point& pt) {return                    true;},
								   [](const vq3::demo2d::Point& pt) {return                      pt;},
								   [](const vq3::demo2d::Point& pt) {return                       1;},
								   [](const vq3::demo2d::Point& pt) {return cv::Scalar(250, 50, 50);},
								   [](const vq3::demo2d::Point& pt) {return                      -1;});

  cv::namedWindow("random",    CV_WINDOW_AUTOSIZE);
  cv::namedWindow("grid",      CV_WINDOW_AUTOSIZE);
  cv::namedWindow("triangles", CV_WINDOW_AUTOSIZE);
  
  int keycode = 0;
 
  std::cout << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << "press ESC to quit." << std::endl
	    << "press <space> to toggle rotation." << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl;
  
  auto sampler_random    = vq3::demo2d::sample::base_sampler::random   (random_device, NB_SAMPLES_PER_M2);
  auto sampler_grid      = vq3::demo2d::sample::base_sampler::grid     (random_device, NB_SAMPLES_PER_M2);
  auto sampler_triangles = vq3::demo2d::sample::base_sampler::triangles(random_device, NB_SAMPLES_PER_M2);

  bool rotate = true;
  while(keycode != 27) {
    
    image = cv::Scalar(255,255,255);
    auto S1 = vq3::demo2d::sample::sample_set(random_device, sampler_random, density);
    std::copy(S1.begin(), S1.end(), dd);
    cv::imshow("random", image);
    
    image = cv::Scalar(255,255,255);
    auto S2 = vq3::demo2d::sample::sample_set(random_device, sampler_grid, density);
    std::copy(S2.begin(), S2.end(), dd);
    cv::imshow("grid", image);
    
    image = cv::Scalar(255,255,255);
    auto S3 = vq3::demo2d::sample::sample_set(random_device, sampler_triangles, density);
    std::copy(S3.begin(), S3.end(), dd);
    cv::imshow("triangles", image);

    
    keycode = cv::waitKey(10) & 0xFF;
    if(keycode == 32)
      rotate = ! rotate;
    if(rotate) ++theta;
  }

  return 0;
}
