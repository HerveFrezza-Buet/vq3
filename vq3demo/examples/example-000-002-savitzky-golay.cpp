#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>
#include <deque>
#include <algorithm>

#include <cmath>

#define INPUT_HISTORY_SIZE 100
#define PERIOD_MS            1
#define OMEGA                5
#define SPEED_SCALE        .05
#define ACCEL_SCALE       .005

demo2d::Point curve(double t) {
  double tt = OMEGA*t;
  return {std::cos(tt), std::sin(2*tt)};
}

int main(int argc, char* argv[]) {
  std::mt19937 random_device(0);

  std::deque<demo2d::Point> history;
  int noise_radius = 0;
  
  auto image = cv::Mat(480, 640, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = demo2d::opencv::direct_orthonormal_frame(image.size(), image.size().width*.3, true);
  auto dd    = demo2d::opencv::dot_drawer<demo2d::Point>(image, frame,
								   [](const demo2d::Point& pt) {return                  true;},
								   [](const demo2d::Point& pt) {return                    pt;},
								   [](const demo2d::Point& pt) {return                     1;},
								   [](const demo2d::Point& pt) {return cv::Scalar(200, 0, 0);},
								   [](const demo2d::Point& pt) {return                    -1;});
  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  cv::createTrackbar("10000*noise_radius", "image", &noise_radius, 1000, nullptr);

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
  double t = 0;
  auto estimator = vq3::utils::savitzky_golay::constant_timestep::estimator<demo2d::Point, 2, 21, 2>();
  estimator.set_timestep(PERIOD_MS*1e-3);
  
  while(keycode != 27) {
    auto sample = demo2d::alter(random_device, curve(t), noise_radius*1e-4);
    history.push_back(sample);
    if(history.size() > INPUT_HISTORY_SIZE)
      history.pop_front();

    // Let us update the savitsky-golay estimator.
    estimator += sample; 

    image = cv::Scalar(255, 255, 255);
    std::copy(history.begin(), history.end(), dd);
    
    // Let us draw the position and speed if they are available.
    if(auto& pos = estimator.get<0>(); pos) { // If the position is available. Pos is an std::optional<Point>.
      auto P = pos.value();
      cv::circle(image, frame(P), 5, cv::Scalar(0, 0, 200), -1);
      if(auto& speed = estimator.get<1>(); speed) {
  	auto dP = speed.value();
  	auto PP = P - dP*SPEED_SCALE; // The inverse of the speed is rather plotted.
  	cv::line(image, frame(P), frame(PP), cv::Scalar(0, 0, 200), 3);
      }
      if(auto& acceleration = estimator.get<2>(); acceleration) {
  	auto ddP = acceleration.value();
  	auto PPP = P + ddP*ACCEL_SCALE;
  	cv::line(image, frame(P), frame(PPP), cv::Scalar(0, 200, 0), 1);
      }
    }
    
    cv::imshow("image", image);
    keycode = cv::waitKey(PERIOD_MS) & 0xFF;
    t += PERIOD_MS*1e-3;
  }

  return 0;
}
