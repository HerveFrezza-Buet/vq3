#include <vq3demo.hpp>
#include <random>

int main(int argc, char* argv[]) {
  std::random_device rd;  
  std::mt19937 random_device(rd());
  
  int N_slider =  5000;
  int T_slider =   127;

  auto video_data = vq3::demo2d::opencv::sample::video_data(0, [&T_slider](const unsigned char* rgb_pixel) {if(rgb_pixel[1] < (unsigned char)(T_slider)) return 1.; else return 0.;});

  auto webcam_distribution = vq3::demo2d::opencv::sample::webcam(video_data);

  auto input_size  = video_data.image.size();
  video_data.frame = vq3::demo2d::opencv::direct_orthonormal_frame(input_size, .5*input_size.width, true);

  cv::namedWindow("webcam", CV_WINDOW_AUTOSIZE);
  
  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  cv::createTrackbar("nb/m^2", "image", &N_slider, 50000, nullptr);
  cv::createTrackbar("threshold", "image", &T_slider, 255, nullptr);
  
  auto image       = cv::Mat(300, 400, CV_8UC3, cv::Scalar(255,255,255));
  auto frame       = vq3::demo2d::opencv::direct_orthonormal_frame(image.size(), .4*image.size().width, true);
  auto dd          = vq3::demo2d::opencv::dot_drawer<vq3::demo2d::Point>(image, frame,
									 [](const vq3::demo2d::Point& pt) {return                      true;},
									 [](const vq3::demo2d::Point& pt) {return                        pt;},
									 [](const vq3::demo2d::Point& pt) {return                         1;},
									 [](const vq3::demo2d::Point& pt) {return cv::Scalar(200, 200, 200);},
									 [](const vq3::demo2d::Point& pt) {return                        -1;});
  
  std::cout << std::endl
	    << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << "click in the image to restart, press ESC to quit." << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl;
  
  int keycode = 0;
  while(keycode != 27) {
    ++video_data; // get next frame.
    auto S = vq3::demo2d::sample::sample_set(random_device, webcam_distribution, N_slider);
    image = cv::Scalar(255, 255, 255);
    std::copy(S.begin(), S.end(), dd);
    cv::imshow("image", image);
    cv::imshow("webcam", video_data.image);
    keycode = cv::waitKey(1) & 0xFF;
  }
  
  return 0;
}

