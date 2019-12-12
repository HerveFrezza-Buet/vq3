#include <vector>
#include <string>
#include <iterator>
#include <cstdlib>

#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

// This examples illustrates the use of frames in order to convert
// mathematical point into image pixel locations.


#define HOUSE_RADIUS .50
#define HOUSE_ROOF   .75

struct Drawing {
  cv::Mat image;
  std::string title;
  vq3::demo2d::opencv::Frame frame;
  std::vector<vq3::demo2d::Point> house;
  vq3::demo2d::Point& click_pos;
  
  Drawing(unsigned int image_width,
	  unsigned int image_height,
	  vq3::demo2d::Point& click_pos,
	  const std::string& title,
	  const vq3::demo2d::opencv::Frame& frame)
    : image(image_height, image_width, CV_8UC3, cv::Scalar(255,255,255)),
      title(title),
      frame(frame),
      click_pos(click_pos) {
    auto out = std::back_inserter(house);
      *(out++) = {+HOUSE_RADIUS, -HOUSE_RADIUS};
      *(out++) = {-HOUSE_RADIUS, +HOUSE_RADIUS};
      *(out++) = {-HOUSE_RADIUS, -HOUSE_RADIUS};
      *(out++) = {+HOUSE_RADIUS, -HOUSE_RADIUS};
      *(out++) = {+HOUSE_RADIUS, +HOUSE_RADIUS};
      *(out++) = {            0,  HOUSE_ROOF  };
      *(out++) = {-HOUSE_RADIUS, +HOUSE_RADIUS};
      *(out++) = {+HOUSE_RADIUS, +HOUSE_RADIUS};
      *(out++) = {-HOUSE_RADIUS, -HOUSE_RADIUS};
  }

  void draw() {
    image = cv::Scalar(255, 255, 255);
    auto prev = house[0];
    for(auto& pt : house) {
      cv::line(image, frame(prev), frame(pt), cv::Scalar(0, 0, 0), 3);
      prev = pt;
    }
    cv::circle(image, frame(vq3::demo2d::Point(+HOUSE_RADIUS, -HOUSE_RADIUS)), 10, cv::Scalar(0, 0, 0), -1);
    cv::circle(image, frame(vq3::demo2d::Point(+HOUSE_RADIUS, +HOUSE_RADIUS)), 10, cv::Scalar(0, 0, 0), -1);

    // Let us draw the reference frame.
    cv::line(image, frame(vq3::demo2d::Point(0,0)), frame(vq3::demo2d::Point(0,1)), cv::Scalar(0, 0, 255), 3);
    cv::line(image, frame(vq3::demo2d::Point(0,0)), frame(vq3::demo2d::Point(1,0)), cv::Scalar(255, 0, 0), 3);

    cv::circle(image, frame(click_pos), 5, cv::Scalar(0, 255, 0), -1);

    cv::imshow(title, image);
  }

  void click(int x, int y) {
    click_pos = frame(cv::Point(x, y));
  }
};

void on_mouse( int event, int x, int y, int, void* user_data) {
  auto& drawing = *(reinterpret_cast<Drawing*>(user_data));
  drawing.click(x, y);
}
    

int main(int argc, char* argv[]) {

  if(argc != 3) {
    std::cout << "Usage : " << argv[0] << " <img-width> <img-height>" << std::endl;
    return 0;
  }

  unsigned int img_width  = std::atoi(argv[1]);
  unsigned int img_height = std::atoi(argv[2]);

  vq3::demo2d::Point click_pos(0,0);
  
  std::vector<Drawing> drawings;
  auto out = std::back_inserter(drawings);

  *(out++) = Drawing(img_width, img_height, click_pos,
		     "direct orthonormal centered",
		     vq3::demo2d::opencv::direct_orthonormal_frame(cv::Size(img_width, img_height),
								   .5*img_height,
								   true));
  *(out++) = Drawing(img_width, img_height, click_pos,
		     "direct orthonormal not centered",
		     vq3::demo2d::opencv::direct_orthonormal_frame(cv::Size(img_width, img_height),
								   .5*img_height,
								   false));
  *(out++) = Drawing(img_width, img_height, click_pos,
		     "direct orthonormal bbox",
		     vq3::demo2d::opencv::direct_orthonormal_frame(cv::Size(img_width, img_height),
								   vq3::demo2d::sample::BBox({-HOUSE_RADIUS, HOUSE_RADIUS}, {0, HOUSE_ROOF}),
								   10));
		     
  *(out++) = Drawing(img_width, img_height, click_pos,
		     "direct orthogonal",
		     vq3::demo2d::opencv::direct_orthogonal_frame(200,
								  50,
								  {200, 400}));
  *(out++) = Drawing(img_width, img_height, click_pos,
		     "Frame",
		     vq3::demo2d::opencv::Frame(vq3::demo2d::Point(200, 100),
						vq3::demo2d::Point(150, -50),
						vq3::demo2d::Point( 30, 120)));
		     

  for(auto& drawing : drawings) {
    cv::namedWindow(drawing.title, CV_WINDOW_AUTOSIZE);
    cv::setMouseCallback(drawing.title, on_mouse, reinterpret_cast<void*>(&drawing));
  }

  int keycode = 0;
  while(keycode != 27) {
    for(auto& drawing : drawings)
      drawing.draw();
    keycode = cv::waitKey(10) & 0xFF;
  }
  
  return 0;
}
