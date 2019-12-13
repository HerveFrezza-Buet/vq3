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
#define PIXEL_MARGIN 10

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
    // The frame transforms mathematical points a pixel positions on
    // the image.
    
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
    // From a pixel position, the frame gives the corresponding point
    // in the mathematical frame.
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


  // We draw a "house" made of a square with a triangle (roof) on the
  // top of it. The square side is 1, its center is (0,0). The drawing
  // are expressed in a "mathematical" frame, i.e. an orthonormal
  // frame where coordinates are floating point values.
  //
  // On an image, the frame coordinate of pixel positions is
  // indirect. The coordinates are integers, (0,0) is the top left
  // corner of the image. The x-axis points to the right and the y
  // axis points downwards. A point expressed in this "image" frame is
  // denoted by an "image coordinate".
  //
  // The mathematical frame (i.e. the x and w unit vectors at the
  // origin) is displayed on the image.
  

  // Here, the frame is such as mathematical (0,0) is at the center of
  // the image. The length 1 in the mathematical frame corresponds to
  // .5*img_height pixels in the image. 
  *(out++) = Drawing(img_width, img_height, click_pos,
		     "direct orthonormal centered",
		     vq3::demo2d::opencv::direct_orthonormal_frame(cv::Size(img_width, img_height),
								   .5*img_height,
								   true));
  
  // This is same as previously, except that mathematical (0,0) is at
  // the bottom left pixel of the image.
  *(out++) = Drawing(img_width, img_height, click_pos,
		     "direct orthonormal not centered",
		     vq3::demo2d::opencv::direct_orthonormal_frame(cv::Size(img_width, img_height),
								   .5*img_height,
								   false));

  // This frame is such as the bounding box in the mathematical frame
  // fits the image (leaving a border margin), keeping the aspect
  // ration.
  *(out++) = Drawing(img_width, img_height, click_pos,
		     "direct orthonormal bbox",
		     vq3::demo2d::opencv::direct_orthonormal_frame(cv::Size(img_width, img_height),
								   vq3::demo2d::sample::BBox({-HOUSE_RADIUS, HOUSE_RADIUS}, {0, HOUSE_ROOF}),
								   PIXEL_MARGIN));

  // The origin of the mathematical frame, here, is at pixel (200,
  // 400). The length of the x unit vector is 200, the length of the y
  // unit vector is 50.
  *(out++) = Drawing(img_width, img_height, click_pos,
		     "direct orthogonal",
		     vq3::demo2d::opencv::direct_orthogonal_frame(200,
								  50,
								  {200, 400}));
  // The origin of the mathematical frame is at pixel (200, 100). The
  // x unit vector of the mathematical frame is transformed into a
  // vector (200, -50) in the image frame. The y unit vector of the
  // mathematical frame is transformed into a vector (30, 120) in the
  // image frame.
  *(out++) = Drawing(img_width, img_height, click_pos,
		     "Frame",
		     vq3::demo2d::opencv::Frame({200, 100},
						{100, -50},
						{ 30, 120}));
		     

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
