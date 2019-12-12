#include <vector>
#include <string>

#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

// This examples illustrates the use of frames in order to convert
// mathematical point into image pixel locations.

#define IMG_WIDTH  640
#define IMG_HEIGHT 480

#define HOUSE_RADIUS .5
#define HOUSE_ROOF    1
void draw_house(const std::string& title,
		const vq3::demo2d::opencv::Frame& frame) {
  std::vector<vq3::demo2d::Point> house {
    {+HOUSE_RADIUS, -HOUSE_RADIUS},
    {-HOUSE_RADIUS, +HOUSE_RADIUS},
    {-HOUSE_RADIUS, -HOUSE_RADIUS},
    {+HOUSE_RADIUS, -HOUSE_RADIUS},
    {+HOUSE_RADIUS, +HOUSE_RADIUS},
    {            0,  HOUSE_ROOF  },
    {-HOUSE_RADIUS, +HOUSE_RADIUS},
    {+HOUSE_RADIUS, +HOUSE_RADIUS},
    {-HOUSE_RADIUS, -HOUSE_RADIUS}
  };
  
  auto image = cv::Mat(IMG_HEIGHT, IMG_WIDTH, CV_8UC3, cv::Scalar(255,255,255));
  auto prev = house[0];
  for(auto& pt : house) {
    cv::line(image, frame(prev), frame(pt), cv::Scalar(0, 0, 0), 3);
    prev = pt;
  }
  cv::circle(image, frame(vq3::demo2d::Point(+HOUSE_RADIUS, -HOUSE_RADIUS)), 10, cv::Scalar(0, 0, 0), -1);
  cv::circle(image, frame(vq3::demo2d::Point(+HOUSE_RADIUS, +HOUSE_RADIUS)), 10, cv::Scalar(0, 0, 0), -1);
  cv::circle(image, frame(vq3::demo2d::Point(            0,             0)), 20, cv::Scalar(0, 0, 0),  3);

  cv::imshow(title, image);
}

int main(int argc, char* argv[]) {

  draw_house("direct orthonormal centered", vq3::demo2d::opencv::direct_orthonormal_frame(cv::Size(IMG_WIDTH, IMG_HEIGHT),
									      .5*IMG_HEIGHT,
									      true));
  draw_house("direct orthonormal not centered", vq3::demo2d::opencv::direct_orthonormal_frame(cv::Size(IMG_WIDTH, IMG_HEIGHT),
										  .5*IMG_HEIGHT,
										  false));
  draw_house("direct orthogonal", vq3::demo2d::opencv::direct_orthogonal_frame(200,
									       50,
									       {200, 400}));
  cv::waitKey(0);
  return 0;
}
