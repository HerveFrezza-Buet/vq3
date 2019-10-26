
/*
 *   Copyright (C) 2018,  CentraleSupelec
 *
 *   Author : Hervé Frezza-Buet
 *
 *   Contributor :
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU General Public
 *   License (GPL) as published by the Free Software Foundation; either
 *   version 3 of the License, or any later version.
 *   
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *   General Public License for more details.
 *   
 *   You should have received a copy of the GNU General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 *   Contact : herve.frezza-buet@centralesupelec.fr
 *
 */


#pragma once

#include <opencv2/opencv.hpp>
#include <utility>
#include <functional>
#include <optional>
#include <algorithm>
#include <map>
#include <random>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <array>
#include <tuple>
#include <string>

#include <vq3.hpp>
#include <vq3demoPoint.hpp>
#include <vq3demoSample.hpp>
#include <algorithm>


namespace vq3 {
  
  namespace demo2d {
    namespace opencv {
      
      class HueSelector {
      private:
	cv::Mat hsv;
	cv::Mat hsv_line;
	int hsv_thickness;
	int stride;
	
	void check_line(int rows, int cols) {
	  if(cols == hsv_line.cols)
	    return;

	  hsv_line.create(1, cols, CV_8UC3);

	  unsigned char* it  = hsv_line.data;
	  unsigned char* end = hsv_line.data + 3*cols;
	  int col = 0;
	  while(it != end) {
	    *(it++) = (unsigned char)(180.0*(col++)/(cols-1)+.5); // H
	    *(it++) = 255;                                        // S
	    *(it++) = 255;                                        // V
	  }
	  cv::cvtColor(hsv_line, hsv_line, CV_HSV2BGR);
	  stride = cols*3;
	  hsv_thickness = std::min(cols, 20);
	}

	bool hue_test(unsigned int h) const {
	  int href = (int)(H_slider*180e-3+.5);
	  int dh   = href - (int)h;
	  dh       = std::min(std::abs(dh), std::abs(dh+180));
	  return dh < T_slider;
	}
	
	bool hsv_test(const unsigned char* hsv_pixel) const {
	  unsigned char s = (unsigned char)(S_slider*255e-3+.5);
	  unsigned char v = (unsigned char)(V_slider*255e-3+.5);

	  return hsv_pixel[1] > s
	    &&   hsv_pixel[2] > v
	    &&   hue_test(hsv_pixel[0]);
	}
	
      public:
	cv::Mat image;
	int H_slider =   0;
	int S_slider = 500;
	int V_slider = 500;
	int T_slider =  20;
	double darken = .2;

	HueSelector() = default;


	auto build_pixel_test() {
	  return [this](const unsigned char* rgb_pixel) {
	    cv::Mat src(1, 1, CV_8UC3, (void*)rgb_pixel);
	    cv::Mat dst;
	    cv::cvtColor(src, dst, CV_BGR2HSV);
	    return this->hsv_test(dst.data);
	  };
	}
	
	void build_sliders(const std::string& window_name) {
	  cv::createTrackbar("Hue",            window_name, &H_slider, 1000, nullptr);
	  cv::createTrackbar("Hue Tolerance",  window_name, &T_slider,  100, nullptr);
	  cv::createTrackbar("Min Saturation", window_name, &S_slider, 1000, nullptr);
	  cv::createTrackbar("Min Value",      window_name, &V_slider, 1000, nullptr);
	}

	void build_image(cv::Mat img) {
	  check_line(img.rows, img.cols);
	  image.create(img.rows, img.cols, CV_8UC3);
	  cv::cvtColor(img, hsv, CV_BGR2HSV);
	  
	  unsigned char* src = img.data;
	  unsigned char* dst = image.data;
	  unsigned char* end = src + img.rows*stride;

	  for(int l = 0; l < hsv_thickness; ++l, dst += stride)
	    std::copy(hsv_line.data, hsv_line.data + stride, dst);
	  auto offset = hsv_thickness*stride;
	  src += offset;

	  double h = H_slider*1e-3*img.cols;
	  cv::line(image, cv::Point(h, 0), cv::Point(h, hsv_thickness), cv::Scalar(0,0,0), 5);
	  
	  for(unsigned char* hsv_it = hsv.data + offset; src != end; hsv_it +=3) 
	    if(hsv_test(hsv_it)) {
	      *(dst++) = *(src++);
	      *(dst++) = *(src++);
	      *(dst++) = *(src++);
	    }
	    else {
	      *(dst++) = *(src++)*darken;
	      *(dst++) = *(src++)*darken;
	      *(dst++) = *(src++)*darken;
	    }
	}
      };

      
      namespace colormap {

	class jet {
	private :

	  std::array<cv::Scalar, 11> colormap = {
	    cv::Scalar(130,   0,   0),  // 0.0
	    cv::Scalar(255,   0,   0),  // 0.1
	    cv::Scalar(255,  75,   0),  // 0.2
	    cv::Scalar(255, 175,   0),  // 0.3
	    cv::Scalar(200, 255,  40),  // 0.4
	    cv::Scalar(155, 255, 155),  // 0.5
	    cv::Scalar( 40, 255, 210),  // 0.6
	    cv::Scalar(  0, 200, 255),  // 0.7
	    cv::Scalar(  0,  90, 255),  // 0.8
	    cv::Scalar(  0,   0, 240),  // 0.9
	    cv::Scalar(  0,   0, 120)}; // 1.0
	  
	  double min = 0;
	  double max = 1;
	  
	public:
	  jet()                      = default;
	  jet(const jet&)            = default;
	  jet& operator=(const jet&) = default;

	  /**
	   * Set (min, max)
	   */
	  void operator=(const std::pair<double, double>& minmax) {
	    std::tie(min,max) = minmax;
	  }

	  /**
	   * Get a color from a value.
	   */
	  cv::Scalar operator()(double value) const {
	    auto c = std::min(std::max((value - min)/(max - min), 0.0), 1-1e-3);
	    auto nb = colormap.size()-1;

	    auto cnb  = c*nb;
	    unsigned int idx_min = (unsigned int)(cnb);
	    auto dc   = cnb - idx_min;
	    auto cmin = colormap[idx_min];
	    auto cmax = colormap[idx_min+1];
	    return cv::Scalar((int)(dc*cmax[0] + (1-dc)*cmin[0] + .5),
			      (int)(dc*cmax[1] + (1-dc)*cmin[1] + .5),
			      (int)(dc*cmax[2] + (1-dc)*cmin[2] + .5));
	  }
	};
	
	template<typename RANDOM_DEVICE>
	class random {
	private:
	  RANDOM_DEVICE& rd;
	  double min_dist;
	  unsigned int nb_retry;
	  std::map<unsigned int, cv::Scalar> colors;
	  

	  cv::Scalar random_color() {
	    std::array<unsigned int, 3> idx = { 0,  1, 2};
	    std::array<double,       3> rgb = {0., 0., 0.};

	    std::shuffle(idx.begin(), idx.end(), rd);
	    rgb[idx[0]] = 0;
	    rgb[idx[1]] = 1;
	    rgb[idx[2]] = std::uniform_real_distribution<double>(0,1)(rd);

	    switch(std::uniform_int_distribution<int>(1,3)(rd)) {
	    case 1:
	      rgb[0] *= .5;
	      rgb[1] *= .5;
	      rgb[2] *= .5;
	      break;
	    case 2:
	      rgb[0] = .5 * (1 + rgb[0]);
	      rgb[1] = .5 * (1 + rgb[1]);
	      rgb[2] = .5 * (1 + rgb[2]);
	      break;
	    default:
	      break;
	    }

	    return cv::Scalar((int)(255*rgb[0]+.5),
			      (int)(255*rgb[1]+.5),
			      (int)(255*rgb[2]+.5));
	  }

	  double distance(const cv::Scalar& c1, const cv::Scalar& c2) {
	    double diff = c1[0] - c2[0];
	    double d2   = diff*diff;
	    diff        = c1[1] - c2[1];
	    d2         += diff*diff;
	    diff        = c1[2] - c2[2];
	    d2         += diff*diff;
	    return std::sqrt(d2);
	  }
	  
	  cv::Scalar random_distinct_color() {
	    if(colors.size() == 0 || nb_retry < 1)
	      return random_color();
	    
	    auto color = cv::Scalar(0, 0, 0);
	    double d = 0;
	    for(unsigned int t = 0; t < nb_retry && d < min_dist; ++t) {
	      color = random_color();
	      d = 10; // big enough
	      for(auto& id_color : colors)
		if(auto dd = distance(color, id_color.second); dd < d) {
		  d = dd;
		  if(d < min_dist)
		    break;
		}
	    }

	    return color;
	  }

	public:
	  /**
	     This color map chooses random colors. Each color is
	     supposed do be very different from the others. To do so,
	     when a new random color is needed, we try at most
	     nb_retry times a random color, until the distance with
	     this color and an existing one is grather that
	     min_dist. The distance is the 3D euclidian distance, with
	     RGB values in [0,1]^3.
	   */
	  random(RANDOM_DEVICE& rd, double min_dist=.01, unsigned int nb_retry=10)
	    : rd(rd), min_dist(min_dist*255.0), nb_retry(nb_retry), colors() {}
	  
	  random()                         = delete;
	  random(const random&)            = default;
	  random(random&&)                 = default;
	  random& operator=(const random&) = delete;
	  random& operator=(random&&)      = delete;

	  cv::Scalar operator()(unsigned int lbl) {
	    auto it = colors.find(lbl);
	    if(it != colors.end())
	      return it->second;
	    auto new_color = random_distinct_color();
	    colors[lbl] = new_color;
	    return new_color;
	  }
	  
	};
      }

      class Frame {
      private:

	demo2d::Point cv_center;
	demo2d::Point cv_Ox;
	demo2d::Point cv_Oy;

	demo2d::Point to_cv_xline;
	demo2d::Point to_cv_yline;
	
	demo2d::Point to_demo_xline;
	demo2d::Point to_demo_yline;
	
      public:

	Frame()                        = default;
	Frame(const Frame&)            = default;
	Frame& operator=(const Frame&) = default;
	/**
	 * @param cv_center The position of the frame origin (in pixel) in the CV image.
	 * @param cv_Ox The abscissa unity vector of the frame (in pixel size).
	 * @param cv_Oy The ordinate unity vector of the frame (in pixel size).
	 */
	Frame(demo2d::Point cv_center,
	      demo2d::Point cv_Ox,
	      demo2d::Point cv_Oy)
	  : cv_center(cv_center),
	    cv_Ox(cv_Ox),
	    cv_Oy(cv_Oy),
	    to_cv_xline(cv_Ox.x, cv_Oy.x),
	    to_cv_yline(cv_Ox.y, cv_Oy.y),
	    to_demo_xline( cv_Oy.y, -cv_Oy.x),
	    to_demo_yline(-cv_Ox.y,  cv_Ox.x) {
	  double coef = 1/(cv_Ox.x * cv_Oy.y - cv_Ox.y * cv_Oy.x);
	  to_demo_xline *= coef;
	  to_demo_yline *= coef;
	}

	cv::Point operator()(const demo2d::Point& p) const {
	  return cv::Point(p*to_cv_xline + cv_center.x,
			   p*to_cv_yline + cv_center.y);
	}

	
	demo2d::Point operator()(const cv::Point& p) const {
	  demo2d::Point x = {double(p.x), double(p.y)};
	  x -= cv_center;
	  return {x*to_demo_xline, x*to_demo_yline};
	}
      };


      /**
       * @param x_unit_size The size of the horizontal unit vector in pixels.
       * @param y_unit_size The size of the vertical unit vector in pixels.
       * @param origin The center of the frame in the CV image.
       */
      inline Frame direct_orthonormal_frame(double x_unit_size,
					    double y_unit_size,
					    demo2d::Point origin) {
	demo2d::Point Ox(x_unit_size, 0);
	demo2d::Point Oy(0, -y_unit_size);
	return Frame(origin, Ox, Oy);
      }
      
      /**
       * @param size The image size in pixels.
       * @param unit_size The size of the vertical and horizontal unit vectors in pixels.
       * @param centered true puts the origin at the center, false puts it at the bottom left corner.
       */
      inline Frame direct_orthonormal_frame(const cv::Size& size,
					    double unit_size,
					    bool centered) {
	demo2d::Point Ox(unit_size, 0);
	demo2d::Point Oy(0, -unit_size);
	if(centered)
	  return Frame(demo2d::Point(size.width*.5, size.height*.5), Ox, Oy);
	else
	  return Frame(demo2d::Point(0, size.height-1), Ox, Oy);
      
      }
      
      /**
       * This draws a bounding box.
       */
      inline void draw(cv::Mat& image, const Frame& frame, const vq3::demo2d::sample::BBox& bbox, const cv::Scalar& color, int thickness) {
	cv::rectangle(image, frame(bbox.bottom_left()), frame(bbox.top_right()), color, thickness);
      }


      template<typename OBJECT>
      class DotDrawer {

      private:
	
	cv::Mat image; // a share pointer.
	Frame frame;
	std::function<bool (const OBJECT&)>          do_draw;
	std::function<demo2d::Point (const OBJECT&)> point_of;
	std::function<int (const OBJECT&)>           radius_of;
	std::function<cv::Scalar (const OBJECT&)>    color_of;
	std::function<int (const OBJECT&)>           thickness_of;
	
	
      public:

        using difference_type   = long;
        using value_type        = OBJECT;
        using pointer           = OBJECT*;
        using reference         = OBJECT&;
        using iterator_category = std::output_iterator_tag;
	
	template<typename DO_DRAW, typename POINT_OF, typename RADIUS_OF, typename COLOR_OF, typename THICKNESS_OF>
	DotDrawer(cv::Mat& image,
		  Frame frame,
		  const DO_DRAW&      do_draw,
		  const POINT_OF&     point_of,
		  const RADIUS_OF&    radius_of,
		  const COLOR_OF&     color_of,
		  const THICKNESS_OF& thickness_of)
	  : image(image),
	    frame(frame),
	    do_draw(do_draw),
	    point_of(point_of),
	    radius_of(radius_of),
	    color_of(color_of),
	    thickness_of(thickness_of) {}

	DotDrawer()                            = delete;
	DotDrawer(const DotDrawer&)            = default;
	DotDrawer& operator=(const DotDrawer&) = default; 

	DotDrawer& operator++()    {return *this;}
	DotDrawer& operator++(int) {return *this;}
	DotDrawer& operator*()     {return *this;}
	DotDrawer& operator=(const OBJECT& o) {
	  if(do_draw(o))
	    cv::circle(image, frame(point_of(o)), radius_of(o), color_of(o), thickness_of(o));
	  return *this;
	}
      };

      template<typename OBJECT, typename DO_DRAW, typename POINT_OF, typename RADIUS_OF, typename COLOR_OF, typename THICKNESS_OF>
      DotDrawer<OBJECT> dot_drawer(cv::Mat& image,
				   Frame frame,
				   const DO_DRAW&      do_draw,
				   const POINT_OF&     point_of,
				   const RADIUS_OF&    radius_of,
				   const COLOR_OF&     color_of,
				   const THICKNESS_OF& thickness_of) {
	return DotDrawer<OBJECT>(image, frame, do_draw, point_of, radius_of, color_of, thickness_of);
      }
      
      template<typename OBJECT>
      class SegmentDrawer {

      private:
	
	cv::Mat image; // a share pointer.
	Frame frame;
	std::function<bool (const OBJECT&)>                                    do_draw;
	std::function<std::pair<demo2d::Point, demo2d::Point> (const OBJECT&)> points_of;
	std::function<cv::Scalar (const OBJECT&)>                              color_of;
	std::function<int (const OBJECT&)>                                     thickness_of;
	
	
      public:
	
        using difference_type   = OBJECT;
        using value_type        = OBJECT;
        using pointer           = OBJECT*;
        using reference         = OBJECT&;
        using iterator_category = std::output_iterator_tag;

	template<typename DO_DRAW, typename POINTS_OF, typename COLOR_OF, typename THICKNESS_OF>
	SegmentDrawer(cv::Mat& image,
		      Frame frame,
		      const DO_DRAW&      do_draw,
		      const POINTS_OF&    points_of,
		      const COLOR_OF&     color_of,
		      const THICKNESS_OF& thickness_of)
	  : image(image),
	    frame(frame),
	    do_draw(do_draw),
	    points_of(points_of),
	    color_of(color_of),
	    thickness_of(thickness_of) {}

	SegmentDrawer()                                = delete;
	SegmentDrawer(const SegmentDrawer&)            = default;
	SegmentDrawer& operator=(const SegmentDrawer&) = default; 

	SegmentDrawer& operator++()    {return *this;}
	SegmentDrawer& operator++(int) {return *this;}
	SegmentDrawer& operator*()     {return *this;}
	SegmentDrawer& operator=(const OBJECT& o) {
	  if(do_draw(o)) {
	    auto pts = points_of(o);
	    cv::line(image, frame(pts.first), frame(pts.second), color_of(o), thickness_of(o));
	  }
	  return *this;
	}
      };

      template<typename OBJECT, typename DO_DRAW, typename POINTS_OF, typename COLOR_OF, typename THICKNESS_OF>
      SegmentDrawer<OBJECT> segment_drawer(cv::Mat& image,
					   Frame frame,
					   const DO_DRAW&      do_draw,
					   const POINTS_OF&    points_of,
					   const COLOR_OF&     color_of,
					   const THICKNESS_OF& thickness_of) {
	return SegmentDrawer<OBJECT>(image, frame, do_draw, points_of, color_of, thickness_of);
      }


      template<typename REF_VERTEX>
      class VertexDrawer {

      private:
	
	using vertex_value_type = typename REF_VERTEX::element_type::value_type;
	
	mutable cv::Mat image; // a share pointer.
	Frame frame;
	std::function<bool (const vertex_value_type&)>          do_draw;
	std::function<demo2d::Point (const vertex_value_type&)> point_of;
	std::function<int (const vertex_value_type&)>           radius_of;
	std::function<cv::Scalar (const vertex_value_type&)>    color_of;
	std::function<int (const vertex_value_type&)>           thickness_of;
	
	
      public:
	
	template<typename DO_DRAW, typename POINT_OF, typename RADIUS_OF, typename COLOR_OF, typename THICKNESS_OF>
	VertexDrawer(cv::Mat& image,
		     Frame frame,
		     const DO_DRAW&      do_draw,
		     const POINT_OF&     point_of,
		     const RADIUS_OF&    radius_of,
		     const COLOR_OF&     color_of,
		     const THICKNESS_OF& thickness_of)
	  : image(image),
	    frame(frame),
	    do_draw(do_draw),
	    point_of(point_of),
	    radius_of(radius_of),
	    color_of(color_of),
	    thickness_of(thickness_of) {}

	VertexDrawer()                               = delete;
	VertexDrawer(const VertexDrawer&)            = default;
	VertexDrawer& operator=(const VertexDrawer&) = default;

	void operator()(REF_VERTEX v) const {
	  auto& o = (*v)();
	  if(do_draw(o))
	    cv::circle(image, frame(point_of(o)), radius_of(o), color_of(o), thickness_of(o));
	}
      };

      template<typename REF_VERTEX, typename DO_DRAW, typename POINT_OF, typename RADIUS_OF, typename COLOR_OF, typename THICKNESS_OF>
      auto vertex_drawer(cv::Mat& image,
			 Frame frame,
			 const DO_DRAW&      do_draw,
			 const POINT_OF&     point_of,
			 const RADIUS_OF&    radius_of,
			 const COLOR_OF&     color_of,
			 const THICKNESS_OF& thickness_of) {
	return VertexDrawer<REF_VERTEX>(image, frame, do_draw, point_of, radius_of, color_of, thickness_of);
      }




      template<typename REF_VERTEX>
      class VertexPrinter {

      private:
	
	using vertex_value_type = typename REF_VERTEX::element_type::value_type;
	
	mutable cv::Mat image; // a share pointer.
	Frame frame;
	std::function<bool (const vertex_value_type&)>                do_draw;
	std::function<std::string (const vertex_value_type&)>         text_of;
	std::function<demo2d::Point (const vertex_value_type&)>       point_of;
	std::function<cv::Scalar (const vertex_value_type&)>          color_of;
	std::function<int (const vertex_value_type&)>                 thickness_of;
	std::function<std::pair<int, int> (const vertex_value_type&)> offset_of;
	std::function<double (const vertex_value_type&)>              scale_of;
	
	
      public:
	
	template<typename DO_DRAW, typename TEXT_OF, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF, typename OFFSET_OF, typename SCALE_OF>
	VertexPrinter(cv::Mat& image,
		      Frame frame,
		      const DO_DRAW&       do_draw,
		      const TEXT_OF&       text_of,
		      const POINT_OF&      point_of,
		      const COLOR_OF&      color_of,
		      const THICKNESS_OF&  thickness_of,
		      const OFFSET_OF&     offset_of,
		      const SCALE_OF&      scale_of)
	  : image(image),
	    frame(frame),
	    do_draw(do_draw),
	    text_of(text_of),
	    point_of(point_of),
	    color_of(color_of),
	    thickness_of(thickness_of),
	    offset_of(offset_of),
	    scale_of(scale_of) {}

	VertexPrinter()                               = delete;
	VertexPrinter(const VertexPrinter&)            = default;
	VertexPrinter& operator=(const VertexPrinter&) = default;

	void operator()(REF_VERTEX v) const {
	  auto& o = (*v)();
	  if(do_draw(o)) {
	    auto pos    = frame(point_of(o));
	    auto offset = offset_of(o);
	    pos.x += offset.first;
	    pos.y += offset.second;
	    cv::putText(image, text_of(o), pos, CV_FONT_HERSHEY_DUPLEX, scale_of(o), color_of(o), thickness_of(o));
	  }
	}
      };

      template<typename REF_VERTEX, typename DO_DRAW, typename TEXT_OF, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF, typename OFFSET_OF, typename SCALE_OF>
      auto vertex_printer(cv::Mat& image,
			  Frame frame,
			  const DO_DRAW&      do_draw,
			  const TEXT_OF&      text_of,
			  const POINT_OF&     point_of,
			  const COLOR_OF&     color_of,
			  const THICKNESS_OF& thickness_of,
			  const OFFSET_OF&    offset_of,
			  const SCALE_OF&     scale_of) {
	return VertexPrinter<REF_VERTEX>(image, frame, do_draw, text_of, point_of, color_of, thickness_of, offset_of, scale_of);
      }

      template<typename REF_EDGE, typename EDGE_VALUE>
      class EdgeDrawer {

      private:
	
      	using edge_value_type   = EDGE_VALUE;
      	using vertex_value_type = typename REF_EDGE::element_type::vertex_val_type;
	
      	mutable cv::Mat image; // a share pointer.
      	Frame frame;
	std::function<bool (const vertex_value_type&,
			    const vertex_value_type&,
			    const edge_value_type&)>            do_draw;
      	std::function<demo2d::Point (const vertex_value_type&)> point_of;
      	std::function<cv::Scalar (const edge_value_type&)>      color_of;
      	std::function<int (const edge_value_type&)>             thickness_of;
	
	
      public:
	
      	template<typename DO_DRAW, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF>
      	EdgeDrawer(cv::Mat& image,
      		   Frame frame,
		   const DO_DRAW&      do_draw,
      		   const POINT_OF&     point_of,
      		   const COLOR_OF&     color_of,
      		   const THICKNESS_OF& thickness_of)
      	  : image(image),
      	    frame(frame),
	    do_draw(do_draw),
      	    point_of(point_of),
      	    color_of(color_of),
      	    thickness_of(thickness_of) {}

      	EdgeDrawer()                            = delete;
      	EdgeDrawer(const EdgeDrawer&)            = default;
      	EdgeDrawer& operator=(const EdgeDrawer&) = default;

      	void operator()(REF_EDGE e) const {
	  auto extr = e->extremities();
	  if(vq3::invalid_extremities(extr)) {
	    e->kill();
	    return;
	  }
      	  auto& o  = (*e)();
	  auto& v1 = (*extr.first)();
	  auto& v2 = (*extr.second)();
	  if(do_draw(v1, v2, o))
	    cv::line(image,
		     frame(point_of(v1)),
		     frame(point_of(v2)),
		     color_of(o),
		     thickness_of(o));
	}
      };

      template<typename REF_EDGE>
      class EdgeDrawer<REF_EDGE, void> {

      private:
	
      	using edge_value_type   = void;
      	using vertex_value_type = typename REF_EDGE::element_type::vertex_val_type;
	
      	mutable cv::Mat image; // a share pointer.
      	Frame frame;
	std::function<bool (const vertex_value_type&,
			    const vertex_value_type&)>          do_draw;
      	std::function<demo2d::Point (const vertex_value_type&)> point_of;
      	std::function<cv::Scalar ()>                            color_of;
      	std::function<int ()>                                   thickness_of;
	
	
      public:
	
      	template<typename DO_DRAW, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF>
      	EdgeDrawer(cv::Mat& image,
      		   Frame frame,
		   const DO_DRAW&      do_draw,
      		   const POINT_OF&     point_of,
      		   const COLOR_OF&     color_of,
      		   const THICKNESS_OF& thickness_of)
      	  : image(image),
      	    frame(frame),
	    do_draw(do_draw),
      	    point_of(point_of),
      	    color_of(color_of),
      	    thickness_of(thickness_of) {}

      	EdgeDrawer()                            = delete;
      	EdgeDrawer(const EdgeDrawer&)            = default;
      	EdgeDrawer& operator=(const EdgeDrawer&) = default;

      	void operator()(REF_EDGE e) const {
	  auto extr = e->extremities();
	  if(vq3::invalid_extremities(extr)) {
	    e->kill();
	    return;
	  }
	  auto& v1 = (*extr.first)();
	  auto& v2 = (*extr.second)();
	  if(do_draw(v1, v2))
	    cv::line(image,
		     frame(point_of(v1)),
		     frame(point_of(v2)),
		     color_of(),
		     thickness_of());
      	}
      };
      
      template<typename REF_EDGE, typename DO_DRAW, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF>
      auto edge_drawer(cv::Mat& image,
		       Frame frame,
		       const DO_DRAW&      do_draw,
		       const POINT_OF&     point_of,
		       const COLOR_OF&     color_of,
		       const THICKNESS_OF& thickness_of) {
      	return EdgeDrawer<REF_EDGE, typename REF_EDGE::element_type::value_type>(image, frame, do_draw, point_of, color_of, thickness_of);
      }


      template<typename REF_EDGE, typename EDGE_VALUE>
      class EdgePrinter {

      private:
	
      	using edge_value_type   = EDGE_VALUE;
      	using vertex_value_type = typename REF_EDGE::element_type::vertex_val_type;
	
      	mutable cv::Mat image; // a share pointer.
      	Frame frame;
	std::function<bool (const vertex_value_type&,
			    const vertex_value_type&,
			    const edge_value_type&)>                do_draw;
	std::function<std::string (const edge_value_type&)>         text_of;
      	std::function<demo2d::Point (const vertex_value_type&)>     point_of;
      	std::function<cv::Scalar (const edge_value_type&)>          color_of;
      	std::function<int (const edge_value_type&)>                 thickness_of;
	std::function<std::pair<int, int> (const edge_value_type&)> offset_of;
	std::function<double (const edge_value_type&)>              scale_of;
	
	
      public:
	
      	template<typename DO_DRAW, typename TEXT_OF, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF, typename OFFSET_OF, typename SCALE_OF>
      	EdgePrinter(cv::Mat& image,
		    Frame frame,
		    const DO_DRAW&       do_draw,
		    const TEXT_OF&       text_of,
		    const POINT_OF&      point_of,
		    const COLOR_OF&      color_of,
		    const THICKNESS_OF&  thickness_of,
		    const OFFSET_OF&     offset_of,
		    const SCALE_OF&      scale_of)
      	  : image(image),
      	    frame(frame),
	    do_draw(do_draw),
	    text_of(text_of),
      	    point_of(point_of),
      	    color_of(color_of),
      	    thickness_of(thickness_of),
	    offset_of(offset_of),
	    scale_of(scale_of) {}

      	EdgePrinter()                            = delete;
      	EdgePrinter(const EdgePrinter&)            = default;
      	EdgePrinter& operator=(const EdgePrinter&) = default;

      	void operator()(REF_EDGE e) const {
	  auto extr = e->extremities();
	  if(vq3::invalid_extremities(extr)) {
	    e->kill();
	    return;
	  }
      	  auto& o  = (*e)();
	  auto& v1 = (*extr.first)();
	  auto& v2 = (*extr.second)();
	  if(do_draw(v1, v2, o)) {
	    auto pos    = frame((point_of(v1) + point_of(v2))*.5);
	    auto offset = offset_of(o);
	    pos.x += offset.first;
	    pos.y += offset.second;
	    cv::putText(image, text_of(o), pos, CV_FONT_HERSHEY_DUPLEX, scale_of(o), color_of(o), thickness_of(o));
	  }
	}
      };

      template<typename REF_EDGE>
      class EdgePrinter<REF_EDGE, void> {

      private:
	
      	using edge_value_type   = void;
      	using vertex_value_type = typename REF_EDGE::element_type::vertex_val_type;
	
      	mutable cv::Mat image; // a share pointer.
      	Frame frame;
	std::function<bool (const vertex_value_type&,
			    const vertex_value_type&)>          do_draw;
	std::function<std::string ()>                           text_of;
      	std::function<demo2d::Point (const vertex_value_type&)> point_of;
      	std::function<cv::Scalar ()>                            color_of;
      	std::function<int ()>                                   thickness_of;
	std::function<std::pair<int, int> ()>                   offset_of;
	std::function<double ()>                                scale_of;
	
	
      public:
	
      	template<typename DO_DRAW, typename TEXT_OF, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF, typename OFFSET_OF, typename SCALE_OF>
      	EdgePrinter(cv::Mat& image,
		    Frame frame,
		    const DO_DRAW&       do_draw,
		    const TEXT_OF&       text_of,
		    const POINT_OF&      point_of,
		    const COLOR_OF&      color_of,
		    const THICKNESS_OF&  thickness_of,
		    const OFFSET_OF&     offset_of,
		    const SCALE_OF&      scale_of)
      	  : image(image),
      	    frame(frame),
	    do_draw(do_draw),
	    text_of(text_of),
      	    point_of(point_of),
      	    color_of(color_of),
      	    thickness_of(thickness_of),
	    offset_of(offset_of),
	    scale_of(scale_of) {}

      	EdgePrinter()                            = delete;
      	EdgePrinter(const EdgePrinter&)            = default;
      	EdgePrinter& operator=(const EdgePrinter&) = default;

      	void operator()(REF_EDGE e) const {
	  auto extr = e->extremities();
	  if(vq3::invalid_extremities(extr)) {
	    e->kill();
	    return;
	  }
	  auto& v1 = (*extr.first)();
	  auto& v2 = (*extr.second)();
	  if(do_draw(v1, v2)) {
	    auto pos    = frame((point_of(v1) + point_of(v2))*.5);
	    auto offset = offset_of();
	    pos.x += offset.first;
	    pos.y += offset.second;
	    cv::putText(image, text_of(), pos, CV_FONT_HERSHEY_DUPLEX, scale_of(), color_of(), thickness_of());
	  }
      	}
      };
      
      template<typename REF_EDGE, typename DO_DRAW, typename TEXT_OF, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF, typename OFFSET_OF, typename SCALE_OF>
      auto edge_printer(cv::Mat& image,
			Frame frame,
			const DO_DRAW&       do_draw,
			const TEXT_OF&       text_of,
			const POINT_OF&      point_of,
			const COLOR_OF&      color_of,
			const THICKNESS_OF&  thickness_of,
			const OFFSET_OF&     offset_of,
			const SCALE_OF&      scale_of) {
      	return EdgePrinter<REF_EDGE, typename REF_EDGE::element_type::value_type>(image, frame, do_draw, text_of, point_of, color_of, thickness_of, offset_of, scale_of);
      }

      template<typename REF_VERTEX>
      class SegmentAtVertexDrawer {

      private:
	
	using vertex_value_type = typename REF_VERTEX::element_type::value_type;
	
	mutable cv::Mat image; // a share pointer.
	Frame frame;
	std::function<bool (const vertex_value_type&)>          do_draw;
	std::function<demo2d::Point (const vertex_value_type&)> point_of;
	std::function<demo2d::Point (const vertex_value_type&)> segment_of;
	std::function<cv::Scalar (const vertex_value_type&)>    color_of;
	std::function<int (const vertex_value_type&)>           thickness_of;
	
	
      public:
	
	template<typename DO_DRAW, typename POINT_OF, typename SEGMENT_OF, typename COLOR_OF, typename THICKNESS_OF>
	SegmentAtVertexDrawer(cv::Mat& image,
			      Frame frame,
			      const DO_DRAW&      do_draw,
			      const POINT_OF&     point_of,
			      const SEGMENT_OF&   segment_of,
			      const COLOR_OF&     color_of,
			      const THICKNESS_OF& thickness_of)
	  : image(image),
	    frame(frame),
	    do_draw(do_draw),
	    point_of(point_of),
	    segment_of(segment_of),
	    color_of(color_of),
	    thickness_of(thickness_of) {}

	SegmentAtVertexDrawer()                                        = delete;
	SegmentAtVertexDrawer(const SegmentAtVertexDrawer&)            = default;
	SegmentAtVertexDrawer& operator=(const SegmentAtVertexDrawer&) = default;

	void operator()(REF_VERTEX v) const {
	  auto& o = (*v)();
	  if(do_draw(o)) {
	    auto A = point_of(o);
	    auto B = A + segment_of(o);
	    cv::line(image, frame(A), frame(B), color_of(o), thickness_of(o));
	  }
	}
      };

      template<typename REF_VERTEX, typename DO_DRAW, typename POINT_OF, typename SEGMENT_OF, typename COLOR_OF, typename THICKNESS_OF>
      auto segment_at_vertex_drawer(cv::Mat& image,
				    Frame frame,
				    const DO_DRAW&      do_draw,
				    const POINT_OF&     point_of,
				    const SEGMENT_OF&   segment_of,
				    const COLOR_OF&     color_of,
				    const THICKNESS_OF& thickness_of) {
	return SegmentAtVertexDrawer<REF_VERTEX>(image, frame, do_draw, point_of, segment_of, color_of, thickness_of);
      }
      
      class histogram : public vq3::stats::histogram {
      public:

	demo2d::Point min, max;
	std::optional<std::string> title;
	std::optional<std::pair<double,double>> range;
	std::optional<std::pair<double,double>> value_bounds;
	std::optional<double> NT;
	std::optional<double> max_histo_display;
	double frame_margin = 0;
	std::optional<std::pair<double, double>> sci_range;
	cv::Scalar frame_foreground = {  0,   0,   0};
	cv::Scalar frame_background = {255, 255, 255};
	cv::Scalar frame_bar        = {255,  85,  85};
	cv::Scalar frame_sci_bar    = {170,  85,  85};
	cv::Scalar NT_color         = {  0,   0, 255};
	cv::Scalar range_color      = {  0,   0,   0};
	int        NT_thickness     = 3;
	int        frame_thickness  = 3;
	int        axis_thickness   = 1;
	int        range_thickness  = 3;
	
	
	histogram(const demo2d::Point& min, const demo2d::Point& max)
	  : vq3::stats::histogram(),
	  min(min), max(max){
	}
	
      private:
	double y_coef;
	demo2d::Point pA, pB; // plot drawing area.
	
	double value_min, value_max;
	double value_coef;

	double bin_coef;

	double maxh;

	bool draw_sci_histo = false;

	double abscissa_of(double value) {
	  return (value - value_min)*value_coef + pA.x;
	}
	
	double value_of(double bin) {
	  return bin_min + bin*bin_coef;
	}
	
	void draw_bar(cv::Mat& image, const Frame& frame, unsigned int bin) {

	  double value = value_of(bin+.5);
	  double a     = abscissa_of(value_of(bin+ 0));
	  double b     = abscissa_of(value_of(bin+ 1));
	  cv::Scalar color = frame_bar;
	  if(sci_range && sci_range.value().first <= value && value < sci_range.value().second)
	    color = frame_sci_bar;

	  auto AA = demo2d::Point(a, pA.y);
	  auto BB = demo2d::Point(b, pA.y + maxh*std::min(1.0, h[bin]*y_coef));
	  cv::rectangle(image, frame(AA), frame(BB), color, -1);
	}

	bool in_drawing_area(const demo2d::Point& pt) {
	  return pA.x <= pt.x && pt.x <= pB.x
	    &&   pA.y <= pt.y && pt.y <= pB.y;
	}
						     
	
      public:

	/**
	 * Toggles the drawing of smallest confidence interval (if sci computation is set).
	 */
	void draw_sci(bool on) {
	  draw_sci_histo = on && sci_conf;
	}

	/**
	 * This draw the histogram. Call this before any other supplementary drawings.
	 */
	void draw(cv::Mat& image, const Frame& frame) {

	  if(draw_sci_histo)
	    sci_range = sci;
	  else
	    sci_range.reset();
	  
	  cv::rectangle(image, frame(min), frame(max), frame_background ,              -1);
	  cv::rectangle(image, frame(min), frame(max), frame_foreground , frame_thickness);
	  if(title) 
	    cv::putText(image, *title, frame(min) + cv::Point(0,20), cv::FONT_HERSHEY_PLAIN, 1., cv::Scalar(0,0,0), 1);

	  pA = min + demo2d::Point( frame_margin,  frame_margin);
	  pB = max + demo2d::Point(-frame_margin, -frame_margin);
	  maxh  = pB.y - pA.y;
	  
	  if(max_histo_display)
	    y_coef = 1/max_histo_display.value();
	  else {
	    unsigned int max = 0;
	    if(h.begin() != h.end())
	      max = *(std::max_element(h.begin(), h.end()));
	    if(max == 0)
	      max = 1;
	    y_coef = 1 / (double)max;
	  }

	  bin_coef = (bin_max - bin_min)/bin_nb;
	  

	  if(value_bounds) 
	    std::tie(value_min, value_max) = value_bounds.value();
	  else if(NT) {
	    double vmin = value_of(0 - .5);
	    double vmax = value_of(bin_nb + .5);

	    if(value_min > NT.value()) {
	      value_max = vmax;
	      value_min = 2*NT.value() - value_max;
	    }
	    else if(vmax < NT.value()) {
	      value_min = vmin;
	      value_max = 2*NT.value() - value_min;
	    }
	    else {
	      auto r = std::max(vmax - NT.value(), NT.value() - vmin);
	      value_min = NT.value() - r;
	      value_max = NT.value() + r;
	    }
	  }
	  else {
	    value_min = value_of(0 - .5);
	    value_max = value_of(bin_nb + .5);
	  }
	  value_coef = (pB.x - pA.x)/(value_max - value_min);

	  for(unsigned int bin = 0; bin < bin_nb; ++bin)
	    draw_bar(image, frame, bin);

	  if(NT) {
	    double v = abscissa_of(NT.value());
	    auto AA = demo2d::Point(v, pA.y);
	    auto BB = demo2d::Point(v, max.y - frame_margin);
	    cv::line(image, frame(AA), frame(BB), NT_color, NT_thickness);
	  }

	  if(range) {
	    
	    auto [rmin, rmax] = range.value();
	    auto rA = demo2d::Point(abscissa_of(rmin), pA.y);
	    auto rB = demo2d::Point(abscissa_of(rmax), pA.y);

	    if(rA.x < pA.x)
	      rA.x = pA.x;
	    else if(rA.x > pB.x)
	      rA.x = pB.x;

	    if(rB.x < pA.x)
	      rB.x = pA.x;
	    else if(rB.x > pB.x)
	      rB.x = pB.x;
	    
	    cv::line(image, frame(rA), frame(rB), range_color, range_thickness);
	  }
	  
	  cv::line(image, frame(pA), frame(demo2d::Point(pB.x,pA.y)), frame_foreground, axis_thickness);
	}

	void interval(cv::Mat& image, const Frame& frame, double vmin, double vmax, int pix_offset, const cv::Scalar& color, int thickness) {
	  auto rA = demo2d::Point(abscissa_of(vmin), pA.y);
	  auto rB = demo2d::Point(abscissa_of(vmax), pA.y);

	  if(rA.x < pA.x)
	    rA.x = pA.x;
	  else if(rA.x > pB.x)
	    rA.x = pB.x;

	  if(rB.x < pA.x)
	    rB.x = pA.x;
	  else if(rB.x > pB.x)
	    rB.x = pB.x;

	  auto A = frame(rA);
	  auto B = frame(rB);
	  A.y += pix_offset;
	  B.y += pix_offset;
	  cv::line(image, A, B, range_color, range_thickness);
	}
	
	/**
	 * This draws a vertical line at abscissa value. Be sure that draw has been called before.
	 */
	void vline(cv::Mat& image, const Frame& frame, double value, const cv::Scalar& color, int thickness) {
	  double x = abscissa_of(value);
	  if(pA.x <= x && x <= pB.x) {
	    auto lA = demo2d::Point(abscissa_of(value), pA.y);
	    auto lB = demo2d::Point(abscissa_of(value), pB.y);
	    cv::line(image, frame(lA), frame(lB), color, thickness);
	  }
	}

	/**
	 * This draws a curve over the histogram. Be sure that draw has been called before.
	 * @param begin, end The range of dato to be plotted.
	 * @param point_of The call point_of(*it) must produce a demo2d::Point.
	 * @param clip true means that the curve is clipped in the drawing area.
	 */
	template<typename Iter, typename PointOf>
	void curve(cv::Mat& image, const Frame& frame, const Iter& begin, const Iter& end, const PointOf& point_of, const cv::Scalar& color, int thickness, bool clip) {

	  if(begin == end) return;
	  auto prev = begin;
	  auto curr = begin;
	  std::advance(curr,1);
	  if(curr == end)  return;
	  while(curr != end) {
	    auto A = point_of(*prev);
	    auto B = point_of(*curr);
	    A.x = abscissa_of(A.x);
	    A.y = pA.y + maxh*A.y*y_coef;
	    B.x = abscissa_of(B.x);
	    B.y = pA.y + maxh*B.y*y_coef;
	    prev = (curr++);
	    if((!clip) || (in_drawing_area(A) && in_drawing_area(B)))
	      cv::line(image, frame(A), frame(B), color, thickness);
	  }
	}

	/**
	 * This draws a gaussian function, using nb_points points.
	 */
	void gaussian_stddev(cv::Mat& image, const Frame& frame,
			     double mean, double std_dev, double amplitude, unsigned nb_points,
			     const cv::Scalar& color, int thickness, bool clip) {
	  gaussian_var(image,frame, mean, std_dev*std_dev, amplitude, nb_points, color, thickness, clip);
	}
	
	/**
	 * This draws a gaussian function, using nb_points points.
	 */
	void gaussian_var(cv::Mat& image, const Frame& frame,
			  double mean, double var, double amplitude, unsigned nb_points,
			  const cv::Scalar& color, int thickness, bool clip) {
	  std::vector<demo2d::Point> points;
	  auto out = std::back_inserter(points);
	  double coef  = (value_max - value_min)/nb_points;
	  double gamma = -.5/var;
	  for(unsigned int i=0; i <= nb_points; ++i) {
	    double x   = value_min + i*coef;
	    double dx = x - mean;
	    *(out++) = {x, amplitude*std::exp(gamma*dx*dx)};
	  }

	  curve(image, frame, points.begin(), points.end(), [](const demo2d::Point& p) {return p;},
		color, thickness, clip);
	}
      };


      namespace sample {

	/**
	 * This handles image-based distributions information. This is
	 * passed to an actual image based distribution for storing all the
	 * required information.
	 */
	struct ImageData {
	  cv::Mat image; //!< a share pointer to the handled image.
	  Frame frame;   //!< The frame for converting pixel positions into math coordinates.
	  std::function<double (const unsigned char*)> pixel_to_density; //!< This converts pixel values (RGB) into density (in [0,1]).
	  ImageData()                            = default;
	  ImageData(const ImageData&)            = default;
	  ImageData& operator=(const ImageData&) = default;
	  ImageData(ImageData&&)                 = default;
	  ImageData& operator=(ImageData&&)      = default;
	  
	  template<typename PIXEL_TO_DENSITY>
	  ImageData(const PIXEL_TO_DENSITY& pixel_to_density) : image(), frame(), pixel_to_density(pixel_to_density) {}
	};
	
	/**
	 * Makes a data that handles image-based distributions
	 * information. This is passed to an actual image distribution
	 * for storing all the required information.
	 */
	template<typename PIXEL_TO_DENSITY>
	ImageData image_data(const PIXEL_TO_DENSITY& pixel_to_density) {
	  return ImageData(pixel_to_density);
	}
	
	/**
	 * This handles video-based distributions information. This is
	 * passed to an actual video-based distribution for storing all the
	 * required information.
	 */

	struct VideoData : ImageData {
	  cv::VideoCapture cap;
	  
	  VideoData()                            = default;
	  VideoData(const VideoData&)            = default;
	  VideoData& operator=(const VideoData&) = default;
	  VideoData(VideoData&&)                 = default;
	  VideoData& operator=(VideoData&&)      = default;
	  
	  template<typename PIXEL_TO_DENSITY>
	  VideoData(int device, const PIXEL_TO_DENSITY& pixel_to_density) : ImageData(pixel_to_density), cap(device) {
	    if(!cap.isOpened()) {
	      std::ostringstream ostr;
	      ostr << "Cannot open video from device " << device << " (/dev/video" << device << ").";
	      throw std::runtime_error(ostr.str());
	    }
	    image = cv::Mat((int)(cap.get(CV_CAP_PROP_FRAME_HEIGHT)),
			    (int)(cap.get(CV_CAP_PROP_FRAME_WIDTH)),
			    CV_8UC3, cv::Scalar(255,255,255));
	  }

	  /**
	   * This grabs the next frame.
	   */
	  void operator++() {
	    cap >> image;
	  }
	};
	
	/**
	 * Makes a data that handles video-based distributions
	 * information. This is passed to an actual video-based
	 * distribution for storing all the required information.
	 */
	template<typename PIXEL_TO_DENSITY>
	VideoData video_data(int device, const PIXEL_TO_DENSITY& pixel_to_density) {
	  return VideoData(device, pixel_to_density);
	}
	

	class Webcam;

	/**
	 * This is a density related to an image (see
	 * vq3::demo2d::opencv::sample::image).
	 */
	class Image : public vq3::demo2d::sample::Density {
	private:

	  ImageData& image_data;
	  
	  Image()                        = delete;
	  Image(const Image&)            = delete;
	  Image& operator=(const Image&) = delete;
	  Image(ImageData& image_data) : image_data(image_data) {}

	  friend vq3::demo2d::sample::density image(ImageData& image_data);
	  friend class Webcam;
	  
	public:
	  
	  virtual vq3::demo2d::sample::BBox bbox() const override {
	    return vq3::demo2d::sample::BBox(image_data.frame(cv::Point(0, image_data.image.size().height - 1)),
					     image_data.frame(cv::Point(image_data.image.size().width - 1, 0)));
	  }
	  
	  virtual double operator()(const Point& p) const override {
	    if(bbox().contains(p)) {
	      auto cv_pt = image_data.frame(p);
	      return image_data.pixel_to_density(reinterpret_cast<const unsigned char*>(&(image_data.image.at<unsigned char[3]>(cv_pt.y, cv_pt.x))));
	    }
	    else
	      return 0;
	  }
	  
	};

	/**
	 * This builds an image-based distribution.
	 * @param image_data The image data for storing image-related information.
	 */
	inline vq3::demo2d::sample::density image(ImageData& image_data) {
	  return vq3::demo2d::sample::density(new Image(image_data));
	}


	/**
	 * This is a density related to images taken from a webcam (see
	 * vq3::demo2d::opencv::sample::webcam).
	 */

	class Webcam : public vq3::demo2d::sample::Density {
	  
	private:
	  
	  VideoData& video_data;
	  Image image;
	  
	  Webcam()                         = delete;
	  Webcam(const Webcam&)            = delete;
	  Webcam& operator=(const Webcam&) = delete;

	  Webcam(VideoData& video_data) : vq3::demo2d::sample::Density(), video_data(video_data), image(video_data) {}

	  friend vq3::demo2d::sample::density webcam(VideoData& video_data);
	
	public:
	  
	  virtual vq3::demo2d::sample::BBox bbox() const override {
	    return image.bbox();
	  }
	  virtual double operator()(const Point& p) const override {
	    return image(p);
	  }

	};
	
	/**
	 * This builds webcam-based distribution.
	 * @param video_data The video data for storing webcam-related information.
	 */
	inline vq3::demo2d::sample::density webcam(VideoData& video_data) {
	  return vq3::demo2d::sample::density(new Webcam(video_data));
	}
      }

    }
  }
}
