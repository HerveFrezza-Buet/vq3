
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
#include <demo2d.hpp>
#include <algorithm>


namespace vq3 {
  
  namespace demo2d {
      
    namespace opencv {
      
      template<typename REF_VERTEX>
      class VertexDrawer {

      private:
	
	using vertex_value_type = typename REF_VERTEX::element_type::value_type;
	
	mutable cv::Mat image; // a share pointer.
	::demo2d::opencv::Frame frame;
	std::function<bool (const vertex_value_type&)>          do_draw;
	std::function<::demo2d::Point (const vertex_value_type&)> point_of;
	std::function<int (const vertex_value_type&)>           radius_of;
	std::function<cv::Scalar (const vertex_value_type&)>    color_of;
	std::function<int (const vertex_value_type&)>           thickness_of;
	
	
      public:
	
	template<typename DO_DRAW, typename POINT_OF, typename RADIUS_OF, typename COLOR_OF, typename THICKNESS_OF>
	VertexDrawer(cv::Mat& image,
		     ::demo2d::opencv::Frame frame,
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
			 ::demo2d::opencv::Frame frame,
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
	::demo2d::opencv::Frame frame;
	std::function<bool (const vertex_value_type&)>                do_draw;
	std::function<std::string (const vertex_value_type&)>         text_of;
	std::function<::demo2d::Point (const vertex_value_type&)>       point_of;
	std::function<cv::Scalar (const vertex_value_type&)>          color_of;
	std::function<int (const vertex_value_type&)>                 thickness_of;
	std::function<std::pair<int, int> (const vertex_value_type&)> offset_of;
	std::function<double (const vertex_value_type&)>              scale_of;
	
	
      public:
	
	template<typename DO_DRAW, typename TEXT_OF, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF, typename OFFSET_OF, typename SCALE_OF>
	VertexPrinter(cv::Mat& image,
		      ::demo2d::opencv::Frame frame,
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
	    cv::putText(image, text_of(o), pos, cv::FONT_HERSHEY_DUPLEX, scale_of(o), color_of(o), thickness_of(o));
	  }
	}
      };

      template<typename REF_VERTEX, typename DO_DRAW, typename TEXT_OF, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF, typename OFFSET_OF, typename SCALE_OF>
      auto vertex_printer(cv::Mat& image,
			  ::demo2d::opencv::Frame frame,
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
      	::demo2d::opencv::Frame frame;
	std::function<bool (const vertex_value_type&,
			    const vertex_value_type&,
			    const edge_value_type&)>            do_draw;
      	std::function<::demo2d::Point (const vertex_value_type&)> point_of;
      	std::function<cv::Scalar (const edge_value_type&)>      color_of;
      	std::function<int (const edge_value_type&)>             thickness_of;
	
	
      public:
	
      	template<typename DO_DRAW, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF>
      	EdgeDrawer(cv::Mat& image,
      		   ::demo2d::opencv::Frame frame,
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
      	::demo2d::opencv::Frame frame;
	std::function<bool (const vertex_value_type&,
			    const vertex_value_type&)>          do_draw;
      	std::function<::demo2d::Point (const vertex_value_type&)> point_of;
      	std::function<cv::Scalar ()>                            color_of;
      	std::function<int ()>                                   thickness_of;
	
	
      public:
	
      	template<typename DO_DRAW, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF>
      	EdgeDrawer(cv::Mat& image,
      		   ::demo2d::opencv::Frame frame,
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
		       ::demo2d::opencv::Frame frame,
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
      	::demo2d::opencv::Frame frame;
	std::function<bool (const vertex_value_type&,
			    const vertex_value_type&,
			    const edge_value_type&)>                do_draw;
	std::function<std::string (const edge_value_type&)>         text_of;
      	std::function<::demo2d::Point (const vertex_value_type&)>     point_of;
      	std::function<cv::Scalar (const edge_value_type&)>          color_of;
      	std::function<int (const edge_value_type&)>                 thickness_of;
	std::function<std::pair<int, int> (const edge_value_type&)> offset_of;
	std::function<double (const edge_value_type&)>              scale_of;
	
	
      public:
	
      	template<typename DO_DRAW, typename TEXT_OF, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF, typename OFFSET_OF, typename SCALE_OF>
      	EdgePrinter(cv::Mat& image,
		    ::demo2d::opencv::Frame frame,
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
	    cv::putText(image, text_of(o), pos, cv::FONT_HERSHEY_DUPLEX, scale_of(o), color_of(o), thickness_of(o));
	  }
	}
      };

      template<typename REF_EDGE>
      class EdgePrinter<REF_EDGE, void> {

      private:
	
      	using edge_value_type   = void;
      	using vertex_value_type = typename REF_EDGE::element_type::vertex_val_type;
	
      	mutable cv::Mat image; // a share pointer.
      	::demo2d::opencv::Frame frame;
	std::function<bool (const vertex_value_type&,
			    const vertex_value_type&)>          do_draw;
	std::function<std::string ()>                           text_of;
      	std::function<::demo2d::Point (const vertex_value_type&)> point_of;
      	std::function<cv::Scalar ()>                            color_of;
      	std::function<int ()>                                   thickness_of;
	std::function<std::pair<int, int> ()>                   offset_of;
	std::function<double ()>                                scale_of;
	
	
      public:
	
      	template<typename DO_DRAW, typename TEXT_OF, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF, typename OFFSET_OF, typename SCALE_OF>
      	EdgePrinter(cv::Mat& image,
		    ::demo2d::opencv::Frame frame,
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
	    cv::putText(image, text_of(), pos, cv::FONT_HERSHEY_DUPLEX, scale_of(), color_of(), thickness_of());
	  }
      	}
      };
      
      template<typename REF_EDGE, typename DO_DRAW, typename TEXT_OF, typename POINT_OF, typename COLOR_OF, typename THICKNESS_OF, typename OFFSET_OF, typename SCALE_OF>
      auto edge_printer(cv::Mat& image,
			::demo2d::opencv::Frame frame,
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
	::demo2d::opencv::Frame frame;
	std::function<bool (const vertex_value_type&)>          do_draw;
	std::function<::demo2d::Point (const vertex_value_type&)> point_of;
	std::function<::demo2d::Point (const vertex_value_type&)> segment_of;
	std::function<cv::Scalar (const vertex_value_type&)>    color_of;
	std::function<int (const vertex_value_type&)>           thickness_of;
	
	
      public:
	
	template<typename DO_DRAW, typename POINT_OF, typename SEGMENT_OF, typename COLOR_OF, typename THICKNESS_OF>
	SegmentAtVertexDrawer(cv::Mat& image,
			      ::demo2d::opencv::Frame frame,
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
				    ::demo2d::opencv::Frame frame,
				    const DO_DRAW&      do_draw,
				    const POINT_OF&     point_of,
				    const SEGMENT_OF&   segment_of,
				    const COLOR_OF&     color_of,
				    const THICKNESS_OF& thickness_of) {
	return SegmentAtVertexDrawer<REF_VERTEX>(image, frame, do_draw, point_of, segment_of, color_of, thickness_of);
      }
      
      class histogram : public stats::histogram {
      public:

	::demo2d::Point min, max;
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
	
	
	histogram(const ::demo2d::Point& min, const ::demo2d::Point& max)
	  : stats::histogram(),
	  min(min), max(max){
	}
	
      private:
	double y_coef;
	::demo2d::Point pA, pB; // plot drawing area.
	
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
	
	void draw_bar(cv::Mat& image, const ::demo2d::opencv::Frame& frame, unsigned int bin) {

	  double value = value_of(bin+.5);
	  double a     = abscissa_of(value_of(bin+ 0));
	  double b     = abscissa_of(value_of(bin+ 1));
	  cv::Scalar color = frame_bar;
	  if(sci_range && sci_range.value().first <= value && value < sci_range.value().second)
	    color = frame_sci_bar;

	  auto AA = ::demo2d::Point(a, pA.y);
	  auto BB = ::demo2d::Point(b, pA.y + maxh*std::min(1.0, h[bin]*y_coef));
	  cv::rectangle(image, frame(AA), frame(BB), color, -1);
	}

	bool in_drawing_area(const ::demo2d::Point& pt) {
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
	void draw(cv::Mat& image, const ::demo2d::opencv::Frame& frame) {

	  if(draw_sci_histo)
	    sci_range = sci;
	  else
	    sci_range.reset();
	  
	  cv::rectangle(image, frame(min), frame(max), frame_background ,              -1);
	  cv::rectangle(image, frame(min), frame(max), frame_foreground , frame_thickness);
	  if(title) 
	    cv::putText(image, *title, frame(min) + cv::Point(0,20), cv::FONT_HERSHEY_PLAIN, 1., cv::Scalar(0,0,0), 1);

	  pA = min + ::demo2d::Point( frame_margin,  frame_margin);
	  pB = max + ::demo2d::Point(-frame_margin, -frame_margin);
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
	    auto AA = ::demo2d::Point(v, pA.y);
	    auto BB = ::demo2d::Point(v, max.y - frame_margin);
	    cv::line(image, frame(AA), frame(BB), NT_color, NT_thickness);
	  }

	  if(range) {
	    
	    auto [rmin, rmax] = range.value();
	    auto rA = ::demo2d::Point(abscissa_of(rmin), pA.y);
	    auto rB = ::demo2d::Point(abscissa_of(rmax), pA.y);

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
	  
	  cv::line(image, frame(pA), frame(::demo2d::Point(pB.x,pA.y)), frame_foreground, axis_thickness);
	}

	void interval(cv::Mat& image, const ::demo2d::opencv::Frame& frame, double vmin, double vmax, int pix_offset) {
	  auto rA = ::demo2d::Point(abscissa_of(vmin), pA.y);
	  auto rB = ::demo2d::Point(abscissa_of(vmax), pA.y);

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
	void vline(cv::Mat& image, const ::demo2d::opencv::Frame& frame, double value, const cv::Scalar& color, int thickness) {
	  double x = abscissa_of(value);
	  if(pA.x <= x && x <= pB.x) {
	    auto lA = ::demo2d::Point(abscissa_of(value), pA.y);
	    auto lB = ::demo2d::Point(abscissa_of(value), pB.y);
	    cv::line(image, frame(lA), frame(lB), color, thickness);
	  }
	}

	/**
	 * This draws a curve over the histogram. Be sure that draw has been called before.
	 * @param begin, end The range of dato to be plotted.
	 * @param point_of The call point_of(*it) must produce a ::demo2d::Point.
	 * @param clip true means that the curve is clipped in the drawing area.
	 */
	template<typename Iter, typename PointOf>
	void curve(cv::Mat& image, const ::demo2d::opencv::Frame& frame, const Iter& begin, const Iter& end, const PointOf& point_of, const cv::Scalar& color, int thickness, bool clip) {

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
	void gaussian_stddev(cv::Mat& image, const ::demo2d::opencv::Frame& frame,
			     double mean, double std_dev, double amplitude, unsigned nb_points,
			     const cv::Scalar& color, int thickness, bool clip) {
	  gaussian_var(image,frame, mean, std_dev*std_dev, amplitude, nb_points, color, thickness, clip);
	}
	
	/**
	 * This draws a gaussian function, using nb_points points.
	 */
	void gaussian_var(cv::Mat& image, const ::demo2d::opencv::Frame& frame,
			  double mean, double var, double amplitude, unsigned nb_points,
			  const cv::Scalar& color, int thickness, bool clip) {
	  std::vector<::demo2d::Point> points;
	  auto out = std::back_inserter(points);
	  double coef  = (value_max - value_min)/nb_points;
	  double gamma = -.5/var;
	  for(unsigned int i=0; i <= nb_points; ++i) {
	    double x   = value_min + i*coef;
	    double dx = x - mean;
	    *(out++) = {x, amplitude*std::exp(gamma*dx*dx)};
	  }

	  curve(image, frame, points.begin(), points.end(), [](const ::demo2d::Point& p) {return p;},
		color, thickness, clip);
	}
      };

      
    }
  }
}
