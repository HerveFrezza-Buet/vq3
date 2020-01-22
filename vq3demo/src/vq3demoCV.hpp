
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
    }
  }
}
