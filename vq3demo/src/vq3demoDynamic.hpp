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

#include <cstdlib>
#include <vq3demoPoint.hpp>

namespace vq3 {
  namespace demo {
    namespace dyn {
      enum class bound : char {saturate = 's', wrap = 'w', bounce = 'b'};

      template<enum bound>
      class linear {
      private:

	double value;
	double min, max;
	double dv;
      
      public:
      
	linear()                         = delete;
	linear(const linear&)            = default;
	linear& operator=(const linear&) = default;

	linear(double init, double min, double max, double dv)
	  : value(init), min(min), max(max), dv(dv) {}

	const double& operator()() const {return value;}

	void operator++(int){++(*this);}
	void operator++();
      };

      template<> 
      inline void linear<bound::saturate>::operator++() {
	value += dv;
	if(value > max)
	  value = max;
	else if (value < min)
	  value = min;
      }
    
      template<> 
      inline void linear<bound::wrap>::operator++() {
	value += dv;
	if(value > max)
	  value -= max-min;
	else if(value < min)
	  value += max - min;
      }
    
      template<> 
      inline void linear<bound::bounce>::operator++() {
	value += dv;
	if(value > max) {
	  value = 2*max-value;
	  dv    = -dv;
	}
	else if(value < min) {
	  value = 2*min-value;
	  dv    = -dv;
	}
      }

      class cos : private linear<bound::wrap> {
      private:

	double min;
	double ampl;
	double cos_value;
	void update() {
	  cos_value = min+ampl*(1+std::cos(this->linear<bound::wrap>::operator()()*3.141592653589793238463/180.0));
	}
      
      public:
	cos()                      = delete;
	cos(const cos&)            = default;
	cos& operator=(const cos&) = default;

	cos(double min, double max, double init_angle_deg, double dangle_deg)
	  : linear<bound::wrap>(init_angle_deg, -180, 180, dangle_deg), min(min), ampl(.5*(max-min)) {
	  update();
	}
      
	const double& operator()() const {return cos_value;}

	void operator++(int){++(*this);}
	void operator++() {
	  this->linear<bound::wrap>::operator++();
	  update();
	}
      };
      
      class sin : private linear<bound::wrap> {
      private:

	double min;
	double ampl;
	double sin_value;
	void update() {
	  sin_value = min+ampl*(1+std::sin(this->linear<bound::wrap>::operator()()*3.141592653589793238463/180.0));
	}
      
      public:
	sin()                      = delete;
	sin(const sin&)            = default;
	sin& operator=(const sin&) = default;

	sin(double min, double max, double init_angle_deg, double dangle_deg)
	  : linear<bound::wrap>(init_angle_deg, -180, 180, dangle_deg), min(min), ampl(.5*(max-min)) {
	  update();
	}
      
	const double& operator()() const {return sin_value;}

	void operator++(int){++(*this);}
	void operator++() {
	  this->linear<bound::wrap>::operator++();
	  update();
	}
      };

    }
  }
  
  namespace demo2d {
    namespace dyn {
      
      class circle : private vq3::demo::dyn::linear<vq3::demo::dyn::bound::wrap> {
      private:
      
	double r;
	Point pt;
      
	void update() {
	  double angle = this->vq3::demo::dyn::linear<vq3::demo::dyn::bound::wrap>::operator()();
	  pt.x = r*std::cos(angle*3.141592653589793238463/180.0);
	  pt.y = r*std::sin(angle*3.141592653589793238463/180.0);
	}
      
      public:
	circle()                      = delete;
	circle(const circle&)            = default;
	circle& operator=(const circle&) = default;

	circle(double r, double init_angle_deg, double dangle_deg)
	  : vq3::demo::dyn::linear<vq3::demo::dyn::bound::wrap>(init_angle_deg, -180, 180, dangle_deg), r(r) {
	  update();
	}
      
	const Point& operator()() const {return pt;}

	void operator++(int){++(*this);}
	void operator++() {
	  this->vq3::demo::dyn::linear<vq3::demo::dyn::bound::wrap>::operator++();
	  update();
	}
      };

      template<typename DYNX, typename DYNY>
      class Point {
      private:

	const DYNX& dynx;
	const DYNY& dyny;

	mutable vq3::demo2d::Point point;

	void update() {
	  point.x = dynx();
	  point.y = dyny();
	}
	
      public:
	Point()                        = delete;
	Point(const Point&)            = default;
	Point& operator=(const Point&) = default;

	Point(const DYNX& dynx, const DYNY& dyny) : dynx(dynx), dyny(dyny), point(dynx(), dyny()) {}
      
	const vq3::demo2d::Point& operator()() const {return point;}
	void operator++(int){++(*this);}
	void operator++() {update();}
	
      };

      template<typename DYNX>
      class Point<DYNX, double> {
      private:

	const DYNX& dynx;

	mutable vq3::demo2d::Point point;
	
	void update() {
	  point.x = dynx();
	}
	
      public:
	Point()                        = delete;
	Point(const Point&)            = default;
	Point& operator=(const Point&) = default;

	Point(const DYNX& dynx, const double& y) : dynx(dynx), point(dynx(),y) {}
      
	const vq3::demo2d::Point& operator()() const {return point;}
	void operator++(int){++(*this);}
	void operator++() {update();}
      };

      template<typename DYNY>
      class Point<double, DYNY> {
      private:

	const DYNY& dyny;

	mutable vq3::demo2d::Point point;
	
	void update() {
	  point.y = dyny();
	}
	
      public:
	Point()                        = delete;
	Point(const Point&)            = default;
	Point& operator=(const Point&) = default;

	Point(const double& x, const DYNY& dyny) : dyny(dyny), point(x,dyny()) {}
      
	const vq3::demo2d::Point& operator()() const {return point;}
	void operator++(int){++(*this);}
	void operator++() {update();}
      };

      
      template<typename DYNX, typename DYNY>
      Point<DYNX, DYNY> point(const DYNX& x, const DYNY& y) {
	return Point<DYNX, DYNY>(x,y);
      }
      
    }
  }
}
  
