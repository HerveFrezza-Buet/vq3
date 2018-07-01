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

#include <random>
#include <optional>
#include <utility>

namespace vq3 {
  namespace demo2d {

    class Point {
    public:
      double x,y;
    
      Point() : x(0), y(0) {}
      Point(double xx, double yy) : x(xx), y(yy) {}
      Point(const Point&)            = default;
      Point& operator=(const Point&) = default;

      Point& operator=(double val) {
	x = val;
	y = val;
	return *this;
      }
  
      bool operator==(const Point& p) const {
	return x == p.x && y == p.y;
      }
  
      bool operator!=(const Point& p) const {
	return x != p.x || y != p.y;
      }

      Point operator+(const Point& p) const {
	return {x+p.x, y+p.y};
      }
  
      Point operator-() const {
	return {-x,-y};
      }

      Point operator+() const {
	return {x,y};
      }

      Point operator-(const Point& p) const {
	return {x-p.x, y-p.y};
      }

      /**
       * The dot product
       */
      double operator*(const Point& p) const {
	return x*p.x + y*p.y;
      }

      /**
       * Component-wise product.
       */
      Point operator&(const Point& p) const {
	return {x*p.x, y*p.y};
      }
      
      Point operator/(const Point& p) const {
	return {x/p.x, y/p.y};
      }

      /**
       * vector cross product
       */
      double operator^(const Point& p) const {
	return x*p.y - y*p.x;
      }

      /**
       * Unitary vector.
       */
      Point operator*() const {
	const Point& me = *this;
	return me/std::sqrt(me*me);
      }

      Point operator*(double a) const {
	return {x*a, y*a};
      }

      Point operator/(double a) const {
	return (*this)*(1/a);
      }

      Point& operator+=(const Point& p) {
	x += p.x;
	y += p.y;
	return *this;
      }

      Point& operator-=(const Point& p) {
	x -= p.x;
	y -= p.y;
	return *this;
      }
  
      Point& operator*=(double a) {
	x *= a;
	y *= a;
	return *this;
      }
  
      Point& operator/=(double a) {
	(*this)*=(1/a);
	return *this;
      }

      Point transpose(const Point& p) {
	return {y, x};
      }
    };

    inline Point operator*(double a, const Point p) {
      return p*a;
    }

    inline std::ostream& operator<<(std::ostream& os, 
				    const Point& p) {
      os << '(' << p.x << ", " << p.y << ')';
      return os;
    }

    inline std::ostream& operator<<(std::ostream& os, 
				    const std::pair<Point, Point>& p) {
      os << '[' << p.first << ", " << p.second << ']';
      return os;
    }

    inline std::istream& operator>>(std::istream& is, 
				    Point& p) {
      char c;
      is >> c >> p.x >> c >> p.y >> c;
      return is;
    }

    template<typename RANDOM_ENGINE>
    Point uniform(RANDOM_ENGINE& rd, const Point& A, const Point& B) {
      std::uniform_real_distribution<double> uniform_dist(0., 1.);
      Point d = {uniform_dist(rd), uniform_dist(rd)};
      return A + (d & (B-A));
    }

    /**
     * Adds a value in [-r, +r[ to each components of A.
     */
    template<typename RANDOM_ENGINE>
    Point alter(RANDOM_ENGINE& rd, const Point& A, double r) {
      return A + uniform(rd, {-r, -r}, {r, r});
    }

    
  
    inline double d2(const Point& A, const Point& B) {
      Point tmp = B-A;
      return tmp*tmp;
    }

    inline double d(const Point& A, const Point& B) {
      return std::sqrt(d2(A,B));
    }

    inline Point min(const Point& A, const Point& B) {
      return {std::min(A.x,B.x),std::min(A.y,B.y)};
    }

    inline Point max(const Point& A, const Point& B) {
      return {std::max(A.x,B.x),std::max(A.y,B.y)};
    }


    /**
     * This returns (if it exists) the single point where the segments do intersect. If the segments overlap (i.e they share a common segment), no intersection is returned. If on segment is a point, no intersection is returned.
     */
    std::optional<Point> operator&&(const std::pair<Point, Point>& seg1, const std::pair<Point, Point>& seg2) {
      auto& P = seg1.first;
      auto  R = seg1.second - P;
      auto& Q = seg2.first;
      auto  S = seg2.second - Q;

      std::optional<Point> res;


      if(R == Point() || S == Point()) 
	return res;
	    
      double rs  = R^S;
      double qpr = (Q-P)^R;
	    
      if(rs != 0) {

	double t = ((Q-P)^S)/rs;
	double u = qpr/rs;
	if(0 <= t && t <= 1 && 0 <= u && u <= 1)
	  res = P + t*R;
      }
      
      return res;
    }
    
  }
}
