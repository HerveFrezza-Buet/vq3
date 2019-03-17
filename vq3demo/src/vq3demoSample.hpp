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

#include <vq3demoPoint.hpp>
#include <memory>
#include <iterator>
#include <cstdlib>
#include <stdexcept>
#include <functional>
#include <algorithm>
#include <iostream>
#include <tuple>

namespace vq3 {
  namespace demo2d {
    namespace sample {

      class bbox_exception : public std::runtime_error {
	using std::runtime_error::runtime_error;
      };
      
      class BBox {
      private:
	Point min,max;
	double a;
	bool empty;

	void check_empty() const {
	  if(empty)
	    throw bbox_exception("vq3::demo2d : empty bbox");
	}

	
	
      public:
	BBox() : a(0), empty(true) {}
	BBox(const BBox&) = default;
	BBox& operator=(const BBox&) = default;

	double area() const {return a;}

	friend std::ostream& operator<<(std::ostream& os, const BBox& bb) {
	  os << '[' << bb.min << ", " << bb.max << " : " << bb.a << ']';
	  return os;
	}
	
	BBox(const Point& min, const Point& max)
	  : min(min), max(max), empty(false)  {
	  auto diff = max - min;
	  a = diff.x*diff.y;
	}
	
	BBox(double xmin, double ymin, double width, double height)
	  : BBox(Point(xmin, ymin), Point(xmin+width, ymin+height)) {}

	BBox(double width, double height)
	  : BBox(-width/2, -height/2, width, height) {}
	
	BBox(double radius)
	  : BBox(-radius, -radius, 2*radius, 2*radius) {}

	BBox(const Point& A, const Point& B, const Point& C, const Point& D)
	  : BBox(demo2d::min(A, demo2d::min(B, demo2d::min(C, D))),
		 demo2d::max(A, demo2d::max(B, demo2d::max(C, D)))) {}

	template<typename ITER, typename POINT_OF>
	BBox(const ITER& begin, const ITER& end, const POINT_OF& point_of) {
	  empty = true;
	  a = 0;
	  if(auto it = begin; it != end) {
	    min = point_of(*(it++));
	    max = min;
	    while(it != end) {
	      auto pt = point_of(*(it++));
	      min = demo2d::min(min, pt);
	      max = demo2d::max(max, pt);
	    }
	    if(min != max) {
	      empty = false;
	      auto diff = max - min;
	      a = diff.x*diff.y;
	    }
	  }
	}

	Point size() const {
	  return max-min;
	}


	Point bottom_left()  const {check_empty(); return min;}
	Point bottom_right() const {check_empty(); return {max.x, min.y};}
	Point top_left()     const {check_empty(); return {min.x, max.y};}
	Point top_right()    const {check_empty(); return max;}
	
	template<typename RANDOM_ENGINE>
	Point uniform(RANDOM_ENGINE& rd) const {
	  check_empty();
	  return demo2d::uniform(rd, min, max);
	}

	BBox operator+(const Point& p) const {
	  if(empty)
	    return BBox();
	  return {min+p, max+p};
	}

	BBox operator*(const Point& p) const {
	  if(empty)
	    return BBox();
	  auto a = min & p;
	  auto b = max & p;
	  return {vq3::demo2d::min(a,b), vq3::demo2d::max(a,b)};
	}

	bool contains(const Point& p) const{
	  if(empty)
	    return false;
	  return p.x >= min.x
	    && p.x < max.x
	    && p.y >= min.y
	    && p.y < max.y;
	}
	
	bool intersects(const BBox& bb) const {
	  if(empty || bb.empty)
	    return false;
	  return min.x <= bb.max.x
	    && max.x >= bb.min.x
	    && min.y <= bb.max.y
	    && max.y >= bb.min.y;
	}

	BBox operator||(const BBox& bb) const {
	  if(empty)
	    return bb;
	  if(bb.empty)
	    return *this;
	  return {vq3::demo2d::min(min,bb.min), vq3::demo2d::max(max,bb.max)};
	}

	BBox operator&&(const BBox& bb) const {
	  if(!intersects(bb))
	    return BBox();
	  return {vq3::demo2d::max(min,bb.min), vq3::demo2d::min(max,bb.max)};
	}
      };

      namespace base_sampler {
	template<typename RANDOM_DEVICE>
	class Random {
	private:
	  RANDOM_DEVICE& rd;
	  double nb_samples_per_m2;
	  
	public:
	  
	  class iterator {
	  private:

	    friend class Random<RANDOM_DEVICE>;
	    
	    RANDOM_DEVICE* rd = nullptr;
	    BBox bbox;
	    unsigned int N;
	    
	    iterator(RANDOM_DEVICE& rd, const BBox& bbox, unsigned int N) : rd(&rd), bbox(bbox), N(N) {}
	    
	  public:
	    
	    using difference_type   = long;
	    using value_type        = vq3::demo2d::Point;
	    using pointer           = vq3::demo2d::Point*;
	    using reference         = vq3::demo2d::Point&;
	    using iterator_category = std::forward_iterator_tag;

	    iterator()                           = default;
	    iterator(const iterator&)            = default;
	    iterator& operator=(const iterator&) = default;
	  
	    iterator&  operator++()      {++N; return *this;}
	    iterator&  operator++(int)   {iterator res = *this; (*this)++; return res;}
	    value_type operator*() const {return bbox.uniform(*rd);}
	    bool       operator==(const iterator& other) const {return N == other.N;}
	    bool       operator!=(const iterator& other) const {return N != other.N;}
	  };

	  Random(RANDOM_DEVICE& rd, double nb_samples_per_m2)
	    : rd(rd), nb_samples_per_m2(nb_samples_per_m2) {}
	  Random()                         = delete;
	  Random(const Random&)            = delete;
	  Random& operator=(const Random&) = delete;


	  void operator=(double nb_samples_per_m2) {this->nb_samples_per_m2 = nb_samples_per_m2;}
	  std::pair<iterator, iterator> operator()(const BBox& bbox) const {
	    return {iterator(rd, bbox, 0), iterator(rd, bbox, (unsigned int)(nb_samples_per_m2*bbox.area()+.5))};
	  }
	};

	/**
	 * This samples an bounding box by tossing random samples uniformly distributed in it.
	 */
	template<typename RANDOM_DEVICE>
	auto random(RANDOM_DEVICE& rd, double nb_samples_per_m2) {return Random<RANDOM_DEVICE>(rd,nb_samples_per_m2);}
      }
      

      class Density {
      public:
	virtual BBox bbox() const = 0;
	virtual double operator()(const Point& p) const = 0;
	Density(const Density&)            = delete;
	Density& operator=(const Density&) = delete;
      protected:
	Density()                          = default;
      };

      using density = std::shared_ptr<const Density>;

      
      template<typename RANDOM_ENGINE, typename BASIC_SAMPLER>
      demo2d::Point get_one_sample(RANDOM_ENGINE& rd, const BASIC_SAMPLER& bs, const density& d) {
	auto bb = d->bbox();
	demo2d::Point value;
	bool notfound = true;
	auto unform   = std::uniform_real_distribution<double>(0,1);
	while(notfound) {
	  auto [begin, end] = bs(bb);
	  for(auto it = begin; notfound && it != end; ++it)
	    if(auto p = *it; uniform(rd) < (*d)(p)) {
	      value = p;
	      notfound = false;
	    }
	}
	return value;
      }
      
      template<typename RANDOM_ENGINE, typename BASE_SAMPLER>
      class SampleSet {
      private:
	density d;
	RANDOM_ENGINE& rd;

	typename BASE_SAMPLER::iterator begin_, end_;
	
      public:
	
	class iterator {
	private :
	  friend class SampleSet<RANDOM_ENGINE, BASE_SAMPLER>;
	  
	  const SampleSet<RANDOM_ENGINE, BASE_SAMPLER>* owner = nullptr;
	  typename BASE_SAMPLER::iterator iter;
	  demo2d::Point value;
	  
	  iterator(const SampleSet<RANDOM_ENGINE, BASE_SAMPLER>* owner, const typename BASE_SAMPLER::iterator& iter) : owner(owner), iter(iter), value() {}

	  void find_data() {
	    auto uniform = std::uniform_real_distribution<double>(0,1);
	    for(; iter != owner->end_; ++iter) 
	      if(auto p = *iter; uniform(owner->rd) < (*(owner->d))(p)) {
		value = p;
		break;
	      }
	  }
	  
	public:
	  
	  using difference_type   = demo2d::Point;
	  using value_type        = demo2d::Point;
	  using pointer           = demo2d::Point*;
	  using reference         = demo2d::Point&;
	  using iterator_category = std::forward_iterator_tag;
	  
	  iterator()                           = default;
	  iterator(const iterator&)            = default;
	  iterator& operator=(const iterator&) = default;
	  
	  bool operator==(const iterator& it) const {
	    return iter == it.iter;
	  }
	  
	  bool operator!=(const iterator& it) const {
	    return iter != it.iter;
	  }
	  
	  const demo2d::Point& operator*() const {
	    return value;
	  }
	  
	  iterator& operator++() {++iter, find_data(); return *this;}
	  iterator operator++(int) {iterator res = *this; ++*this; return res;}
	};
	
	SampleSet()                            = delete;
	SampleSet(const SampleSet&)            = default;
	SampleSet& operator=(const SampleSet&) = delete; // due to random engine reference.
	SampleSet(RANDOM_ENGINE& rd, BASE_SAMPLER& bs, density d) : d(d), rd(rd), begin_(), end_() {
	  std::tie(begin_, end_) = bs(d->bbox());
	}

	iterator begin() const {return iterator(this, begin_);}
	iterator end()   const {return iterator(this, end_);}
	
      };

      template<typename RANDOM_ENGINE, typename BASE_SAMPLER>
      SampleSet<RANDOM_ENGINE, BASE_SAMPLER> sample_set(RANDOM_ENGINE& rd, BASE_SAMPLER& bs, density d) {
	return SampleSet<RANDOM_ENGINE, BASE_SAMPLER>(rd, bs, d);
      }

      class Rectangle1 : public Density {
      private:

	const double& width;
	const double& height;	
	const double& value;
	
	Rectangle1(const Rectangle1&)            = delete;
	Rectangle1& operator=(const Rectangle1&) = delete;
	Rectangle1()                             = delete;

	Rectangle1(const double& width, const double& height, const double& value) : Density(), width(width), height(height), value(value) {}

	friend density rectangle(const double& width, const double& height, const double& value);
	
      public:
	virtual BBox bbox() const override {return BBox(width, height);}
	virtual double operator()(const Point& p) const override {
	  if(bbox().contains(p))
	    return value;
	  else
	    return 0;
	}
      };

      class Rectangle2 : public Density {
      private:
	const demo2d::Point& min;
	const demo2d::Point& max;
	const double&        value;
	
	Rectangle2(const Rectangle2&)            = delete;
	Rectangle2& operator=(const Rectangle2&) = delete;
	Rectangle2()                             = delete;

	Rectangle2(const demo2d::Point& min, const demo2d::Point& max, const double& value) : Density(), min(min), max(max), value(value) {}

	friend density rectangle(const demo2d::Point& min, const demo2d::Point& max, const double& value);
	
      public:
	virtual BBox bbox() const override {return BBox(min, max);}
	virtual double operator()(const Point& p) const override {
	  if(bbox().contains(p))
	    return value;
	  else
	    return 0;
	}
      };

      inline density rectangle(const double& width, const double& height, const double& value) {
	return density(new Rectangle1(width, height, value));
      }

      inline density rectangle(const demo2d::Point& min, const demo2d::Point& max, const double& value) {
	return density(new Rectangle2(min, max, value));
      }

      class Disk : public Density {
      private:

	const double& radius;
	const double& value;
	
	Disk(const Disk&)            = delete;
	Disk& operator=(const Disk&) = delete;
	Disk()                       = delete;

	Disk(const double& radius, const double& value) : Density(), radius(radius), value(value) {}

	friend density disk(const double& radius, const double& value);
	
      public:
	virtual BBox bbox() const override {return BBox(radius);}
	virtual double operator()(const Point& p) const override {
	  if(p*p <= radius*radius)
	    return value;
	  else
	    return 0;
	}
      };

      inline density disk(const double& radius, const double& value) {
	return density(new Disk(radius, value));
      }

      class Custom : public Density {
      private:

	BBox bb;
	std::function<double (const Point&)> f;
	
	Custom(const Custom&)            = delete;
	Custom& operator=(const Custom&) = delete;
	Custom()                         = delete;


	
      public:
	template<typename Func>
	Custom(const BBox& bb, const Func& f) : Density(), bb(bb), f(f) {}
	
	virtual BBox bbox() const override {return bb;}
	virtual double operator()(const Point& p) const override {
	  if(bb.contains(p))
	    return f(p);
	  else
	    return 0;
	}
      };

      template<typename Func>
      inline density custom(const BBox& bb, const Func& f) {
	return density(new Custom(bb, f));
      }

      class Translate : public Density {
      private:

	density op;
	const demo2d::Point& t;
	
	Translate(const Translate&)            = delete;
	Translate& operator=(const Translate&) = delete;
	Translate()                            = delete;

	Translate(density op, const demo2d::Point& translation) : op(op), t(translation) {}

	friend density operator+(density op, const demo2d::Point& translation);
	
      public:
	virtual BBox bbox() const override {return op->bbox() + t;}
	virtual double operator()(const Point& p) const override {
	  return (*op)(p - t);
	}
      };

      inline density operator+(density op, const demo2d::Point& translation) {
	return density(new Translate(op, translation));
      }

      class Scale : public Density {
      private:

	density op;
	const demo2d::Point& s;
	
	Scale(const Scale&)            = delete;
	Scale& operator=(const Scale&) = delete;
	Scale()                        = delete;

	Scale(density op, const demo2d::Point& scale) : op(op), s(scale) {}

	friend density operator*(density op, const demo2d::Point& scale);
	
      public:
	virtual BBox bbox() const override {return op->bbox() * s;}
	virtual double operator()(const Point& p) const override {
	  return (*op)(p/s);
	}
      };

      inline density operator*(density op, const demo2d::Point& scale) {
	return density(new Scale(op, scale));
      }
      
      class Rotate : public Density {
      private:

	density op;
	const double& theta;
	mutable double t;
	mutable double ct;
	mutable double st;
	
	Rotate(const Rotate&)            = delete;
	Rotate& operator=(const Rotate&) = delete;
	Rotate()                         = delete;

	void check_angle() const{
	  if(t != theta) {
	    t = theta;
	    double t_rad = t*3.141592653589793238463/180.0;
	    ct = cos(t_rad);
	    st = sin(t_rad);
	  }
	}
	
	Rotate(density op, const double& angle_deg)
	  : op(op), theta(angle_deg), t(0),
	    ct(1), st(0) {}

	Point positive_rotation(const Point& p) const {
	  return {ct*p.x - st*p.y, st*p.x + ct*p.y};
	}
	
	Point negative_rotation(const Point& p) const {
	  return {ct*p.x + st*p.y, -st*p.x + ct*p.y};
	}
	  
	friend density operator%(density op, const double& angle_deg);
	
      public:
	
	virtual BBox bbox() const override {
	  check_angle();
	  auto bb = op->bbox();
	  return  {positive_rotation(bb.bottom_left()),
	      positive_rotation(bb.bottom_right()),
	      positive_rotation(bb.top_left()),
	      positive_rotation(bb.top_right())};
	}
	
	virtual double operator()(const Point& p) const override {
	  check_angle();
	  return (*op)(negative_rotation(p));
	}
      };

      inline density operator%(density op, const double& angle_deg) {
	return density(new Rotate(op, angle_deg));
      }


      class Extrude : public Density {
      private:

	density op1;
	density op2;
	
	Extrude(const Extrude&)            = delete;
	Extrude& operator=(const Extrude&) = delete;
	Extrude()                          = delete;

	
	Extrude(density op1, density op2)
	  : op1(op1), op2(op2){}

	  
	friend density operator-(density op1, density op2);
	
      public:
	
	virtual BBox bbox() const override {
	  return op1->bbox();
	}
	
	virtual double operator()(const Point& p) const override {
	  return std::min((*op1)(p), 1 - (*op2)(p));
	}
      };

      inline density operator-(density op1, density op2) {
	return density(new Extrude(op1, op2));
      }
      
      class Intersection : public Density {
      private:

	density op1;
	density op2;
	
	Intersection(const Intersection&)            = delete;
	Intersection& operator=(const Intersection&) = delete;
	Intersection()                               = delete;

	
	Intersection(density op1, density op2)
	  : op1(op1), op2(op2){}

	  
	friend density operator&&(density op1, density op2);
	
      public:
	
	virtual BBox bbox() const override {
	  return op1->bbox() && op2->bbox();
	}
	
	virtual double operator()(const Point& p) const override {
	  if(bbox().contains(p))
	    return std::min((*op1)(p),(*op2)(p));
	  else
	    return 0;
	}
      };

      inline density operator&&(density op1, density op2) {
	return density(new Intersection(op1, op2));
      }
      
      class Union : public Density {
      private:

	density op1;
	density op2;
	
	Union(const Union&)            = delete;
	Union& operator=(const Union&) = delete;
	Union()                               = delete;

	
	Union(density op1, density op2)
	  : op1(op1), op2(op2){}

	  
	friend density operator||(density op1, density op2);
	
      public:
	
	virtual BBox bbox() const override {
	  return op1->bbox() || op2->bbox();
	}
	
	virtual double operator()(const Point& p) const override {
	    return std::max((*op1)(p),(*op2)(p));
	}
      };

      inline density operator||(density op1, density op2) {
	return density(new Union(op1, op2));
      }
    }
  }
}
