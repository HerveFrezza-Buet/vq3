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

#include <sstream>
#include <iomanip>
#include <string>

namespace vq3 {
  namespace demo {
    class irange {
    private:
      int first, last, step;
    public:
      class iterator {
      private :
	friend class irange;
	int i;
	iterator(int i) : i(i) {}
	
      public:
	
        using difference_type   = int;
        using value_type        = int;
        using pointer           = int*;
        using reference         = int&;
        using iterator_category = std::forward_iterator_tag;
	
	
	iterator() : iterator(0) {}
	iterator(const iterator&)             = default;
	iterator& operator= (const iterator&) = default;
	
	bool      operator==(const iterator& it) const {return i == it.i;}
	bool      operator!=(const iterator& it) const {return i != it.i;}
	int       operator*()                    const {return i;}
	iterator& operator++()                         {++i; return *this;}
	iterator  operator++(int)                      {iterator res = *this; ++*this; return res;}
      };
      
      irange(int first, int last, int step) : first(first), last(last), step(step) {}
      irange(int first, int last)           : irange(first, last, 1)                {}
      irange(int last)                      : irange(0, last, 1)                    {}
      irange()                              : irange(0, 0, 1)                       {}
      
      irange(const irange&)            = default;
      irange& operator=(const irange&) = default;

      iterator begin() const {return iterator(first);}
      iterator end()   const {return iterator(last);}
    };

    template<typename VALUE>
    class Range {
    private:

      irange rg;
      std::function<VALUE (int)> f;
      
    public:
      class iterator {
      private :
	friend class Range;
	irange::iterator i;
	std::function<VALUE (int)> f;
	
	iterator(irange::iterator i, std::function<VALUE (int)> f) : i(i), f(f) {}
	
      public:
        using difference_type   = int;
        using value_type        = VALUE;
        using pointer           = VALUE*;
        using reference         = VALUE&;
        using iterator_category = std::forward_iterator_tag;
	
	iterator()                            = default;
	iterator(const iterator&)             = default;
	iterator& operator= (const iterator&) = default;
	
	bool      operator==(const iterator& it) const {return i == it.i;}
	bool      operator!=(const iterator& it) const {return i != it.i;}
	VALUE     operator*()                    const {return f(*i);}
	iterator& operator++()                         {++i; return *this;}
	iterator  operator++(int)                      {iterator res = *this; ++*this; return res;}
      };
      
      Range()                        = default;
      Range(const Range&)            = default;
      Range& operator=(const Range&) = default;
      template<typename FUN> Range(const irange rg, const FUN& f) : rg(rg), f(f) {}

      iterator begin() const {return iterator(rg.begin(), f);}
      iterator end()   const {return iterator(rg.end(),   f);}
    };

    template<typename FUN>
    auto range(irange rg, const FUN& f)  {return Range<decltype(f(0))>(rg, f);}

    inline Range<double> range(double first, double last, unsigned int nb_intervals) {
      return range(irange(0, (int)(nb_intervals+1)),
		   [first, coef = (last - first)/nb_intervals](int i) {return first + i*coef;});
    }
    
    class videoframe_name {
    public:
      unsigned int num = 0;
      std::string name = "frame";
      std::string suffix = "png";
      
      videoframe_name()                                  = default;
      videoframe_name(const videoframe_name&)            = default;
      videoframe_name& operator=(const videoframe_name&) = default;
      videoframe_name(const std::string& name, const std::string& suffix)
	: num(0), name(name), suffix(suffix) {}

      std::string operator()(std::string frame_label) {
	std::ostringstream str;
	str << name << '-' << frame_label << '.' << suffix;
	return str.str();
      }
      
      std::string operator()(unsigned int num_frame) {
	num = num_frame;
	return (*this)();
      }
      
      std::string operator()() {
	std::ostringstream str;
	str << name << '-' << std::setw(6) << std::setfill('0') << num++ << '.' << suffix;
	return str.str();
      }
    };
  }
}
