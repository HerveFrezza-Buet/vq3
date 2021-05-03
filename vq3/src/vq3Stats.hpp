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

#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <limits>
#include <cmath>
#include <utility>
#include <tuple>


namespace vq3 {
  namespace stats {

    /**
     * This allows for computing the mean of a scalar collection.
     */
    class Mean {
    protected:
      double K;
      unsigned int n;
      double s;
      double tmp;
	
      friend std::ostream& operator<<(std::ostream& os, const Mean& ms) {
	os << "[mu = " << ms.mean() << ", n = " << ms.n << "]";
	return os;
      }
      
    public:

      class iterator {
      private:
	Mean& ms;
	friend class Mean;
	iterator(Mean& ms) : ms(ms) {}
      public:
	using difference_type   = long;
	using value_type        = double;
	using pointer           = double*;
	using reference         = double&;
	using iterator_category = std::output_iterator_tag;
	
	iterator()                           = delete;
	iterator(const iterator&)            = default;
	iterator& operator=(const iterator&) = default;
	iterator& operator++()    {return *this;}
	iterator& operator++(int) {return *this;}
	iterator& operator*()     {return *this;}
	iterator& operator=(double x) {ms = x; return *this;}
      };
	
      Mean(): K(0), n(0), s(0), tmp(0) {}
      Mean(const Mean&)            = default;
      Mean& operator=(const Mean&) = default;
      Mean(Mean&&)                 = default;
      Mean& operator=(Mean&&)      = default;

      /** 
       * This provides an output iterator for submitting values.
       */
      iterator output_iterator() {return iterator(*this);}
      

      void clear() {
	K=0;
	n=0;
	s=0;
	tmp=0;
      }
      
      /**
       * This notifies a value of the scalar collection.
       */
      Mean& operator=(double x) {
	if (n == 0) K = x;
	n++;
	tmp = x - K;
	s  += tmp;
	return *this;
      }

      double sum()  const {return n * K + s;}
      double mean() const {return K + s / n;}
	

      unsigned int nb_samples() {return n;}

      /**
       * @returns the mean.
       */
      double operator()() const {
	return mean();
      }
    };

    
    inline Mean mean() {return Mean();}

    /**
     * This allows for computing the mean and the variance or standard deviation of a scalar collection.
     */
    class MeanStd : public Mean {
    private:
      double s2;
	
      friend std::ostream& operator<<(std::ostream& os, const MeanStd& ms) {
	os << "[mu = " << ms.mean() << ", sigma^2 = " << std::sqrt(ms.variance()) << ", n = " << ms.n << "]";
	return os;
      }
      
    public:

      class iterator {
      private:
	MeanStd& ms;
	friend class MeanStd;
	iterator(MeanStd& ms) : ms(ms) {}
      public:
	using difference_type   = long;
	using value_type        = double;
	using pointer           = double*;
	using reference         = double&;
	using iterator_category = std::output_iterator_tag;
	
	iterator()                           = delete;
	iterator(const iterator&)            = default;
	iterator& operator=(const iterator&) = default;
	iterator& operator++()    {return *this;}
	iterator& operator++(int) {return *this;}
	iterator& operator*()     {return *this;}
	iterator& operator=(double x) {ms = x; return *this;}
      };
	
      MeanStd(): Mean(), s2(0) {}
      MeanStd(const MeanStd&)            = default;
      MeanStd& operator=(const MeanStd&) = default;
      MeanStd(MeanStd&&)                 = default;
      MeanStd& operator=(MeanStd&&)      = default;

      /** 
       * This provides an output iterator for submitting values.
       */
      iterator output_iterator() {return iterator(*this);}
      

      void clear() {
	Mean::clear();
	s2=0;
      }
      
      /**
       * This notifies a value of the scalar collection.
       */
      MeanStd& operator=(double x) {
	Mean::operator=(x);
	s2 += tmp*tmp;
	return *this;
      }
	
      double variance() const {
	if(n > 1)
	  return (s2 - (s*s)/double(n)) / (n-1.0);
	else
	  return 0;
      }

      /**
       * @returns the (mean, standard_deviation) pair.
       */
      std::pair<double,double> operator()() const {
	return {mean(), std::sqrt(variance()) };
      }
    };

    
    inline MeanStd mean_std() {return MeanStd();}
  }

  namespace spec {
    struct OnlineParam {
      OnlineParam();                    //!< the default constructor isrequired.
      OnlineParam(const OnlineParam&);  //!< the copy constructor is required.
      double alpha();                   //!< returns the update coefficient. This can be static.
      unsigned int min_updates();       //!< tells how many updates are required at least to be able to provide some trustable value. This can be static.
    };
  }
  
  namespace stats {
    namespace online {

      template<typename VALUE, typename ONLINE_PARAM>
      class MeanStd;
      
      /**
       * This computes the mean approximation thanks to a low-pass
       * recursive filter.
       *
       * \f${\overline x}_{n+1} = (1-\alpha){\overline x}_n + \alpha x_n \f$
       *
       * The value of \f$\alpha \in [0,1] \f$ is given by ONLINE_PARAM, that must fit vq3::spec::OnlineParam
       */
      template<typename VALUE, typename ONLINE_PARAM>
      class Mean {
      private:
	ONLINE_PARAM param;
	VALUE mu;
	bool ok = false;
	unsigned int nb = 0;

	friend class MeanStd<VALUE, ONLINE_PARAM>;
    
      public:
    
	Mean() : param(), mu(), ok(false), nb(0) {}
	Mean(const ONLINE_PARAM& p) : param(p), mu(), ok(false), nb(0) {} 
	Mean(const Mean&) = default;
	Mean& operator=(const Mean&) = default;

	/**
	 * @returns an estimation of the mean.
	 */
	VALUE mean() const {return mu;}

	/**
	 * Clears and forget the previous samples.
	 * @param v The initial value of the mean.
	 */
	void clear(const VALUE& v = VALUE()) {
	  ok = false;
	  nb = 0;
	  mu = v;
	}

	/**
	 * This enables to use an instance as a boolean, in order to test wether it is ready for providing a trustable value.
	 */
	operator bool() const {
	  return ok;
	}

	/**
	 * This considers a new sample.
	 */ 
	Mean& operator+=(const VALUE& v) {
	  mu += param.alpha()*(v-mu);
	  if(!(*this))
	    ok = ++nb >= param.min_updates();
	  return *this;
	}
      };
  
      template<typename VALUE, typename ONLINE_PARAM>
      Mean<VALUE, ONLINE_PARAM> mean(const ONLINE_PARAM& p) {return  Mean<VALUE, ONLINE_PARAM>(p);}
  

      /**
       * This computes mean and variance approximations thanks to a low-pass
       * recursive filter.
       *
       * \f${\overline x}_{n+1} = (1-\alpha){\overline x}_n + \alpha x_n \f$
       *
       * \f${\mathrm{var}}_{n+1}(x) = (1-\alpha){\mathrm{var}}_n(x) + \alpha (x_n - {\overline x}_{n+1})^2\f$
       *
       * The value of \f$\alpha \in [0,1] \f$ is given by ONLINE_PARAM, that must fit vq3::spec::OnlineParam
       */
      template<typename VALUE, typename ONLINE_PARAM>
      class MeanStd : public Mean<VALUE, ONLINE_PARAM> {
      private:
	Mean<VALUE, ONLINE_PARAM> var;
    
      public:
    
	MeanStd() = default;
	MeanStd(const ONLINE_PARAM& p) :  Mean<VALUE, ONLINE_PARAM>(p), var(p) {} 
	MeanStd(const MeanStd&) = default;
	MeanStd& operator=(const MeanStd&) = default;

	/**
	 * @returns an estimation of the variance.
	 */
	VALUE variance() const {return var.mean();}
    
	/**
	 * @returns the (mean, standard_deviation) pair.
	 */
	std::pair<double,double> operator()() const {
	  return {this->mean(), std::sqrt(variance())};
	}

	/**
	 * This considers a new sample.
	 */ 
	MeanStd& operator+=(const VALUE& v) {
	  this->Mean<VALUE, ONLINE_PARAM>::operator+=(v);
	  auto tmp = v - this->mean();
	  var += tmp*tmp;
	  return *this;
	}
      };
  
      template<typename VALUE, typename ONLINE_PARAM>
      MeanStd<VALUE, ONLINE_PARAM> mean_std(const ONLINE_PARAM& p) {return  MeanStd<VALUE, ONLINE_PARAM>(p);}
    }
        
    
    /**
     * This computes the d-shortest confidence interval from a sorted collection of double.
     * @param d belongs to [.5,1], it is the confidence level.
     * @param begin,end The sorted collection of double.
     * @return (min,max), the interval bounds.
     */
    template<typename ITER>
    std::pair<double, double> shortest_confidence_interval(double d,
							   const ITER& begin,
							   const ITER& end) {
      auto size = std::distance(begin, end);

      if(size == 0)
	throw std::runtime_error("vq3::utils::shortest_confidence_interval : empty collection");
      
      if(size <= 2) {
	auto v = *begin;
	return {v, v};
      }

      
      unsigned int window_size = std::max((unsigned int)(size * d + .5), (unsigned int)1);

      std::pair<double, double> res;
      double delta_min = std::numeric_limits<double>::max();

      ITER wit = begin;
      std::advance(wit,window_size);
      if(wit == end)
	throw std::runtime_error("vq3::stats::shortest_confidence_interval : the whole collection is in the interval.");
      for(auto it = begin; wit != end; ++it, ++wit)
	if(double delta = *wit - *it; delta < delta_min) {
	  delta_min = delta;
	  res = {*it, *wit};
	}
      
      return res;
    }
    
    /**
     * This computes the d-shortest confidence interval.
     * @param d belongs to [.5,1], it is the confidence level.
     * @param begin,end 
     * @param value_of A functor that extracts a double from *iter
     * @return (min,max), the interval bounds.
     */
    template<typename ITER, typename VALUE_OF>
    std::pair<double, double> shortest_confidence_interval(double d,
							   const ITER& begin,
							   const ITER& end,
							   const VALUE_OF& value_of) {
      auto size = std::distance(begin, end);

      if(size == 0)
	throw std::runtime_error("vq3::utils::shortest_confidence_interval : empty collection");
      
      if(size <= 2) {
	auto v = value_of(*begin);
	return {v, v};
      }

      std::vector<double> values(size);
      auto it = begin;
      for(auto& v : values) v = value_of(*(it++));
	
      std::sort(values.begin(), values.end());
      return shortest_confidence_interval(d, values.begin(), values.end());
    }

    
    class histogram {
    private:
	
      std::vector<double> values;


      unsigned int bin_quantile = 0;
      
      double nfit(double std, double nb) {
	return nb*bin_width*0.3989422804014327/std; // 0.39 = 1/sqrt(2*pi)
      }
      
    protected:

      std::optional<double> sci_conf;
      double bin_min           =  0;
      double bin_max           =  1;
      unsigned int bin_nb      = 10;
      unsigned int nb_hits     =  0;
      unsigned int nb_hits_sci =  0;
      double bin_width         =  1;

      std::vector<unsigned int> h;
      std::pair<double, double> sci;
      
      double mean;
      double std_dev;
      
      double smean;
      double sstd_dev;
	
    public:

      histogram()                 = default;
      histogram(const histogram&) = default;
      histogram(histogram&&)      = default;

      /**
       * Use this output iterator to insert values in the histogram.
       */
      auto output_iterator() {return std::back_inserter(values);}

      /**
       * Set sci related stuff. This triggers the SCI computation.
       * @param confidence in [.5, 1[. The confidence of the shortest confidence interval.
       */
      void set_sci(double confidence) {
	sci_conf = confidence;
      }

      /**
       * @param min, max [min, max] is the range of values used for the hostogram.
       * @param nb The number of bins in [min, max]
       */
      void set_bins(double min, double max, unsigned int nb) {
	bin_min      = min;
	bin_max      = max;
	bin_nb       = nb;
	bin_quantile = 0;
      }

      void set_bins(unsigned int quantile, unsigned int nb) {
	if(quantile < 50 || quantile > 100 || nb < 2) {
	  std::ostringstream ostr;
	  ostr << "vq3::demo::gnuplot::histogram::set_bins(quant=" << quantile << ", nb=" << nb << ") : quant must be in [50, 100] and nb must be >= 2.";
	  throw std::runtime_error(ostr.str());
	}
	bin_nb       = nb;
	bin_quantile = quantile;
      }
      
      void set_bins(unsigned int nb) {
	bin_nb = nb;
	if(values.begin() == values.end())
	  throw std::runtime_error("vq3::stats::histogram::set_bins(n) : empty values.");
	  
	auto [min_iter, max_iter] = std::minmax_element(values.begin(), values.end());
	bin_min = *min_iter;
	bin_max = (*max_iter)*1.001; // In order to include the maximal value.
      }

      /**
       * This builds the histogram and clears the stored values (for a next use of output_iterator()).
       */
      void make() {

	nb_hits = values.size();
	if(nb_hits < 2)
	  throw std::runtime_error("vq3::stats::histogram::make : at least 2 samples are required.");

	if(sci_conf) {
	  std::sort(values.begin(), values.end());
	  sci = shortest_confidence_interval(*sci_conf, values.begin(), values.end());
	}
	
	auto   ms = mean_std();
	auto  sms = mean_std();
	auto  out = ms.output_iterator();
	auto sout = sms.output_iterator();

	if(sci_conf) {
	  nb_hits_sci = 0;
	  for(auto it = values.begin(); it != values.end(); ++it) {
	    auto v = *it;
	    *(out++) =  v;
	    if(sci.first <= v && v <= sci.second) {
	      *(sout++) =  v;
	      ++nb_hits_sci;
	    }
	  }
	}
	else
	  std::copy(values.begin(), values.end(), out);
	
	std::tie(mean,  std_dev ) = ms ();
	if(sci_conf) 
	  std::tie(smean, sstd_dev) = sms(); 
	
	if(bin_quantile != 0) {
	  auto bin_sci = shortest_confidence_interval(.01*bin_quantile, values.begin(), values.end());
	  if(bin_sci.first == bin_sci.second)
	    throw std::runtime_error("vq3::demo::gnuplot::histogram::make : bin quantile leads to an empty value range.");
	      
	  bin_min = *(values.begin());
	  bin_max = *(values.end()-1);
	  bin_nb  = (unsigned int)((bin_max - bin_min) * bin_nb / (sci.second - sci.first) + .5);
	}

	  
	if(bin_max == bin_min)
	  throw std::runtime_error("vq3::::histogram::make : 0-sized value range.");

	bin_width = (bin_max-bin_min)/bin_nb;
	
	h = std::vector<unsigned int>(bin_nb, 0);
	double coef = bin_nb/(bin_max-bin_min);
	for(auto v : values) 
	  if(bin_min <= v && v < bin_max) {
	    unsigned int idx = (unsigned int)((v - bin_min)*coef);
	    h[idx]++;
	  }

	values.clear();
      }

      /**
       * @returns an vector of pairs (v,nb)  (ordered according to v) where v is the center of the value corresponding to the bin and nb the number of data that felt in that bin.
       */
      std::vector<std::pair<double, unsigned int>> bins() {
	unsigned int i = 0;
	std::vector<std::pair<double, unsigned int>> res;
	auto out =  std::back_inserter(res);
	double coef = 0;
	if(h.size() != 0)
	  coef = (bin_max-bin_min)/h.size();
	for(auto hit : h)
	  *(out++) = {bin_min + ((i++)+.5)*coef, hit};
	return res;
      }

      /**
       * @return (mu, sigma^2, ampl) such as ampl*exp(-.5*((x-mu)/sigma)^2) is the gaussian that approximates the histogram (i.e. same mean and variance).
       */
      std::tuple<double, double, double> gaussian_fit() {
	return {mean, std_dev*std_dev, nfit(std_dev, nb_hits)};
      }
      
      /**
       * @return (mu, sigma^2, ampl) such as ampl*exp(-.5*((x-mu)/sigma)^2) is the gaussian that approximates the histogram (i.e. same mean and variance). Here, the histogram is only considered in the SCI interval for the approximation.
       */
      std::tuple<double, double, double> gaussian_fit_sci() {

	  
	return {smean, sstd_dev*sstd_dev, nfit(sstd_dev, nb_hits_sci)};
      }
      
    };
  }
}
