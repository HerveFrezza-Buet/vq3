
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

#include <limits>
#include <utility>
#include <vector>
#include <list>
#include <map>
#include <stack>
#include <type_traits>
#include <iterator>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <vq3Graph.hpp>

namespace vq3 {

  namespace error {
    struct empty_accumulation : public std::runtime_error {
      using std::runtime_error::runtime_error;
    };
  }
  
  namespace utils {
    
    template<typename INCREMENTABLE, typename NB_TYPE>
    struct accum {
      NB_TYPE nb;
      INCREMENTABLE     value;

      accum(const accum&)             = default;
      accum& operator=(const accum&)  = default;
      accum(accum&&)                  = default;
      accum& operator=(accum&&)       = default;

      accum() : nb(0), value() {}

      /**
       * Returns the average. The template argument is the type for casting nb befor the division.
       * @return value/(NB_AVG_TYPE)nb;
       */
      template<typename NB_AVG_TYPE = NB_TYPE>
      auto average() const {return value/static_cast<NB_AVG_TYPE>(nb);}
	
      void clear() {
	nb = 0;
	value = INCREMENTABLE();
      }

      /** nb = 1, value = v */
      /** nb += 1; value += v; */
      accum& operator=(const INCREMENTABLE& v) {
	nb    = 1;
	value = v;
	return *this;
      }

      /** nb += coef; value += v; */
      accum& increment(NB_TYPE coef, const INCREMENTABLE& v) {
	nb    += coef;
	value += coef * v;
	return *this;
      }
      
      /** nb += 1; value += v; */
      accum& operator+=(const INCREMENTABLE& v) {
	++nb;
	value += v;
	return *this;
      }

      /** nb += a.nb; value += a.value; */
      accum& operator+=(const accum& a) {
	nb    += a.nb;
	value += a.value;
	return *this;
      }
    };


    /**
     * @param g The graph
     * @param value_of value_of(ref_vertex) gives the value to be stored in the map.
     * @return A map where keys are the vertex references and the values the one calculated by value_of.
     */
    template<typename GRAPH, typename VALUE_OF_REF_VERTEX>
    auto make_vertex_table(GRAPH& g, const VALUE_OF_REF_VERTEX& value_of) {
      std::map<typename GRAPH::ref_vertex, typename std::result_of<VALUE_OF_REF_VERTEX (typename GRAPH::ref_vertex)>::type> res;
      g.foreach_vertex([&res, &value_of](const typename GRAPH::ref_vertex& ref_v) {res[ref_v] =  value_of(ref_v);});
      return res;
    }

    /**
     * Collects all the vertex references in a collection.
     * @param g The graph
     * @param outit An output iterator for storing ref_vertex values.
     */
    template<typename GRAPH, typename OUTPUT_IT>
    void collect_vertices(GRAPH& g, const OUTPUT_IT& outit) {
      auto out = outit;
      g.foreach_vertex([&out](const typename GRAPH::ref_vertex& ref_v) {*(out++) = ref_v;});
    }

    /**
     * Collects all the edge references in a collection.
     * @param g The graph
     * @param outit An output iterator for storing ref_edge values.
     */
    template<typename GRAPH, typename OUTPUT_IT>
    void collect_edges(GRAPH& g, const OUTPUT_IT& outit) {
      auto out = outit;
      g.foreach_edge([&out](const typename GRAPH::ref_edge& ref_e) {
	  auto extr = ref_e->extremities();
	  if(invalid_extremities(extr))
	    ref_e->kill();
	  else
	    *(out++) = ref_e;});
    }
 

    /**
     * Finds the closest vertex.
     * @param g the graph.
     * @param sample We want the vertex closest to this sample.
     * @param distance computes the distance as distance(vertex_value, sample).
     * @param closest_distance_value returns by reference the closest distance value.
     * @return The closest vertex reference. It can be nullptr if the graph is empty or if all distances to the sample are +infty.
     */
    template<typename GRAPH, typename SAMPLE, typename DISTANCE>
    typename GRAPH::ref_vertex closest(GRAPH& g, const SAMPLE& sample, const DISTANCE& distance, double& closest_distance_value) {
      typename GRAPH::ref_vertex res = nullptr;
      double dist = std::numeric_limits<double>::max();
      g.foreach_vertex([&dist, &res, &sample, &distance](const typename GRAPH::ref_vertex& ref_v) {
	  if(double d = distance((*ref_v)(), sample); d < dist) {
	    dist = d;
	    res  = ref_v;
	  }

	});
      closest_distance_value = dist;
      return res;
    }

    /**
     * Finds the closest vertex.
     * @param g the graph.
     * @param sample We want the vertex closest to this sample.
     * @param distance computes the distance as distance(vertex_value, sample).
     * @return The closest vertex reference. It can be nullptr if the graph is empty or if all distances to the sample are +infty.
     */
    template<typename GRAPH, typename SAMPLE, typename DISTANCE>
    typename GRAPH::ref_vertex closest(GRAPH& g, const SAMPLE& sample, const DISTANCE& distance) {
      double d;
      return closest(g, sample, distance, d);
    }

    /**
     * Finds the two closest vertices.
     * @param g the graph.
     * @param sample We want the vertex closest to this sample.
     * @param distance computes the distance as distance(vertex_value, sample).
     * @param closest_distance_values returns by reference the two closest distance values.
     * @return The closest vertices references as a pair. They can be nullptr if the graph is empty or if all distances to the sample are +infty.
     */
    template<typename GRAPH, typename SAMPLE, typename DISTANCE>
    typename std::pair<typename GRAPH::ref_vertex, typename GRAPH::ref_vertex> two_closest(GRAPH& g, const SAMPLE& sample, const DISTANCE& distance, std::pair<double, double>& closest_distance_values) {
      std::pair<typename GRAPH::ref_vertex, typename GRAPH::ref_vertex> res = {nullptr, nullptr};
      double dist1 = std::numeric_limits<double>::max();
      double dist2 = std::numeric_limits<double>::max();
      g.foreach_vertex([&dist1, &dist2, &res, &sample, &distance](const typename GRAPH::ref_vertex& ref_v) {
	  double d = distance((*ref_v)(), sample);
	  if(d < dist1) {
	    dist2      = dist1;
	    res.second = res.first;
	    dist1      = d;
	    res.first  = ref_v;
	  }
	  else if(d < dist2) {
	    dist2      = d;
	    res.second = ref_v;
	  }
	});
      closest_distance_values = {dist1, dist2};
      return res;
    }


    /**
     * Finds the two closest vertices.
     * @param g the graph.
     * @param sample We want the vertex closest to this sample.
     * @param distance computes the distance as distance(vertex_value, sample).
     * @return The closest vertices references as a pair. They can be nullptr if the graph is empty or if all distances to the sample are +infty.
     */
    template<typename GRAPH, typename SAMPLE, typename DISTANCE>
    typename std::pair<typename GRAPH::ref_vertex, typename GRAPH::ref_vertex> two_closest(GRAPH& g, const SAMPLE& sample, const DISTANCE& distance) {
      std::pair<double, double> dists;
      return two_closest(g, sample, distance, dists);
    }

    /**
     * Builds a graph as a grid (indeed, it add vertices and edges, so the graph is usually empty when this function is called).
     * @param g the graph.
     * @param width, height The grid dimensions.
     * @param v_of A function such as v_of(w, h) is the vertex value at position (w,h).
     * @param e_of A function such as e_of(w, h, ww, hh) is the edge value for the edge linking (w,h) to (ww, hh).
     * @param loop_w True means that left and right are connected.
     * @param loop_h True means that top and bottom are connected.
     * @return a height-sized vector of width-sized vector of vertex references.
     */
    template<typename GRAPH, typename VERTEX_VALUE_OF, typename EDGE_VALUE_OF>
    auto make_grid(GRAPH& g, unsigned int width, unsigned int height,
		   const VERTEX_VALUE_OF& v_of, const EDGE_VALUE_OF& e_of,
		   bool loop_w = false, bool loop_h = false) {
      std::vector<std::vector<typename GRAPH::ref_vertex>> lines;
      std::vector<typename GRAPH::ref_vertex> line(width);
      auto lout = std::back_inserter(lines);

      unsigned int w, h;
      unsigned int width_  = width  - 1;
      unsigned int height_ = height - 1;
      
      for(h = 0; h < height; ++h) {
	for(w = 0; w < width; ++w)
	  line[w] = (g += v_of(w, h));
	*(lout++) = line;
      }
      
      for(h = 0; h < height_; ++h) {
	for(w = 0; w < width_; ++w) {
	  g.connect(lines[h][w], lines[h][w + 1], e_of(w, h, w + 1, h));
	  g.connect(lines[h][w], lines[h + 1][w], e_of(w, h, w, h + 1));
	}
	g.connect(lines[h][w], lines[h + 1][w], e_of(w, h, w, h + 1));
      }
      for(w = 0; w < width_; ++w)
       	g.connect(lines[h][w], lines[h][w + 1], e_of(w, h, w + 1, h));

      if(loop_w && width > 1) 
	for(h = 0; h < height; ++h)
	  g.connect(lines[h][0], lines[h][width_], e_of(0, h, width_, h));

      if(loop_h && height > 1) 
	for(w = 0; w < width; ++w)
	  g.connect(lines[0][w], lines[height_][w], e_of(w, 0, w, height_));
	  
      
      return lines;
    }


    /**
     * Builds a graph as a grid (indeed, it add vertices and edges, so the graph is usually empty when this function is called). This is dedicated for graphs where edges have no values.
     * @param g the graph.
     * @param width, height The grid dimensions.
     * @param v_of A function such as v_of(w, h) is the vertex value at position (w,h).
     * @param loop_w True means that left and right are connected.
     * @param loop_h True means that top and bottom are connected.
     * @return a height-sized vector of width-sized vector of vertex references.
     */
    template<typename GRAPH, typename VERTEX_VALUE_OF>
    auto make_grid(GRAPH& g, unsigned int width, unsigned int height,
		   const VERTEX_VALUE_OF& v_of,
		   bool loop_w = false, bool loop_h = false) {
      std::vector<std::vector<typename GRAPH::ref_vertex>> lines;
      std::vector<typename GRAPH::ref_vertex> line(width);
      auto lout = std::back_inserter(lines);

      unsigned int w, h;
      unsigned int width_  = width  - 1;
      unsigned int height_ = height - 1;
      
      for(h = 0; h < height; ++h) {
	for(w = 0; w < width; ++w)
	  line[w] = (g += v_of(w, h));
	*(lout++) = line;
      }
      
      for(h = 0; h < height_; ++h) {
	for(w = 0; w < width_; ++w) {
	  g.connect(lines[h][w], lines[h][w + 1]);
	  g.connect(lines[h][w], lines[h + 1][w]);
	}
	g.connect(lines[h][w], lines[h + 1][w]);
      }
      for(w = 0; w < width_; ++w)
       	g.connect(lines[h][w], lines[h][w + 1]);

      if(loop_w && width > 1) 
	for(h = 0; h < height; ++h)
	  g.connect(lines[h][0], lines[h][width_]);

      if(loop_h && height > 1) 
	for(w = 0; w < width; ++w)
	  g.connect(lines[0][w], lines[height_][w]);
      
      return lines;
    }

    /**
     * Set a value to each vertex tag.
     */
    template<typename GRAPH>
    void clear_vertex_tags(GRAPH& g, bool value) {
      g.foreach_vertex([value](const typename GRAPH::ref_vertex& ref_v) {(*ref_v)().vq3_tag = value;});
    }

    /**
     * Set a value to each edge tag.
     */
    template<typename GRAPH>
    void clear_edge_tags(GRAPH& g, bool value) {
      g.foreach_edge([value](const typename GRAPH::ref_edge& ref_e) {
	  auto extr = ref_e->extremities();
	  if(invalid_extremities(extr))
	    ref_e->kill();
	  else
	    (*ref_e)().vq3_tag = value;
	});
    }

    /**
     * Sets a value to each edge and vertice tag.
     */
    template<typename GRAPH>
    void clear_all_tags(GRAPH& g, bool value) {
      clear_vertex_tags(g, value);
      clear_edge_tags(g, value);
    }


    /**
     * Set a value to each vertex label.
     */
    template<typename GRAPH>
    void clear_vertex_labels(GRAPH& g, unsigned int value) {
      g.foreach_vertex([value](const typename GRAPH::ref_vertex& ref_v) {(*ref_v)().vq3_label = value;});
    }

    /**
     * Set a value to each edge label.
     */
    template<typename GRAPH>
    void clear_edge_labels(GRAPH& g, unsigned int value) {
      g.foreach_edge([value](const typename GRAPH::ref_edge& ref_e) {
	  auto extr = ref_e->extremities();
	  if(invalid_extremities(extr))
	    ref_e->kill();
	  else
	    (*ref_e)().vq3_label = value;
	});
    }

    /**
     * Sets a value to each edge and vertice label.
     */
    template<typename GRAPH>
    void clear_all_labels(GRAPH& g, unsigned int value) {
      clear_vertex_labels(g, value);
      clear_edge_labels(g, value);
    }
    
    /**
     * Set a value to each vertex efficiency.
     */
    template<typename GRAPH>
    void clear_vertex_efficiencies(GRAPH& g, bool value) {
      g.foreach_vertex([value](const typename GRAPH::ref_vertex& ref_v) {(*ref_v)().vq3_efficient = value;});
    }

    /**
     * Set a value to each edge efficiency.
     */
    template<typename GRAPH>
    void clear_edge_efficiencies(GRAPH& g, bool value) {
      g.foreach_edge([value](const typename GRAPH::ref_edge& ref_e) {
	  auto extr = ref_e->extremities();
	  if(invalid_extremities(extr))
	    ref_e->kill();
	  else
	    (*ref_e)().vq3_efficient = value;
	});
    }

    /**
     * Sets a value to each edge and vertice efficiency.
     */
    template<typename GRAPH>
    void clear_all_efficiencies(GRAPH& g, bool value) {
      clear_vertex_efficiencies(g, value);
      clear_edge_efficiencies(g, value);
    }

    /**
     * This splits a collection into nb parts.
     * @return a vector of {begin, end} pairs for each subpart.
     */
    template<typename IT>
    auto split(const IT& begin, const IT& end, unsigned int nb) {
      std::vector< std::pair<IT, IT> > pairs;
      pairs.reserve(nb);
      auto out       = std::back_inserter(pairs);
      auto coef      = std::distance(begin, end)/double(nb);
      unsigned int q = 0;
      IT last        = begin;
      IT it          = begin;
      for(unsigned int i = 1; i <= nb; ++i) {
	auto first = last;
	auto p     = q;
	q = (unsigned int)(i*coef+.5);
	std::advance(it, q-p);
	last = it;
	*(out++) = {first, last};
      }
      return pairs;
    }


    /**
     * Iterates on efficient edges of a node.
     */
    template<typename REF_VERTEX, typename EDGE_FUN>
    void foreach_efficient_edge(const REF_VERTEX& ref_v, const EDGE_FUN& fun) {
      ref_v->foreach_edge([fun](const typename REF_VERTEX::element_type::ref_edge_type& ref_e) {
	  if((*ref_e)().vq3_efficient)
	    fun(ref_e);
	});
    }


    namespace savitzky_golay {

      /**
       * This is a base class for Savitsky-Golay signal estimators.
       */
      template<typename VALUE, unsigned int ORDER, unsigned int WINDOW_SIZE>
      class Estimator {
      protected:
	
	std::deque<VALUE> window;
	mutable std::array<std::optional<VALUE>, ORDER+1> value; // value[i] is the ith order derivative... optional since it may not be available.
	
      public:
	
	Estimator()                            = default;
	Estimator(const Estimator&)            = default;
	Estimator(Estimator&&)                 = default;
	Estimator& operator=(const Estimator&) = default;
	Estimator& operator=(Estimator&&)      = default;

	/**
	 * This registers a new sample.
	 */
	void operator+=(const VALUE& v) {
	  window.emplace_back(v);
	  if(window.size() > WINDOW_SIZE)
	    window.pop_front();
	  for(auto& v : value) v.reset();
	}

	/**
	 * Returns the estimate of the Oth derivative. It is an optional since it may not be available.
	 */
	template<unsigned int O>
	const std::optional<VALUE>& get() const {
	  static_assert(O <= ORDER, "vq3::utils::savitzky_golay::get<O> : O is too high.");
	  return value[O];
	}
      };

      namespace constant_timestep {
	
	template<unsigned int ORDER, unsigned int WINDOW_SIZE, unsigned int DEGREE>
	struct sg_norm {
	  // We use "ORDER < 0" instead of "false" to make the static
	  // assert dependent on T. Otherwise, the compilation failure
	  // occurs systematically.
	  static_assert(ORDER < 0, "The parameters passed to vq3::utils::savitzky_golay::constant_timestep::estimator<...> leads to an unimplemented case (sg_norm fails)."); 
	};
	
	template<typename VALUE, unsigned int ORDER, unsigned int WINDOW_SIZE, unsigned int DEGREE>
	struct sg_compute {
	  // We use "ORDER < 0" instead of "false" to make the static
	  // assert dependent on T. Otherwise, the compilation failure
	  // occurs systematically.
	  static_assert(ORDER < 0, "The parameters passed to vq3::utils::savitzky_golay::constant_timestep::estimator<...> leads to an unimplemented case (sg_compute fails).");
	};

	// The coefficients are obtained from Robert Mellet's web page.
	// http://robert.mellet.pagesperso-orange.fr/rgrs_pol/regrs_06.htm
	
	//
	// win = 9, degree = 2.
	//
	
	template<> struct sg_norm<0, 9, 2> {double value; sg_norm(double step) {value = 1/231.0;}};
	template<typename VALUE> struct sg_compute<VALUE, 0, 9, 2> {
	  VALUE value;
	  template<typename ITER>
	  sg_compute(ITER iter) {
	    value  = (*(iter++)) *  -21;
	    value += (*(iter++)) *   14;
	    value += (*(iter++)) *   39;
	    value += (*(iter++)) *   54;
	    value += (*(iter++)) *   59;
	    value += (*(iter++)) *   54;
	    value += (*(iter++)) *   39;
	    value += (*(iter++)) *   14;
	    value += (*(iter++)) *  -21;
	  }
	};
	
	template<> struct sg_norm<1, 9, 2> {double value; sg_norm(double step) {value = 1/(60*step);}};
	template<typename VALUE> struct sg_compute<VALUE, 1, 9, 2> {
	  VALUE value;
	  template<typename ITER>
	  sg_compute(ITER iter) {
	    value  = (*(iter++)) *    -4;
	    value += (*(iter++)) *    -3;
	    value += (*(iter++)) *    -2;
	    value += (*(iter++)) *    -1;
	    ++iter;
	    value += (*(iter++)) *     1;
	    value += (*(iter++)) *     2;
	    value += (*(iter++)) *     3;
	    value += (*(iter++)) *     4;
	  }
	};
	
	template<> struct sg_norm<2, 9, 2> {double value; sg_norm(double step) {value = 1/(462*step*step);}};
	template<typename VALUE> struct sg_compute<VALUE, 2, 9, 2> {
	  VALUE value;
	  template<typename ITER>
	  sg_compute(ITER iter) {
	    value  = (*(iter++)) *   28;
	    value += (*(iter++)) *    7;
	    value += (*(iter++)) *   -8;
	    value += (*(iter++)) *  -17;
	    value += (*(iter++)) *  -20;
	    value += (*(iter++)) *  -17;
	    value += (*(iter++)) *   -8;
	    value += (*(iter++)) *    7;
	    value += (*(iter++)) *   28;
	  }
	};
	

	//
	// win = 15, degree = 2.
	//
	
	template<> struct sg_norm<0, 15, 2> {double value; sg_norm(double step) {value = 1/1105.0;}};
	template<typename VALUE> struct sg_compute<VALUE, 0, 15, 2> {
	  VALUE value;
	  template<typename ITER>
	  sg_compute(ITER iter) {
	    value  = (*(iter++)) *  -78;
	    value += (*(iter++)) *  -13;
	    value += (*(iter++)) *   42;
	    value += (*(iter++)) *   87;
	    value += (*(iter++)) *  122;
	    value += (*(iter++)) *  147;
	    value += (*(iter++)) *  162;
	    value += (*(iter++)) *  167;
	    value += (*(iter++)) *  162;
	    value += (*(iter++)) *  147;
	    value += (*(iter++)) *  122;
	    value += (*(iter++)) *   87;
	    value += (*(iter++)) *   42;
	    value += (*(iter++)) *  -13;
	    value += (*(iter++)) *  -78;
	  }
	};
	
	template<> struct sg_norm<1, 15, 2> {double value; sg_norm(double step) {value = 1/(280*step);}};
	template<typename VALUE> struct sg_compute<VALUE, 1, 15, 2> {
	  VALUE value;
	  template<typename ITER>
	  sg_compute(ITER iter) {
	    value  = (*(iter++)) *  -7;
	    value += (*(iter++)) *  -6;
	    value += (*(iter++)) *  -5;
	    value += (*(iter++)) *  -4;
	    value += (*(iter++)) *  -3;
	    value += (*(iter++)) *  -2;
	    value += (*(iter++)) *  -1;
	    ++iter;
	    value += (*(iter++)) *   1;
	    value += (*(iter++)) *   2;
	    value += (*(iter++)) *   3;
	    value += (*(iter++)) *   4;
	    value += (*(iter++)) *   5;
	    value += (*(iter++)) *   6;
	    value += (*(iter++)) *   7;
	  }
	};
	
	template<> struct sg_norm<2, 15, 2> {double value; sg_norm(double step) {value = 1/(6188*step*step);}};
	template<typename VALUE> struct sg_compute<VALUE, 2, 15, 2> {
	  VALUE value;
	  template<typename ITER>
	  sg_compute(ITER iter) {
	    value  = (*(iter++)) *   91;
	    value += (*(iter++)) *   52;
	    value += (*(iter++)) *   19;
	    value += (*(iter++)) *   -8;
	    value += (*(iter++)) *  -29;
	    value += (*(iter++)) *  -44;
	    value += (*(iter++)) *  -53;
	    value += (*(iter++)) *  -56;
	    value += (*(iter++)) *  -53;
	    value += (*(iter++)) *  -44;
	    value += (*(iter++)) *  -29;
	    value += (*(iter++)) *   -8;
	    value += (*(iter++)) *   19;
	    value += (*(iter++)) *   52;
	    value += (*(iter++)) *   91;
	  }
	};

	//
	// win = 21, degree = 2.
	//
	
	template<> struct sg_norm<0, 21, 2> {double value; sg_norm(double step) {value = 1/3059.0;}};
	template<typename VALUE> struct sg_compute<VALUE, 0, 21, 2> {
	  VALUE value;
	  template<typename ITER>
	  sg_compute(ITER iter) {
	    value  = (*(iter++)) *  -171;
	    value += (*(iter++)) *   -76;
	    value += (*(iter++)) *     9;
	    value += (*(iter++)) *    84;
	    value += (*(iter++)) *   149;
	    value += (*(iter++)) *   204;
	    value += (*(iter++)) *   249;
	    value += (*(iter++)) *   284;
	    value += (*(iter++)) *   309;
	    value += (*(iter++)) *   324;
	    value += (*(iter++)) *   329;
	    value += (*(iter++)) *   324;
	    value += (*(iter++)) *   309;
	    value += (*(iter++)) *   284;
	    value += (*(iter++)) *   249;
	    value += (*(iter++)) *   204;
	    value += (*(iter++)) *   149;
	    value += (*(iter++)) *    84;
	    value += (*(iter++)) *     9;
	    value += (*(iter++)) *   -76;
	    value += (*(iter++)) *  -171;
	  }
	};
	
	template<> struct sg_norm<1, 21, 2> {double value; sg_norm(double step) {value = 1/(770*step);}};
	template<typename VALUE> struct sg_compute<VALUE, 1, 21, 2> {
	  VALUE value;
	  template<typename ITER>
	  sg_compute(ITER iter) {
	    value  = (*(iter++)) * -11;
	    value += (*(iter++)) * -10;
	    value += (*(iter++)) *  -9;
	    value += (*(iter++)) *  -7;
	    value += (*(iter++)) *  -6;
	    value += (*(iter++)) *  -5;
	    value += (*(iter++)) *  -4;
	    value += (*(iter++)) *  -3;
	    value += (*(iter++)) *  -2;
	    value += (*(iter++)) *  -1;
	    ++iter;
	    value += (*(iter++)) *   1;
	    value += (*(iter++)) *   2;
	    value += (*(iter++)) *   3;
	    value += (*(iter++)) *   4;
	    value += (*(iter++)) *   5;
	    value += (*(iter++)) *   6;
	    value += (*(iter++)) *   7;
	    value += (*(iter++)) *   9;
	    value += (*(iter++)) *  10;
	    value += (*(iter++)) *  11;
	  }
	};

	
	  template<> struct sg_norm<2, 21, 2> {double value; sg_norm(double step) {value = 1/(33649*step*step);}};
	template<typename VALUE> struct sg_compute<VALUE, 2, 21, 2> {
	  VALUE value;
	  template<typename ITER>
	  sg_compute(ITER iter) {
	    value  = (*(iter++)) *   190;
	    value += (*(iter++)) *   133;
	    value += (*(iter++)) *    82;
	    value += (*(iter++)) *    37;
	    value += (*(iter++)) *    -2;
	    value += (*(iter++)) *   -35;
	    value += (*(iter++)) *   -62;
	    value += (*(iter++)) *   -83;
	    value += (*(iter++)) *   -98;
	    value += (*(iter++)) *  -107;
	    value += (*(iter++)) *  -110;
	    value += (*(iter++)) *  -107;
	    value += (*(iter++)) *   -98;
	    value += (*(iter++)) *   -83;
	    value += (*(iter++)) *   -62;
	    value += (*(iter++)) *   -35;
	    value += (*(iter++)) *    -2;
	    value += (*(iter++)) *    37;
	    value += (*(iter++)) *    82;
	    value += (*(iter++)) *   133;
	    value += (*(iter++)) *   190;
	  }
	};
	  
	//
	// win = 21, degree = 3.
	//
	
	template<> struct sg_norm<0, 21, 3> {double value; sg_norm(double step) {value = sg_norm<0, 21, 2>(step).value;}};
	template<typename VALUE> struct sg_compute<VALUE, 0, 21, 3> {
	  VALUE value;
	  template<typename ITER>
	  sg_compute(ITER iter) {
	    value = sg_compute<VALUE, 0, 21, 2>(iter).value;
	  }
	};
	
	template<> struct sg_norm<1, 21, 3> {double value; sg_norm(double step) {value = 1/(3634092*step);}};
	template<typename VALUE> struct sg_compute<VALUE, 1, 21, 3> {
	  VALUE value;
	  template<typename ITER>
	  sg_compute(ITER iter) {
	    value  = (*(iter++)) *   84075;
	    value += (*(iter++)) *   10032;
	    value += (*(iter++)) *  -43284;
	    value += (*(iter++)) *  -78176;
	    value += (*(iter++)) *  -96947;
	    value += (*(iter++)) * -101900;
	    value += (*(iter++)) *  -95338;
	    value += (*(iter++)) *  -79654;
	    value += (*(iter++)) *  -56881;
	    value += (*(iter++)) *  -29592;
	    ++iter;
	    value += (*(iter++)) *   29592;
	    value += (*(iter++)) *   56881;
	    value += (*(iter++)) *   79654;
	    value += (*(iter++)) *   95338;
	    value += (*(iter++)) *  101900;
	    value += (*(iter++)) *   96947;
	    value += (*(iter++)) *   78176;
	    value += (*(iter++)) *   43284;
	    value += (*(iter++)) *  -10032;
	    value += (*(iter++)) *  -84075;
	  }
	};
	
	template<> struct sg_norm<2, 21, 3> {double value; sg_norm(double step) {value = sg_norm<2, 21, 2>(step).value;}};
	template<typename VALUE> struct sg_compute<VALUE, 2, 21, 3> {
	  VALUE value;
	  template<typename ITER>
	  sg_compute(ITER iter) {
	    value = sg_compute<VALUE, 2, 21, 2>(iter).value;
	  }
	};

	
	
	//
	// Estimator
	//
	
	template<typename VALUE, unsigned int ORDER, unsigned int WINDOW_SIZE, unsigned int DEGREE>
	class estimator : public savitzky_golay::Estimator<VALUE, ORDER, WINDOW_SIZE> {
	private:
	  std::array<double, ORDER+1> hcoef;

	  template<unsigned int O>
	  void __set_timestep(double step) {
	    if constexpr (O == 0) 
	      hcoef[0] = sg_norm<0, WINDOW_SIZE, DEGREE>(step).value;
	    else {
	      hcoef[O] = sg_norm<O, WINDOW_SIZE, DEGREE>(step).value;
	      __set_timestep<O-1>(step);
	    }
	  }
	  
	public:
	
	  using savitzky_golay::Estimator<VALUE, ORDER, WINDOW_SIZE>::Estimator;

	  /**
	   * This sets the time step.
	   */
	  void set_timestep(double step) {
	    __set_timestep<ORDER>(step);
	  }
	  
	  /**
	   * Returns the estimate of the Oth derivative. It is an optional since it may not be available.
	   */
	  template<unsigned int O>
	  const std::optional<VALUE>& get() const {
	    static_assert(O <= ORDER, "vq3::utils::savitzky_golay::get<O> : O is too high.");
	    auto& v = this->value[O];
	    if(!v && this->window.size() == WINDOW_SIZE)
	      v = sg_compute<VALUE, O, WINDOW_SIZE, DEGREE>(this->window.begin()).value * hcoef[O];
	    return v;
	  }
	};
      }
    }


    
  }
}
