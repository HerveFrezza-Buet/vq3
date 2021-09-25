#pragma once


#include <utility>
#include <functional>
#include <iostream>
#include <limits>

#include <vq3Utils.hpp>
#include <vq3ShortestPath.hpp>

namespace vq3 {
  namespace spec {

    struct GIT_Traits {
      /**
       * This is the value type (int is fake here).
       */
      using value_type = int;

      /**
       * This is the underlying graph type (int is fake here).
       */
      using graph_type = int;

      
      /**
       * This is the underlying graph (can be a reference).
       */
      graph_type g;

      /**
       * This function gives the distance between a vertex value and a value (int is fake here). 
       */
      double D(int vertex_value, const value_type& v) const;

      /**
       * This function gives the distance between two values. 
       */
      double d(const value_type& v1, const value_type& v2) const;

      /**
       * This function interpolates, returning lambda*a + (1-lambda)*b, lambda in [0,1].
       */
      value_type interpolate(const value_type& a, const value_type& b, double lambda) const;

      /**
       * This computes the shortest path from a to b. Int is fake here, it should be graph_type::ref_vertex. If a is nullptr, all distances to b are computed at each vertex.
       */ 
      void shortest_path(int a,  int b) const;
    };
  }
  
  namespace topology {
    namespace gi { // graph induced

      template<typename VALUE,
	       typename GRAPH,
	       typename VAL_DIST,
	       typename GRAPH_DIST,
	       typename INTERPOLATE,
	       typename SHORTEST_PATH>
      struct Traits {
	using value_type = VALUE;
	using graph_type = GRAPH;
	graph_type&                                                                              g;
	std::function<double (const value_type&, const value_type&)>                             d_func;
	std::function<double (const typename graph_type::vertex_value_type&, const value_type&)> D_func;
	std::function<value_type (const value_type&, const value_type&, double)>                 i_func;
	std::function<void (GRAPH&, typename GRAPH::ref_vertex, typename GRAPH::ref_vertex)>     sp_func;

	double d(const value_type& a, const value_type& b) const  {return d_func(a, b);}
	
	double D(const typename graph_type::vertex_value_type& a, const value_type& b) const  {return D_func(a, b);}
	
	value_type interpolate(const value_type& a, const value_type& b, double lambda) const {return i_func(a, b, lambda);}

	void shortest_path(typename GRAPH::ref_vertex a, typename GRAPH::ref_vertex b) const {
	  sp_func(g, a, b);
	}

	Traits(graph_type& g,
	       const VAL_DIST& d,
	       const GRAPH_DIST& D,
	       const INTERPOLATE& i,
	       const SHORTEST_PATH& sp)
	  : g(g), d_func(d), D_func(D), i_func(i), sp_func(sp) {}
	Traits()                         = delete;
	Traits(const Traits&)            = default;
	Traits& operator=(const Traits&) = delete;
      };

      template<typename VALUE,
	       typename GRAPH,
	       typename VAL_DIST,
	       typename GRAPH_DIST,
	       typename INTERPOLATE,
	       typename SHORTEST_PATH>
      auto traits(GRAPH& g, const VAL_DIST& d, const GRAPH_DIST& D, const INTERPOLATE& i, const SHORTEST_PATH& sp) {
	return Traits<VALUE, GRAPH, VAL_DIST, GRAPH_DIST, INTERPOLATE, SHORTEST_PATH>(g, d, D, i, sp);
      }

      /**
       * This is intended to be used in non executed code. Typically:
       * @code{.cpp} using traits = decltype(vq3::topology::gi::traits_val<sample, graph>(d2, D2, interpolate, shortest_path)); @endcode
       */
      template<typename VALUE,
	       typename GRAPH,
	       typename VAL_DIST,
	       typename GRAPH_DIST,
	       typename INTERPOLATE,
	       typename SHORTEST_PATH>
      auto traits_val(const VAL_DIST& d, const GRAPH_DIST& D, const INTERPOLATE& i, const SHORTEST_PATH& sp) {
	GRAPH g;
	return Traits<VALUE, GRAPH, VAL_DIST, GRAPH_DIST, INTERPOLATE, SHORTEST_PATH>(g, d, D, i, sp);
      }

      template<typename GIT_TRAITS> class Value;

      template<typename GIT_TRAITS> 
      typename GIT_TRAITS::value_type operator*(double alpha, const std::pair<const Value<GIT_TRAITS>&, const Value<GIT_TRAITS>&>& diff);
      
      template<typename GIT_TRAITS>
      std::ostream& operator<<(std::ostream& os, const Value<GIT_TRAITS>& val);
      
      template<typename GIT_TRAITS>
      class Value {
	
      public:

	using traits_type = GIT_TRAITS;
	
      private:
	
	typename traits_type::value_type value;
	
	friend typename traits_type::value_type operator*<>(double alpha, const std::pair<const Value<GIT_TRAITS>&, const Value<GIT_TRAITS>&>& diff);

	friend std::ostream& operator<<<>(std::ostream& os, const Value<GIT_TRAITS>& val);
	
	traits_type&                                         traits;
	mutable typename traits_type::graph_type::ref_vertex closest;
	mutable double                                       dist_to_closest;

	void update_closest() const {
	  if(!closest)
	    closest = vq3::utils::closest(traits.g, value,
					  [this](const typename traits_type::graph_type::vertex_value_type& a, const typename traits_type::value_type& b) {
					    return traits.D(a, b); 
					  },
					  dist_to_closest);
	}
	
      public:

	Value()                        = delete;
	Value(const Value&)            = default;
	Value(traits_type& traits, const typename traits_type::value_type& value) : value(value), traits(traits) {}
	Value& operator=(const Value&) = default;

	const typename traits_type::value_type& operator()() const {
	  return value;
	}

	/**
	 * The value is not linked to a vertex in the auxiliary graph
	 * anymore. The link will be made the next time it is
	 * required. This is usefull when the auxiliary graph has
	 * changed.
	 */
	auto reset() {closest = nullptr;}

	auto closest_vertex() const {
	  update_closest();
	  return closest;
	}

	double distance_to_closest_vertex() const {
	  update_closest();
	  return dist_to_closest;
	};
	
	Value& operator=(const typename traits_type::value_type& v) {
	  closest.reset();
	  value = v;
	  return *this;
	}

	// w += alpha*(xi - w) : alpha*(xi - w) returns the new value
	Value& operator+=(const typename traits_type::value_type& v) {
	  *this = v;
	  return *this;
	}
	

	// w += alpha*(xi - w) : (xi - w) sets the shortest path and returns (w, xi)
	std::pair<const Value&, const Value&> operator-(const Value& other) const {
	  traits.shortest_path(other.closest_vertex(), closest_vertex()); 
	  return {other, *this};
	}
      };

      
      template<typename GIT_TRAITS>
      auto value(GIT_TRAITS& traits, const typename GIT_TRAITS::value_type& value) {return Value<GIT_TRAITS>(traits, value);}

      template<typename GIT_TRAITS>
      std::ostream& operator<<(std::ostream& os, const Value<GIT_TRAITS>& v) {
	os << v.val << '[';
	if(v.closest) {
	  os << v.closest << ", " << v.dist_to_closest;
	}
	else
	  os << "not computed" << std::endl;
	os << ']';
	return os;	  
      }

      // Previous implementation...
      /*

      template<typename GIT_TRAITS>
      auto operator*(double alpha, const std::pair<const Value<GIT_TRAITS>&, const Value<GIT_TRAITS>&>& diff) -> typename GIT_TRAITS::value_type {
	auto& [w, xi] = diff; 
	
	auto l_graph  = (*(w.closest))().vq3_shortest_path.cost;
	if(l_graph == std::numeric_limits<double>::max())
	  return w.value;
	
	auto l1       = w.dist_to_closest;
	auto l2       = xi.dist_to_closest;
	double l      = l1 + l_graph + l2;
	auto alpha_l  = alpha * l;

	l = l1;
	if(alpha_l <= l) {
	  double lambda = 0;
	  if(l > 0) lambda = alpha_l/l;
	  return w.traits.interpolate(w.value, (*(w.closest))().vq3_value, lambda);
	}
	
	l += l_graph;
	if(alpha_l <= l) {
	  double lambda = 0;
	  if(l_graph > 0)
	    lambda = (alpha_l - l1)/l_graph;
	  if(auto val = vq3::path::travel(w.closest, xi.closest, lambda,
					  [&w](const typename GIT_TRAITS::graph_type::ref_vertex& a, const typename GIT_TRAITS::graph_type::ref_vertex& b, double lambda){
					    return w.traits.interpolate((*a)().vq3_value, (*b)().vq3_value, lambda);
					  }); val)
	    return *val;
	  else
	    return w.value; // Trying to reach a sample in another connected component has no effect.
	}
	
	double lambda = 0;
	if(l2 > 0)
	  lambda = (alpha_l - l)/l2;
	return w.traits.interpolate((*(xi.closest))().vq3_value, xi.value, lambda);
      }

      */


#define vq3DEBUG_GIT
      
      // w += alpha*(xi - w) : alpha, (w, xi) --> the interpolated value.
#define vq3_GIT_DLAMBDA .01
      template<typename GIT_TRAITS>
      auto operator*(double alpha, const std::pair<const Value<GIT_TRAITS>&, const Value<GIT_TRAITS>&>& diff) -> typename GIT_TRAITS::value_type {	
	auto& [w, xi] = diff;
	// Here, internal values (closest, dist_to_closest, ...) are
	// up to date, since operator- did the update.

#ifdef vq3DEBUG_GIT
	std::cout << std::endl
		  << "####" << std::endl
		  << std::endl
		  << "w = " << w.value << ", xi = " << xi.value << std::endl
		  << "w* = " << (*(w.closest))().vq3_value << ", xi* = " << (*(xi.closest))().vq3_value << std::endl;
#endif
	
	auto l_graph  = (*(w.closest))().vq3_shortest_path.cost;

	// If no path exists to reach xi, we do not modify w.
	if(l_graph == std::numeric_limits<double>::max())
	  return w.value;


	// If both xi and w have the same closest vertex, we move
	// directly.
	if(xi.closest == w.closest)
	  return w.traits.interpolate(w.value, xi.value, alpha);

	// The path as at least two vertices...

	// We have to determine the two anchors to the graph, one for
	// xi and one for w. Using xi.closest and w.closest as anchors
	// it naive, we well try to do better.

	// Searching for w and xi anchors.
	auto path_it  = path::begin(w.closest);
	auto path_end = path::end_by_vertex(w.closest);
	auto wA = *(path_it++);    // Immediate closest node
	auto wB = *(path_it++);    // Next node in the path toward xi.
	auto xA = wB;
	auto xB = wA;
	while(path_it != path_end) {
	  xB = xA;
	  xA = *(path_it++);
	}
	
#ifdef vq3DEBUG_GIT
	std::cout << "w  anchors : A = " << ((*wA)()).vq3_value << " --> B = " << ((*wB)()).vq3_value << std::endl
	 	  << "xi anchors : A = " << ((*xA)()).vq3_value << " --> B = " << ((*xB)()).vq3_value << std::endl;
#endif

	// Do we anchor w to wA or wB ?
	
#ifdef vq3DEBUG_GIT
	std::cout << std::endl
		  << "-- Anchoring w :" << std::endl;
#endif
	auto L = w.traits.interpolate((*wA)().vq3_value, (*wB)().vq3_value, vq3_GIT_DLAMBDA);
#ifdef vq3DEBUG_GIT
	std::cout << "   L    = " << L << std::endl
		  << "   |wL| = " << w.traits.d(w.value, L) << std::endl
		  << "   |wA| = " << w.dist_to_closest << std::endl;
#endif
	auto wK = wA;
	if(w.traits.d(w.value, L) < w.dist_to_closest)
	  wK = wB;
	
#ifdef vq3DEBUG_GIT
	if(wK == wA) std::cout << "   anchor = A" << std::endl;
	else         std::cout << "   anchor = B" << std::endl;
		       
#endif

	
	// Do we anchor xi to xA or xB ?
	
#ifdef vq3DEBUG_GIT
	std::cout << std::endl
		  << "-- Anchoring xi :" << std::endl;
#endif
	L = w.traits.interpolate((*xA)().vq3_value, (*xB)().vq3_value, vq3_GIT_DLAMBDA);
#ifdef vq3DEBUG_GIT
	std::cout << "   L     = " << L << std::endl
		  << "   |xiL| = " << w.traits.d(xi.value, L) << std::endl
		  << "   |xiA| = " << xi.dist_to_closest << std::endl;
#endif
	auto xK = xA;
	if(w.traits.d(xi.value, L) < xi.dist_to_closest)
	  xK = xB;
	
#ifdef vq3DEBUG_GIT
	if(xK == xA) std::cout << "   anchor = A" << std::endl;
	else         std::cout << "   anchor = B" << std::endl;
		       
#endif
	  
	  

	return w.value;
	
	auto l1       = w.dist_to_closest;
	auto l2       = xi.dist_to_closest;
	double l      = l1 + l_graph + l2;
	auto alpha_l  = alpha * l;

	l = l1;
	if(alpha_l <= l) {
	  double lambda = 0;
	  if(l > 0) lambda = alpha_l/l;
	  return w.traits.interpolate(w.value, (*(w.closest))().vq3_value, lambda);
	}
	
	l += l_graph;
	if(alpha_l <= l) {
	  double lambda = 0;
	  if(l_graph > 0)
	    lambda = (alpha_l - l1)/l_graph;
	  if(auto val = vq3::path::travel(w.closest, xi.closest, lambda,
					  [&w](const typename GIT_TRAITS::graph_type::ref_vertex& a, const typename GIT_TRAITS::graph_type::ref_vertex& b, double lambda){
					    return w.traits.interpolate((*a)().vq3_value, (*b)().vq3_value, lambda);
					  }); val)
	    return *val;
	  else
	    return w.value; // Trying to reach a sample in another connected component has no effect.
	}
	
	double lambda = 0;
	if(l2 > 0)
	  lambda = (alpha_l - l)/l2;
	return w.traits.interpolate((*(xi.closest))().vq3_value, xi.value, lambda);
      }

      template<typename VERTEX_VALUE,
	       typename GI_TRAITS>
      class Distance {

	mutable bool compute_paths = true;
	GI_TRAITS& traits;
	std::function<const Value<GI_TRAITS>& (const VERTEX_VALUE&)> prototype_of;

      public:

	template<typename PROTOTYPE_OF_VERTEX_VALUE>
	Distance(GI_TRAITS& traits, const PROTOTYPE_OF_VERTEX_VALUE& prototype_of)
	  : compute_paths(true), traits(traits), prototype_of(prototype_of) {}
	Distance(const Distance&)            = default;
	Distance& operator=(const Distance&) = delete;
	Distance()                           = delete;
	 

	void vq3_closest_init() const {compute_paths = true;}
	
	double operator()(const VERTEX_VALUE& vertex_value, const Value<GI_TRAITS>& xi) const {
	  if(compute_paths) {
	    compute_paths = false;
	    traits.shortest_path(nullptr, xi.closest_vertex());
	  }

	  auto& proto = prototype_of(vertex_value);
	  return proto.distance_to_closest_vertex()
	    +    xi.distance_to_closest_vertex()
	    +    (*(proto.closest_vertex()))().vq3_shortest_path.cost;
	}
      };

      
      template<typename VERTEX_VALUE, typename GI_TRAITS, typename PROTOTYPE_OF_VERTEX_VALUE>
      Distance<VERTEX_VALUE, GI_TRAITS> distance(GI_TRAITS& traits, const PROTOTYPE_OF_VERTEX_VALUE& prototype_of) {
	return Distance<VERTEX_VALUE, GI_TRAITS>(traits, prototype_of);
      }
    }
  }
}
