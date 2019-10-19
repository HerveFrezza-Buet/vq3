#pragma once

#include <vq3Utils.hpp>

namespace vq3 {
  namespace online {
    namespace wtm {
      /**
       * Applies WTM learning rule. If no vertex is closest to the prototype (it may append if the distance returns +infty), this function does nothing.
       * @return The closest vertex (or nullptr is it doesn't exist).
       */
      template<typename TABLE, typename SAMPLE, typename PROTOTYPE_OF_VERTEX_VALUE, typename DISTANCE>
      auto learn(TABLE& table,
		 const typename TABLE::neighborhood_key_type& topology_key, 
		 const PROTOTYPE_OF_VERTEX_VALUE& prototype_of,
		 const DISTANCE& distance,
		 const SAMPLE& xi, double alpha) {
	auto closest = vq3::utils::closest(table.g, xi, distance);
	if(closest) {
	  auto& n = table.neighborhood(closest, topology_key);
	  for(auto& info : n) {
	    auto& w = prototype_of((*(table(info.index)))());
	    w += (alpha*info.value)*(xi-w);
	  }
	}
	
	return closest;
      }
    }
    
    /**
     * Applies WTM learning rule. If no vertex is closest to the prototype (it may append if the distance returns +infty), this function does nothing.
     * @return The closest vertex (or nullptr is it doesn't exist).
     */
    namespace wta {
      template<typename GRAPH, typename SAMPLE, typename PROTOTYPE_OF_VERTEX_VALUE, typename DISTANCE>
      auto learn(GRAPH& g, 
		 const PROTOTYPE_OF_VERTEX_VALUE& prototype_of,
		 const DISTANCE& distance,
		 const SAMPLE& xi, double alpha) {
	auto closest = vq3::utils::closest(g, xi, distance);
	if(closest) {
	  auto& w = prototype_of((*(closest))());
	  w += alpha*(xi-w);
	}
	
	return closest;
      }
    }
  }
}
