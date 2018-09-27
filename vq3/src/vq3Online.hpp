#pragma once

#include <vq3Utils.hpp>

namespace vq3 {
  namespace online {
    namespace wtm {
      template<typename GRAPH, typename VERTICES, typename SAMPLE, typename DISTANCE, typename NTABLE>
      auto learn(GRAPH& g, VERTICES& vertices, const DISTANCE& dist, const NTABLE& neighborhood_table, const SAMPLE& xi, double alpha) {
	auto ref_v = vq3::utils::closest(g, xi , [&dist](const typename GRAPH::vertex_type::value_type& v, const SAMPLE& p) {return dist(v.vq3_value, p);});
	vq3::utils::clear_vertex_tags(g, false);
	auto& n = neighborhood_table.find(ref_v)->second;
	for(auto& info : n) {
	  auto& w = (*(vertices(info.index)))().vq3_value;
	  w += alpha*info.value*(xi-w);
	}
      }
    }
    
    namespace wta {
      template<typename GRAPH, typename VERTICES, typename SAMPLE, typename DISTANCE>
      auto learn(GRAPH& g, VERTICES& vertices, const DISTANCE& dist, const SAMPLE& xi, double alpha) {
	auto ref_v = vq3::utils::closest(g, xi , [&dist ](const typename GRAPH::vertex_type::value_type& v, const SAMPLE& p) {return dist(v.vq3_value, p);});
	auto& w = (*(ref_v))().vq3_value;
	w += alpha*(xi-w);
      }
    }
  }
}
