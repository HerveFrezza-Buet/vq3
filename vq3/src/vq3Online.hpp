#pragma once

#include <vq3Utils.hpp>

namespace vq3 {
  namespace online {
    namespace wtm {
      template<typename TABLE, typename SAMPLE, typename DISTANCE>
      auto learn(TABLE& table, const DISTANCE& dist, const SAMPLE& xi, double alpha) {
	auto ref_v = vq3::utils::closest(table.g, xi , [&dist](const typename TABLE::graph_type::vertex_type::value_type& v, const SAMPLE& p) {return dist(v.vq3_value, p);});
	auto& n = table[ref_v];
	for(auto& info : n) {
	  auto& w = (*(table(info.index)))().vq3_value;
	  w += alpha*info.value*(xi-w);
	}

	return ref_v;
      }
    }
    
    namespace wta {
      template<typename GRAPH, typename SAMPLE, typename DISTANCE>
      auto learn(GRAPH& g, const DISTANCE& dist, const SAMPLE& xi, double alpha) {
	auto ref_v = vq3::utils::closest(g, xi , [&dist ](const typename GRAPH::vertex_type::value_type& v, const SAMPLE& p) {return dist(v.vq3_value, p);});
	auto& w = (*(ref_v))().vq3_value;
	w += alpha*(xi-w);
	
	return ref_v;
      }
    }
  }
}
