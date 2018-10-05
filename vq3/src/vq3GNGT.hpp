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

#include <functional>
#include <vector>
#include <iterator>
#include <algorithm>

#include <vq3Graph.hpp>
#include <vq3Epoch.hpp>
#include <vq3Utils.hpp>
#include <vq3Stats.hpp>
#include <vq3Online.hpp>
#include <vq3Topology.hpp>

namespace vq3 {
  namespace algo {
    namespace gngt {
      
      namespace by_default {
	struct Evolution {
	  double T          =   0;
	  double density    =   0;
	  double sigma_coef = 1.5;
	  
	  Evolution()                 = default;
	  Evolution(const Evolution&) = default;
	  
	  template<typename TABLE, typename BMU_RESULT, typename CLONE_PROTOTYPE>
	  void operator()(TABLE& table,
			  const BMU_RESULT& bmu_epoch_result,
			  const CLONE_PROTOTYPE& clone_prototype) {

	    std::vector<typename TABLE::graph_type::ref_vertex> above;
	    std::vector<typename TABLE::graph_type::ref_vertex> below;
	    auto out_above = std::back_inserter(above);
	    auto out_below = std::back_inserter(below);
    
	    double NT = density*T;
    
	    vq3::stats::MeanStd mean_std;
    
	    for(auto& res : bmu_epoch_result)
	      if(res.vq3_bmu_accum.nb != 0)
		mean_std  = res.vq3_bmu_accum.value;
	    double spatial_dmean = std::sqrt(mean_std.variance());
    
	    unsigned int idx = 0;
	    for(auto& res : bmu_epoch_result) {
	      auto& ref_v = table(idx++);
	      if(res.vq3_bmu_accum.nb == 0)
		ref_v->kill(); // We kill a vertex which has never won the competition.
	      else if(auto& omstd = (*ref_v)().vq3_online_mean_std; omstd) {
		auto [m, std] = omstd(); // We get the vertex distortion statistics.
	  
		double radius = (std + spatial_dmean)*sigma_coef;
	  
		if(NT < m - radius) // There are not enough nodes (cold color).
		  *(out_above++) = ref_v;
		else if(m + radius < NT)  // There are too many nodes (hot color).
		  *(out_below++) = ref_v;
	      }
	    }

	    if(above.size() > 0)
	      for(auto it = above.begin(); it != above.end(); ++it)
		table.g += clone_prototype((*(*it))().vq3_value);
    
	    if(below.size() > 0)
	      for(auto it = below.begin(); it != below.end(); ++it)
		(*it)->kill();

	  }
	};

	inline Evolution evolution() {
	  return Evolution();
	}
      }
      
      template<typename PROTOTYPE, typename SAMPLE, typename TABLE>
      class Processor {
      public:

	using topology_table_type = TABLE;
	using graph_type = typename topology_table_type::graph_type;
	using ref_vertex = typename graph_type::ref_vertex;
	using ref_edge   = typename graph_type::ref_edge;
	using vertex     = typename graph_type::vertex_value_type;
	using edge       = typename graph_type::edge_value_type;
      
	
	using epoch_bmu = vq3::epoch::data::online::bmu_mean_std<vq3::epoch::data::none<SAMPLE, vertex, PROTOTYPE> >;
	using epoch_wta = vq3::epoch::data::wta<vq3::epoch::data::none<SAMPLE, vertex, PROTOTYPE> >;
	       

      public:
	topology_table_type& table;
	vq3::epoch::wta::Processor<topology_table_type> wta;
	vq3::epoch::wta::Processor<topology_table_type> bmu;
	vq3::epoch::chl::Processor<topology_table_type> chl;

	Processor(topology_table_type& table)
	  : table(table), wta(table), bmu(table), chl(table.g) {
	}
	
	Processor()                            = delete;
	Processor(const Processor&)            = default;
	Processor(Processor&&)                 = default;
	Processor& operator=(const Processor&) = default;
	Processor& operator=(Processor&&)      = default;

	/**
	 * This updates the prototypes and the graph topology (i.e vertex and/or edge modification).
	 * @param begin, end The samples
	 * @param sample_of The samples are obtained from sample_of(*it).
	 * @param ref_prototype_of_vertex Returns a reference to the prototype from the vertex value.
	 * @param clone_prototype Computes a prototype value that is close to (*ref_v)().vq3_value.
	 * @param distance Compares the vertex value to a sample.
	 * @param evolution Modifies the graph. See vq3::algo::gngt::by_default::evolution for an example.
	 */
	template<typename ITER, typename PROTOTYPE_OF_VERTEX, typename SAMPLE_OF, typename EVOLUTION, typename CLONE_PROTOTYPE, typename DISTANCE>
	void process(unsigned int nb_threads,
		     const ITER& begin, const ITER& end,
		     const SAMPLE_OF& sample_of,
		     const PROTOTYPE_OF_VERTEX& ref_prototype_of_vertex,
		     const CLONE_PROTOTYPE& clone_prototype,
		     const DISTANCE& distance,
		     EVOLUTION& evolution) {
 	  if(begin == end) {
	    table.g.foreach_vertex([](const ref_vertex& ref_v) {ref_v->kill();});
	    table();
	    return;
	  }

	  if(table.nb_vertices() == 0) {
	    // empty graph, we create one vertex, and do one wta pass.
	    table.g += sample_of(*begin);
	    table();
	    wta.template process<epoch_wta>(nb_threads, begin, end, sample_of, ref_prototype_of_vertex, distance);
	    return;
	  }

	  auto bmu_results = bmu.template process<epoch_bmu>(nb_threads, begin, end, sample_of, ref_prototype_of_vertex, distance);
	  
	  evolution(table, bmu_results, clone_prototype);
	  table();
	  
	  chl.process(nb_threads, begin, end, sample_of, ref_prototype_of_vertex, distance, edge());
	}
	
      };

      template<typename PROTOTYPE, typename SAMPLE, typename TABLE>
      Processor<PROTOTYPE, SAMPLE, TABLE> processor(TABLE& table) {
	return Processor<PROTOTYPE, SAMPLE, TABLE>(table);
      }
      
    }
  }
}
