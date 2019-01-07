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

	
	// inline Evolution evolution() {
	//   return Evolution();
	// }
      }


      
      template<typename SAMPLE, typename TABLE>
      class Processor {
      public:
	
	using topology_table_type = TABLE;

	using graph_type = typename topology_table_type::graph_type;
	using ref_vertex = typename graph_type::ref_vertex;
	using ref_edge   = typename graph_type::ref_edge;
	using vertex     = typename graph_type::vertex_value_type;
	using edge       = typename graph_type::edge_value_type;

	using prototype_type = typename vertex::decorated_type;
	using sample_type    = SAMPLE;
      
	
	using epoch_bmu = vq3::epoch::data::bmu<vq3::epoch::data::none<sample_type, vertex, prototype_type> >;
	using epoch_wta = vq3::epoch::data::wta<vq3::epoch::data::none<sample_type, vertex, prototype_type> >;
	using epoch_wtm = vq3::epoch::data::wtm<vq3::epoch::data::none<sample_type, vertex, prototype_type> >;
	       

      private:

	topology_table_type& table;
	topology_table_type  avg_table;
	std::optional<unsigned int> old_avg_radius;
	
      public:

	double alpha             = 0.05; //!< The learning rate for the on-line SOM update performed at the beginning of an epoch.
	double sample_per_vertex = 10;   //!< The number of samples used for the on-line SOM update is at most sample_per_vertex * nb_vertices.
	
	vq3::epoch::wta::Processor<topology_table_type> wta;
	vq3::epoch::wtm::Processor<topology_table_type> wtm;
	vq3::epoch::chl::Processor<graph_type>          chl;

	std::vector<epoch_bmu> bmu_results;

	Processor(topology_table_type& table)
	  : table(table), avg_table(table.g),
	    wta(table), wtm(table), chl(table.g) {}
	
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
	 * @param average_radius The error is averaged from the vertices around. The neihgborhood is made of the vertices with an edge distance d <= average_radius. Provide 0 for averaging with all vertices in the connected component. It is an optional value, so if the value is unset, no averaging is performed.
	 * @param evolution Modifies the graph. See vq3::algo::gngt::by_default::evolution for an example.
	 */
	template<typename ITER, typename PROTOTYPE_OF_VERTEX, typename SAMPLE_OF, typename EVOLUTION, typename CLONE_PROTOTYPE, typename DISTANCE, typename VALUE_OF_EDGE_DISTANCE>
	void epoch(unsigned int nb_threads,
		   const ITER& begin, const ITER& end,
		   const SAMPLE_OF& sample_of,
		   const PROTOTYPE_OF_VERTEX& ref_prototype_of_vertex,
		   const CLONE_PROTOTYPE& clone_prototype,
		   const DISTANCE& distance,
		   const VALUE_OF_EDGE_DISTANCE& voed, unsigned int max_dist, double min_val,
		   const std::optional<unsigned int>& average_radius,
		   EVOLUTION& evolution) {
 	  if(begin == end) {
	    // No samples. We kill all nodes and keep the two topological tables updated.
	    table.g.foreach_vertex([](const ref_vertex& ref_v) {ref_v->kill();});
	    table(voed, max_dist, min_val);
	    if(average_radius)
	      avg_table([](unsigned int) {return 1;}, *average_radius, 0);
	    return;
	  }

	  auto nb_vertices = table.g.nb_vertices();

	  if(nb_vertices == 0) {
	    // Empty graph. We create one vertex, keep the two
	    // topological tables updated, and do one wta pass.
	    table.g += sample_of(*begin);
	    table(voed, max_dist, min_val);
	    if(average_radius)
	      avg_table([](unsigned int) {return 1;}, *average_radius, 0);
	    wta.template process<epoch_wta>(nb_threads, begin, end, sample_of, ref_prototype_of_vertex, distance);
	    return;
	  }

	  if(average_radius != old_avg_radius) {
	    // If the average radius is changed from last time, the
	    // related topological table has to be recomputed.
	    old_avg_radius = average_radius;
	    if(average_radius)
	      avg_table([](unsigned int) {return 1;}, *average_radius, 0);
	  }

	  // This is an online SOM update in order to quickly move the
	  // graph so that it fits the samples, while topological
	  // evolution is freezed.
	  auto sample_it = begin;
	  auto sample_end = begin + std::min(std::distance(begin, end), (decltype(std::distance(begin, end)))(nb_vertices*sample_per_vertex));
	  while(sample_it != sample_end)
	    vq3::online::wtm::learn(table, distance, sample_of(*(sample_it++)), alpha);

	  // We compute the error cumulated values for all the vertices.
	  bmu_results = wta.template process<epoch_bmu>(nb_threads, begin, end, sample_of, ref_prototype_of_vertex, distance);

	  // bmu_results_ptr refers the collection of values actually
	  // used in the evolution algorithms. It can be directly
	  // bmu_results or another set of values obtained by spatial
	  // averaging, according to average_radius argument.
	  std::vector<epoch_bmu>* bmu_results_ptr = &bmu_results;
	  
	  std::vector<epoch_bmu> avg_bmu_results;
	  if(average_radius) {
	    // If average radius is set, we average the bmu results.
	    avg_bmu_results.resize(bmu_results.size());
	    typename topology_table_type::index_type idx = 0;
	    for(auto& data : avg_bmu_results) {
	      auto& neighborhood = avg_table[idx++];
	      for(auto& info : neighborhood) 
		if(auto& acc = bmu_results[info.index].vq3_bmu_accum; acc.nb > 0)
		  data.vq3_bmu_accum += acc.value;
	    }
	    
	    for(auto& data : avg_bmu_results)
	      if(data.vq3_bmu_accum.nb > 0)
		data.vq3_bmu_accum = data.vq3_bmu_accum.average();
	    
	    bmu_results_ptr = &avg_bmu_results;
	  }

	  // We call the user evolution method
	  if(evolution(table, *bmu_results_ptr, clone_prototype)) {
	    // If the tolopogy has changed, we update both tables.
	    table(voed, max_dist, min_val);
	    if(average_radius)
	      avg_table([](unsigned int) {return 1;}, *average_radius, 0);
	  }

	  // We adjust the vertices position with a single batch update.
	  wtm.template process<epoch_wtm>(nb_threads, begin, end, sample_of, ref_prototype_of_vertex, distance);

	  // We update the edges thanks to Competitive Hebbian learning.
	  if(chl.process(nb_threads, begin, end, sample_of, ref_prototype_of_vertex, distance, edge())) {
	    // If the tolopogy has changed, we update both tables.
	    table(voed, max_dist, min_val);
	    if(average_radius)
	      avg_table([](unsigned int) {return 1;}, *average_radius, 0);
	  }
	}
	
      };

      template<typename SAMPLE, typename TABLE>
      Processor<SAMPLE, TABLE> processor(TABLE& table) {
	return Processor<SAMPLE, TABLE>(table);
      }
      
      
    }
  }
}
