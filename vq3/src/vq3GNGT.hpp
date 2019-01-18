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

	/**
	 * This is the default evolution algorithm for the GNG-T algorithm.
	 */
	class Evolution {

	public:
  
	  double T;                   //!< The target.
	  double density;             //!< The N value. 
	  double margin_above = .30;  //!< If bmu_error > T*density*(1+margin_above), the vertex is cloned.
	  double margin_below = .15;  //!< If bmu_error < T*density*(1-margin_below), the vertex is killed.
	  double topo_ratio   = .30;  //!< Add/remove topo_ratio*nb_vertices_out_of_margin.

	private:
	  // Collect the (vertex_idx, bmu_error) whose BMU error is too high and too_low.
  
	  std::vector<std::pair<std::size_t, double>> above;
	  std::vector<std::pair<std::size_t, double>> below;

	public:

	  	  
	  Evolution()                            = default;
	  Evolution(const Evolution&)            = default;
	  Evolution& operator=(const Evolution&) = default;

	  // GNG-T calls this method when it considers to perform a
	  // modification of the number of vertices.
	  template<typename TABLE, typename BMU_RESULT, typename CLONE_PROTOTYPE>
	  bool operator()(TABLE&                 topology,
			  const BMU_RESULT&      bmu_epoch_result,
			  const CLONE_PROTOTYPE& clone_prototype) {
	    double topology_changed = false;

	    above.clear();
	    auto above_out = std::back_inserter(above);

	    below.clear();
	    auto below_out = std::back_inserter(below);
    
	    double NT = density*T;

	    double above_bound = NT*(1+margin_above);
	    double below_bound = NT*(1-margin_below);

	    // Let us consider all the errors.
	    std::size_t vertex_idx = 0;
	    for(auto& res : bmu_epoch_result) {
	      if(res.vq3_bmu_accum.nb != 0) { // We consider only the vertices which have been a bmu at least once.
		double error  = res.vq3_bmu_accum.value;
		if     (error > above_bound) *(above_out++) = {vertex_idx, error};
		else if(error < below_bound) *(below_out++) = {vertex_idx, error};
	      }
	      else {
		// We remove from the graph the vertices that have never been selected as a BMU.
		topology(vertex_idx)->kill();
		topology_changed = true;
	      }
	      ++vertex_idx;
	    }

	    // We sort out of bounds vertices (small error first for below, big error first for above).
	    std::sort(above.begin(), above.end(),
		      [](const std::pair<std::size_t, double>& p1, const std::pair<std::size_t, double>& p2) {
			return p1.second > p2.second;});
	    std::sort(below.begin(), below.end(),
		      [](const std::pair<std::size_t, double>& p1, const std::pair<std::size_t, double>& p2) {
			return p1.second < p2.second;});


	    // We clone a topo_ration fraction of the above vertices.
	    auto above_end = above.begin();
	    if(above_end != above.end()) {// if not empty
	      std::advance(above_end, std::max((std::size_t)(above.size()*topo_ratio), (std::size_t)1));
	      for(auto it = above.begin(); it != above_end; ++it) topology.g += clone_prototype((*(topology(it->first)))().vq3_value);
	    }

	    // we delete a topo_ration fraction of the below vertices.
	    auto below_end = below.begin();
	    if(below_end != below.end()) {// if not empty
	      std::advance(below_end, std::max((std::size_t)(below.size()*topo_ratio), (std::size_t)1));
	      for(auto it = below.begin(); it != below_end; ++it) topology(it->first)->kill();
	    }
      
	    topology_changed = (above.begin() != above_end) || (below.begin() != below_end);

	    // We tell GNG-T if any change in the vertices (add/remove) has occurred.
	    return topology_changed;
	  }
	};

	
	inline Evolution evolution() {
	  return Evolution();
	}
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

	double alpha              = 0.05; //!< The learning rate for the on-line SOM update performed at the beginning of an epoch.
	double samples_per_vertex = 10;   //!< The number of samples used for the on-line SOM update is at most samples_per_vertex * nb_vertices.
	
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
	 * @param nb_post_evolution_steps After the topology update, several batch WTM are performed.
	 * @param evolution Modifies the graph. See vq3::algo::gngt::by_default::evolution for an example.
	 */
	template<typename ITER, typename PROTOTYPE_OF_VERTEX, typename SAMPLE_OF, typename EVOLUTION, typename CLONE_PROTOTYPE, typename DISTANCE, typename VALUE_OF_EDGE_DISTANCE>
	void process(unsigned int nb_threads,
		     const ITER& begin, const ITER& end,
		     const SAMPLE_OF& sample_of,
		     const PROTOTYPE_OF_VERTEX& ref_prototype_of_vertex,
		     const CLONE_PROTOTYPE& clone_prototype,
		     const DISTANCE& distance,
		     const VALUE_OF_EDGE_DISTANCE& voed, unsigned int max_dist, double min_val,
		     const std::optional<unsigned int>& average_radius,
		     unsigned int nb_post_evolution_steps,
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
	  decltype(std::distance(begin, end)) nb_remaining_samples = (int)(nb_vertices*samples_per_vertex);
	  auto sample_it = begin;
	  while(nb_remaining_samples != 0) {
	    auto to_do = std::min(std::distance(sample_it, end), nb_remaining_samples);
	    auto sample_end = sample_it + to_do;
	    while(sample_it != sample_end) vq3::online::wtm::learn(table, distance, sample_of(*(sample_it++)), alpha);
	    nb_remaining_samples -= to_do;
	    if(sample_it == end) sample_it = begin;
	  }
	    

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

	  // We performe more batch updates.
	  if(nb_post_evolution_steps > 1) {
	    --nb_post_evolution_steps;
	    for(unsigned int i=0; i < nb_post_evolution_steps; ++i)
	      wtm.template process<epoch_wtm>(nb_threads, begin, end, sample_of, ref_prototype_of_vertex, distance);
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
