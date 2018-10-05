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

#include <vq3Epoch.hpp>
#include <vq3Graph.hpp>
#include <vq3Utils.hpp>
#include <vq3Topology.hpp>

#include <iterator>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <iomanip>

namespace vq3 {
  namespace algo {
    /**
     * This is the Linde-Buzo-Gray algorithm.
     * @param rd the random engine.
     * @param nb_threads the number of threads used for computing LBG.
     * @param g the graph, it is cleared, and then k vertices are added.
     * @param k the number of vertices required.
     * @param begin, end The samples
     * @param sample_of The samples are obtained from sample_of(*it).
     * @param prototype_of prototype_of(vertex_value) is a **reference** to the prototype.
     * @param distance distance(prototype, sample) is used internally.
     * @param nearly nearly(p) produces a prototype that is slightly different from p.
     * @param check The result of check(previous_vertex_value, current_vertex_value) should be false for all vertex once convergence is considered to be reached.
     * @param verbose Toggles verbosity.
     */
    template<typename PROTOTYPE, typename RANDOM_ENGINE, typename GRAPH, typename ITERATOR, typename SAMPLE_OF, typename PROTOTYPE_OF_VERTEX_VALUE, typename DISTANCE, typename NEARLY, typename CHECK>
    void lbg(RANDOM_ENGINE& rd,
	     unsigned int nb_threads, GRAPH& g, unsigned int k,
	     const ITERATOR& begin, const ITERATOR& end, const SAMPLE_OF& sample_of,
	     const PROTOTYPE_OF_VERTEX_VALUE& prototype_of,
	     const DISTANCE& distance,
	     const NEARLY& nearly,
	     const CHECK& check,
	     bool verbose) {
      if(begin == end)
	throw std::runtime_error("vq3::algo::lbg : empty dataset.");
      if(k == 0)
	throw std::runtime_error("vq3::algo::lbg : k==0.");
      auto d = (unsigned int)(std::distance(begin, end));
      if(d < k) {
	std::ostringstream msg;
	msg << "vq3::algo::lbg : not enough data : only " << d << " samples while k = " << k << ".";
	throw std::runtime_error(msg.str());
      }
      
      g.foreach_vertex([](const typename GRAPH::ref_vertex& ref_v){ref_v->kill();});

      auto table = vq3::topo::table(g);
      auto wta   = vq3::epoch::wta::processor(table);
      
      g += sample_of(*begin);
      unsigned int nb_nodes = 1;
      table();

      if(verbose)
	std::cout << std::endl
		  << "Starting Linde-Buzo-Gray with K =" << std::setw(4) << k << "." << std::endl
		  << "--------------------------------------" << std::endl;

      wta.template processs<typename vq3::epoch::data::delta<typename vq3::epoch::data::wta<typename epoch::data::none<PROTOTYPE,
														       typename GRAPH::vertex_value_type,
														       PROTOTYPE> > > >(nb_threads,
																	begin, end,
																	sample_of,
																	prototype_of,
																	distance);
      
      while(nb_nodes < k) {
	unsigned int new_nb_nodes = std::min(k, 2*nb_nodes);
	unsigned int added        = new_nb_nodes - nb_nodes;
	nb_nodes                  = new_nb_nodes;
	if(verbose)
	  std::cout << "  Adding "<< std::setw(4) << added << " nodes => " << std::setw(4) << nb_nodes << std::endl;

	std::vector<typename GRAPH::ref_vertex> V;
	utils::collect_vertices(g, std::back_inserter(V));
	std::shuffle(V.begin(), V.end(), rd);
	auto vend = V.begin() + added;
	for(auto it = V.begin(); it != vend; ++it)
	  g += nearly(prototype_of((*(*it))()));
	
	table();

	bool stop = false;
	while(!stop) {
	  auto res = wta.template processs<typename vq3::epoch::data::delta<typename vq3::epoch::data::wta<typename epoch::data::none<PROTOTYPE,
																      typename GRAPH::vertex_value_type,
																      PROTOTYPE> > > >(nb_threads,
																		       begin, end,
																		       sample_of,
																		       prototype_of,
																		       distance);
	  stop = true;
	  for(auto& d : res)
	    if(check(d.vq3_previous_prototype, d.vq3_current_prototype)) {
	      stop = false;
	      break;
	    }
	}
      }

      if(verbose)
	std::cout << "Done." << std::endl
		  << std::endl;
    }
  }
}
