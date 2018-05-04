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


#include <list>
#include <utility>
#include <deque>
#include <vector>
#include <map>
#include <iterator>
#include <stdexcept>


#include <vq3Graph.hpp>
#include <vq3Utils.hpp>

namespace vq3 {

  namespace topo {


    /**
     * This structure stores neihgborhood information.
     */
    template<typename INDEX, typename DIST_VALUE>
    struct Info {
      DIST_VALUE value;
      INDEX      index;

      Info(const DIST_VALUE& v, const INDEX& r) : value(v), index(r) {}
      Info()                       = default;
      Info(const Info&)            = default;
      Info(Info&&)                 = default;
      Info& operator=(const Info&) = default;
      Info& operator=(Info&&)      = default;
    };
    
    /**
     * The graph vertices must be tag-decorated. WARNING ! Clear the tags (false) before calling that function.
     * @param vertex_table A vertex table (see topo::vertices) of the current vertices of the graph.
     * @param ref_v the origin vertex.
     * @param voed A function providing a value (double >= 0) according to the number of edges (unsigned int) separating a vertex in the neighborhood from the central vertex.
     * @param max_dist The maximal distance considered. 0 means "no limit".
     * @param min_val if voed(dist) < min_val, the node is not included in the neighborhood.
     * @return The list of (value, idx) pairs corresponding to the neighborhood. idx is the index of the vertex in a vertices structure. The origin vertex index is in the list (at first position).
     */
    template<typename VERTICES, typename VALUE_OF_EDGE_DISTANCE>
    auto edge_based_neighborhood(const VERTICES& vertex_table, const typename VERTICES::ref_vertex_type& ref_v, const VALUE_OF_EDGE_DISTANCE& voed, unsigned int max_dist, double min_val) {
      std::list<Info<typename VERTICES::index_type, decltype(voed(0))> > res;
      std::deque<std::pair<unsigned int, typename VERTICES::ref_vertex_type> > to_do;

      auto res_out = std::back_inserter(res);
      auto job_out = std::back_inserter(to_do);

      *(res_out++) = {voed(0), vertex_table(ref_v)};
      (*ref_v)().vq3_tag = true;
      ref_v->foreach_edge([&job_out, &ref_v](const typename VERTICES::ref_vertex_type::element_type::ref_edge_type ref_e) {
	  auto extr = ref_e->extremities();
	  if(invalid_extremities(extr)) {ref_e->kill(); return;}
	  auto& other = other_extremity(extr, ref_v);
	  (*other)().vq3_tag = true;
	  *(job_out++) = {1, other};
	});

      while(!(to_do.empty())) {
	auto d_v = to_do.front();
	to_do.pop_front();
	double val = voed(d_v.first);
	if(val > min_val) {
	  *(res_out++) = {val, vertex_table(d_v.second)};
	  if(d_v.first != max_dist)
	    d_v.second->foreach_edge([&job_out, &v = d_v.second, dist = d_v.first + 1](const typename VERTICES::ref_vertex_type::element_type::ref_edge_type ref_e) {
		auto extr = ref_e->extremities();
		if(invalid_extremities(extr)) {ref_e->kill(); return;}
		auto& other = other_extremity(extr, v);
		auto& tag = (*other)().vq3_tag;
		if(!tag) {
		  tag = true;
		  *(job_out++) = {dist, other};
		}
	      });
	}
      }
      return res;
    }


    
  }
  
}
