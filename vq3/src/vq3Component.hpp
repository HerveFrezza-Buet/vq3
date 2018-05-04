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
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <limits>
#include <vq3Utils.hpp>

namespace vq3 {
  namespace connected_components {
    
    template<typename GRAPH>
    class component;
    
    template<typename GRAPH>
    std::list<vq3::connected_components::component<GRAPH>> make(GRAPH& g);
    
    /**
     * This is a connected component.
     */
    template<typename GRAPH>
    class component {
    private:

      typename GRAPH::weak_vertices V;
    
      friend std::list<vq3::connected_components::component<GRAPH>> make<GRAPH>(GRAPH& g);
    
    public:

      using ref_vertex = typename GRAPH::ref_vertex;
      
      component()                            = default;
      component(const component&)            = default;
      component(component&&)                 = default;
      component& operator=(const component&) = default;
      component& operator=(component&&)      = default;
    
      template<typename VERTEX_FUN>
      void foreach_vertex(const VERTEX_FUN& fun) {
	auto it =  V.begin();
	while(it != V.end()) {
	  auto v = (*it).lock();
	  if(v == nullptr || v->is_killed())
	    it = V.erase(it);
	  else {
	    fun(v);
	    if(v->is_killed())
	      it = V.erase(it);
	    else
	      ++it;
	  }
	}
      }

    };

    
    /**
     * Computes the connected component in a list. The vertex value
     * must be decorated with vq3::decorator::efficiency and
     * vq3::decorator::tagged. Edge values only need
     * vq3::decorator::efficiency.
     */
    template<typename GRAPH>
    std::list<vq3::connected_components::component<GRAPH>> make(GRAPH& g) {
      std::list<vq3::connected_components::component<GRAPH>> components;
      std::vector<typename GRAPH::ref_vertex> vertices;
      std::stack<typename GRAPH::ref_vertex>  untagged;

      vq3::utils::clear_vertex_tags(g,true);
      vq3::utils::collect_vertices(g, std::back_inserter(vertices));
      
      for(auto ref_v : vertices)
	if((*ref_v)().vq3_tag && (*ref_v)().vq3_efficient) {
	  component<GRAPH> comp;
	  (*ref_v)().vq3_tag = false;
	  untagged.push(ref_v);
	  while(!untagged.empty()) {
	    auto ref_vv = untagged.top();
	    untagged.pop();
	    comp.V.push_back(ref_vv);
	    vq3::utils::foreach_efficient_edge(ref_vv,[&untagged](const typename GRAPH::ref_edge& ref_e) {
		auto extr_pair = ref_e->extremities();           
		if(vq3::invalid_extremities(extr_pair)) {
		  ref_e->kill();
		  return;
		}
		auto& v1 = (*(extr_pair.first))();
		auto& v2 = (*(extr_pair.second))();
		if(v1.vq3_tag && v1.vq3_efficient) {
		  v1.vq3_tag = false;
		  untagged.push(extr_pair.first);
		}
		if(v2.vq3_tag && v2.vq3_efficient) {
		  v2.vq3_tag = false;
		  untagged.push(extr_pair.second);
		}
	      });
	  }
	  components.push_back(std::move(comp));
	}
      return components;
    }
  }

  namespace labelling {

    /**
     * This labels connected components without previous labelling
     * consideration.
     */
    template<typename COMPONENT_ITERATOR>
    void basic(const COMPONENT_ITERATOR& begin,
	       const COMPONENT_ITERATOR& end) {
      unsigned int label = 1;
      for(auto it = begin; it != end; ++it)
	vq3::utils::clear_vertex_labels(*it, label++);
    }

    /**
     * This tries to preserve the previous labelling.
     */
    template<typename COMPONENT_ITERATOR>
    void conservative(const COMPONENT_ITERATOR& begin,
		      const COMPONENT_ITERATOR& end) {
      
      // label_count[nb] = (label, comp_it) 
      std::multimap<unsigned int, std::pair<unsigned int, COMPONENT_ITERATOR>, std::greater<unsigned int> > label_count;
      // labels[comp_id] = label
      std::vector<unsigned int> labels(std::distance(begin,end), 0);

      // Let us count....
      for(auto it = begin; it != end; ++it) {
	std::map<unsigned int, unsigned int> l; // l[label] = nb
	it->foreach_vertex([&l](typename COMPONENT_ITERATOR::value_type::ref_vertex& ref_v) {
	    auto lbl = (*ref_v)().vq3_label;
	    auto lit = l.find(lbl);
	    if(lit == l.end())
	      l[lbl] = 1;
	    else
	      ++(lit->second);
	  });
	for(auto& ln : l)
	  label_count.insert({ln.second, {ln.first, it}});
      }


      // ... and then label the components.
      std::set<unsigned int> used_labels;
      unsigned int min_label = std::numeric_limits<unsigned int>::max();
      unsigned int max_label = 0;

      
      for(auto& count : label_count) {
	auto lbl = count.second.first;
	auto it  = count.second.second;
	auto idx = std::distance(begin, it);
	
	unsigned int the_label = 0;
	if(labels[idx] == 0) {// The component has no label yet
	  if(lbl != 0 && used_labels.find(lbl) == used_labels.end()) {
	    the_label = lbl;  // The label is free
	    used_labels.insert(the_label);
	    min_label = std::min(min_label, the_label);
	    max_label = std::max(max_label, the_label);
	  }
	  else 
	    the_label = std::numeric_limits<unsigned int>::max(); // request for new label.
	  labels[idx] = the_label;
	}
      }
      auto lit = labels.begin();
      for(auto it = begin; it != end; ++it, ++lit) {
	auto the_label = *lit;
	if(the_label == std::numeric_limits<unsigned int>::max()) {
	  if(min_label > 1) 
	    the_label = --min_label;
	  else
	    the_label = ++max_label;
	}
	vq3::utils::clear_vertex_labels(*it, the_label);
      }

      

    }

    template<typename GRAPH>
    void edges_from_vertices(GRAPH& g) {
      vq3::utils::clear_edge_labels(g, 0);
      g.foreach_edge([](typename GRAPH::ref_edge ref_e) {
	  auto extr_pair = ref_e->extremities();           
	  if(vq3::invalid_extremities(extr_pair)) {
	    ref_e->kill();
	    return;
	  }

	  auto l1 = (*(extr_pair.first))().vq3_label;
	  auto l2 = (*(extr_pair.first))().vq3_label;
	  if(l1 != l2)
	    (*ref_e)().vq3_label = 0;
	  else
	    (*ref_e)().vq3_label = l1;
	});
    }
	
  }
}
