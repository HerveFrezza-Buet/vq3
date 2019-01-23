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
#include <functional>


#include <vq3Graph.hpp>
#include <vq3Utils.hpp>

namespace vq3 {

  namespace topology {

    /**
     * This builds a array of the vertex currently in the graph,
     * associating to each vertex an integer idf. Access to vertices
     * can be done from the index, and recipocally, the index of a
     * vertex can be retrieved (log(n) complexity).
     *
     * Table also hosts the computation of neighborhoods for each vertex.
     *
     * The graph vertices must be tag-decorated. 
     */
    template<typename GRAPH, typename NEIGHBORHOOD_KEY>
    class Table {
    public:
      using graph_type            = GRAPH;
      using index_type            = typename std::vector<typename GRAPH::ref_vertex>::size_type;
      using neighborhood_key_type = NEIGHBORHOOD_KEY;
      
      /**
       * This structure stores neihgborhood information.
       */
      struct Info {
	double     value;
	index_type index;

	Info(double v, index_type r) : value(v), index(r) {}
	Info()                       = default;
	Info(const Info&)            = default;
	Info(Info&&)                 = default;
	Info& operator=(const Info&) = default;
	Info& operator=(Info&&)      = default;
      };

      struct EdgeDistanceInfo {
	std::function<double (unsigned int)> voed;     //!< A function providing a value (double >= 0) according to the number of edges (unsigned int) separating a vertex in the neighborhood from the central vertex.
	unsigned int                         max_dist; //!< The maximal distance considered. 0 means "no limit".
	double                               min_val;  //!< if voed(dist) < min_val, the node is not included in the neighborhood.

	template<typename VALUE_OF_EDGE_DISTANCE>
	EdgeDistanceInfo(const VALUE_OF_EDGE_DISTANCE& voed, unsigned int max_dist, double min_val)
	  : voed(voed), max_dist(max_dist), min_val(min_val) {}

	EdgeDistanceInfo()                                   = default;
	EdgeDistanceInfo(const EdgeDistanceInfo&)            = default;
	EdgeDistanceInfo& operator=(const EdgeDistanceInfo&) = default;
      };
      
      using neighborhood_type      = std::list<Info>;
      using neighborhoods_type     = std::map<NEIGHBORHOOD_KEY, neighborhood_type>;
      using neighborhoods_map_type = std::map<typename GRAPH::ref_vertex, neighborhoods_type>;
      using edge_distances_type    = std::map<neighborhood_key_type, EdgeDistanceInfo>;
      

    private:
      
      std::vector<typename graph_type::ref_vertex> idx2vertex;
      std::map<typename graph_type::vertex_type*, index_type> vertex2idx;
      neighborhoods_map_type neighborhood_tables;
      
      

      friend std::ostream& operator<<(std::ostream& os, Table<graph_type, neighborhood_key_type>& v) {
	os << "Vertex map : " << std::endl;
	unsigned int idx = 0;
	for(auto& ref_v : v.idx2vertex)
	  os << "  " << std::setw(3) << idx++ << " : " << ref_v.get() << std::endl;
	return os;
      }
      
      /**
       * Clears the table. Consider the fill method for a refill.
       */
      void clear_vertices() {
	idx2vertex.clear();
	vertex2idx.clear();
      }

      /**
       * Fills or refills the table from the graph. It is called in the constructor.
       */
      void fill_vertices() {
	g.foreach_vertex([this](const typename graph_type::ref_vertex& ref_v) {
	    auto idx = this->idx2vertex.size();
	    this->idx2vertex.push_back(ref_v);
	    this->vertex2idx[ref_v.get()] = idx;
	  });
      }

      /**
       * WARNING ! Clear the tags (false) before calling that function.
       * @param vertex_index the origin vertex.
       * @param ref_v the origin vertex.
       * @param begin, end Iterators such as it->first is a key, and it->second is an EdgeDistanceInfo. Do not provide duplicate keys.
       * @return a map associating keys with the list of (value, idx) pairs corresponding to the neighborhood. idx is the index of the vertex in a vertices structure. The origin vertex index is in the list (at first position).
       */
      template<typename EDGE_DISTANCE_INFO_IT>
      auto edge_based_neighborhood(index_type vertex_index, const typename graph_type::ref_vertex& ref_v,
				   const EDGE_DISTANCE_INFO_IT& begin, const EDGE_DISTANCE_INFO_IT& end) {
	std::map<decltype(begin->first),std::list<Info>> res;
	std::deque<std::pair<unsigned int, typename graph_type::ref_vertex> > to_do;

	// Let us init info lists with the first vertex (one list for each key).
	for(auto it =  begin; it != end; ++it) {
	  std::list<Info> l;
	  l.push_back({(double)(it->second.voed(0)),vertex_index});
	  res[it->first] = std::move(l);
	}
	(*ref_v)().vq3_tag = true;

	auto job_out = std::back_inserter(to_do);

	// Let us put direct neighbors in the job queue.
	ref_v->foreach_edge([&job_out, &ref_v](const typename graph_type::ref_edge ref_e) {
	    auto extr = ref_e->extremities();
	    if(invalid_extremities(extr)) {ref_e->kill(); return;}
	    auto& other = other_extremity(extr, ref_v);
	    (*other)().vq3_tag = true;
	    *(job_out++) = {1, other};
	  });

	// Let us flush the to do list.
	while(!(to_do.empty())) {
	  auto d_v = to_do.front();
	  to_do.pop_front();
	  bool push_others = false;
	  
	  for(auto it =  begin; it != end; ++it) {
	    auto& edi = it->second;
	    if(double val = edi.voed(d_v.first); val > edi.min_val) {
	      res[it->first].push_back({val, (*this)(d_v.second)});
	      if(edi.max_dist == 0 || d_v.first < edi.max_dist)
		push_others = true;
	    }
	  }
	  
	  if(push_others) 
	    d_v.second->foreach_edge([&job_out, &v = d_v.second, dist = d_v.first + 1](const typename graph_type::ref_edge ref_e) {
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
	
	return res;
      }
      
      /**
       * This builds the neighborhoods.
       */
      void make_neighborhood_tables() {
	neighborhood_tables =  utils::make_vertex_table(g,
							[this](const typename graph_type::ref_vertex& ref_v) {
							  utils::clear_vertex_tags(g, false);
							  return edge_based_neighborhood((*this)(ref_v), ref_v,
											 edge_distances.begin(), edge_distances.end());
							});
      }
    
      
    public:
      
      graph_type& g;
      edge_distances_type edge_distances; //!< This contains different distances. It can be managed directly, or by the use of declare_distance.
      

      Table(graph_type& g) : g(g) {}
      Table()                        = delete;
      Table(const Table&)            = delete;
      Table(Table&&)                 = default;
      Table& operator=(const Table&) = delete;
      Table& operator=(Table&&)      = delete;

      /**
       * This adds a (voed, max_dist, min_val) distance for computing neighborhoods in the edge_distances map.
       * @param key The key for further reference to this distance.
       * @param voed A function providing a value (double >= 0) according to the number of edges (unsigned int) separating a vertex in the neighborhood from the central vertex.
       * @param max_dist The maximal distance considered. 0 means "no limit".
       * @param min_val if voed(dist) < min_val, the node is not included in the neighborhood.
       */
      template<typename VALUE_OF_EDGE_DISTANCE>
      void declare_distance(const neighborhood_key_type& key,
			    const VALUE_OF_EDGE_DISTANCE& voed, unsigned int max_dist, double min_val) {
	edge_distances[key] = EdgeDistanceInfo(voed, max_dist, min_val);
      }
			    
			    

      /**
       * @return the number of vertices in the table.
       */
      const index_type size() const {return idx2vertex.size();}

      /**
       * Updates the vertices only (and not the neighborhoods) (typically after the adding or removal of vertices in the graph).
       */
      void update() {
	clear_vertices();
	fill_vertices();
      }
      
      /**
       * Updates the vertices and neighborhoods (typically after the adding or removal of vertices in the graph).
       */
      void update_full() {
	update();
	make_neighborhood_tables();
      }
      

    
      
      /**
       * WARNING ! Clear the tags (false) before calling that function.
       * @param ref_v the origin vertex.
       * @param voed A function providing a value (double >= 0) according to the number of edges (unsigned int) separating a vertex in the neighborhood from the central vertex.
       * @param max_dist The maximal distance considered. 0 means "no limit".
       * @param min_val if voed(dist) < min_val, the node is not included in the neighborhood.
       * @return The list of (value, idx) pairs corresponding to the neighborhood. idx is the index of the vertex in a vertices structure. The origin vertex index is in the list (at first position).
       */
      template<typename VALUE_OF_EDGE_DISTANCE>
      auto neighborhood(const typename graph_type::ref_vertex& ref_v, const VALUE_OF_EDGE_DISTANCE& voed, unsigned int max_dist, double min_val) {
	EdgeDistanceInfo edi(voed, max_dist, min_val);
	auto minimap = std::make_pair(0, edi);
	return edge_based_neighborhood((*this)(ref_v), ref_v, &minimap, &minimap+1)[0];
      }
      
      /**
       * WARNING ! Clear the tags (false) before calling that function.
       * @param vertex_index the index of the origin vertex.
       * @param voed A function providing a value (double >= 0) according to the number of edges (unsigned int) separating a vertex in the neighborhood from the central vertex.
       * @param max_dist The maximal distance considered. 0 means "no limit".
       * @param min_val if voed(dist) < min_val, the node is not included in the neighborhood.
       * @return The list of (value, idx) pairs corresponding to the neighborhood. idx is the index of the vertex in a vertices structure. The origin vertex index is in the list (at first position).
       */
      template<typename VALUE_OF_EDGE_DISTANCE>
      auto neighborhood(index_type vertex_index, const VALUE_OF_EDGE_DISTANCE& voed, unsigned int max_dist, double min_val) {
	EdgeDistanceInfo edi(voed, max_dist, min_val);
	auto minimap = std::make_pair(0, edi);
	return edge_based_neighborhood(vertex_index, (*this)(vertex_index), &minimap, &minimap+1)[0];
      }
	  

      /**
       * @return the vertex (reference) whose index is idx. Complexiy is contant.
       */
      const typename graph_type::ref_vertex& operator()(index_type idx) const {return idx2vertex[idx];}
      
      /**
       * @return the index of the vertex (reference). Complexiy is logarithmic.
       */
      const index_type operator()(const typename graph_type::ref_vertex& ref_v) const {
	auto pt = ref_v.get();
	auto it = vertex2idx.find(pt);
	if(it == vertex2idx.end()) {
	  std::ostringstream ostr;
	  ostr << "vq3::topology::Table::operator(" << pt << ") : bad vertex reference";
	  throw std::runtime_error(ostr.str());
	}
	return it->second;
      }

      /**
       * @returns the neighborhood named key of vertex #idx. The method update() should be called first in order to update the neigborhoods of all the vertices.
       */
      auto& neighborhood(index_type idx, const neighborhood_key_type& key) const {
	auto ref_v = idx2vertex[idx];
	return neighborhood(ref_v, key);
      }

      /**
       * @returns the neighborhood named key of vertex ref_v. The method update() should be called first in order to update the neigborhood of all the vertices.
       */
      auto& neighborhood(const typename graph_type::ref_vertex& ref_v, const neighborhood_key_type& key) const {
	auto it = neighborhood_tables.find(ref_v);
	return it->second[key];
      }
    };

    template<typename KEY, typename GRAPH>
    auto table(GRAPH& g) {
      return Table<GRAPH, KEY>(g);
    }
    
  }
  
}
