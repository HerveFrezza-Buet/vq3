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
    template<typename GRAPH>
    class Table {
    public:
      using graph_type      = GRAPH;
      using index_type      = typename std::vector<typename GRAPH::ref_vertex>::size_type;
      
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

      using neighborhood_table_type = std::map<typename graph_type::ref_vertex, std::list<Info> >;

    private:

      
      std::vector<typename graph_type::ref_vertex> idx2vertex;
      std::map<typename graph_type::vertex_type*, index_type> vertex2idx;
      neighborhood_table_type neighborhood_table;




      

      friend std::ostream& operator<<(std::ostream& os, Table<graph_type>& v) {
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
       * @param ref_v the origin vertex.
       * @param voed A function providing a value (double >= 0) according to the number of edges (unsigned int) separating a vertex in the neighborhood from the central vertex.
       * @param max_dist The maximal distance considered. 0 means "no limit".
       * @param min_val if voed(dist) < min_val, the node is not included in the neighborhood.
       * @return The list of (value, idx) pairs corresponding to the neighborhood. idx is the index of the vertex in a vertices structure. The origin vertex index is in the list (at first position).
       */
      template<typename VALUE_OF_EDGE_DISTANCE>
      auto edge_based_neighborhood(index_type vertex_index, const typename graph_type::ref_vertex& ref_v, const VALUE_OF_EDGE_DISTANCE& voed, unsigned int max_dist, double min_val) {
	std::list<Info> res;
	std::deque<std::pair<unsigned int, typename graph_type::ref_vertex> > to_do;

	auto res_out = std::back_inserter(res);
	auto job_out = std::back_inserter(to_do);
	*(res_out++) = {(double)(voed(0)), vertex_index};
	(*ref_v)().vq3_tag = true;
	ref_v->foreach_edge([&job_out, &ref_v](const typename graph_type::ref_edge ref_e) {
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
	    *(res_out++) = {val, (*this)(d_v.second)};
	    if(d_v.first != max_dist)
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
	}
	return res;
      }
      
      /**
       * This builds a std::map, keys are ref_vertex type, values are
       * the neighborhood returned by the neighborhood method.
       */
      template<typename VALUE_OF_EDGE_DISTANCE>
      void make_neighborhood_table(const VALUE_OF_EDGE_DISTANCE& voed, unsigned int max_dist, double min_val) {
	neighborhood_table =  utils::make_vertex_table(g,
						       [this, &voed, max_dist, min_val](const typename graph_type::ref_vertex& ref_v) {
							 utils::clear_vertex_tags(g, false); 
							 return neighborhood(ref_v, voed, max_dist, min_val); 
						       });
      }
    
      
    public:
      
      graph_type& g;

      Table(graph_type& g) : g(g) {}
      Table()                        = delete;
      Table(const Table&)            = delete;
      Table(Table&&)                 = default;
      Table& operator=(const Table&) = delete;
      Table& operator=(Table&&)      = delete;

      /**
       * @return the number of vertices in the table.
       */
      const index_type size() const {return idx2vertex.size();}

      /**
       * Updates the vertices only (typically after the adding or removal of vertices in the graph).
       */
      void operator()() {
	clear_vertices();
	fill_vertices();
      }
      
      /**
       * Updates the vertices and neighbouts (typically after a topology change in terms of vertices and/or edges of the graph).
       */
      template<typename VALUE_OF_EDGE_DISTANCE>
      void operator()(const VALUE_OF_EDGE_DISTANCE& voed, unsigned int max_dist, double min_val) {
	  clear_vertices();
	  fill_vertices();
	  make_neighborhood_table(voed, max_dist, min_val);
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
	return edge_based_neighborhood((*this)(ref_v), ref_v, voed, max_dist, min_val);
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
	return edge_based_neighborhood(vertex_index, (*this)(vertex_index),voed, max_dist, min_val);
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
       * @returns the neighborhood of vertex #idx.
       */
      auto& operator[](index_type idx) const {
	return neighborhood_table[idx];
      }

      /**
       * @returns the neighborhood of vertex ref_v.
       */
      auto& operator[](const typename graph_type::ref_vertex& ref_v) const {
	auto it = neighborhood_table.find(ref_v);
	return it->second;
      }
    };

    template<typename GRAPH>
    auto table(GRAPH& g) {
      return Table<GRAPH>(g);
    }
    
  }
  
}
