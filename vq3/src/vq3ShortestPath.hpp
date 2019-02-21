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

#include <memory>
#include <queue>

namespace vq3 {
  namespace path {

    enum class status : char {
      unprocessed, //!< not processed yet.
      processing,  //!< considered, but best path not found yet.
      done         //!< best path found.
    };
    
    /**
     * This is what shortest path calculation algorithms need to store at each node.
     */
    struct Info {
      status state;                 //!< Tells whether best path is found for this node.
      std::shared_ptr<void> to_src; //!< The ref edge leading to the source.
      double cost;                  //!< The cumulated cost.
      
      void raz() {
	to_src = nullptr;
	cost   = std::numeric_limits<double>::max();
	state  = status::unprocessed;
      }

      template<typename REF_VERTEX>
      void set(double new_val, const REF_VERTEX& parent_edge) {
	to_src     = parent_edge;
	cost = new_val;
	state      = status::processing;
      }

      void set(double new_val) {
	to_src     = nullptr;
	cost = new_val;
	state      = status::processing;
      }

      void ended() {
	state = status::done;
      }

      bool operator<(const Info& other) {return accum < other.accum;}
      
      Info() {raz();}
      Info(const Info&)            = default;
      Info& operator=(const Info&) = default;
    };
  }

 

  namespace decorator {
    namespace path {


      /* ################# */
      /* #               # */
      /* # Shortest path # */
      /* #               # */
      /* ################# */

      
      template<typename MOTHER, typename KIND>
      struct ShortestPath: public MOTHER {
	using decorated_type = typename MOTHER::decorated_type;
	vq3::path::Info vq3_shortest_path;
	ShortestPath(const decorated_type& val) : MOTHER(val), vq3_shortest_path() {}
	ShortestPath& operator=(const decorated_type& val) {this->vq3_value = val;}
      };

      //When we decorate a non decorated value
      template<typename MOTHER>
      struct ShortestPath<MOTHER, not_decorated> {
	using decorated_type = MOTHER;
	MOTHER vq3_value;
	vq3::path::Info vq3_shortest_path;
	ShortestPath(const decorated_type& val) : vq3_value(val), vq3_shortest_path() {}
	/** Affectation from a value has to work */
	ShortestPath& operator=(const decorated_type& val) {vq3_value = val;}
      };

      // When we decorate a decorated type with no value.
      template<typename MOTHER>
      struct ShortestPath<MOTHER, unvalued_decoration> : public MOTHER {
	using decorated_type = MOTHER;
	vq3::path::Info vq3_shortest_path;
	ShortestPath() : MOTHER(), vq3_shortest_path() {}
      };

      // When we decorate void
      template<>
      struct ShortestPath<void, not_decorated> {
	using decorated_type = void;
	vq3::path::Info vq3_shortest_path;
	ShortestPath() : vq3_shortest_path() {}
      };

      // User-friendly decorator  
      template<typename MOTHER>
      using shortest_path = ShortestPath<MOTHER, typename decoration<MOTHER>::value_type>;
    }

    
  }


  namespace path {

    /**
     * This function fills the vq3_shortest_path values at each vertex
     * according to the shortest path linking each vertex to the destination vertex.
     * @param g the graph
     * @param start the start vertex. Provide nullptr for computing shortest path from all the vertices to dest.
     * @param dest the destination vertex
     * @param edge_cost a function such as edge_cost(ref_edge) is the cost of the edge.
     */
    template<typename GRAPH, typename EDGE_COST>
    void dijkstra(GRAPH& g, typename GRAPH::ref_vertex& start, typename GRAPH::ref_vertex& dest, const EDGE_COST& edge_cost) {
      // Init
      g.foreach_vertex([](typename GRAPH::ref_vertex ref_v){(*ref_v)().vq3_shortest_path.raz();});

      auto accum_comp = [](const typename GRAPH::ref_vertex& ref_v1,
			   const typename GRAPH::ref_vertex& ref_v2) {return (*ref_v1)().vq3_shortest_path < (*ref_v2)().vq3_shortest_path;};
      std::set<typename GRAPH::ref_vertex, decltype(accum_comp)> q(accum_comp);

      (*dest)().vq3_shortest_path.set(0);
      q.insert(dest);
      while(!q.empty()) {
	auto curr = q.top();
	q.pop();
	curr.ended();
	if(curr == start) break;
	
	curr->foreach_edge([curr, &edge_cost, &q](typename GRAPH::ref_ref_edge ref_e) {
	    auto extr_pair = ref_e->extremities();           
	    if(vq3::invalid_extremities(extr_pair)) {
	      ref_e->kill();
	      return;
	    }
	    double cost = edge_cost(ref_e);
	    auto& other = vq3::other_extremity(extr_pair, curr);
	    auto& other_path_info = (*other)().vq3_shortest_path;
	    switch(other_path_info.state) {
	    case status::done :
	      break;
	    case status::unprocessed :
	      other_path_info.set(cost + (*curr)().vq3_shortest_path.cost, ref_e);
	      q.insert(other);
	      break;
	    case status::processing :
	      if(double cost_candidate = cost + (*curr)().vq3_shortest_path.cost; cost_candidate < other_path_info.cost) {
		q.erase(other);
		other_path_info.set(cost_candidate, ref_e);
		q.insert(other);
	      }
	      break;
	    }
	  });
      }
    }
  }

  
}
