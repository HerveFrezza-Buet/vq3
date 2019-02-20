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

namespace vq3 {
  namespace path {

    /**
     * This is what shortest path calculation algorithms need to store at each node.
     */
    struct Info {
      std::shared_ptr<void> to_src; //!< The ref edge leading to the source.
      double cost_accum;            //!< The cumulated cost.
      
      void raz(double accum = std::numeric_limits<double>::max()) {
	to_src     = nullptr;
	cost_accum = accum; 
      }
      
      Info()                       = default;
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
     * @param dest the destination vertex
     * @param edge_cost a function such as edge_cost(ref_edge) is the cost of the edge.
     */
    template<typename GRAPH, typename EDGE_COST>
    void dijkstra(GRAPH& g, typename GRAPH::ref_vertex& dest, const EDGE_COST& edge_cost) {
    }
  }

  
}
