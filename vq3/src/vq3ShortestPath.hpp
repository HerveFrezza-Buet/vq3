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
#include <set>
#include <iterator>

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
      double hcost;                 //!< The cumulated cost plus the heuristics.
      
      void raz() {
	to_src = nullptr;
	cost   = std::numeric_limits<double>::max();
	hcost  = std::numeric_limits<double>::max();
	state  = status::unprocessed;
      }

      template<typename REF_VERTEX>
      void set(double new_val, const REF_VERTEX& parent_edge) {
	to_src = parent_edge;
	cost   = new_val;
	state  = status::processing;
      }

      void set(double new_val) {
	to_src = nullptr;
	cost   = new_val;
	state  = status::processing;
      }

      template<typename REF_VERTEX>
      void set(double new_val, double estimation_to_dest, const REF_VERTEX& parent_edge) {
	to_src = parent_edge;
	cost   = new_val;
	hcost  = new_val + estimation_to_dest;
	state  = status::processing;
      }

      template<typename REF_VERTEX>
      void set_estimated(double new_val, const REF_VERTEX& parent_edge) {
	double h = hcost - cost;
	to_src = parent_edge;
	cost   = new_val;
	hcost  = new_val + h;
	// This is used for vertices which are already in status::processing state. 
      }

      void set(double new_val, double estimation_to_dest) {
	to_src = nullptr;
	cost   = new_val;
	hcost  = new_val + estimation_to_dest;
	state  = status::processing;
      }

      void ended() {
	state = status::done;
      }

      Info() {raz();}
      Info(const Info&)            = default;
      Info& operator=(const Info&) = default;
    };

    inline std::ostream& operator<<(std::ostream& os, const Info& i) {
      if(i.cost == std::numeric_limits<double>::max())
	os << "{Inf, ";
      else
	os << '{' << i.cost << ", ";
      switch(i.state) {
      case status::done : os << "done"; break;
      case status::processing : os << "processing"; break;
      case status::unprocessed : os << "unprocessed"; break;
      }
      os << '}';
      return os;
    }


    template<typename REF_VERTEX>
    class iterator {
    private:
      REF_VERTEX curr;

      
    public:
      using difference_type   = long;
      using value_type        = REF_VERTEX;
      using pointer           = REF_VERTEX*;
      using reference         = REF_VERTEX&;
      using iterator_category = std::input_iterator_tag;

      using edge_type = typename REF_VERTEX::element_type::ref_edge_type::element_type;
      
      /**
       * This is the end iterator.
       */
      iterator()                           = default;
      iterator(const iterator&)            = default;
      iterator& operator=(const iterator&) = default;
      iterator(const REF_VERTEX& curr) : curr(curr) {}

      bool operator==(const iterator& other) {return curr == other.curr;}
      bool operator!=(const iterator& other) {return curr != other.curr;}
      REF_VERTEX operator*() const           {return curr;}

      /**
       * This returns the edge from the vertex pointed by the iterator to the next vertex.
       */
      auto get_edge() {
	return std::reinterpret_pointer_cast<edge_type>((*curr)().vq3_shortest_path.to_src);
      }
      
      auto& operator++() {
	if(auto& sp = (*curr)().vq3_shortest_path.to_src; sp) {
	  auto to_src = std::reinterpret_pointer_cast<edge_type>(sp);
	  auto extr   = to_src->extremities();
	  curr        = vq3::other_extremity(extr, curr);
	}
	else
	  curr = nullptr;
	return *this;
      }
    
      auto operator++(int) {
	auto res = *this;
	(*this)++;
	return res;
      }
    };

  /**
   * Get an iterator on the past starting from a vertex.
   */
  template<typename REF_VERTEX>
  auto begin(const REF_VERTEX& ref_v) {
    return iterator(ref_v);
  }

  /**
   * Get a end interator. The graph is provided for type inference only.
   */
  template<typename GRAPH>
  auto end(const GRAPH& g) {return iterator<typename GRAPH::ref_vertex>();}
  
}

 

  namespace decorator {
    namespace path {


      /* ################# */
      /* #               # */
      /* # Shortest path # */
      /* #               # */
      /* ################# */

      
      template<typename MOTHER, typename KIND>
      struct Shortest: public MOTHER {
	using decorated_type = typename MOTHER::decorated_type;
	vq3::path::Info vq3_shortest_path;
	Shortest(const decorated_type& val) : MOTHER(val), vq3_shortest_path() {}
	Shortest& operator=(const decorated_type& val) {this->vq3_value = val;}
      };

      //When we decorate a non decorated value
      template<typename MOTHER>
      struct Shortest<MOTHER, not_decorated> {
	using decorated_type = MOTHER;
	MOTHER vq3_value;
	vq3::path::Info vq3_shortest_path;
	Shortest(const decorated_type& val) : vq3_value(val), vq3_shortest_path() {}
	/** Affectation from a value has to work */
	Shortest& operator=(const decorated_type& val) {vq3_value = val;}
      };

      // When we decorate a decorated type with no value.
      template<typename MOTHER>
      struct Shortest<MOTHER, unvalued_decoration> : public MOTHER {
	using decorated_type = MOTHER;
	vq3::path::Info vq3_shortest_path;
	Shortest() : MOTHER(), vq3_shortest_path() {}
      };

      // When we decorate void
      template<>
      struct Shortest<void, not_decorated> {
	using decorated_type = void;
	vq3::path::Info vq3_shortest_path;
	Shortest() : vq3_shortest_path() {}
      };

      // User-friendly decorator  
      template<typename MOTHER>
      using shortest = Shortest<MOTHER, typename decoration<MOTHER>::value_type>;
    }

    
  }


  namespace path {

    /**
     * This function fills the vq3_shortest_path values at each vertex
     * according to the shortest path linking each vertex to the
     * destination vertex. The boolean template parameters
     * VERTEX_EFFICIENCY and EDGE_EFFICIENCY toggle the efficiency
     * test on vertices and edges (see the efficiency decorator).
     * @param g the graph
     * @param start the start vertex. Provide nullptr for computing shortest path from all the vertices to dest.
     * @param dest the destination vertex
     * @param edge_cost a function such as edge_cost(ref_edge) is the cost of the edge.
     */
    template<bool VERTEX_EFFICIENCY, bool EDGE_EFFICIENCY, typename GRAPH, typename EDGE_COST>
    void dijkstra(GRAPH& g, typename GRAPH::ref_vertex& start, typename GRAPH::ref_vertex& dest, const EDGE_COST& edge_cost) {
      // Init
      g.foreach_vertex([](typename GRAPH::ref_vertex ref_v){(*ref_v)().vq3_shortest_path.raz();});

      auto accum_comp = [](const typename GRAPH::ref_vertex& ref_v1,
			   const typename GRAPH::ref_vertex& ref_v2) {return (*ref_v1)().vq3_shortest_path.cost < (*ref_v2)().vq3_shortest_path.cost;};
      std::set<typename GRAPH::ref_vertex, decltype(accum_comp)> q(accum_comp);

      if constexpr(VERTEX_EFFICIENCY) {
	  if((*dest)().vq3_efficient) {
	    (*dest)().vq3_shortest_path.set(0);
	    q.insert(dest);
	  }
	  else
	    return;
	}
      else {
	(*dest)().vq3_shortest_path.set(0);
	q.insert(dest);
      }
      
      while(!q.empty()) {
	auto curr = *(q.begin());
	auto& curr_path_info = (*curr)().vq3_shortest_path;
	curr_path_info.ended();
	if(curr == start) break;
	
	q.erase(q.begin());
	curr->foreach_edge([curr, &edge_cost, &q, &curr_path_info](typename GRAPH::ref_edge ref_e) {
	    auto extr_pair = ref_e->extremities();           
	    if(vq3::invalid_extremities(extr_pair)) {
	      ref_e->kill();
	      return;
	    }

	    if constexpr(EDGE_EFFICIENCY) {
		if(!(*ref_e)().vq3_efficient)
		  return;
	      }
			  
	    double cost           = edge_cost(ref_e);
	    auto& other           = vq3::other_extremity(extr_pair, curr);
	    auto& other_path_info = (*other)().vq3_shortest_path;

	    if constexpr(VERTEX_EFFICIENCY) {
		if(!(*other)().vq3_efficient) {
		  other_path_info.ended();
		  return;
		}
	      }
	    
	    switch(other_path_info.state) {
	    case status::done :
	      break;
	    case status::unprocessed :
	      other_path_info.set(cost + curr_path_info.cost, ref_e);
	      q.insert(other);
	      break;
	    case status::processing :
	      if(double cost_candidate = cost + curr_path_info.cost; cost_candidate < other_path_info.cost) {
		q.erase(other);
		other_path_info.set(cost_candidate, ref_e);
		q.insert(other);
	      }
	      break;
	    }
	  });
      }
    }

    /**
     * This function fills the vq3_shortest_path values at each vertex
     * according to the shortest path linking each vertex to the
     * destination vertex. The boolean template parameters
     * VERTEX_EFFICIENCY and EDGE_EFFICIENCY toggle the efficiency
     * test on vertices and edges (see the efficiency decorator).
     * @param g the graph
     * @param start the start vertex. Provide nullptr for computing shortest path from all the vertices to dest.
     * @param dest the destination vertex
     * @param edge_cost a function such as edge_cost(ref_edge) is the cost of the edge.
     * @param to_start_estimation a function such as to_dest_estimation(ref_vertex) gives a pessimistic optimation of the cumulated cost from ref_vertex to start.
     */
    template<bool VERTEX_EFFICIENCY, bool EDGE_EFFICIENCY, typename GRAPH, typename EDGE_COST, typename TO_START>
    void a_star(GRAPH& g, typename GRAPH::ref_vertex& start, typename GRAPH::ref_vertex& dest,
		const EDGE_COST& edge_cost,
		const TO_START& to_start_estimation) {
      // Init
      g.foreach_vertex([](typename GRAPH::ref_vertex ref_v){(*ref_v)().vq3_shortest_path.raz();});

      auto accum_comp = [](const typename GRAPH::ref_vertex& ref_v1,
			   const typename GRAPH::ref_vertex& ref_v2) {return (*ref_v1)().vq3_shortest_path.hcost < (*ref_v2)().vq3_shortest_path.hcost;};
      std::set<typename GRAPH::ref_vertex, decltype(accum_comp)> q(accum_comp);

      if constexpr(VERTEX_EFFICIENCY) {
	  if((*dest)().vq3_efficient) {
	    (*dest)().vq3_shortest_path.set(0, to_start_estimation(dest));
	    q.insert(dest);
	  }
	  else
	    return;
	}
      else {
	(*dest)().vq3_shortest_path.set(0, to_start_estimation(dest));
	q.insert(dest);
      }
      
      while(!q.empty()) {
	auto curr = *(q.begin());
	auto& curr_path_info = (*curr)().vq3_shortest_path;
	curr_path_info.ended();
	if(curr == start) break;
	
	q.erase(q.begin());
	curr->foreach_edge([curr, &edge_cost, &to_start_estimation, &q, &curr_path_info](typename GRAPH::ref_edge ref_e) {
	    auto extr_pair = ref_e->extremities();           
	    if(vq3::invalid_extremities(extr_pair)) {
	      ref_e->kill();
	      return;
	    }

	    if constexpr(EDGE_EFFICIENCY) {
		if(!(*ref_e)().vq3_efficient)
		  return;
	      }
			  
	    double cost           = edge_cost(ref_e);
	    auto& other           = vq3::other_extremity(extr_pair, curr);
	    auto& other_path_info = (*other)().vq3_shortest_path;

	    if constexpr(VERTEX_EFFICIENCY) {
		if(!(*other)().vq3_efficient) {
		  other_path_info.ended();
		  return;
		}
	      }
	    
	    switch(other_path_info.state) {
	    case status::done :
	      break;
	    case status::unprocessed :
	      other_path_info.set(cost + curr_path_info.cost, to_start_estimation(other), ref_e);
	      q.insert(other);
	      break;
	    case status::processing :
	      if(double cost_candidate = cost + curr_path_info.cost; cost_candidate < other_path_info.cost) {
		q.erase(other);
		other_path_info.set_estimated(cost_candidate, ref_e);
		q.insert(other);
	      }
	      break;
	    }
	  });
      }
    }





    
  }

  
}
