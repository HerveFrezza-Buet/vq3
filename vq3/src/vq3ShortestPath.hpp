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
#include <optional>
#include <limits>

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
      unsigned int qpos;            //!< The position in ht priority queue (temporary information).
      
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

      Info() : qpos(0) {raz();}
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
     * Get an iterator on the path starting from a vertex.
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
	Shortest() = default;
	Shortest(const decorated_type& val) : MOTHER(val), vq3_shortest_path() {}
	Shortest& operator=(const decorated_type& val) {this->vq3_value = val;}
      };

      //When we decorate a non decorated value
      template<typename MOTHER>
      struct Shortest<MOTHER, not_decorated> {
	using decorated_type = MOTHER;
	MOTHER vq3_value;
	vq3::path::Info vq3_shortest_path;
	Shortest() = default;
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

    namespace priority_queue {

      namespace internal {
      
	template<bool USE_HCOST, typename REF_VERTEX>
	inline unsigned int swap_if_1_is_higher(std::vector<REF_VERTEX>& q, unsigned int pos1, unsigned int pos2) {

	  // if(pos1 == 0 || pos1 >= q.size() || pos2  >= q.size()) {
	  //   std::cout << "Error : swap_if_1_is_higher(q(" << q.size() << "), " << pos1 << ", " << pos2 << ')' << std::endl;
	  //   ::exit(0);
	  // }
	  
	  auto it1  = q.begin() + pos1;
	  auto it2  = q.begin() + pos2;
	  auto& sp1 = (*(*it1))().vq3_shortest_path;
	  auto& sp2 = (*(*it2))().vq3_shortest_path;

	  if constexpr(USE_HCOST) {
	      if(sp1.hcost < sp2.hcost)
		return 0;
	    }
	  else {
	    if(sp1.cost < sp2.cost)
	      return 0;
	  }
	
	  sp1.qpos = pos2;
	  sp2.qpos = pos1;
	  std::iter_swap(it1, it2);

	  return pos1;
	}
      
	template<bool USE_HCOST, typename REF_VERTEX>
	inline unsigned int swap_if_1_is_higher_than_min_2_3(std::vector<REF_VERTEX>& q, unsigned int pos1, unsigned int pos2, unsigned int pos3) {
	  // if(pos1 >= q.size() || pos2  >= q.size() || pos3  >= q.size()) {
	  //   std::cout << "Error : swap_if_1_is_higher_than_min_2_3(q(" << q.size() << "), " << pos1 << ", " << pos2 << ", " << pos3 <<')' << std::endl;
	  //   ::exit(0);
	  // }
	  
	  auto it1  = q.begin() + pos1;
	  auto it2  = q.begin() + pos2;
	  auto it3  = q.begin();

	  Info* sp1_ptr = &((*(*it1))().vq3_shortest_path);
	  Info* sp2_ptr = &((*(*it2))().vq3_shortest_path);
	  Info* sp3_ptr = nullptr;

	  if(pos3 != 0) {
	    it3 += pos3;
	    sp3_ptr = &((*(*it3))().vq3_shortest_path);
	  }

	  if constexpr(USE_HCOST) {
	    if(pos3 == 0) {                       // min(sp2, 0) = sp2
	      if(sp1_ptr->hcost > sp2_ptr->hcost) { //   if sp1 > min, we swap 1 and 2
		sp1_ptr->qpos = pos2;
		sp2_ptr->qpos = pos1;
		std::iter_swap(it1, it2);
		return pos2;
	      }
	      else
		return 0;                        //  else, no swap. done.
	    }

	    if(sp2_ptr->hcost < sp3_ptr->hcost) {   // min(sp2, sp3) = sp2
	      if(sp1_ptr->hcost > sp2_ptr->hcost) { //   if sp1 > min, we swap 1 and 2
		sp1_ptr->qpos = pos2;
		sp2_ptr->qpos = pos1;
		std::iter_swap(it1, it2);
		return pos2;
	      }
	      else
		return 0;                         //  else, no swap. done.
	    }

	    // min(sp2, sp3) = sp3
	    if(sp1_ptr->hcost > sp3_ptr->hcost) { //   if sp1 > min, we swap 1 and 3
		sp1_ptr->qpos = pos3;
		sp3_ptr->qpos = pos1;
		std::iter_swap(it1, it3);
		return pos3;
	      }
	      else
		return 0;  
	    }
	  else {
	    if(pos3 == 0) {                       // min(sp2, 0) = sp2
	      if(sp1_ptr->cost > sp2_ptr->cost) { //   if sp1 > min, we swap 1 and 2
		sp1_ptr->qpos = pos2;
		sp2_ptr->qpos = pos1;
		std::iter_swap(it1, it2);
		return pos2;
	      }
	      else
		return 0;                        //  else, no swap. done.
	    }

	    if(sp2_ptr->cost < sp3_ptr->cost) {   // min(sp2, sp3) = sp2
	      if(sp1_ptr->cost > sp2_ptr->cost) { //   if sp1 > min, we swap 1 and 2
		sp1_ptr->qpos = pos2;
		sp2_ptr->qpos = pos1;
		std::iter_swap(it1, it2);
		return pos2;
	      }
	      else
		return 0;                         //  else, no swap. done.
	    }

	    // min(sp2, sp3) = sp3
	    if(sp1_ptr->cost > sp3_ptr->cost) { //   if sp1 > min, we swap 1 and 3
		sp1_ptr->qpos = pos3;
		sp3_ptr->qpos = pos1;
		std::iter_swap(it1, it3);
		return pos3;
	      }
	      else
		return 0;  
	  }
	}
      
	template<bool USE_HCOST, typename REF_VERTEX>
	inline void percolate_up(std::vector<REF_VERTEX>& q, unsigned int pos) {
	  do
	    if(pos == 1) return;
	  while((pos = swap_if_1_is_higher<USE_HCOST, REF_VERTEX>(q, pos / 2, pos)) != 0);
	}
      
	template<bool USE_HCOST, typename REF_VERTEX>
	inline void percolate_down(std::vector<REF_VERTEX>& q, unsigned int pos) {
	  unsigned int son_left  = 2 * pos;
	  unsigned int son_right;

	  do {
	    if(son_left = 2 * pos; son_left >= q.size()) return;
	    if(son_right = son_left + 1; son_right >= q.size()) son_right = 0;
	  }
	  while((pos = swap_if_1_is_higher_than_min_2_3<USE_HCOST>(q, pos, son_left, son_right)) != 0);
	}
      }
      
      template<typename REF_VERTEX>
      inline void make_empty(std::vector<REF_VERTEX>& q) {
	q.clear();
	q.push_back(nullptr);
      }
      
      template<typename REF_VERTEX>
      inline bool is_empty(std::vector<REF_VERTEX>& q) {
	return q.size() == 1;
      }

      template<bool USE_HCOST, typename REF_VERTEX>
      inline void push(std::vector<REF_VERTEX>& q, const REF_VERTEX& ref_v) {
	// std::cout << "Push(q(" << q.size() << "))" << std::endl;
	q.push_back(ref_v);
	(*ref_v)().vq3_shortest_path.qpos = q.size() - 1;
	internal::percolate_up<USE_HCOST, REF_VERTEX>(q, q.size() - 1);
      }
      
      template<bool USE_HCOST, typename REF_VERTEX>
      inline REF_VERTEX pop(std::vector<REF_VERTEX>& q) {
	// std::cout << "Pop(q(" << q.size() << "))" << std::endl;
	auto top  = q.begin() + 1;
	auto last = q.end()   - 1;
	auto res  = *top;
	*top = *last;
	(*(*top))().vq3_shortest_path.qpos = 1;
	q.pop_back();
	internal::percolate_down<USE_HCOST, REF_VERTEX>(q, 1);
	return res;
      }
      
      template<bool USE_HCOST, typename REF_VERTEX>
      inline void notify_decrease(std::vector<REF_VERTEX>& q, unsigned int pos) {
	// std::cout << "Decrease(q(" << q.size() << "), " << pos << ")" << std::endl;
	internal::percolate_up<USE_HCOST, REF_VERTEX>(q, pos);
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
     */
    template<bool VERTEX_EFFICIENCY, bool EDGE_EFFICIENCY, typename GRAPH, typename EDGE_COST>
    void dijkstra(GRAPH& g, typename GRAPH::ref_vertex start, typename GRAPH::ref_vertex dest, const EDGE_COST& edge_cost) {

      
      // Init
      priority_queue::make_empty(g.heap);
      g.foreach_vertex([](typename GRAPH::ref_vertex ref_v){(*ref_v)().vq3_shortest_path.raz();});

      if constexpr(VERTEX_EFFICIENCY) {
	  if((*dest)().vq3_efficient) {
	    (*dest)().vq3_shortest_path.set(0);
	    priority_queue::push<false>(g.heap, dest);
	  }
	  else
	    return;
	}
      else {
	(*dest)().vq3_shortest_path.set(0);
	priority_queue::push<false>(g.heap, dest);
      }

      while(!priority_queue::is_empty(g.heap)) {
	auto curr = priority_queue::pop<false>(g.heap);
	auto& curr_path_info = (*curr)().vq3_shortest_path;
	curr_path_info.ended();
	if(curr == start) break;
	
	curr->foreach_edge([curr, &g, &edge_cost, &curr_path_info](typename GRAPH::ref_edge ref_e) {
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
	      priority_queue::push<false>(g.heap, other);
	      break;
	    case status::processing :
	      if(double cost_candidate = cost + curr_path_info.cost; cost_candidate < other_path_info.cost) {
		other_path_info.set(cost_candidate, ref_e);
		priority_queue::notify_decrease<false>(g.heap, other_path_info.qpos);
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
     * @param out It is an output iterator for storing the vertices ib the order they are extracted from the inner queue.
     */
    template<bool VERTEX_EFFICIENCY, bool EDGE_EFFICIENCY, typename GRAPH, typename EDGE_COST, typename OUTPUT_ITERATOR>
    void dijkstra(GRAPH& g, typename GRAPH::ref_vertex start, typename GRAPH::ref_vertex dest, const EDGE_COST& edge_cost, OUTPUT_ITERATOR out) {

      
      // Init
      priority_queue::make_empty(g.heap);
      g.foreach_vertex([](typename GRAPH::ref_vertex ref_v){(*ref_v)().vq3_shortest_path.raz();});

      if constexpr(VERTEX_EFFICIENCY) {
	  if((*dest)().vq3_efficient) {
	    (*dest)().vq3_shortest_path.set(0);
	    priority_queue::push<false>(g.heap, dest);
	  }
	  else
	    return;
	}
      else {
	(*dest)().vq3_shortest_path.set(0);
	priority_queue::push<false>(g.heap, dest);
      }

      while(!priority_queue::is_empty(g.heap)) {
	auto curr = priority_queue::pop<false>(g.heap);
	*(out++) = curr;
	auto& curr_path_info = (*curr)().vq3_shortest_path;
	curr_path_info.ended();
	if(curr == start) break;
	
	curr->foreach_edge([curr, &g, &edge_cost, &curr_path_info](typename GRAPH::ref_edge ref_e) {
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
	      priority_queue::push<false>(g.heap, other);
	      break;
	    case status::processing :
	      if(double cost_candidate = cost + curr_path_info.cost; cost_candidate < other_path_info.cost) {
		other_path_info.set(cost_candidate, ref_e);
		priority_queue::notify_decrease<false>(g.heap, other_path_info.qpos);
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
    void a_star(GRAPH& g, typename GRAPH::ref_vertex start, typename GRAPH::ref_vertex dest,
		const EDGE_COST& edge_cost,
		const TO_START& to_start_estimation) {      

      if(start == nullptr) {
	dijkstra<VERTEX_EFFICIENCY, EDGE_EFFICIENCY>(g, nullptr, dest, edge_cost);
	return;
      }
      
      // Init
      priority_queue::make_empty(g.heap);
      g.foreach_vertex([](typename GRAPH::ref_vertex ref_v){(*ref_v)().vq3_shortest_path.raz();});

      if constexpr(VERTEX_EFFICIENCY) {
	  if((*dest)().vq3_efficient) {
	    (*dest)().vq3_shortest_path.set(0, to_start_estimation(dest));
	    priority_queue::push<true>(g.heap, dest);
	  }
	  else
	    return;
	}
      else {
	(*dest)().vq3_shortest_path.set(0, to_start_estimation(dest));
	priority_queue::push<true>(g.heap, dest);
      }
      
      while(!priority_queue::is_empty(g.heap)) {
	auto curr = priority_queue::pop<true>(g.heap);
	auto& curr_path_info = (*curr)().vq3_shortest_path;
	curr_path_info.ended();
	if(curr == start) break;
	
	curr->foreach_edge([curr, &edge_cost, &to_start_estimation, &g, &curr_path_info](typename GRAPH::ref_edge ref_e) {
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
	      priority_queue::push<true>(g.heap, other);
	      break;
	    case status::processing :
	      if(double cost_candidate = cost + curr_path_info.cost; cost_candidate < other_path_info.cost) {
		other_path_info.set_estimated(cost_candidate, ref_e);
		priority_queue::notify_decrease<true>(g.heap, other_path_info.qpos);
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
     * @param out It is an output iterator for storing the vertices ib the order they are extracted from the inner queue.
     */
    template<bool VERTEX_EFFICIENCY, bool EDGE_EFFICIENCY, typename GRAPH, typename EDGE_COST, typename TO_START, typename OUTPUT_ITERATOR>
    void a_star(GRAPH& g, typename GRAPH::ref_vertex start, typename GRAPH::ref_vertex dest,
		const EDGE_COST& edge_cost,
		const TO_START& to_start_estimation,
		OUTPUT_ITERATOR out) {      

      if(start == nullptr) {
	dijkstra<VERTEX_EFFICIENCY, EDGE_EFFICIENCY>(g, nullptr, dest, edge_cost);
	return;
      }
      
      // Init
      priority_queue::make_empty(g.heap);
      g.foreach_vertex([](typename GRAPH::ref_vertex ref_v){(*ref_v)().vq3_shortest_path.raz();});

      if constexpr(VERTEX_EFFICIENCY) {
	  if((*dest)().vq3_efficient) {
	    (*dest)().vq3_shortest_path.set(0, to_start_estimation(dest));
	    priority_queue::push<true>(g.heap, dest);
	  }
	  else
	    return;
	}
      else {
	(*dest)().vq3_shortest_path.set(0, to_start_estimation(dest));
	priority_queue::push<true>(g.heap, dest);
      }
      
      while(!priority_queue::is_empty(g.heap)) {
	auto curr = priority_queue::pop<true>(g.heap);
	*(out++) = curr;
	auto& curr_path_info = (*curr)().vq3_shortest_path;
	curr_path_info.ended();
	if(curr == start) break;
	
	curr->foreach_edge([curr, &edge_cost, &to_start_estimation, &g, &curr_path_info](typename GRAPH::ref_edge ref_e) {
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
	      priority_queue::push<true>(g.heap, other);
	      break;
	    case status::processing :
	      if(double cost_candidate = cost + curr_path_info.cost; cost_candidate < other_path_info.cost) {
		other_path_info.set_estimated(cost_candidate, ref_e);
		priority_queue::notify_decrease<true>(g.heap, other_path_info.qpos);
	      }
	      break;
	    }
	  });
      }
    }



    

    /**
     * Computes a position along the path. Some shortest path
     * algorithm allowing to compute the path from start to dest has
     * to have been performed beforehand. The vertex dest is the one with
     * a 0 cummulated cost.
     * @param start, dest The extremities of the path
     * @param lambda 0 means start, 1 means dest, inbetween value leads to an interpolation.
     * @param interpolate A function such as interpolate(ref_v1, ref_v2, lambda) returnes something inbetween ref_v1 (lambda = 0) and ref_v2 (lambda = 1). If ref_v1 is nullptr, the result corresponds to ref_v2, whatever lambda. The same stands is ref_v2 is nullptr.
     * @return An interpolated value. It is optional, since the path may not extist.
     */
    template<typename REF_VERTEX, typename INTERPOLATE>
    auto travel(const REF_VERTEX& start, const REF_VERTEX& dest, double lambda, const INTERPOLATE& interpolate) -> std::optional<decltype(interpolate(start, nullptr, 0))> {
      if(start == dest)
	return interpolate(start, nullptr, 0);

      if((*start)().vq3_shortest_path.state == status::unprocessed)
	return std::optional<decltype(interpolate(start, nullptr, 0))>();
      
      double l = 1.0 - std::max(0.0, std::min(lambda, 1.0));
      if(l == 1)
	return interpolate(start, nullptr, 0);
      if(l == 0)
	return interpolate(dest, nullptr, 0);

      auto prev_cost    = (*start)().vq3_shortest_path.cost;
      auto to_dest_cost = l * prev_cost;
      REF_VERTEX prev   = start;
      auto it           = begin(start);
      auto cur          = *(++it);
      auto cur_cost     = (*cur)().vq3_shortest_path.cost;
      while(cur_cost > to_dest_cost) {
	prev      = cur;
	prev_cost = cur_cost;
	cur       = *(++it);
	cur_cost  = (*cur)().vq3_shortest_path.cost;
      }
      
      double lbd = .5;
      if(double dif = prev_cost - cur_cost; dif > 0)
	lbd = (to_dest_cost - cur_cost) / dif;
      return interpolate(cur, prev, lbd);
    }
  }
}
