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
#include <list> 
#include <utility>
#include <atomic>
#include <iostream>
#include <map>
#include <functional>


namespace vq3 {

  template<typename PAIR>
  bool invalid_extremities(const PAIR& p) {return p.first == nullptr || p.second == nullptr;}
  
  template<typename PAIR>
  auto& other_extremity(const PAIR& p, const typename PAIR::first_type& me) {
    if(me == p.second)
      return p.first;
    return p.second;
  }

  template<typename VERTEX_VALUE, typename EDGE_VALUE> class graph;
  template<typename VERTEX_VALUE, typename EDGE_VALUE> class graph_;
  template<typename VERTEX_VALUE, typename EDGE_VALUE> class vertex;
  template<typename VERTEX_VALUE, typename EDGE_VALUE> class edge;

  template<typename VERTEX_VALUE, typename EDGE_VALUE>
  class graph_element {
  protected:

    friend class graph_<VERTEX_VALUE, EDGE_VALUE>;
    
    std::atomic_bool killed;
      
    graph_element() : killed(false) {}
      
  public:
      
    graph_element(const graph_element&)            = delete;
    graph_element& operator=(const graph_element&) = delete;

    /**
     * This is a self-destruction request. This object will be
     * ignored in further iterations and freed soon.
     */
    void kill() {killed = true;}

    /**
     * Tells wether a kill request is pending on this element.
     */
    bool is_killed() {return killed;}
  };
  
  template<typename VALUE, typename VERTEX_VALUE, typename EDGE_VALUE>
  class valued_graph_element : public graph_element<VERTEX_VALUE, EDGE_VALUE> {
  protected:

    VALUE v;
      
    valued_graph_element() : graph_element<VERTEX_VALUE, EDGE_VALUE>(), v() {}
    valued_graph_element(const VALUE& v) : graph_element<VERTEX_VALUE, EDGE_VALUE>(), v(v) {}
      
  public:

    using value_type = VALUE;
      
    valued_graph_element(const valued_graph_element&)            = delete;
    valued_graph_element& operator=(const valued_graph_element&) = delete;
    
    VALUE&       operator()()       {return v;}
    const VALUE& operator()() const {return v;}
  };
  
  template<typename VERTEX_VALUE, typename EDGE_VALUE>
  class valued_graph_element<void, VERTEX_VALUE, EDGE_VALUE> : public graph_element<VERTEX_VALUE, EDGE_VALUE> {
  protected:
      
    valued_graph_element() : graph_element<VERTEX_VALUE, EDGE_VALUE>() {}
      
  public:

    using value_type = void;
      
    valued_graph_element(const valued_graph_element&)            = delete;
    valued_graph_element& operator=(const valued_graph_element&) = delete;
  };


  template<typename VERTEX_VALUE, typename EDGE_VALUE>
  class vertex : public valued_graph_element<VERTEX_VALUE, VERTEX_VALUE, EDGE_VALUE> {
  private:

    friend class graph_<VERTEX_VALUE, EDGE_VALUE>;
    friend class graph<VERTEX_VALUE, EDGE_VALUE>;
    
    std::list<std::weak_ptr<edge<VERTEX_VALUE, EDGE_VALUE> > > E;
    bool garbaging = true;

    vertex(const VERTEX_VALUE& v) : valued_graph_element<VERTEX_VALUE, VERTEX_VALUE, EDGE_VALUE>(v), E() {}

    
    template<typename EDGE_FUN>
    void foreach_edge_garbaging_on(const EDGE_FUN& fun) {
      auto it =  E.begin();
      while(it != E.end()) {
	auto e = it->lock();
	if(e == nullptr || e->is_killed())
	  it = E.erase(it);
	else {
	  fun(e);
	  if(e->is_killed())
	    it = E.erase(it);
	  else
	    ++it;
	}
      }
    }
    
    template<typename EDGE_FUN>
    void foreach_edge_garbaging_off(const EDGE_FUN& fun) {
      for(auto& wp : E) 
	if(auto e = wp.lock(); e && !(e->is_killed()))
	  fun(e);
    }
      
  public:
    
    using edge_val_type = EDGE_VALUE;
    using ref_edge_type = std::shared_ptr<edge<VERTEX_VALUE, EDGE_VALUE> >;

    vertex()                         = delete;
    vertex(const vertex&)            = delete;
    vertex& operator=(const vertex&) = delete;

    template<typename EDGE_FUN>
    void foreach_edge(const EDGE_FUN& fun) {
      if(garbaging) foreach_edge_garbaging_on(fun);
      else          foreach_edge_garbaging_off(fun);
    }
  };

  template<typename VERTEX_VALUE, typename EDGE_VALUE>
  class edge : public valued_graph_element<EDGE_VALUE, VERTEX_VALUE, EDGE_VALUE> {
  private:

    friend class graph<VERTEX_VALUE, EDGE_VALUE>;
    friend class graph_<VERTEX_VALUE, EDGE_VALUE>;
    
    std::weak_ptr<vertex<VERTEX_VALUE, EDGE_VALUE> > v1, v2;
    
    edge(const std::shared_ptr<vertex<VERTEX_VALUE, EDGE_VALUE> >& v1, const std::shared_ptr<vertex<VERTEX_VALUE, EDGE_VALUE> >& v2, const EDGE_VALUE& v) : valued_graph_element<EDGE_VALUE, VERTEX_VALUE, EDGE_VALUE>(v), v1(v1), v2(v2) {}
      
  public:
    
    using vertex_val_type = VERTEX_VALUE;
      
    edge()                       = delete;
    edge(const edge&)            = delete;
    edge& operator=(const edge&) = delete;
    
    auto extremities() {
      auto res = std::make_pair(v1.lock(), v2.lock());
      if(res.first != nullptr && res.first->is_killed())
	res.first = nullptr;
      if(res.second != nullptr && res.second->is_killed())
	res.second = nullptr;
      return res;
    }
  };

  template<typename VERTEX_VALUE>
  class edge<VERTEX_VALUE, void> : public valued_graph_element<void, VERTEX_VALUE, void> {
  private:

    friend class graph_<VERTEX_VALUE, void>;
    friend class graph<VERTEX_VALUE, void>;
    
    std::weak_ptr<vertex<VERTEX_VALUE, void> > v1, v2;
    
    edge(const std::shared_ptr<vertex<VERTEX_VALUE, void> >& v1, const std::shared_ptr<vertex<VERTEX_VALUE, void> >& v2) : valued_graph_element<void, VERTEX_VALUE, void>(), v1(v1), v2(v2) {}
      
  public:
    
    using vertex_val_type = VERTEX_VALUE;
      
    edge()                       = delete;
    edge(const edge&)            = delete;
    edge& operator=(const edge&) = delete;
    
    auto extremities() {
      auto res = std::make_pair(v1.lock(), v2.lock());
      if(res.first != nullptr && res.first->is_killed())
	res.first = nullptr;
      if(res.second != nullptr && res.second->is_killed())
	res.second = nullptr;
      return res;
    }
  };
  
  
  template<typename VERTEX_VALUE, typename EDGE_VALUE>
  class graph_ {
  public:

    using vertex_type       = vertex<VERTEX_VALUE, EDGE_VALUE>;
    using vertex_value_type = VERTEX_VALUE;
    using ref_vertex        = std::shared_ptr<vertex_type>;

    using edge_type         = edge<VERTEX_VALUE, EDGE_VALUE>;
    using edge_value_type   = EDGE_VALUE;
    using ref_edge          = std::shared_ptr<edge_type>;

    using wref_vertex       = std::weak_ptr<vertex_type>;
    using wref_edge         = std::weak_ptr<edge_type>;
    using vertices          = std::list<ref_vertex>;
    using edges             = std::list<ref_edge>;
    using weak_vertices     = std::list<wref_vertex>;
    using weak_edges        = std::list<wref_edge>;
    
  private:

    vertices V;

  protected :
    
    edges    E;

  private :

    bool garbaging = true;
    
    template<typename VERTEX_FUN>
    void foreach_vertex_garbaging_on(const VERTEX_FUN& fun) {
      auto it =  V.begin();
      while(it != V.end()) {
	auto& v = *it;
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
    
    template<typename VERTEX_FUN>
    void foreach_vertex_garbaging_off(const VERTEX_FUN& fun) {
      for(auto& v : V)
	if(v && !(v->is_killed()))
	  fun(v);
    }

    template<typename EDGE_FUN>
    void foreach_edge_garbaging_on(const EDGE_FUN& fun) {
      auto it =  E.begin();
      while(it != E.end()) {
	auto& e = *it;
	if(e == nullptr || e->is_killed())
	  it = E.erase(it);
	else {
	  fun(e);
	  if(e->is_killed())
	    it = E.erase(it);
	  else
	    ++it;
	}
      }
    }

    template<typename EDGE_FUN>
    void foreach_edge_garbaging_off(const EDGE_FUN& fun) {
      for(auto& e : E)
	if(e && !(e->is_killed()))
	  fun(e);
    }
    

  public:
    /**
     * This is the support of minimal heaps for shortest path computation.
     */
    std::vector<ref_vertex> heap;

    graph_()                         = default;
    graph_(const graph_&)            = delete;
    graph_& operator=(const graph_&) = delete;

    /**
     * Counts the number of vertices (with internal foreach)
     */
    unsigned int nb_vertices() {
      unsigned int res = 0;
      this->foreach_vertex([&res](const ref_vertex&) {++res;});
      return res;
    }

    /**
     * Counts the number of edges (with internal foreach)
     */
    unsigned int nb_edges() {
      unsigned int res = 0;
      this->foreach_edge([&res](const ref_edge& ref_e) {
	  auto extr_pair = ref_e->extremities();           
	  if(vq3::invalid_extremities(extr_pair)) {
	    ref_e->kill();
	    return;
	  }
	  ++res;
	});
      return res;
    }
    
    /**
     * Creates a new vertex in the graph
     */
    ref_vertex operator+=(const vertex_value_type& v) {
      auto res = ref_vertex(new vertex_type(v)); // no std::make_shared since constructor is private.
      V.push_back(res);
      return res;
      return nullptr;
    }

    /**
     * This function do not modify the graph, so it is thread-safe.
     */
    ref_edge get_edge(const ref_vertex& v1, const ref_vertex& v2) const {
      ref_vertex v, vv;
      
      if(v1->E.size() > v2->E.size()) {
	v  = v2;
	vv = v1;
      }
      else {
	v  = v1;
	vv = v2;
      }

      for(auto& we : v->E) {
	auto e = we.lock();
	if(e != nullptr) {
	  auto ref_v1 = e->v1.lock();
	  if(ref_v1 != nullptr) {
	    auto ref_v2 = e->v2.lock();
	    if(ref_v2 != nullptr && (ref_v1 == vv || ref_v2 == vv))
	      return e;
	  }
	}
      }

      return nullptr;
    }
    
    template<typename VERTEX_FUN>
    void foreach_vertex(const VERTEX_FUN& fun) {
      if(garbaging) foreach_vertex_garbaging_on(fun);
      else          foreach_vertex_garbaging_off(fun);
    }

    template<typename EDGE_FUN>
    void foreach_edge(const EDGE_FUN& fun) {
      if(garbaging) foreach_edge_garbaging_on(fun);
      else          foreach_edge_garbaging_off(fun);
    }

    /**
     * Memory garbaging is not thread safe. You have to lock it before
     * splitting your graph computations into threads.
     */
    void garbaging_lock() {
      garbaging = false;
      foreach_vertex_garbaging_on([](const ref_vertex& ref_v) {ref_v->garbaging = false;});
    }
    
    /**
     * Memory garbaging is not thread safe. This unlocks the garbaging
     * for single-thread graph computations. Indeed, this is usefull
     * when garbaging has been locked previously before splitting the
     * computation into threads. In this case, call this method when
     * you're back to sequential processing.
     */
    void garbaging_unlock() {
      garbaging = true;
      foreach_vertex_garbaging_on([](const ref_vertex& ref_v) {ref_v->garbaging = true;});
    }

    /**
     * This cleans up pending kills. It may be useless since garbaging
     * is performed all along graph computations (e.g. if vertices and
     * edges are regularly accessed for display in a sequential
     * section). Nevertheless, if garbaging has been suspended for
     * parallel computing (see garbaging_lock and garbaging_unlock),
     * calling garbaging_force() may be usefull. The method
     * garbaging_force() should be called outside parallel sections
     * (i.e. it is not thread safe).
     */
    void garbaging_force() {
      foreach_edge_garbaging_on([](const ref_edge& ref_e) {
	  if(auto extr_pair = ref_e->extremities(); vq3::invalid_extremities(extr_pair)) 
	    ref_e->kill();
	});
      
      foreach_vertex_garbaging_on([](const ref_vertex& ref_v) {
	  ref_v->foreach_edge_garbaging_on([](const ref_edge&) {});
	});
    }

    
    
  };

   
  template<typename VERTEX_VALUE,
	   typename EDGE_VALUE>
  class graph : public graph_<VERTEX_VALUE, EDGE_VALUE> {

  public:

    // I use this since I have ADL issues....
    std::function<void (std::ostream&, const VERTEX_VALUE&)> vertex_to_stream;
    std::function<void (std::ostream&, const EDGE_VALUE&)>   edge_to_stream;
    std::function<void (std::istream&, VERTEX_VALUE&)>       vertex_from_stream;
    std::function<void (std::istream&, EDGE_VALUE&)>         edge_from_stream;

    graph()                        = default;
    graph(const graph&)            = delete;
    graph& operator=(const graph&) = delete;
    
    /**
     * Add an edge. The vertex references have to be vertices that
     * actually belongs to the graph.
     */
    typename graph_<VERTEX_VALUE, EDGE_VALUE>::ref_edge connect(const typename graph_<VERTEX_VALUE, EDGE_VALUE>::ref_vertex& v1, const typename graph_<VERTEX_VALUE, EDGE_VALUE>::ref_vertex& v2, const typename graph_<VERTEX_VALUE, EDGE_VALUE>::edge_value_type& v) {
      auto res = typename graph_<VERTEX_VALUE, EDGE_VALUE>::ref_edge(new typename graph_<VERTEX_VALUE, EDGE_VALUE>::edge_type(v1, v2, v)); // no std::make_shared since constructor is private.
      v1->E.push_back(res);
      v2->E.push_back(res);
      this->E.push_back(res);
      return res;
    }
    
    /**
     * Add a default-valued edge. The vertex references have to be vertices that
     * actually belongs to the graph.
     */
    typename graph_<VERTEX_VALUE, EDGE_VALUE>::ref_edge connect(const typename graph_<VERTEX_VALUE, EDGE_VALUE>::ref_vertex& v1, const typename graph_<VERTEX_VALUE, EDGE_VALUE>::ref_vertex& v2) {
      auto res = typename graph_<VERTEX_VALUE, EDGE_VALUE>::ref_edge(new typename graph_<VERTEX_VALUE, EDGE_VALUE>::edge_type(v1, v2, typename graph_<VERTEX_VALUE, EDGE_VALUE>::edge_value_type())); // no std::make_shared since constructor is private.
      v1->E.push_back(res);
      v2->E.push_back(res);
      this->E.push_back(res);
      return res;
    }

    /**
     * This writes the graph into a stream. You can use the << operator as well.
     */
    void write(std::ostream& os) {
      std::map<typename graph_<VERTEX_VALUE, EDGE_VALUE>::ref_vertex, unsigned int> vertex_idf;
      unsigned int idf = 0;
      
      os << this->nb_vertices() << std::endl;
      this->foreach_vertex([&idf, &vertex_idf, &os, this](auto& ref_v){
	  vertex_idf[ref_v] = idf++;
	  vertex_to_stream(os, (*ref_v)());
	  os << std::endl;
	});

      os << this->nb_edges() << std::endl;
      this->foreach_edge([&vertex_idf, &os, this](auto& ref_e){
	  auto extr_pair = ref_e->extremities();
	  os << vertex_idf[extr_pair.first] << ' ' <<  vertex_idf[extr_pair.second] << std::endl;
	  edge_to_stream(os, (*ref_e)());
	  os << std::endl;
	});
    }
     

    /**
     * This reads the graph from a stream. You can use the >> operator as well.
     */
    void read(std::istream& is) {
      this->foreach_vertex([](auto ref_v){ref_v->kill();});
      this->foreach_edge  ([](auto ref_e){ref_e->kill();});

      unsigned int nb;
      char c;
      typename graph_<VERTEX_VALUE, EDGE_VALUE>::vertex_value_type v {};
      typename graph_<VERTEX_VALUE, EDGE_VALUE>::edge_value_type e {};
      
      is >> nb; is.get(c);
      std::vector<typename graph_<VERTEX_VALUE, EDGE_VALUE>::ref_vertex> vtx(nb);
      for(auto& ref_v : vtx) {
	vertex_from_stream(is, v);
	is.get(c);
	ref_v = ((*this) += v);
      }
      
      is >> nb; is.get(c);
      for(unsigned int i=0; i < nb; ++i) {
      	unsigned int v1, v2;
      	is >> v1 >> v2; is.get(c);
	edge_from_stream(is, e);
      	is.get(c);
	
      	this->connect(vtx[v1], vtx[v2], e);
      }
    }
  };
  
  template<typename VERTEX_VALUE>
  class graph<VERTEX_VALUE, void> : public graph_<VERTEX_VALUE, void> {

  private:


    // I use this since I have ADL issues....
    std::function<void (std::ostream&, const VERTEX_VALUE&)> vertex_to_stream;
    std::function<void (std::istream&, VERTEX_VALUE&)>       vertex_from_stream;

  public:
    

    graph()                        = default;
    graph(const graph&)            = delete;
    graph& operator=(const graph&) = delete;
    
    /**
     * Add an edge. The vertex references have to be vertices that
     * actually belongs to the graph.
     */
    typename graph_<VERTEX_VALUE, void>::ref_edge connect(const typename graph_<VERTEX_VALUE, void>::ref_vertex& v1, const typename graph_<VERTEX_VALUE, void>::ref_vertex& v2) {
      auto res = typename graph_<VERTEX_VALUE, void>::ref_edge(new typename graph_<VERTEX_VALUE, void>::edge_type(v1, v2)); // no std::make_shared since constructor is private.
      v1->E.push_back(res);
      v2->E.push_back(res);
      this->E.push_back(res);
      return res;
    }
    
    /**
     * This writes the graph into a stream. You can use the << operator as well.
     */
    void write(std::ostream& os) {
      std::map<typename graph_<VERTEX_VALUE, void>::ref_vertex, unsigned int> vertex_idf;
      unsigned int idf = 0;
      
      os << this->nb_vertices() << std::endl;
      this->foreach_vertex([&idf, &vertex_idf, &os, this](auto& ref_v){
	  vertex_idf[ref_v] = idf++;
	  vertex_to_stream(os, (*ref_v)());
	  os << std::endl;
	});

      os << this->nb_edges() << std::endl;
      this->foreach_edge([&vertex_idf, &os](auto& ref_e){
	  auto extr_pair = ref_e->extremities();
	  os << vertex_idf[extr_pair.first] << ' ' <<  vertex_idf[extr_pair.second] << std::endl;
	});
    }
     

    /**
     * This reads the graph from a stream. You can use the >> operator as well.
     */
    void read(std::istream& is) {
      this->foreach_vertex([](auto ref_v){ref_v->kill();});
      this->foreach_edge  ([](auto ref_e){ref_e->kill();});

      unsigned int nb;
      char c;
      typename graph_<VERTEX_VALUE, void>::vertex_value_type v {};
      
      is >> nb; is.get(c);
      std::vector<typename graph_<VERTEX_VALUE, void>::ref_vertex> vtx(nb);
      for(auto& ref_v : vtx) {
	vertex_from_stream(is, v);
	is.get(c);
	ref_v = ((*this) += v);
      }
      
      is >> nb; is.get(c);
      for(unsigned int i=0; i < nb; ++i) {
      	unsigned int v1, v2;
      	is >> v1 >> v2; is.get(c);
      	this->connect(vtx[v1], vtx[v2]);
      }
    }
  };

  
  template<typename VERTEX_VALUE, typename EDGE_VALUE>
  std::ostream& operator<<(std::ostream& os, graph<VERTEX_VALUE, EDGE_VALUE>& g) {
    g.write(os);
    return os;
  }
  
  template<typename VERTEX_VALUE, typename EDGE_VALUE>
  std::istream& operator>>(std::istream& is, graph<VERTEX_VALUE, EDGE_VALUE>& g) {
    g.read(is);
    return is;
  }
}
