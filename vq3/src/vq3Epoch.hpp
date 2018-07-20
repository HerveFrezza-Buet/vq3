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

#include <future>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <iterator>

#include <vq3Topology.hpp>
#include <vq3Utils.hpp>

namespace vq3 {
  
  namespace concept {

    /**
     * This stores data collected about a vertex (handling a prototype) during an epoch.
     */
    template<typename SAMPLE_TYPE, typename VERTEX_VALUE_TYPE, typename PROTOTYPE_TYPE>
    struct EpochData {

      using sample_type       = SAMPLE_TYPE;
      using vertex_value_type = VERTEX_VALUE_TYPE;
      using prototype_type    = PROTOTYPE_TYPE;

      EpochData();

      /**
       * This notifies that the prototype handled by the vertex is the closest.
       * @param sample The current sample that is submitted.
       * @param closest_distance The distance from the vertex prototype and the sample (which is the cmallest among all vertices in the graph).
       */
      void notify_closest(const sample_type& sample, double closest_distance);

      /**
       * This notifies that the prototype should be updated.
       * @param sample The current sample that is submitted.
       * @param topo_coef This is the WTM coefficient related to edge-distance between the BMU and the prototype.
       */
      void notify_wtm_update(const sample_type& sample, double topo_coef);

      /**
       * This notifies that the prototype should be updated. No coefficient is provides, this is used in a WTA context.
       * @param sample The current sample that is submitted.
       */
      void notify_wta_update(const sample_type& sample);

      /**
       * This is learning, at the end of an epoch.
       * @prototype This is the prototype of the vertex, passed by reference in order do be modified by the call.
       */
      void set_prototype(prototype_type& prototype);

      /**
       * This is prototype value updating, at the end of an epoch.
       * @vertex_vale This is the vertex value of the vertex, passed by reference in order do be modified by the call.
       */
      void set_content(vertex_value_type& vertex_value);
 
      /**
       * This combine the current data with another one (computed in another thread) for the same vertex.
       */
      EpochData& operator+=(const EpochData& other);
    };
  }
  
  namespace epoch {

    /**
     * Computes the distortion for several epoch data. The epoch data given by iterators may implement the vq3::epoch::data::bmu behavior.
     */
    template<typename EPOCH_DATA_ITER>
    auto distortion(const EPOCH_DATA_ITER& begin, const EPOCH_DATA_ITER& end) {
      if(begin == end)
	throw std::runtime_error("vq3::epoch::distortion : empty collection");
      auto it = begin;
      auto res = it->vq3_bmu_accum;
      for(++it; it != end; ++it)
	res += it->vq3_bmu_accum;
      return res;
    }
    
    namespace data {

      /**
       * Root data type, it fits vq3::concept::EpochData. It only notifies the types.
       */
      template<typename SAMPLE_TYPE, typename VERTEX_VALUE_TYPE, typename PROTOTYPE_TYPE>
      struct none {
	using sample_type       = SAMPLE_TYPE;
	using vertex_value_type = VERTEX_VALUE_TYPE;
	using prototype_type    = PROTOTYPE_TYPE;
	
	none() = default;

	void  notify_closest    (const sample_type&, double) {}
	//!< nop. 
	
	void  notify_wtm_update (const sample_type&, double) {}
	//!< nop. 
	
	void  notify_wta_update (const sample_type&)         {}
	//!< nop.

	void set_prototype(prototype_type& prototype) {}
	//!< nop.
	
	void set_content(vertex_value_type& vertex_value) {}
	//!< nop.
	
	none& operator+=(const none& other)                  {return *this;}
      };

      /**
       * This collects winner-take-all data and computes the
       * prototype. There is no check of stability (rather use
       * wta_ckeck for that).
       */ 
      template<typename MOTHER>
      struct wta : MOTHER {
	using sample_type = typename MOTHER::sample_type;
	using vertex_value_type = typename MOTHER::vertex_value_type;
	using prototype_type = typename MOTHER::prototype_type;
	
	utils::accum<sample_type, double> vq3_wta_accum; //!< The sum of samples in the Voronoï cell
	
	wta() = default;
	
	void notify_closest(const sample_type& sample, double dist) {
	  this->MOTHER::notify_closest(sample, dist);
	}
	//!< nop. 
	
	void notify_wtm_update(const sample_type& sample, double coef) {
	  this->MOTHER::notify_wtm_update(sample, coef);
	}
	//!< nop. 
	
	void notify_wta_update(const sample_type& sample) {
	  this->MOTHER::notify_wta_update(sample);
	  vq3_wta_accum += sample;
	}
	//!< accum += sample

	void set_prototype(prototype_type& prototype) {
	  this->MOTHER::set_prototype(prototype);
	  if(this->vq3_wta_accum.nb != 0)
	    prototype = this->vq3_wta_accum.template average<double>();
	}
	//!< prototype = accum.average();
	
	void set_content(vertex_value_type& vertex_value) {
	  this->MOTHER::set_content(vertex_value);
	}
	//!< nop.

	wta<MOTHER>& operator+=(const wta<MOTHER>& arg) {
	  this->MOTHER::operator+=(arg);
	  vq3_wta_accum += arg.vq3_wta_accum;
	  return *this;
	}
      };


      /**
       * This collects winner-take-most data and computes the prototype.
       */ 
      template<typename MOTHER>
      struct wtm : MOTHER {
	using sample_type = typename MOTHER::sample_type;
	using vertex_value_type = typename MOTHER::vertex_value_type;
	using prototype_type = typename MOTHER::prototype_type;
	
	utils::accum<sample_type, double> vq3_wtm_accum; //!< The weighted sum of samples in the Voronoï cell
	
	
	wtm() = default;
	
	void notify_closest(const sample_type& sample, double dist) {
	  this->MOTHER::notify_closest(sample, dist);
	}
	//!< nop. 
	
	
	void notify_wtm_update(const sample_type& sample, double coef) {
	  this->MOTHER::notify_wtm_update(sample, coef);
	  vq3_wtm_accum.increment(coef, sample);
	}
	//!< accum += coef * sample; 
	
	
	void notify_wta_update(const sample_type& sample) {
	  this->MOTHER::notify_wta_update(sample);
	}
	//!< nop. 

	void set_prototype(prototype_type& prototype) {
	  this->MOTHER::set_prototype(prototype);
	  if(vq3_wtm_accum.nb != 0)
	    prototype = vq3_wtm_accum.template average<double>();
	}
	//!< prototype = accum.average(); 
	
	void set_content(vertex_value_type& vertex_value) {
	  this->MOTHER::set_content(vertex_value);
	}
	//!< nop.

	wtm<MOTHER>& operator+=(const wtm<MOTHER>& arg) {
	  this->MOTHER::operator+=(arg);
	  vq3_wtm_accum += arg.vq3_wtm_accum;
	  return *this;
	}
      };

      
      /**
       * This collects the computation related to the bmu (best matching unit).
       */ 
      template<typename MOTHER>
      struct bmu : MOTHER {
	using sample_type = typename MOTHER::sample_type;
	using vertex_value_type = typename MOTHER::vertex_value_type;
	using prototype_type = typename MOTHER::prototype_type;
	
	utils::accum<double, double> vq3_bmu_accum; //!< The sum of distortions when the nodes is the BMU. 

	bmu() = default;

	void notify_closest(const sample_type& sample, double dist) {
	  this->MOTHER::notify_closest(sample, dist);
	  vq3_bmu_accum += dist;
	}
	//!< acum += distance;
	
	void notify_wtm_update(const sample_type& sample, double coef) {
	  this->MOTHER::notify_wtm_update(sample, coef);
	}
	//!< nop. 
	
	void notify_wta_update(const sample_type& sample) {
	  this->MOTHER::notify_wta_update(sample);
	}
	//!< nop. 

	void set_prototype(prototype_type& prototype) {
	  this->MOTHER::set_prototype(prototype);
	}
	//!< nop.
	
	void set_content(vertex_value_type& vertex_value) {
	  this->MOTHER::set_content(vertex_value);
	}
	//!< nop.

	bmu<MOTHER>& operator+=(const bmu<MOTHER>& arg) {
	  this->MOTHER::operator+=(arg);
	  vq3_bmu_accum += arg.vq3_bmu_accum;
	  return *this;
	}
      };

      /**
       * This data stores both the current and the
       * previous value that were used for updating the prototype,
       * for the sake of convergence checking.
       */ 
      template<typename MOTHER>
      struct delta : MOTHER {
	
	using sample_type = typename MOTHER::sample_type;
	using vertex_value_type = typename MOTHER::vertex_value_type;
	using prototype_type = typename MOTHER::prototype_type;
	
	prototype_type wq3_wta_previous_prototype; //!< The vertex value computed at the previous step.
	prototype_type wq3_wta_current_prototype;  //!< The vertex value computed at the current step.
	
	delta() = default;
	
	void notify_closest(const sample_type& sample, double dist) {
	  this->MOTHER::notify_closest(sample, dist);
	}
	//!< nop. 
	
	void notify_wtm_update(const sample_type& sample, double coef) {
	  this->MOTHER::notify_wtm_update(sample, coef);
	}
	//!< nop. 
	
	void notify_wta_update(const sample_type& sample) {
	  this->MOTHER::notify_wta_update(sample);
	}
	//!< accum += sample
	
	void set_prototype(prototype_type& prototype) {
	  wq3_wta_previous_prototype = prototype;
	  this->MOTHER::set_prototype(prototype);
	  wq3_wta_current_prototype = prototype;
	}
	//!< previous_proto = prototype; update prototype; current_proto = prototype;
	
	void set_content(vertex_value_type& vertex_value) {
	  this->MOTHER::set_content(vertex_value);
	}
	//!< nop.

	delta<MOTHER>& operator+=(const delta<MOTHER>& arg) {
	  this->MOTHER::operator+=(arg);
	  return *this;
	}
      };
      
      namespace online {
	/**
	 * This collects the computation related to the bmu (best
	 * matching unit), as vq3::epoch::data::bmu does. Moreover, it
	 * supposes that the vertex value has a vq3_online_mean_std
	 * attribute (see vq3::decorator::online::mean_std) and uses it
	 * to estimate the mean and standard deviation of the sum of
	 * distortions obtained at each pass.
	 */ 
	template<typename MOTHER>
	struct bmu_mean_std : MOTHER {
	  using sample_type = typename MOTHER::sample_type;
	  using vertex_value_type = typename MOTHER::vertex_value_type;
	  using prototype_type = typename MOTHER::prototype_type;
	
	  vq3::utils::accum<double, double> vq3_bmu_accum; //!< The sum of distortions when the nodes is the BMU. 

	  bmu_mean_std() = default;

	  void notify_closest(const sample_type& sample, double dist) {
	    this->MOTHER::notify_closest(sample, dist);
	    vq3_bmu_accum += dist;
	  }
	  //!< acum += distance;
	
	  void notify_wtm_update(const sample_type& sample, double coef) {
	    this->MOTHER::notify_wtm_update(sample, coef);
	  }
	  //!< nop. 
	
	  void notify_wta_update(const sample_type& sample) {
	    this->MOTHER::notify_wta_update(sample);
	  }
	  //!< nop. 

	  void set_prototype(prototype_type& prototype) {
	    this->MOTHER::set_prototype(prototype);
	  }
	  //!< nop.

	  void set_content(vertex_value_type& vertex_value) {
	    this->MOTHER::set_content(vertex_value);
	    vertex_value.vq3_online_mean_std += vq3_bmu_accum.value;
	  }
	  //!< vertex_value.vq3_online_mean_std += vq3_bmu_accum.value

	  bmu_mean_std<MOTHER>& operator+=(const bmu_mean_std<MOTHER>& arg) {
	    this->MOTHER::operator+=(arg);
	    vq3_bmu_accum += arg.vq3_bmu_accum;
	    return *this;
	  }
	};
      }
    }

    namespace chl {
      template<typename GRAPH>
      class Processor {
      private:

	using ref_vertex = typename GRAPH::ref_vertex;
	using ref_edge   = typename GRAPH::ref_edge;
	using edge       = typename GRAPH::edge_value_type;
      
	GRAPH& g;

	struct refpair {
	  ref_vertex first, second;
	  refpair()                          = default;
	  refpair(const refpair&)            = default;
	  refpair& operator=(const refpair&) = default;
	  refpair(refpair&&)                 = default;
	  refpair& operator=(refpair&&)      = default;
	  
	  refpair(const std::pair<ref_vertex, ref_vertex>& p)
	    : first(p.first), second(p.second) {
	    if(second < first)
	      std::swap(first, second);
	  }

	  bool operator<(const refpair& other) const {
	    return (first < other.first) || ((first == other.first) && (second < other.second));
	  }
	};

	struct data {
	  std::set<ref_edge> survivors;
	  std::set<refpair>  newedges;
	  data()                       = default;
	  data(const data&)            = default;
	  data& operator=(const data&) = default;
	  data(data&&)                 = default;
	  data& operator=(data&&)      = default;
	};
	  
      public:
      
	Processor(GRAPH& g) : g(g) {}
	Processor()                            = delete;
	Processor(const Processor&)            = default;
	Processor(Processor&&)                 = default;
	Processor& operator=(const Processor&) = default;
	Processor& operator=(Processor&&)      = default;

	/**
	 * This processes Competitive Hebbian learning, adding or removing edges in the graph.
	 * @return true if the process has modified the graph topology. 
	 */
	template<typename ITERATOR, typename SAMPLE_OF, typename PROTOTYPE_OF_VERTEX_VALUE, typename DISTANCE>
	bool update_edges(unsigned int nb_threads,
			  const ITERATOR& samples_begin, const ITERATOR& samples_end, const SAMPLE_OF& sample_of,
			  const PROTOTYPE_OF_VERTEX_VALUE& prototype_of, const DISTANCE& distance,
			  const edge& value_for_new_edges) {
	  if(g.nb_vertices() < 2)
	    return false;
	  auto iters = utils::split(samples_begin, samples_end, nb_threads);
	  std::vector<std::future<data> > futures;
	  auto out = std::back_inserter(futures);

	  for(auto& begin_end : iters) 
	    *(out++) = std::async(std::launch::async,
				  [begin_end, this, &sample_of, &distance]() {
				    data res;
				    for(auto it = begin_end.first; it != begin_end.second; ++it) {
				      auto two = utils::two_closest(g, sample_of(*it), distance);
				      auto ref_e = g.get_edge(two.first, two.second);
				      if(ref_e == nullptr)
					res.newedges.emplace(two);
				      else
					res.survivors.emplace(ref_e);
				    }
				    return res;
				  });
	  
	  std::vector<ref_edge> survivors;
	  std::vector<refpair>  newedges;
	  for(auto& f : futures) {
	    auto d = f.get();
	    
	    std::vector<ref_edge> s;
	    auto outs = std::back_inserter(s);
	    std::merge(survivors.begin(),   survivors.end(),
		       d.survivors.begin(), d.survivors.end(),
		       outs);
	    
	    std::vector<refpair>  n;
	    auto outn = std::back_inserter(n);
	    std::merge(newedges.begin(),   newedges.end(),
		       d.newedges.begin(), d.newedges.end(),
		       outn);
	    std::swap(s, survivors);
	    std::swap(n, newedges);
	  }

	  bool one_kill = false;

	  // Let us remove non surviving edges.
	  utils::clear_edge_tags(g, true);
	  for(auto& ref_e : survivors) (*ref_e)().vq3_tag = false;
	  g.foreach_edge([&one_kill](const ref_edge& ref_e) {
	      if((*ref_e)().vq3_tag) {
		ref_e->kill();
		one_kill = true;
	      }
	    });

	  // Let us add the new edges.
	  for(auto& p : newedges)
	    g.connect(p.first, p.second, value_for_new_edges);

	  return newedges.size() != 0 || one_kill;
	}
      };
    
      template<typename GRAPH>
      auto processor(GRAPH& g) {return Processor<GRAPH>(g);}
    }
    
    namespace wta {

      template<typename GRAPH>
      class Processor {
      public:
	
	using ref_vertex    = typename GRAPH::ref_vertex;
	using vertices_type = utils::Vertices<ref_vertex>;
	using index_type    = typename vertices_type::index_type;
	
      private:

      
	GRAPH& g;
	vertices_type& vertices;
      
      public:
      
	Processor(GRAPH& g, vertices_type& vertices) : g(g), vertices(vertices) {}
	Processor()                            = delete;
	Processor(const Processor&)            = default;
	Processor(Processor&&)                 = default;
	Processor& operator=(const Processor&) = default;
	Processor& operator=(Processor&&)      = default;


	/**
	 * @return A vector, for each prototype index, of the epoch data.
	 */
	template<typename EPOCH_DATA, typename ITERATOR, typename SAMPLE_OF, typename PROTOTYPE_OF_VERTEX_VALUE, typename DISTANCE>
	auto update_prototypes(unsigned int nb_threads, const ITERATOR& samples_begin, const ITERATOR& samples_end, const SAMPLE_OF& sample_of, const PROTOTYPE_OF_VERTEX_VALUE& prototype_of, const DISTANCE& distance) {
	  auto iters = utils::split(samples_begin, samples_end, nb_threads);
	  std::vector<std::future<std::vector<EPOCH_DATA> > > futures;
	  auto out = std::back_inserter(futures);

	  for(auto& begin_end : iters) 
	    *(out++) = std::async(std::launch::async,
				  [begin_end, this, &sample_of, &distance]() {
				    std::vector<EPOCH_DATA> data(vertices.size());
				    for(auto it = begin_end.first; it != begin_end.second; ++it) {
				      double min_dist;
				      const auto&  sample = sample_of(*it);
				      auto        closest = utils::closest(g, sample, distance, min_dist);
				      if(closest != nullptr) {
					auto&             d = data[vertices(closest)];
					d.notify_closest(sample, min_dist);
					d.notify_wta_update(sample);
				      }
				    }
				    return data;
				  });

	  auto fit = futures.begin();
	  auto data0 = (fit++)->get();
	  auto b0 = data0.begin();
	  auto e0 = data0.end();
	  for(unsigned int thread_id = 1; thread_id < nb_threads; ++thread_id) {
	    auto datai = (fit++)->get();
	    auto bi = datai.begin();
	    for(auto b = b0; b != e0; ++b, ++bi)
	      (*b) += *bi;
	  }

	  unsigned int idx = 0;
	  for(auto&  d : data0) {
	    auto& value = (*(vertices(idx++)))();
	    d.set_prototype(prototype_of(value));
	    d.set_content(value);
	  }

	  return data0;
	}
      };
    
      template<typename GRAPH>
      auto processor(GRAPH& g, utils::Vertices<typename GRAPH::ref_vertex>& vertices) {return Processor<GRAPH>(g, vertices);}
    }
    
    namespace wtm {

      template<typename GRAPH>
      class Processor {
      public:
	

	using ref_vertex    = typename GRAPH::ref_vertex;
	using vertices_type = utils::Vertices<ref_vertex>;
	using index_type    = typename vertices_type::index_type;
	
      private:
      
	GRAPH& g;
	vertices_type& vertices;
	std::map<ref_vertex, std::list<topo::Info<index_type, double> > > neighborhood_table;
      
      public:
      
	
	Processor(GRAPH& g, vertices_type& vertices) : g(g), vertices(vertices), neighborhood_table() {}
	Processor()                            = delete;
	Processor(const Processor&)            = default;
	Processor(Processor&&)                 = default;
	Processor& operator=(const Processor&) = default;
	Processor& operator=(Processor&&)      = default;

	/**
	 * Do not forget to update the external vertices topology first.
	 * @param voed A function providing a value (double >= 0) according to the number of edges (unsigned int) separating a vertex in the neighborhood from the central vertex.
	 * @param max_dist The maximal distance considered. 0 means "no limit".
	 * @param min_val if voed(dist) < min_val, the node is not included in the neighborhood.
	 */
	template<typename VALUE_OF_EDGE_DISTANCE>
	void update_topology(const VALUE_OF_EDGE_DISTANCE& voed, unsigned int max_dist, double min_val) {
	  neighborhood_table = utils::make_vertex_table(g,
							[this, &voed, max_dist, min_val](const typename GRAPH::ref_vertex& ref_v) {
							  utils::clear_vertex_tags(g, false); 
							  return topo::edge_based_neighborhood(vertices, ref_v, voed, max_dist, min_val); 
							});
	}


	/**
	 * @return A vector, for each prototype index, of the epoch data.
	 */
	template<typename EPOCH_DATA, typename ITERATOR, typename SAMPLE_OF, typename PROTOTYPE_OF_VERTEX_VALUE, typename DISTANCE>
	auto update_prototypes(unsigned int nb_threads, const ITERATOR& samples_begin, const ITERATOR& samples_end, const SAMPLE_OF& sample_of, const PROTOTYPE_OF_VERTEX_VALUE& prototype_of, const DISTANCE& distance) {
	  auto iters = utils::split(samples_begin, samples_end, nb_threads);
	  std::vector<std::future<std::vector<EPOCH_DATA> > > futures;
	  auto out = std::back_inserter(futures);

	  for(auto& begin_end : iters) 
	    *(out++) = std::async(std::launch::async,
				  [begin_end, this, &sample_of, &distance, size = vertices.size()]() {
				    std::vector<EPOCH_DATA> data(size);
				    for(auto it = begin_end.first; it != begin_end.second; ++it) {
				      double min_dist;
				      const auto&  sample = sample_of(*it);
				      auto        closest = utils::closest(g, sample, distance, min_dist);
				      if(closest != nullptr) {
					auto&  neighborhood = neighborhood_table[closest];
					data[neighborhood.begin()->index].notify_closest(sample, min_dist);
					for(auto& info : neighborhood) data[info.index].notify_wtm_update(sample, info.value);
				      }
				    }
				    return data;
				  });
	  
	  auto fit = futures.begin();
	  auto data0 = (fit++)->get();
	  auto b0 = data0.begin();
	  auto e0 = data0.end();
	  for(unsigned int thread_id = 1; thread_id < nb_threads; ++thread_id) {
	    auto datai = (fit++)->get();
	    auto bi = datai.begin();
	    for(auto b = b0; b != e0; ++b, ++bi)
	      (*b) += *bi;
	  }

	  unsigned int idx = 0;
	  for(auto&  d : data0) {
	    auto& value = (*(vertices(idx++)))();
	    d.set_prototype(prototype_of(value));
	    d.set_content(value);
	  }
	  return data0;
	}
      };
    
      template<typename GRAPH>
      auto processor(GRAPH& g, utils::Vertices<typename GRAPH::ref_vertex>& vertices) {return Processor<GRAPH>(g, vertices);}
    }
  }
}
