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

#include <vq3Component.hpp>
#include <vq3Decorator.hpp>
#include <vq3Epoch.hpp>
#include <vq3GIT.hpp>
#include <vq3Graph.hpp>
#include <vq3GNGT.hpp>
#include <vq3LBG.hpp>
#include <vq3Online.hpp>
#include <vq3ShortestPath.hpp>
#include <vq3SOM.hpp>
#include <vq3ShortestPath.hpp>
#include <vq3Stats.hpp>
#include <vq3Temporal.hpp>
#include <vq3Topology.hpp>
#include <vq3Utils.hpp>

/**

   @mainpage
   
   @section overview Overview

   The vq3 library is a result of the <a href="http://interreg-grone.eu">GRONE project</a>, supported by the Interreg "Grande Région" program of the European Union's European Regional Development Fund.
   
   The vq3 library is a massive rewrinting of a former vq2 library for vector
   quantization. It complies with C++11 and later requirements, and
   integrates parallélism (multi-thread).
   
   Code examples are not in included this package, but rather in the
   vq3demo package. Indeed, the vq3 features are illustrated in that
   vq3demo package with 2D points that can be displayed thanks to
   openCV.
   
   The following describes the overall 'philosophy' underlying the
   design of vq3, for the sake of preparing the reading of vq3demo
   examples.
   
   @section graphs Graphs

   @subsection ve Vertices and edges

   The VQ algorithms handle a graph of prototypes. In vq3, you may
   start with defining your graph type as follows:
   @code
using node_value = std::string;
using edge_value = double; 
using graph      = vq3::graph<node_value, edge_value>;
   @endcode
   If no edge type is needed, you can write something like:
   @code
using graph = vq3::graph<std::string, void>;
   @endcode
   Adding a node in the graph uses the += operator. It returns a reference (it is a shared pointer).
   @code
using graph = vq3::graph<std::string, double>;
graph g; 

// first vertex.
g += "first";  

// second vertex, accessible from ref_v2.
auto ref_v2 = g += "second"; 

// the same, where auto is explicited
graph::ref_vertex ref_v3 = g += "third";

// We get the value from a ref_vertex:
std::string value_2 = (*ref_v2)();
   @endcode

   Edges are handled in a very similar manner, using shared pointers as well.
   
   @code
using graph = vq3::graph<std::string, double>;
graph g; 

auto ref_v1 = g += "first";  
auto ref_v2 = g += "second"; 

// Let us create an edge (the graph is not oriented) hosting a value (double here).
graph::ref_edge ref_e = g.connect(ref_v1, ref_v2, 3.141592654);

// We can retrieve the value as we did for the vertices.
double pi = (*ref_e)();

// We can get the extremities of the edge... they are vertex references.
auto [ref_vv1, ref_vv2] = ref_e->extremities();
   @endcode

   Last, during execution, one can ask for the removing of an edge or a vertex. This is done as follows:
   @code
ref_v->kill(); // kill a vertex from a reference.
ref_e->kill(); // kill an edge from a reference.
   @endcode

   @subsection graphit Graph iterations

   The library vq3 offers "foreach" functions in order to iterate on
   the nodes of a graph, on the edges of a graph, but also on the
   edges of a node. The iteration consists in calling a function for
   each element (node or edge), taking dynamical removals into
   account. This is why a direct acces to the collection of edges or
   nodes is not provided by vq3.

   @code
using graph = vq3::graph<std::string, double>;
graph g;

auto ref_A = g += "A";
auto ref_B = g += "B";
auto ref_C = g += "C";
auto ref_D = g += "D";

auto ref_AB = g.connect(ref_A, ref_B, 1.0);
auto ref_AC = g.connect(ref_A, ref_C, 2.0);
auto ref_AD = g.connect(ref_A, ref_D, 3.0);
auto ref_BC = g.connect(ref_B, ref_C, 4.0);
auto ref_BD = g.connect(ref_B, ref_D, 5.0);
auto ref_CD = g.connect(ref_C, ref_D, 6.0);

// Let us iterate on vertices.
std::string s1;
g.foreach_vertex([&s1](graph::ref_vertex ref_v) {s1 += (*ref_v)();});

// Let us iterate on edges. This is not mandatory, but if the graph is
// dynamical, some vertices may have been removed while the edge is
// still there. In this case, this has to be tested. This test is done
// here for the sake of illustration.
std::string s2;
double sum = 0;
g.foreach_edge([&s2, &sum](graph::ref_edge ref_e) {
  auto extr_pair = ref_e->extremities();           
  if(vq3::invalid_extremities(extr_pair)) {
    ref_e->kill();
    return;
  }
  std::string v1_value = (*(extr_pair.first))();
  std::string v2_value = (*(extr_pair.second))();
  s2  += std::string("(") + v1_value + v2_value + ")";
  sum += (*ref_e)();
});

// The edges of a vertex can be iterated as well.
sum = 0;
ref_A->foreach_edge([&sum](graph::ref_edge ref_e) {sum += (*ref_e)();});
   @endcode

   @subsection graphthread Graph, threads and garbage collector

   Many operation on graphs done by vq3 can be parallelized by using
   multi-threading (thanks to std::async mainly). Nevertheless, the
   graph is designed to be highly dynamic, so having several threads
   modifying the topology is to be avoided. 

   The rule of thumb is to keep the graph "read only" during multi-thread
   parallel sections.

   While inspecting the graph thanks to foreach functions, one can
   call the kill() method on vertices or edges. This stands for a kill
   request, not an immediate kill. Once the request is made, the graph
   behaves as if the killed edge or vertex do not exist, but the
   memory is not necessarily freed. This can happen even when reading
   the graph, for example in foreach_edge loops when bad vertices are
   detected for some edge.

   By default, the actual memory cleaning is done progressively, as
   the foreach calls inspect the graph (even when the graph is only
   read). This memory cleaning (i.e. garbaging) can be explicitly made by calling

   @code
g.garbaging_force()
   @endcode

   As the garbaging runs implicitly during graph exploration, a
   multi-thread computation, even if it only reads the graph, may
   alter the inner graph structure for elements which are waiting to
   be actually killed. To prevent thread memory issue, one should
   suspend the garbaging process during parallel computation.
   
   @code
// sequential code section
g.garbaging_lock();
// Split the computation into threads here.
...
// Threads are done, we are back in the main thread.
g.garbaging_unlock();
g.garbaging_force(); // Not mandatory, further non locked graph explorations may end up doing the job. 
   @endcode

   Last, let us stress here that vq3 processes which allow for multi-threading (they have a nb_threads parameter) may lock and unlock garbaging if they actually involve more than one thread.

   @subsection topo Topology

   We call the "topology" the knowledge of the graph structure, i.e the vertices and for each one its neighborhood. Algorithm may need to retrieve vertices from an index value. The topology table object provides this. 
   @code
graph g;
auto topology = vq3::topology::table(g);
...
topology(); // This updates the topology according to the current graph structure (required when the graph changes).
            // Only knowledge about vertices is updated.
auto ref_vertex = topology(3);          // Gets the 4th vertex from the table (constant time)
auto idx        = topology(ref_vertex); // Should return 3... (logarithmic time).
   @endcode

   Topology tables also provides the computation of neighborhoods. The neighborhood of a vertex v is a list of (value, index) pairs. Each (value, index) is the index of one neighbour of v, value is the distance related coefficient associated to it. To compute this value, we apply a function h(e) where e is the number of edges from v to vertex #index. If e>Emax or if h(e)<Hmin, vertices are not considered as neighbours. This can be computed as follows using the topology table.
   @code
graph g;
auto topology = vq3::topology::table(g);
...
auto ref_vertex   = topology(3); // we get some vertex.
auto neighborhood = topology.neighborhood(3, // ref_vertex works as well
                                          h, Emax, Hmin);
for(auto& info : neighborhood) {
    auto& ref_v = topology(info.index);
    double coef = info.value;
}
   @endcode

   When some algorithm requires the knowledge of all neigborhoods, they can be computed at once.

   @code
graph g;
auto topology = vq3::topology::table(g);
...
topology(h, Emax, Hmin); // Both vertices and their neighborhoods are updated.
auto ref_vertex   = topology(3);          // we get some vertex.
neighborhood      = topology[3];          // Gets the precomputed neighborhood of the 4th vertex
neighborhood      = topology[ref_vertex]; // Does the same.
   @endcode
   
   @subsection graphutils Utilities

   There are several utilities associated with the graph class, see the vq3::utils namespace.

   @subsection decor Decorated values

   Vertex and edge values, i.e. the type arguments provided to the
   vq3::graph template, can be decorated in order to enable the use of
   some specific algorithms. The decoration consists of adding extra
   attributes to some existing values... typically, adding a boolean
   tag for algorithms that needs to mark some nodes or edges as
   visited is a decoration. These addings are performed by inheritance, but
   inheritance is kept simple thanks to pre-defined decorators.
   
   @code
using layer_0 = std::string;
using layer_1 = vq3::decorator::tagged<layer_0>;      // we add a tag for topology computation.
using layer_2 = vq3::decorator::sum<layer_1, double>; // we add the ability to hande a sum of floating values.
using layer_3 = vq3::decorator::grid_pos<layer_2>;    // we add the ability to register a grid position.
using vertex  = layer_3;

using edge    = vq3::decorator::tagged<void>;         // The edge has nothing but a tag.

using graph   = vq3::graph<vertex, edge>;             // we define a graph with the decorated values.

// Let us illustrate the attributes of the values.

graph::ref_vertex ref_v = ...;

std::string                           val = ref_v->vq3_value;
bool                                  tag = ref_v->vq3_tag;
vq3::utils::accum<double, double>     sum = ref_v->vq3_sum;
std::pair<unsigned int, unsigned int> pos = ref_v->vq3_gridpos;

graph::ref_edge ref_e = ...;

bool tag = ref_e->vq3_tag;
   @endcode

   @section proc Epochs and parallel processors

   
   @subsection epoch The need for computing epochs

   In vq3, the concept of epoch is highlighted, even if many VQ
   algorithms are usually online. Indeed, for the sake of multi-thread
   computing, it is more efficient to consider a bunch of data whose
   computation can be split to feed several threads. A computation pass on 
   this bunch is what is called an 'epoch' in vq3.

   When computing an epoch, i.e. presenting a bunch of data to a
   graph, the trick is to bind temporarily each node of the graph and
   a dedicated data, called an 'epoch data'. This temporary binding is
   done internally, within each thread. At the end of the parallel
   computation, and for each vertex, the epoch data values (one for
   each thread) correponding to that vertex are merged.

   Epoch data types are user defined, they have to fit the
   vq3::concept::EpochData. Nevertheless, the basic use of vq3 for
   defining epochs data consists of using pre-defined building blocks
   that are stacked in order to form the epoch data the user needs
   (this is similar to the decoration of vertex and edge values).

   The meaning of epoch data lifecycle is summarized in this piece of <b>pseudo-code</b>.
   @code
// This is done in parallel by splitting that loop between several threads.
for(sample : current_epoch) {

  // We get the vertex that is closest to the sample, as well as the distance
  // between the two. Then, we notify this information to the data associated 
  // with that vertex.
  auto [min_dist, ref_min_v] = closest(sample);
  min_data = epoch_data_associated_with(ref_min_v);
  min_data.notify_closest(sample, min_dist);

  // In the winner-take-all case, the data associated with
  // the current vertex is the only one to be updated. We notify
  // it with the concerned sample for the next coming update.
  if(winner_take_all) 
    min_data.notify_wta_update(sample);
  
  // In the winner-take-most case, all vertices have to be considered for an
  // update. We notify all of them for the next coming update, providing them 
  // with both the sample and their weight h(d) according to their distance d 
  // to the BMU (i.e. ref_min_v). 
  if(winner_take_most) 
    for(ref_v : vertices) {
      data = epoch_data_associated_with(ref_vv);
      data.notify_wtm_update(sample, h(topological_distance(ref_min_v, ref_v)));
    }
}

// Once the epoch is done, all nodes may have received an update
// notifications (i.e. the corresponding epoch data has stored internally
// data for actual prototype update).
// The vertex value can thus be updated now.
for(ref_vv : vertices) { 
  data           = epoch_data_associated_with(ref_vv); 
  vertex_content = (*ref_vv)();
  prototype      = prototype_of(vertex_content);

  data.set_prototype(prototype);    // This setting is specific to prototype update
  data.set_content(vertex_content); // This setting is more general, operating on the values handled by the vertex.
}
   @endcode

   @subsection epochdata Designing custom epoch data

   The pseudo-code mentionned previously is performed in parallel, and vq3 triggers its execution by using processors. So the customization of what happens during the loop summarized in the psuedo-code is done by defining an epoch data class. Epoch data classes are designed by inheritance, which is hidden by the following stack approach.

   Let us consider the following example:
   @code
using sample       = ...;                                            // sample is the type of our samples.
using epoch_data_0 = vq3::epoch::data::none<sample, vertex, sample>; // This is the root of the stack, the sample type vertex value type and prototype type 
                                                                     // have to be provided. Usually, sample and prototypes have the same type.
using epoch_data_1 = vq3::epoch::data::bmu<epoch_data_0>;            // We collect the sum of the distances for each best-matching vertex.
using epoch_data   = epoch_data_1;
...
some_processor p;
p.process_something<epoch_data>(nb_threads, ...);
   @endcode
   
   The purpose of the type stack is to customize the type epoch_data
   used by the processor. Each stack element (epoch_data_0,
   epoch_data_1, ...) provide the functions so that the final
   epoch_data fits the vq3::concept::EpochData concept. For example,
   the required method notify_wta_update can be defines several times,
   by several stack elements. All the definitions of notify_wta_update
   will be called successively when overall notify_wta_update is
   called.

   Here, epoch_data_0 is the ground of the stack, where the type
   sample is notified. Then, epoch_data_1 adds an accumulator of numbers
   (attribute vq3_bmu_accum) that accumulates the distances between the
   sample and the vertex prototype, when the vertex associated to the
   data is the actual best matching unit. 

   See the vq3::epoch::data namespace.

   @subsection processors Processors

   As previously said, processors enable a parallel computing of epochs, thanks to the epoch data associated to each node. 

   @subsubsection chl Competitive Hebbian Learning processor

   The Competitive Hebbian Learning (CHL) processor applies the
   Martinez and Schulten algorithm to a graph, thus updating its
   edges.
   @code
using vertex = ...;
using edge   = vq3::decorator::tagged<void>; // We need tags on the edges.
using graph  = vq3::graph<vertex, edge>; 

graph g;  
auto  processor = vq3::epoch::chl::processor(g);
auto  S         = ...; // a sample set

auto sample_of(const set_content& content) {
  // returns the sample from *it : sample = sample_of(*it)
}

auto prototype_of(const vertex& v) {
  // returns the prototype from the vertex value.
}

double dist(const prototype& p, const sample& s) {
  // compares the sample to the prototype.
}

auto new_edge_value = edge(); // Value for initializing new edges.

bool topology_modified = processor.update_edges(nb_threads, S.begin(), S.end(),
                                                sample_of, prototype_of, dist,
						new_edge_value);
   @endcode

   @subsubsection wta Winner-take-all processor

   This processor applies a k-means update, i.e. each prototype is set to the centroid of its Voronoi cell.
   @code
using sample       = ...;
using vertex       = ...;

using epoch_data_0 = vq3::epoch::data::none<sample, vertex, sample>;        
using epoch_data_1 = vq3::epoch::data::wta<epoch_data_0>; 
using epoch_data   = epoch_data_1;

using graph  = vq3::graph<vertex, void>; 

graph g;  
auto  topology  = vq3::topology::table(g);
auto  processor = vq3::epoch::wta::processor(table);
auto  S         = ...; // a sample set

auto sample_of(const set_content& content) {
  // returns the sample from *it : sample = sample_of(*it)
}

auto& prototype_of(vertex& v) {
  // returns a reference on the prototype from the vertex value.
}

double dist(const prototype& p, const sample& s) {
  // compares the sample to the prototype.
}

auto new_edge_value = edge(); // Value for initializing new edges.

topology(); // Notifies the graph structure, redo it at each topology change.    
auto epoch_result = processor.process<epoch_data>(nb_threads, S.begin(), S.end(),
                                                  sample_of, prototype_of, dist);
for(epoch_data& data : epoch_result) {
  // data is the epoch data computed for each vertex.
  std::cout << data.vq3_wta_accum.average<double>() << std::endl;
}

@endcode

   @subsubsection wta Winner-take-most processor

   This processor applies a SOM update, i.e. each prototype is set the weighted sum of some samples, the weight depending on the topological proximity of the vertex with the best matching unit.
   @code
using sample       = ...;

using prototype = ...;
using vertex    = vq3::decorator::tagged<prototype>; // We need tags on the vertices, for topology computation.
using graph     = vq3::graph<vertex, void>; 

using epoch_data_0 = vq3::epoch::data::none<sample, vertex, sample>;        
using epoch_data_1 = vq3::epoch::data::wtm<epoch_data_0>; 
using epoch_data   = epoch_data_1;

graph g;  
auto  topology  = vq3::topology::table(g);
auto  processor = vq3::epoch::wtm::processor(g);
auto  S         = ...; // a sample set

auto sample_of(const set_content& content) {
  // returns the sample from *it : sample = sample_of(*it)
}

auto& prototype_of(vertex& v) {
  // returns a reference on the prototype from the vertex value.
}

double dist(const prototype& p, const sample& s) {
  // compares the sample to the prototype.
}

auto new_edge_value = edge(); // Value for initializing new edges.

// Let us the graph structure, redo it at each topology change. 
// We have to specify the way the learning coefficient decreases
// with the number of edges to the best matching unit.
// 5, 1e-3 are the maximal edge distance considered, and 1e-3 
// a minimal threshold for considering a coefficient as non null.
topology([](unsigned int edge_distance) {return std::max(0., 1 - edge_distance/5.0;},
         5, 1e-3));   
auto epoch_result = processor.process<epoch_data>(nb_threads, S.begin(), S.end(),
                                                  sample_of, prototype_of, dist);
for(epoch_data& data : epoch_result) {
  // data is the epoch data computed for each vertex.
  std::cout << data.vq3_wtm_accum.average<double>() << std::endl;
}

  @endcode

  @section algo Amgorithms

  The algorithms provided by vq3 are based on the "processor"
  paradigm. They are illustrated in the documentation of the vq3demo
  package. Read the examples provided by the vq3demo package in the
  order suggested by the numbers to gat a practical view of the vq3
  library.

*/
