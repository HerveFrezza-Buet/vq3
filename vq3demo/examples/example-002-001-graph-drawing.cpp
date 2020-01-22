#include <vq3.hpp>
#include <vq3demo.hpp>
#include <opencv2/opencv.hpp>

#include <vector>
#include <iterator>
#include <random>



using intensity = int;

struct edge_value_type {
  static int nb;
  intensity i;
  edge_value_type(intensity i) : i(i) {++nb;}
  edge_value_type() = delete;
  edge_value_type(const edge_value_type& e) : i(e.i) {++nb;}
  ~edge_value_type() {--nb;}
};

struct vertex_value_type {
  static int nb;
  intensity          i;
  demo2d::Point pos;
  
  vertex_value_type(intensity i, double x, double y) : i(i), pos(x,y) {++nb;}
  vertex_value_type() = delete;
  vertex_value_type(const vertex_value_type& v) : i(v.i), pos(v.pos) {++nb;}
  ~vertex_value_type() {--nb;}
};

using graph = vq3::graph<vertex_value_type, edge_value_type>;

int edge_value_type::nb   = 0;
int vertex_value_type::nb = 0;

int main(int argc, char* argv[]) {
  std::random_device rd;  
  std::mt19937 random_device(rd());

  // Let us build a graph. Vertex positions are in [0,1]
  
  graph g;

  // We set up the graph, and keep a temporary access to the created
  // vertices for connections.
  {
    unsigned int mesh_side      = 20;
    unsigned int nb_vertex_side = mesh_side + 1;
    
    // vertices
    
    std::vector<graph::ref_vertex> temporary_vertex_references;
    auto out = std::back_inserter(temporary_vertex_references);
    auto r = demo::range(.1, .9, mesh_side);
    for(auto y : r)
      for(auto x : r) {
	// we create a vertex : g += vertex_value_type
	auto ref_vertex = g += {intensity(x*255), x, y};
	// we store the reference returned for further connections.
	*(out++) = ref_vertex;
      }

    // edges

    auto it  = temporary_vertex_references.begin();
    auto yit = r.begin();
    for(unsigned int h = 0; h < nb_vertex_side - 1; ++h, ++yit, it += nb_vertex_side) {
      auto iit = it;
      for(unsigned int w = 0; w < nb_vertex_side - 1; ++w, ++iit) {
    	g.connect(*iit, *(iit + 1), intensity(*(yit)*255));
    	g.connect(*iit, *(iit + nb_vertex_side), intensity(*(yit)*255));
      }
      g.connect(*iit, *(iit + nb_vertex_side), intensity(*(yit)*255));
    }
    for(unsigned int w = 0; w < nb_vertex_side - 1; ++w, ++it)
      g.connect(*it, *(it + 1), intensity(*(yit)*255));
    
  } // The temporary references to vertices are released here.

  // Let us display the graph

  cv::namedWindow("image", CV_WINDOW_AUTOSIZE);
  int delay = 500;
  auto image = cv::Mat(480, 640, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = demo2d::opencv::direct_orthogonal_frame(image.size().width,
						       image.size().height,
						       {0., image.size().height - 1.0});
  auto draw_edge = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
  								     [](const vertex_value_type&, const vertex_value_type&, const edge_value_type&) {return true;},  // always draw
  								     [](const vertex_value_type& v)                {return                                 v.pos;},  // position
  								     [](const edge_value_type& e)                  {return cv::Scalar(255 - e.i, 255 - e.i, 255);},  // color
  								     [](const edge_value_type&)                    {return                                     3;}); // thickness

  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex_value_type& v) {return                                  true;},  // always draw
									   [](const vertex_value_type& v) {return                                 v.pos;},  // position
									   [](const vertex_value_type& v) {return                                     5;},  // radius
									   [](const vertex_value_type& v) {return cv::Scalar(255, 255 - v.i, 255 - v.i);},  // color
									   [](const vertex_value_type& v) {return                                    -1;}); // thickness

  // This is the drawing....
  g.foreach_edge(draw_edge); 
  g.foreach_vertex(draw_vertex);
  
  cv::imshow ("image", image);
  cv::waitKey(delay);

  // from now, let us remove the nodes which have less than 4
  // neigbours with probability p and systematically nodes with no
  // edge, until the graph is empty.

  double p = .2;
  
  bool empty = false;
  while(!empty) {
    std::cout << vertex_value_type::nb << " vertices, " << edge_value_type::nb << " edges." << std::endl;
    empty = true;
    std::uniform_real_distribution<double> toss(0,1);
    g.foreach_vertex([&random_device, &toss, &empty, p](const graph::ref_vertex& ref_v) {
	empty = false; // ... since we enter an actual iteration.
	unsigned int edge_count = 0;
	ref_v->foreach_edge([&edge_count](const graph::ref_edge& ref_e) {
	    // Let us get the edge extremities
	    auto extremity_pair = ref_e->extremities();
	    if(vq3::invalid_extremities(extremity_pair)) {// if some extremity is a dead node.
	      ref_e->kill(); // The edge is invalid, we kill it. 
	      return;
	    }
	    ++edge_count; // the edge is ok, we count it.
	  });
	if(edge_count == 0 || (edge_count < 4 && toss(random_device) < p)) // The node has to be removed from the graph.
	  ref_v->kill();
      });

    image = cv::Scalar(255, 255, 255);
    g.foreach_edge(draw_edge); 
    g.foreach_vertex(draw_vertex);
    cv::imshow ("image", image);
    cv::waitKey(delay);
  }
  std::cout << "Empty graph : " << vertex_value_type::nb << " vertices, " << edge_value_type::nb << " edges." << std::endl;
  
  return 0;
}
