#include <cmath>
#include <utility>
#include <vector>
#include <iterator>
#include <random>

#include <opencv2/opencv.hpp>
#include <vq3.hpp>
#include <vq3demo.hpp>
#include <demo2d.hpp>

/*

  The purpose of this example is to imlement GNG-T (online here) with
  a custom type of data. Indeed, vq3 is generic and it can apply to
  any kind of data type, as soon as they fit some basic concept. This
  is what we detail here.

*/

// Our custom type is the set of rotations in 2D. It can be
// represented as a complex number, with norm 1 (as queternions
// represent rotations in 3D). We do not reuse std::complex here, for
// the sake of illustration, so that we have to implement very
// necessary functions (and not more).
class Rotation {
private:
  double re;
  double im;

  void normalize() {
    double norm_2 = re*re + im*im;
    if(norm_2 == 0) {
      re = 1;
      im = 0;
    }
    else {
      double inv_n = 1/std::sqrt(norm_2);
      re *= inv_n;
      im *= inv_n;
    }
  }
  
  Rotation(double re, double im) : re(re), im(im) {
    normalize();
  }
  
public:
  static double to_rad(double degree) {return degree * 0.017453292519943295;}

  Rotation() : re(1), im(0) {}
  Rotation(double angle) : re(std::cos(angle)), im(std::sin(angle)) {}

  Rotation(const Rotation&)            = default;
  Rotation& operator=(const Rotation&) = default;

  // This returns a Rotation that is very close to this.
  Rotation very_close() const {
    return {re + 1e-5, im + 1e-5}; 
  }

  // For the sake of plotting into figures, we enable a casting to 2D points.
  operator demo2d::Point () const {return {re, im};}

  // Here are the operations that are required for GNG-T (or online vq
  // computation). Indeed, the learning consists of computing
  //   w += alpha * (xi - w)
  // So we have to make this work. This can be tricky...

  // Here, we want to compute everything as if rotations were complex
  // numbers, so that only the last operation, +=, does the
  // normalization. Intermediate computation have to handle
  // un-normalized (re, im) pairs. Let us store them in actual std::pair.

  using intermediate_type = std::pair<double, double>;

  intermediate_type operator-(const Rotation& arg) const {
    return {re - arg.re, im - arg.im};
  }

  void operator+=(const intermediate_type& alpha_times_diff) {
    re += alpha_times_diff.first;
    im += alpha_times_diff.second;
    normalize(); // Here, at the end of the update process, we normalize.
  }
};

inline Rotation::intermediate_type operator*(double alpha, const Rotation::intermediate_type& difference) {
  return {alpha * difference.first, alpha * difference.second};
}


// That's it, let us use it with an online GNG-T.


// Graph definition

using sample    = Rotation;
using prototype = Rotation;

//                                                 ## Node properties :
using vlayer_0 = prototype;                        // prototypes are 2D points (this is the "user defined" value).
using vlayer_1 = vq3::decorator::tagged<vlayer_0>; // we add a tag for topology computation.
using vertex   = vlayer_1;

//                                                 ## Edge properties :
using elayer_0 = vq3::decorator::tagged<void>;     // we add a tag for CHL computation.
using edge     = elayer_0;

using graph  = vq3::graph<vertex, edge>;

using neighbour_key_type = int;

// Distance

// This is the distance used by closest-like algorithms. We need to
// compare actual vertex values with samples. Here, demo2d::d2 works
// because v.vq3_value and p are Rotation's, that can be casted into
// demo2d::Points.
double dist(const vertex& v, const sample& p) {return demo2d::d2(v.vq3_value, p);}

struct Evolution {
  // here, evolution consists of adding or removing prototypes so as
  // to have a determined number of prototypes.
  unsigned int T = 0;
  
  // GNG-T calls this method when it considers to perform a
  // modification of the number of vertices.
  template<typename TABLE, typename BMU_RESULT, typename CLONE_PROTOTYPE, typename FCT_ERROR_OF_ACCUM>
  bool operator()(TABLE&                    topology,
		  const BMU_RESULT&         neighboring_bmu_epoch_result,
		  const CLONE_PROTOTYPE&    clone_prototype,
		  const FCT_ERROR_OF_ACCUM& error_of_accum) {
    
    // Let us check is some node is isolated and useless.
    if(topology.size() > 1) { 
      bool remove = false;
      std::size_t vertex_idx = 0;
      for(auto& res : neighboring_bmu_epoch_result) {
	if(res.vq3_bmu_accum.nb == 0) { // The node has never won. We have to check if it has edges.
	  // The following counts the edges of the node.
	  auto ref_v = topology(vertex_idx);
	  unsigned int nb_edges = 0;
	  ref_v->foreach_edge([&nb_edges](auto ref_e) {
				auto extr_pair = ref_e->extremities();           
				if(vq3::invalid_extremities(extr_pair))
				  ref_e->kill();
				else
				  ++nb_edges;
			      });
	  if(nb_edges == 0) {
	    ref_v->kill();
	    remove = true;
	  }
	}
	++vertex_idx;
      }
      if(remove)
	return true;
    }
    
    if(topology.size() < T) {
      auto max_iter = std::max_element(neighboring_bmu_epoch_result.begin(), neighboring_bmu_epoch_result.end(),
				       [&error_of_accum](auto& content1, auto& content2){
					 return error_of_accum(content1.vq3_bmu_accum) < error_of_accum(content2.vq3_bmu_accum);});
      auto ref_vertex = topology(std::distance(neighboring_bmu_epoch_result.begin(), max_iter));
      topology.g += clone_prototype((*ref_vertex)().vq3_value);
      return true;
    }

    if(topology.size() > T) {
      auto min_iter = std::min_element(neighboring_bmu_epoch_result.begin(), neighboring_bmu_epoch_result.end(),
				       [&error_of_accum](auto& content1, auto& content2){
					 return error_of_accum(content1.vq3_bmu_accum) < error_of_accum(content2.vq3_bmu_accum);});
      auto ref_vertex = topology(std::distance(neighboring_bmu_epoch_result.begin(), min_iter));
      ref_vertex->kill();
      return true;
    }
    
    return false;
  }
};

inline Evolution make_evolution() {return Evolution();}

int main(int argc, char* argv[]) {
  
  unsigned int nb_threads = std::thread::hardware_concurrency();
  if(nb_threads == 0) {
    nb_threads = 1;
    std::cout << "The harware multi-threading capabilities cannot be queried. I use a single thread." << std::endl;
  }
  else
    std::cout << std::endl
	      << "I use " << nb_threads << " thread(s)." << std::endl;
  
  std::random_device rd;  
  std::mt19937 random_device(rd());
  
  graph g;
  
  int slider_T = 10;
  
  cv::namedWindow("algorithm", cv::WINDOW_AUTOSIZE);
  cv::createTrackbar("T", "algorithm", &slider_T, 100, nullptr);
  
  auto image = cv::Mat(600, 600, CV_8UC3, cv::Scalar(255,255,255));
  auto frame = demo2d::opencv::direct_orthonormal_frame(image.size(), .45*image.size().width, true);
  
  auto dd          = demo2d::opencv::dot_drawer<demo2d::Point>(image, frame,
							       [](const demo2d::Point& pt) {return                      true;},
							       [](const demo2d::Point& pt) {return                        pt;},
							       [](const demo2d::Point& pt) {return                         1;},
							       [](const demo2d::Point& pt) {return cv::Scalar(200, 200, 200);},
							       [](const demo2d::Point& pt) {return                        -1;});
  auto draw_edge   = vq3::demo2d::opencv::edge_drawer<graph::ref_edge>(image, frame,
								       [](const vertex& v1, const vertex& v2, const edge& e) {return true;}, // always draw
								       [](const vertex& v)              {return               v.vq3_value;}, // position
								       [](const edge&   e)              {return cv::Scalar(255, 150, 150);}, // color
								       [](const edge&   e)              {return                        3;}); // thickness
  auto draw_vertex = vq3::demo2d::opencv::vertex_drawer<graph::ref_vertex>(image, frame,
									   [](const vertex& v) {return        true;},                // always draw
									   [](const vertex& v) {return v.vq3_value;},                // position
									   [](const vertex& v) {return           5;},                // radius
									   [](const vertex& v) {return cv::Scalar(200, 000, 000);},  // color
									   [](const vertex& v) {return          -1;});               // thickness
  

  // Some initializations
 
  std::cout << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl
	    << "press ESC to quit." << std::endl
	    << std::endl
	    << "##################" << std::endl
	    << std::endl;
  int keycode = 0;

  // Let us build up a dataset.
  std::vector<Rotation> S;
  auto out = std::back_inserter(S);

#define NB_PER_SLICE 1000
  {
    auto uniform = std::uniform_real_distribution<double>(0, Rotation::to_rad(100));
    for(unsigned int i = 0; i < NB_PER_SLICE; ++i) *(out++) = uniform(random_device);
  }

  {
    auto uniform = std::uniform_real_distribution<double>(Rotation::to_rad(120), Rotation::to_rad(200));
    for(unsigned int i = 0; i < NB_PER_SLICE; ++i) *(out++) = uniform(random_device);
  }

  {
    auto uniform = std::uniform_real_distribution<double>(Rotation::to_rad(-120), Rotation::to_rad(-30));
    for(unsigned int i = 0; i < NB_PER_SLICE; ++i) *(out++) = uniform(random_device);
  }

  std::shuffle(S.begin(), S.end(), random_device);

  
  // This keeps up to date information about the graph topology.
  auto topology = vq3::topology::table<neighbour_key_type>(g);
  topology.declare_distance(0,
			    [](unsigned int edge_distance) {return edge_distance == 0 ? 1 : .2;},
			    1,  0.0);

  
  // This processes the topology evolution (number of vertices and edges)
  auto gngt  = vq3::algo::online::gngt::processor<sample>(topology);
  gngt.alpha = .05;
  
  // This is how the default evolution would have been obtained.
  //   auto evolution = vq3::algo::gngt::size_control::evolution();
  // But we use here the one that we have handcrafted.
  auto evolution = make_evolution();

  
  while(keycode != 27) {       // no ESC key pressed.
    
    evolution.T = slider_T;
    std::shuffle(S.begin(), S.end(), random_device);
    gngt.process(nb_threads,
		 S.begin(), S.end(),                                 // The sample set. Shuffle if the dataset is not sampled randomly.
		 [](const sample& s) {return s;},                    // get sample from *iter (identity here).
		 [](vertex& v) -> prototype& {return v.vq3_value;},  // get a prototype reference from the vertex value.
		 [](const prototype& p) {return p.very_close();},    // get a point close to a prototype.
		 dist,                                               // The distance used for bmu-related stuff.
		 0,                                                  // Neighborhood key.
		 evolution);  
      
    // Display    
    image = cv::Scalar(255, 255, 255);
    std::copy(S.begin(), S.end(), dd);
    g.foreach_edge(draw_edge); 
    g.foreach_vertex(draw_vertex);
      
    cv::imshow("algorithm", image);
    keycode = cv::waitKey(10) & 0xFF;
  }

  
  return 0;
}

