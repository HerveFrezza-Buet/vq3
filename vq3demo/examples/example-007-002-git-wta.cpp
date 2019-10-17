#include <vector>
#include <iterator>
#include <algorithm>

#include <opencv2/opencv.hpp>
#include <vq3.hpp>
#include <vq3demo.hpp>

/*

  This examples show that processors cannot really be exploited with
  GIT values. Indeed, each distance computation involves a A*
  computation on the support graph. On that graph, vertices are tagged
  during A*, so many A* instances cannot run in parallel since they
  cannot benefit from a dedicated tagging process (tags are shared).

  This example shows how to perform with such constraints.

*/



