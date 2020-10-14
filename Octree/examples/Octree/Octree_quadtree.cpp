#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Random.h>

// Type Declarations
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Point_2 Point_2;
typedef std::vector<Point_3> Point_3_vector;
typedef std::vector<Point_2> Point_2_vector;

typedef CGAL::Octree::Octree<Point_3_vector> Octree;
typedef CGAL::Octree::Octree<Point_2_vector> Quadtree;

int main(int argc, char **argv)
{
  Point_3_vector points3;
  for (std::size_t i = 0; i < 5; ++ i)
    points3.emplace_back(r.get_double(-1., 1.),
                         r.get_double(-1., 1.),
                         r.get_double(-1., 1.));

  Octree octree(points3);
  octree.refine(10, 1);
  std::cerr << "Octree = " << std::endl
            << octree << std::endl;

  Point_2_vector points2;
  for (std::size_t i = 0; i < 5; ++ i)
    points2.emplace_back(r.get_double(-1., 1.),
                         r.get_double(-1., 1.));

  Quadtree quadtree(points2);
  quadtree.refine(10, 1);
  std::cerr << "Quadtree = " << std::endl
            << quadtree << std::endl;

  return EXIT_SUCCESS;
}
