
#include <CGAL/Bbox_3.h>
#include <CGAL/Intersections_3/Bbox_3_Bbox_3.h>

int main() {

  // Bounding boxes that don't intersect

  CGAL::Bbox_3 a{1, 2, 3,
                 4, 5, 6};
  CGAL::Bbox_3 b{7, 8, 9,
                 10, 11, 12};

  std::cout << CGAL::do_intersect(a, b);
  CGAL_assertion_msg(!CGAL::do_intersect(a, b), "Bbox a and b should not intersect");

  // Bounding boxes that do intersect

  CGAL::Bbox_3 c{1, 2, 3,
                 4, 5, 6};
  CGAL::Bbox_3 d{1, 2.5, 3,
                 3, 3, 3};

  CGAL_assertion_msg(CGAL::do_intersect(c, d), "Bbox c and d should not intersect");

  // Self-intersection
  CGAL_assertion_msg(CGAL::do_intersect(c, c), "Bbox c should intersect with itself");
}