//
// Created by jcampola on 7/19/21.
//

#include <iostream>

#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<float> K;
typedef K::Ray_3 Ray_3;
typedef CGAL::Bbox_3 Bbox_3;

bool __attribute__((noinline)) single(const Ray_3 &ray, const Bbox_3 &bbox) {
  return CGAL::do_intersect(ray, bbox);
}

std::array<bool, 4> __attribute__((noinline)) __attribute__((flatten)) quad(const Ray_3 &ray, const std::array<Bbox_3, 4> &boxes) {
  return {CGAL::do_intersect(ray, boxes[0]),
          CGAL::do_intersect(ray, boxes[1]),
          CGAL::do_intersect(ray, boxes[2]),
          CGAL::do_intersect(ray, boxes[3])};
}

int main() {

  Bbox_3 box{0, 0, 0, 1, 1, 1};
  Ray_3 ray{};

  std::array<Bbox_3, 4> boxes{box, box, box, box};

  std::cout << single(ray, box) << std::endl;
  for (auto r : quad(ray, boxes)) std::cout << r;
  std::cout << std::endl;

  for (auto b : boxes) std::cout << b << std::endl;
}