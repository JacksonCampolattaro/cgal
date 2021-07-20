//
// Created by jcampola on 7/19/21.
//

#include <iostream>

#include <CGAL/Simple_cartesian.h>

#include <xsimd/xsimd.hpp>

typedef CGAL::Simple_cartesian<float> K;
typedef K::Ray_3 Ray_3;
typedef CGAL::Bbox_3 Bbox_3;

bool __attribute__((noinline)) single_intersection(const Ray_3 &ray, const Bbox_3 &bbox) {
  return CGAL::do_intersect(ray, bbox);
}

template<std::size_t N>
std::array<bool, N> __attribute__((noinline)) __attribute__((flatten)) loop_intersection(const Ray_3 &ray, const std::array<Bbox_3, N> &boxes) {
  std::array<bool, N> results{false};
  for (std::size_t i = 0; i < N; ++i) {
    results[i] = CGAL::do_intersect(ray, boxes[i]);
  }
  return results;
}

template<std::size_t N>
std::array<bool, N> __attribute__((noinline)) xsimd_intersection(const Ray_3 &ray, const std::array<Bbox_3, N> &bbox) {

  // Load values into simd registers
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  xsimd::batch<double, N> ray_source_x(ray.source().x());
  xsimd::batch<double, N> ray_source_y(ray.source().y());
  xsimd::batch<double, N> ray_source_z(ray.source().z());
  xsimd::batch<double, N> ray_target_x(ray.second_point().x());
  xsimd::batch<double, N> ray_target_y(ray.second_point().y());
  xsimd::batch<double, N> ray_target_z(ray.second_point().z());

  xsimd::batch<double, N> box_min_x, box_min_y, box_min_z, box_max_x, box_max_y, box_max_z;
  for (std::size_t i = 0; i < N; ++i) {
    box_min_x[i] = bbox[i].xmin();
    box_min_y[i] = bbox[i].ymin();
    box_min_z[i] = bbox[i].zmin();
    box_max_x[i] = bbox[i].xmax();
    box_max_y[i] = bbox[i].ymax();
    box_max_z[i] = bbox[i].zmax();
  }

  // Use a standard intersection algorithm
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // The ray's direction is the difference between its target and origin
  xsimd::batch<double, N> ray_direction_x = ray_target_x - ray_source_x;
  xsimd::batch<double, N> ray_direction_y = ray_target_y - ray_source_y;
  xsimd::batch<double, N> ray_direction_z = ray_target_z - ray_source_z;

  // Determine signs of the ray's directions
  xsimd::batch_bool<double, N> negative_x = ray_direction_x < 0;
  xsimd::batch_bool<double, N> negative_y = ray_direction_y < 0;
  xsimd::batch_bool<double, N> negative_z = ray_direction_z < 0;

  // Determine box starts and ends, based on the ray's direction
  xsimd::batch<double, N> box_start_x = xsimd::select(negative_x, box_max_x, box_min_x);
  xsimd::batch<double, N> box_start_y = xsimd::select(negative_y, box_max_y, box_min_y);
  xsimd::batch<double, N> box_start_z = xsimd::select(negative_z, box_max_z, box_min_z);
  xsimd::batch<double, N> box_end_x = xsimd::select(!negative_x, box_max_x, box_min_x);
  xsimd::batch<double, N> box_end_y = xsimd::select(!negative_y, box_max_y, box_min_y);
  xsimd::batch<double, N> box_end_z = xsimd::select(!negative_z, box_max_z, box_min_z);

  // Determine bounds of overlap in each direction
  xsimd::batch<double, N> min_x = (box_start_x - ray_source_x) / ray_direction_x;
  xsimd::batch<double, N> min_y = (box_start_y - ray_source_y) / ray_direction_y;
  xsimd::batch<double, N> min_z = (box_start_z - ray_source_z) / ray_direction_z;
  xsimd::batch<double, N> max_x = (box_end_x - ray_source_x) / ray_direction_x;
  xsimd::batch<double, N> max_y = (box_end_y - ray_source_y) / ray_direction_y;
  xsimd::batch<double, N> max_z = (box_end_z - ray_source_z) / ray_direction_z;

  // Determine combined region of overlap
  xsimd::batch<double, N> min = xsimd::max(min_x, xsimd::max(min_y, min_z));
  xsimd::batch<double, N> max = xsimd::max(max_x, xsimd::max(max_y, max_z));

  // If there is positive overlap, then there is an intersection
  xsimd::batch_bool<double, N> result = max > min;

  // Write the result to an array that can be returned
  std::array<bool, N> result_array{};
  result.store_unaligned(result_array.data());
  return result_array;
}

int main() {

  Bbox_3 box{0, 0, 0, 1, 1, 1};
  Ray_3 ray{};

  std::array<Bbox_3, 2> boxes{box, box};

  std::cout << single_intersection(ray, box) << std::endl;

  for (auto r : loop_intersection(ray, boxes)) std::cout << r;
  std::cout << std::endl;

  for (auto r : xsimd_intersection(ray, boxes)) std::cout << r;
  std::cout << std::endl;

  for (auto b : boxes) std::cout << b << std::endl;
}