#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>
#include <fstream>

#include "vector3.h"
#include "ray.h"
#include "bbox.h"

#include "intersection_strategies/smits_method.h"
#include "intersection_strategies/improved.h"
#include "intersection_strategies/clarified.h"
#include "intersection_strategies/branchless.h"
#include "intersection_strategies/xsimd.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

std::vector<std::pair<Ray, std::vector<BBox>>> load_scenarios(std::ifstream &file, long N) {
  std::vector<std::pair<Ray, std::vector<BBox>>> scenarios;

  double px, py, pz, qx, qy, qz,
          bxmin, bymin, bzmin, bxmax, bymax, bzmax;

  for (int i = 0; i < N; ++i) {

    // Read ray data
    file >> px >> py >> pz
         >> qx >> qy >> qz;
    Vector3 origin = {px, py, pz};
    Vector3 direction = {qx, qy, qz};
    Ray ray = {origin, direction};

    // Read box data
    file >> bxmin >> bymin >> bzmin
         >> bxmax >> bymax >> bzmax;
    Vector3 min = {bxmin, bymin, bzmin};
    Vector3 max = {bxmax, bymax, bzmax};
    BBox box = {min, max};

    // Only create a new scenario when the query ray has changed
    if (scenarios.empty() || !(ray == scenarios.back().first))
      scenarios.emplace_back(ray, std::vector<BBox>());

    scenarios.back().second.push_back(box);
  }

  long count = 0;
  for (const auto &scenario : scenarios) count += scenario.second.size();

  return scenarios;
}

double time(const std::function<void(void)> &f) {

  auto start = high_resolution_clock::now();
  {
    f();
  }
  auto end = high_resolution_clock::now();

  return (end - start).count();
}

int main() {

  long N = 3742217;
  long R = 10;

  // Load test data
  auto file = std::ifstream("data/remeshing_intersections_3742217.txt");
  if (!file.is_open()) exit(1);
  auto scenarios = load_scenarios(file, N);

  std::vector<double> smits_method_times, improved_times, clarified_times, branchless_times, xsimd_times;

  for (int i = 0; i < R; ++i) {

    long sum = 0;

    for (const auto &scenario : scenarios) {

      smits_method_times.push_back(time([&] {
        const auto &ray = scenario.first;
        for (const auto &bbox : scenario.second)
          sum += smits_method::intersect(bbox, ray);
      }));

      improved_times.push_back(time([&] {
        const auto &ray = scenario.first;
        for (const auto &bbox : scenario.second)
          sum += improved::intersect(bbox, ray);
      }));

      clarified_times.push_back(time([&] {
        const auto &ray = scenario.first;
        for (const auto &bbox : scenario.second)
          sum += clarified::intersect(bbox, ray);
      }));

      branchless_times.push_back(time([&] {
        const auto &ray = scenario.first;
        for (const auto &bbox : scenario.second)
          sum += branchless::intersect(bbox, ray);
      }));

      xsimd_times.push_back(time([&] {
        const auto &ray = scenario.first;
        for (const auto &bbox : scenario.second)
          sum += xsimd::intersect(bbox, ray);
      }));

    }
    std::cout << i << ":\t" << (float) sum / (float) (N * 5) << std::endl;
  }

  std::cout << "{| class=\"wikitable\"" << std::endl;
  std::cout << "|+ Time to Complete " << N << " Intersection Tests" << std::endl;
  std::cout << "| Smits' Method || "
            << std::accumulate(smits_method_times.begin(), smits_method_times.end(), 0.0) / (double) R
            << " ms"
            << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "| Improved || "
            << std::accumulate(improved_times.begin(), improved_times.end(), 0.0) / (double) R
            << " ms"
            << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "| Clarified || "
            << std::accumulate(clarified_times.begin(), clarified_times.end(), 0.0) / (double) R
            << " ms"
            << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "| Branchless || "
            << std::accumulate(branchless_times.begin(), branchless_times.end(), 0.0) / (double) R
            << " ms"
            << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "| XSimd || "
            << std::accumulate(xsimd_times.begin(), xsimd_times.end(), 0.0) / (double) R
            << " ms"
            << std::endl;
  std::cout << "|-" << std::endl;
  std::cout << "|}";
}