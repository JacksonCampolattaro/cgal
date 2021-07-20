// Copyright (c) 2012  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau


#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_DO_INTERSECT_H

#include <CGAL/disable_warnings.h>

#include <CGAL/double.h>
#include <CGAL/number_utils.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Kernel/Same_uncertainty.h>
#include <CGAL/assertions.h>
#include <CGAL/Coercion_traits.h>
#include <boost/type_traits/is_same.hpp>

// inspired from http://cag.csail.mit.edu/~amy/papers/box-jgt.pdf

// This algorithm intersects the line with the x-, y-, and z-slabs of the
// bounding box, and computes the interval [t1, t2], in the
// parameterization of the line given by the segment (for t=0, that is the
// source of the segment, and for t=1 that is its target), where the line
// intersects the three slabs of the bounding box.

// For a segment, the intersection is non-empty iff
//    [t1, t2] intersects [0, 1].

namespace CGAL {

  namespace Intersections {

    namespace internal {

      template<typename FT, bool bounded_0, bool use_static_filters = false>
      struct Do_intersect_bbox_segment_aux_is_greater {
        typedef typename Same_uncertainty<bool, FT>::type result_type;

        void register_new_input_values(const FT &, const FT &) {}

        void compute_new_error_bound() {}

        bool bound_overflow() { return false; }

        bool value_might_underflow() { return false; }

        static result_type uncertain() {
          return true;
        }

        result_type operator()(const FT &a, const FT &b) const {
          return a > b;
        }
      }; // end struct template Do_intersect_bbox_segment_aux_is_greater

      template<typename FT, bool bounded_0>
      class Do_intersect_bbox_segment_aux_is_greater<FT, bounded_0, true> {
        double error;
        double tmax;
        double dmax;

      public:
        CGAL_static_assertion((boost::is_same<FT, double>::value));

        Do_intersect_bbox_segment_aux_is_greater() : error(0.), tmax(0.), dmax(0.) {}

        void register_new_input_values(const double &t, const double &d) {
          if (bounded_0) {
            if (t > tmax) tmax = t;
            if (d > dmax) dmax = d;
          } else {
            const double at = CGAL::abs(t);
            const double ad = CGAL::abs(d);
            if (at > tmax) tmax = at;
            if (ad > dmax) dmax = ad;
          }
        }

        void compute_new_error_bound() {
          const double EPS = 8.8872057372592798e-16;
          error = tmax * dmax * EPS;
        }

        bool bound_overflow() {
          const double OVERF = 1e153;
          return dmax > OVERF || tmax > OVERF;
        }

        bool value_might_underflow() {
          const double UNDERF = 1e-146;
          return dmax < UNDERF || tmax < UNDERF;
        }

        typedef Uncertain<bool> result_type;

        static result_type uncertain() {
          return result_type::indeterminate();
        }

        result_type operator()(const FT &a, const FT &b) const {
          const FT x = a - b;
          if (x > error) return true;
          else if (x < -error) return false;
          else return uncertain();
        }

      }; // end specialization Do_intersect_bbox_segment_aux_is_greater<FT, true>

      template<typename FT,
              typename BFT,
              bool bounded_0,
              bool bounded_1,
              bool use_static_filters>
      inline
      typename Do_intersect_bbox_segment_aux_is_greater<
              FT,
              bounded_0,
              use_static_filters
      >::result_type
      do_intersect_bbox_segment_aux_old(
              const FT &px, const FT &py, const FT &pz,
              const FT &qx, const FT &qy, const FT &qz,
              const BFT &bxmin, const BFT &bymin, const BFT &bzmin,
              const BFT &bxmax, const BFT &bymax, const BFT &bzmax) {

        // If the entire line segment is contained by the bbox, it intersects
        if (((px >= bxmin) && (px <= bxmax) &&
             (py >= bymin) && (py <= bymax) &&
             (pz >= bzmin) && (pz <= bzmax)) ||
            ((qx >= bxmin) && (qx <= bxmax) &&
             (qy >= bymin) && (qy <= bymax) &&
             (qz >= bzmin) && (qz <= bzmax))) {
          return true;
        }

        // The following code encode t1 and t2 by:
        //    t1 = tmin/dmin
        //    t2 = tmax/dmax
        // For the first lines, dmax==dmin and is not explicitly defined.

        // -----------------------------------
        // treat x coord
        // -----------------------------------
        typedef typename Coercion_traits<double, FT>::Type CFT;
        CFT dmin, tmin, tmax, dmax;

        // If the ray points to the right
        if (qx >= px) {
          // There is no intersection if the ray starts on the right of the box (and doesn't extend past its origin)
          if (bounded_0 && px > bxmax) return false; // segment on the right of bbox

          // There is no intersection if the ray ends on the left of the box (and doesn't extend past its terminus)
          if (bounded_1 && qx < bxmin) return false; // segment on the left of bbox

          // If the ray does not extend past the right of the box
          if (bounded_1 && bxmax > qx) {

            // Set the maximum t bounds to the full extension of the ray
            tmax = 1;
            // todo what is dmax used for?
            dmax = 1;

          } else {

            // Otherwise, shorten the maximum t bounds to the point where it exits the box
            tmax = bxmax - px;
            // todo what is dmax used for?
            dmax = qx - px;
          }

          // Set the minimum t bounds to be the length along the ray where it first encounters the box
          tmin = bxmin - px;
          // todo what is dmin used for?
          dmin = qx - px;
        } else {
          // Otherwise, we know the ray points to the left

          // There is no intersection if the ray ends on the right of the box (and doesn't extend past its terminus)
          if (bounded_1 && qx > bxmax) return false; // segment on the right of bbox

          // There is no intersection if the ray starts on the left of the box (and doesn't extend past its origin)
          if (bounded_0 && px < bxmin) return false; // segment on the left of bbox

          // If the ray does not extend past the left of the box
          if (bounded_1 && bxmin < qx) {

            // Set the maximum t bounds to the full extension of the ray
            tmax = 1;
            // todo what is dmax used for?
            dmax = 1;

          } else {

            // Otherwise, shorten the maximum t bounds to the point where it exits the box
            tmax = px - bxmin;
            // todo what is dmax used for?
            dmax = px - qx;
          }

          // Set the minimum t bounds to be the length along the ray where it first encounters the box
          tmin = px - bxmax;
          // todo what is dmin used for?
          dmin = px - qx;
        }

        // Enforce 0 bounds by making sure tmin >= 0
        if (bounded_0) tmin = (CGAL::max)(CFT(0), tmin);

        // Special case: if we have an infinite ray and it has no x-angle, it's easy to check that it doesn't intersect
        if ((px == qx) && (!(bounded_0 && bounded_1))) {
          // If the ray starts with an x-value outside the box, then it must not intersect with the box
          if (px > bxmax || px < bxmin) return false;
          // Note: for a segment the condition has already been tested by the two
          // previous tests tmax<0 || tmin>dmin (with dmin==0).
        }

        // If dmin == 0, at this point, [t1, t2] == ]-inf, +inf[, or t1 or t2
        // is a NaN. But the case with NaNs is treated as if the interval
        // [t1, t2] was ]-inf, +inf[.

        // d bounds should never be negative
        CGAL_assertion(dmin >= 0);
        CGAL_assertion(dmax >= 0);

        // If we're bounded by the start of the ray, t bounds should never be negative
        if (bounded_0) {
          CGAL_assertion(tmin >= 0);
          CGAL_assertion(tmax >= 0);
        }


        // -----------------------------------
        // treat y coord
        // -----------------------------------
        // (Logic is the same as X)
        CFT dymin, tymin, tymax, dymax;
        if (qy >= py) {
          if (bounded_0 && py > bymax) return false; // segment on the right of bbox
          if (bounded_1 && qy < bymin) return false; // segment on the left of bbox

          if (bounded_1 && bymax > qy) {
            tymax = 1;
            dymax = 1;
          } else {
            tymax = bymax - py;
            dymax = qy - py;
          }

          tymin = bymin - py;
          dymin = qy - py;
        } else {
          if (bounded_1 && qy > bymax) return false; // segment on the right of bbox
          if (bounded_0 && py < bymin) return false; // segment on the left of bbox

          if (bounded_1 && bymin < qy) {
            tymax = 1;
            dymax = 1;
          } else {
            tymax = py - bymin;
            dymax = py - qy;
          }

          tymin = py - bymax;
          dymin = py - qy;
        }

        if (bounded_0) tymin = (CGAL::max)(CFT(0), tymin);

        // If the query is vertical for y, then check its y-coordinate is in
        // the y-slab.
        if ((py == qy) &&  // <=> (dmin == 0)
            (!(bounded_0 && bounded_1))) // do not check for a segment
        {
          if (py > bymax || py < bymin) return false;
        }

        // If dmin == 0, at this point, [t1, t2] == ]-inf, +inf[, or t1 or t2
        // is a NaN. But the case with NaNs is treated as if the interval
        // [t1, t2] was ]-inf, +inf[.

        CGAL_assertion(dymin >= 0);
        CGAL_assertion(dymax >= 0);
        if (bounded_0) {
          CGAL_assertion(tymin >= 0);
          CGAL_assertion(tymax >= 0);
        }


        // -----------------------------------
        // treat z coord
        // -----------------------------------
        // (Logic is the same as X)
        CFT dzmin, tzmin, tzmax, dzmax;
        if (qz >= pz) {
          if (bounded_0 && pz > bzmax) return false; // segment on the right of bbox
          if (bounded_1 && qz < bzmin) return false; // segment on the left of bbox

          if (bounded_1 && bzmax > qz) {
            tzmax = 1;
            dzmax = 1;
          } else {
            tzmax = bzmax - pz;
            dzmax = qz - pz;
          }

          tzmin = bzmin - pz;
          dzmin = qz - pz;
        } else {
          if (bounded_1 && qz > bzmax) return false; // segment on the right of bbox
          if (bounded_0 && pz < bzmin) return false; // segment on the left of bbox

          if (bounded_1 && bzmin < qz) {
            tzmax = 1;
            dzmax = 1;
          } else {
            tzmax = pz - bzmin;
            dzmax = pz - qz;
          }

          tzmin = pz - bzmax;
          dzmin = pz - qz;
        }

        if (bounded_0) tzmin = (CGAL::max)(CFT(0), tzmin);

        // If the query is vertical for z, then check its z-coordinate is in
        // the z-slab.
        if ((pz == qz) &&  // <=> (dmin == 0)
            (!(bounded_0 && bounded_1))) // do not check for a segment
        {
          if (pz > bzmax || pz < bzmin) return false;
        }

        // If dmin == 0, at this point, [t1, t2] == ]-inf, +inf[, or t1 or t2
        // is a NaN. But the case with NaNs is treated as if the interval
        // [t1, t2] was ]-inf, +inf[.

        CGAL_assertion(dzmin >= 0);
        CGAL_assertion(dzmax >= 0);
        if (bounded_0) {
          CGAL_assertion(tzmin >= 0);
          CGAL_assertion(tzmax >= 0);
        }


        // Use this custom type to determine our error bounds
        typedef Do_intersect_bbox_segment_aux_is_greater
                <CFT,
                        bounded_0,
                        use_static_filters
                > Is_greater;
        typedef typename Is_greater::result_type Is_greater_value;
        Is_greater is_greater;

        // Find the largest input values that will be used (each of these updates the internal values if larger)
        is_greater.register_new_input_values(tmin, dmin);
        is_greater.register_new_input_values(tymin, dymin);
        is_greater.register_new_input_values(tmax, dmax);
        is_greater.register_new_input_values(tymax, dymax);

        // Calculate the appropriate margin of error given the previous inputs
        is_greater.compute_new_error_bound();

        // If the values are extremely large or small, the intersection cannot be determined with exactness guarantees
        if (is_greater.bound_overflow() || is_greater.value_might_underflow())
          return Is_greater::uncertain();

        // If the ray has nonzero x and y angle
        if (py != qy && px != qx) { // dmin > 0, dymax >0, dmax > 0, dymin > 0

          // If the ymin could be greater than the ymax, the ray can't be proved to intersect
          const Is_greater_value b1 = is_greater(dymax * tmin, dmin * tymax);
          if (possibly(b1)) return !b1; // if(is_greater) return false; // or uncertain
          const Is_greater_value b2 = is_greater(dmax * tymin, dymin * tmax);
          if (possibly(b2)) return !b2;
        }

        // If the ray has nonzero y angle, but no x angle, and we're sure y's min is greater than the current min bound
        Is_greater_value b = Is_greater_value();
        if ((px == qx) || ((py != qy) &&
                           certainly(b = is_greater(dmin * tymin, dymin * tmin)))) {
          // Update bounds to y's min
          tmin = tymin;
          dmin = dymin;
        }

        // If we're uncertain about whether our bounds should shrink, we can't give a certain answer
        if (is_indeterminate(b)) return b; // Note that the default-constructed
        // Is_greater_value cannot be
        // indeterminate.

        // If tymax/dymax < t2, set t2 = tymax/dymax.
        // If the ray has nonzero y angle, but no x angle, and we're sure y's max is less than the current max bound
        if ((px == qx) || ((py != qy) &&
                           certainly(b = is_greater(dymax * tmax, dmax * tymax)))) {
          // Update bounds to y's max
          tmax = tymax;
          dmax = dymax;
        }

        // If we're uncertain about whether our bounds should shrink, we can't give a certain answer
        if (is_indeterminate(b)) return b;

        // Our bounds should not be negative
        CGAL_assertion(dmin >= 0);
        CGAL_assertion(dmax >= 0);

        // If we have a nonzero angle in the x or y directions as well as z
        if ((px != qx || py != qy) &&
            (pz != qz)) {
          // Update our margin of error
          is_greater.register_new_input_values(tzmin, dzmin);
          is_greater.register_new_input_values(tzmax, dzmax);
          is_greater.compute_new_error_bound();

          // If our z bounds are too large or small, we can't give a certain answer
          if (is_greater.bound_overflow() || is_greater.value_might_underflow())
            return Is_greater::uncertain();

          // If the zmin could be greater than the zmax, the ray can't be proved to intersect
          const Is_greater_value b1 = is_greater(dzmax * tmin, dmin * tzmax);
          if (possibly(b1)) return !b1; // if(is_greater) return false; // or uncertain
          const Is_greater_value b2 = is_greater(dmax * tzmin, dzmin * tmax);
          if (possibly(b2)) return !b2; // if(is_greater) return false; // or uncertain
        }

        // If none of our previous checks eliminated the possibility, than the ray must intersect
        return true;
      }

      template<
              typename FT,
              typename BFT,
              bool bounded_0,
              bool bounded_1,
              bool use_static_filters
      >
      typename Do_intersect_bbox_segment_aux_is_greater<
              FT,
              bounded_0,
              use_static_filters
      >::result_type
      do_intersect_bbox_segment_aux(
              const FT &px, const FT &py, const FT &pz,
              const FT &qx, const FT &qy, const FT &qz,
              const BFT &bxmin, const BFT &bymin, const BFT &bzmin,
              const BFT &bxmax, const BFT &bymax, const BFT &bzmax) {

        const BFT bx[] = {bxmin, bxmax};
        const BFT by[] = {bymin, bymax};
        const BFT bz[] = {bzmin, bzmax};

        const int sx = qx < 0;
        const int sy = qy < 0;
        const int sz = qz < 0;

        // Determine bounds x, y, and z
        FT xmin = (bx[sx] - px) / (qx - px);
        FT xmax = (bx[1 - sx] - px) / (qx - px);
        FT ymin = (by[sy] - py) / (qy - py);
        FT ymax = (by[1 - sy] - py) / (qy - py);
        FT zmin = (bz[sz] - pz) / (qz - pz);
        FT zmax = (bz[1 - sz] - pz) / (qz - pz);

        // Determine the bounds of the overlapping region
        FT min = std::max({xmin, ymin, zmin});
        FT max = std::min({xmax, ymax, zmax});

        // The ray intercepts if this region overlaps with the interval provided
        Uncertain<bool> new_result =
                (max >= min) &            // Check if there is any overlap at all
                !(bounded_0 && max < 0) & // Check if the overlap is outside the ray (before the start)
                !(bounded_1 && min > 1);  // Check if the overlap is outside the ray (after the end)

//        // Calculate the result using the old approach
//        auto old_result = do_intersect_bbox_segment_aux_old<FT, BFT, bounded_0, bounded_1, use_static_filters>(
//                px, py, pz,
//                qx, qy, qz,
//                bxmin, bymin, bzmin,
//                bxmax, bymax, bzmax
//        );
//
//        std::cout << (is_certain(old_result) ? (int) old_result : 2) << " "
//                  << (is_certain(new_result) ? (int) new_result : 2) << std::endl;
//
//        // Check that the two results are the same
//        CGAL_assertion_msg(
//                possibly(old_result) == possibly(new_result),
//                "The new ray-bbox intersection function produced an incorrect result!"
//        );
//
//        CGAL_assertion_msg(
//                is_certain(old_result) == is_certain(new_result),
//                "The new ray-bbox intersection function produced a result with the wrong level of certainty!"
//        );

        return new_result;
      }


      template<typename FT,
              bool bounded_0,
              bool bounded_1,
              bool use_static_filters>
      inline
      typename Do_intersect_bbox_segment_aux_is_greater<
              FT,
              bounded_0,
              use_static_filters
      >::result_type
      do_intersect_bbox_segment_aux(
              const FT &px, const FT &py, const FT &pz,
              const FT &qx, const FT &qy, const FT &qz,
              const Bbox_3 &bb) {
        return do_intersect_bbox_segment_aux<FT, double, bounded_0, bounded_1, use_static_filters>(px, py, pz,
                                                                                                   qx, qy, qz,
                                                                                                   bb.xmin(), bb.ymin(),
                                                                                                   bb.zmin(),
                                                                                                   bb.xmax(), bb.ymax(),
                                                                                                   bb.zmax());
      }


      template<class K>
      bool do_intersect(const typename K::Segment_3 &segment,
                        const CGAL::Bbox_3 &bbox,
                        const K &) {
        typedef typename K::FT FT;
        typedef typename K::Point_3 Point_3;

        const Point_3 &source = segment.source();
        const Point_3 &target = segment.target();

        return do_intersect_bbox_segment_aux<FT, true, true, false>(
                source.x(), source.y(), source.z(),
                target.x(), target.y(), target.z(),
                bbox);
      }

      template<class K>
      bool do_intersect(const CGAL::Bbox_3 &bbox,
                        const typename K::Segment_3 &segment,
                        const K &k) {
        return do_intersect(segment, bbox, k);
      }

    } // namespace internal
  } // namespace Intersections
} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif  // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_DO_INTERSECT_H
