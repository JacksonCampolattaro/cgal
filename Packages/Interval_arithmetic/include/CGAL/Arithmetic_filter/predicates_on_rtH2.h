// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of the
// Computational Geometry Algorithms Library (CGAL).
//
// Every use of CGAL requires a license. Licenses come in three kinds:
//
// - For academic research and teaching purposes, permission to use and
//   copy the software and its documentation is hereby granted free of  
//   charge, provided that
//   (1) it is not a component of a commercial product, and
//   (2) this notice appears in all copies of the software and
//       related documentation.
// - Development licenses grant access to the source code of the library 
//   to develop programs. These programs may be sold to other parties as 
//   executable code. To obtain a development license, please contact
//   the CGAL Consortium (at cgal@cs.uu.nl).
// - Commercialization licenses grant access to the source code and the
//   right to sell development licenses. To obtain a commercialization 
//   license, please contact the CGAL Consortium (at cgal@cs.uu.nl).
//
// This software and documentation is provided "as-is" and without
// warranty of any kind. In no event shall the CGAL Consortium be
// liable for any damage of any kind.
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Free University of Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany) Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : 
// file          : include/CGAL/Arithmetic_filter/predicates_on_rtH2.h
// package       : Interval_arithmetic
// revision      : 1.3.1
// revision_date :
// author(s)     : Sylvain.Pion@sophia.inria.fr
//
// coordinator   : MPI, Saarbruecken
//
// email         : cgal@cs.uu.nl
//
// ======================================================================


#ifndef CGAL_ARITHMETIC_FILTER_PREDICATES_ON_RTH2_H
#define CGAL_ARITHMETIC_FILTER_PREDICATES_ON_RTH2_H

// This file is automatically generated with the script for filtering
// predicates, using Interval arithmetic.

#include <CGAL/Interval_arithmetic.h>

template < class CT, class ET >
// CGAL_KERNEL_INLINE
CGAL_Orientation
CGAL_orientationH2( const CGAL_Filtered_exact<CT,ET>& phx,
	const CGAL_Filtered_exact<CT,ET>& phy,
	const CGAL_Filtered_exact<CT,ET>& phw,
	                    const CGAL_Filtered_exact<CT,ET>& qhx,
	const CGAL_Filtered_exact<CT,ET>& qhy,
	const CGAL_Filtered_exact<CT,ET>& qhw,
	                    const CGAL_Filtered_exact<CT,ET>& rhx,
	const CGAL_Filtered_exact<CT,ET>& rhy,
	const CGAL_Filtered_exact<CT,ET>& rhw )
{
  CGAL_Orientation result;
  CGAL_FPU_set_rounding_to_infinity();
  try
  {
    result = CGAL_orientationH2(
		phx.interval(),
		phy.interval(),
		phw.interval(),
		qhx.interval(),
		qhy.interval(),
		qhw.interval(),
		rhx.interval(),
		rhy.interval(),
		rhw.interval());
    CGAL_FPU_set_rounding_to_nearest();
  } 
  catch (CGAL_Interval_nt_advanced::unsafe_comparison)
  {
    CGAL_FPU_set_rounding_to_nearest();
    result = CGAL_orientationH2(
		phx.exact(),
		phy.exact(),
		phw.exact(),
		qhx.exact(),
		qhy.exact(),
		qhw.exact(),
		rhx.exact(),
		rhy.exact(),
		rhw.exact());
  }
  return result;
}

template < class CT, class ET >
// CGAL_KERNEL_INLINE
bool
CGAL_leftturnH2( const CGAL_Filtered_exact<CT,ET>& phx,
	const CGAL_Filtered_exact<CT,ET>& phy,
	const CGAL_Filtered_exact<CT,ET>& phw,
	                 const CGAL_Filtered_exact<CT,ET>& qhx,
	const CGAL_Filtered_exact<CT,ET>& qhy,
	const CGAL_Filtered_exact<CT,ET>& qhw,
	                 const CGAL_Filtered_exact<CT,ET>& rhx,
	const CGAL_Filtered_exact<CT,ET>& rhy,
	const CGAL_Filtered_exact<CT,ET>& rhw )
{
  bool result;
  CGAL_FPU_set_rounding_to_infinity();
  try
  {
    result = CGAL_leftturnH2(
		phx.interval(),
		phy.interval(),
		phw.interval(),
		qhx.interval(),
		qhy.interval(),
		qhw.interval(),
		rhx.interval(),
		rhy.interval(),
		rhw.interval());
    CGAL_FPU_set_rounding_to_nearest();
  } 
  catch (CGAL_Interval_nt_advanced::unsafe_comparison)
  {
    CGAL_FPU_set_rounding_to_nearest();
    result = CGAL_leftturnH2(
		phx.exact(),
		phy.exact(),
		phw.exact(),
		qhx.exact(),
		qhy.exact(),
		qhw.exact(),
		rhx.exact(),
		rhy.exact(),
		rhw.exact());
  }
  return result;
}

template < class CT, class ET >
// CGAL_KERNEL_INLINE
bool
CGAL_rightturnH2(const CGAL_Filtered_exact<CT,ET>& phx,
	const CGAL_Filtered_exact<CT,ET>& phy,
	const CGAL_Filtered_exact<CT,ET>& phw,
	                 const CGAL_Filtered_exact<CT,ET>& qhx,
	const CGAL_Filtered_exact<CT,ET>& qhy,
	const CGAL_Filtered_exact<CT,ET>& qhw,
	                 const CGAL_Filtered_exact<CT,ET>& rhx,
	const CGAL_Filtered_exact<CT,ET>& rhy,
	const CGAL_Filtered_exact<CT,ET>& rhw )
{
  bool result;
  CGAL_FPU_set_rounding_to_infinity();
  try
  {
    result = CGAL_rightturnH2(
		phx.interval(),
		phy.interval(),
		phw.interval(),
		qhx.interval(),
		qhy.interval(),
		qhw.interval(),
		rhx.interval(),
		rhy.interval(),
		rhw.interval());
    CGAL_FPU_set_rounding_to_nearest();
  } 
  catch (CGAL_Interval_nt_advanced::unsafe_comparison)
  {
    CGAL_FPU_set_rounding_to_nearest();
    result = CGAL_rightturnH2(
		phx.exact(),
		phy.exact(),
		phw.exact(),
		qhx.exact(),
		qhy.exact(),
		qhw.exact(),
		rhx.exact(),
		rhy.exact(),
		rhw.exact());
  }
  return result;
}

template < class CT, class ET >
// CGAL_KERNEL_INLINE
bool
CGAL_collinearH2(const CGAL_Filtered_exact<CT,ET>& phx,
	const CGAL_Filtered_exact<CT,ET>& phy,
	const CGAL_Filtered_exact<CT,ET>& phw,
	                 const CGAL_Filtered_exact<CT,ET>& qhx,
	const CGAL_Filtered_exact<CT,ET>& qhy,
	const CGAL_Filtered_exact<CT,ET>& qhw,
	                 const CGAL_Filtered_exact<CT,ET>& rhx,
	const CGAL_Filtered_exact<CT,ET>& rhy,
	const CGAL_Filtered_exact<CT,ET>& rhw )
{
  bool result;
  CGAL_FPU_set_rounding_to_infinity();
  try
  {
    result = CGAL_collinearH2(
		phx.interval(),
		phy.interval(),
		phw.interval(),
		qhx.interval(),
		qhy.interval(),
		qhw.interval(),
		rhx.interval(),
		rhy.interval(),
		rhw.interval());
    CGAL_FPU_set_rounding_to_nearest();
  } 
  catch (CGAL_Interval_nt_advanced::unsafe_comparison)
  {
    CGAL_FPU_set_rounding_to_nearest();
    result = CGAL_collinearH2(
		phx.exact(),
		phy.exact(),
		phw.exact(),
		qhx.exact(),
		qhy.exact(),
		qhw.exact(),
		rhx.exact(),
		rhy.exact(),
		rhw.exact());
  }
  return result;
}
template < class CT, class ET >
// CGAL_KERNEL_INLINE
CGAL_Bounded_side
CGAL_side_of_bounded_circleH2( const CGAL_Filtered_exact<CT,ET>& qhx,
	const CGAL_Filtered_exact<CT,ET>& qhy,
	const CGAL_Filtered_exact<CT,ET>& qhw,
	                               const CGAL_Filtered_exact<CT,ET>& rhx,
	const CGAL_Filtered_exact<CT,ET>& rhy,
	const CGAL_Filtered_exact<CT,ET>& rhw,
	                               const CGAL_Filtered_exact<CT,ET>& shx,
	const CGAL_Filtered_exact<CT,ET>& shy,
	const CGAL_Filtered_exact<CT,ET>& shw,
	                               const CGAL_Filtered_exact<CT,ET>& thx,
	const CGAL_Filtered_exact<CT,ET>& thy,
	const CGAL_Filtered_exact<CT,ET>& thw )
{
  CGAL_Bounded_side result;
  CGAL_FPU_set_rounding_to_infinity();
  try
  {
    result = CGAL_side_of_bounded_circleH2(
		qhx.interval(),
		qhy.interval(),
		qhw.interval(),
		rhx.interval(),
		rhy.interval(),
		rhw.interval(),
		shx.interval(),
		shy.interval(),
		shw.interval(),
		thx.interval(),
		thy.interval(),
		thw.interval());
    CGAL_FPU_set_rounding_to_nearest();
  } 
  catch (CGAL_Interval_nt_advanced::unsafe_comparison)
  {
    CGAL_FPU_set_rounding_to_nearest();
    result = CGAL_side_of_bounded_circleH2(
		qhx.exact(),
		qhy.exact(),
		qhw.exact(),
		rhx.exact(),
		rhy.exact(),
		rhw.exact(),
		shx.exact(),
		shy.exact(),
		shw.exact(),
		thx.exact(),
		thy.exact(),
		thw.exact());
  }
  return result;
}

template < class CT, class ET >
// CGAL_KERNEL_INLINE
CGAL_Oriented_side
CGAL_side_of_oriented_circleH2(const CGAL_Filtered_exact<CT,ET>& qhx,
	const CGAL_Filtered_exact<CT,ET>& qhy,
	const CGAL_Filtered_exact<CT,ET>& qhw,
	                               const CGAL_Filtered_exact<CT,ET>& rhx,
	const CGAL_Filtered_exact<CT,ET>& rhy,
	const CGAL_Filtered_exact<CT,ET>& rhw,
	                               const CGAL_Filtered_exact<CT,ET>& shx,
	const CGAL_Filtered_exact<CT,ET>& shy,
	const CGAL_Filtered_exact<CT,ET>& shw,
	                               const CGAL_Filtered_exact<CT,ET>& thx,
	const CGAL_Filtered_exact<CT,ET>& thy,
	const CGAL_Filtered_exact<CT,ET>& thw )
{
  CGAL_Oriented_side result;
  CGAL_FPU_set_rounding_to_infinity();
  try
  {
    result = CGAL_side_of_oriented_circleH2(
		qhx.interval(),
		qhy.interval(),
		qhw.interval(),
		rhx.interval(),
		rhy.interval(),
		rhw.interval(),
		shx.interval(),
		shy.interval(),
		shw.interval(),
		thx.interval(),
		thy.interval(),
		thw.interval());
    CGAL_FPU_set_rounding_to_nearest();
  } 
  catch (CGAL_Interval_nt_advanced::unsafe_comparison)
  {
    CGAL_FPU_set_rounding_to_nearest();
    result = CGAL_side_of_oriented_circleH2(
		qhx.exact(),
		qhy.exact(),
		qhw.exact(),
		rhx.exact(),
		rhy.exact(),
		rhw.exact(),
		shx.exact(),
		shy.exact(),
		shw.exact(),
		thx.exact(),
		thy.exact(),
		thw.exact());
  }
  return result;
}
template < class CT, class ET >
// CGAL_KERNEL_INLINE
CGAL_Comparison_result
CGAL_compare_lexicographically_xyH2(const CGAL_Filtered_exact<CT,ET>& phx,
	const CGAL_Filtered_exact<CT,ET>& phy,
	const CGAL_Filtered_exact<CT,ET>& phw,
	                                    const CGAL_Filtered_exact<CT,ET>& qhx,
	const CGAL_Filtered_exact<CT,ET>& qhy,
	const CGAL_Filtered_exact<CT,ET>& qhw)
{
  CGAL_Comparison_result result;
  CGAL_FPU_set_rounding_to_infinity();
  try
  {
    result = CGAL_compare_lexicographically_xyH2(
		phx.interval(),
		phy.interval(),
		phw.interval(),
		qhx.interval(),
		qhy.interval(),
		qhw.interval());
    CGAL_FPU_set_rounding_to_nearest();
  } 
  catch (CGAL_Interval_nt_advanced::unsafe_comparison)
  {
    CGAL_FPU_set_rounding_to_nearest();
    result = CGAL_compare_lexicographically_xyH2(
		phx.exact(),
		phy.exact(),
		phw.exact(),
		qhx.exact(),
		qhy.exact(),
		qhw.exact());
  }
  return result;
}

template < class CT, class ET >
// CGAL_KERNEL_INLINE
CGAL_Comparison_result
CGAL_compare_xH2( const CGAL_Filtered_exact<CT,ET>& phx,
	const CGAL_Filtered_exact<CT,ET>& phw,
	                  const CGAL_Filtered_exact<CT,ET>& qhx,
	const CGAL_Filtered_exact<CT,ET>& qhw )
{
  CGAL_Comparison_result result;
  CGAL_FPU_set_rounding_to_infinity();
  try
  {
    result = CGAL_compare_xH2(
		phx.interval(),
		phw.interval(),
		qhx.interval(),
		qhw.interval());
    CGAL_FPU_set_rounding_to_nearest();
  } 
  catch (CGAL_Interval_nt_advanced::unsafe_comparison)
  {
    CGAL_FPU_set_rounding_to_nearest();
    result = CGAL_compare_xH2(
		phx.exact(),
		phw.exact(),
		qhx.exact(),
		qhw.exact());
  }
  return result;
}

// No CGAL_compare_yH2; use CGAL_compare_xH2( py, pw, qy, qw)

template < class CT, class ET >
// CGAL_KERNEL_INLINE
CGAL_Comparison_result
CGAL_compare_deltax_deltayH2(const CGAL_Filtered_exact<CT,ET>& phx,
	const CGAL_Filtered_exact<CT,ET>& phw,
	                             const CGAL_Filtered_exact<CT,ET>& qhx,
	const CGAL_Filtered_exact<CT,ET>& qhw,
	                             const CGAL_Filtered_exact<CT,ET>& rhy,
	const CGAL_Filtered_exact<CT,ET>& rhw,
	                             const CGAL_Filtered_exact<CT,ET>& shy,
	const CGAL_Filtered_exact<CT,ET>& shw )
  const CGAL_Filtered_exact<CT,ET>  tbc1 = CGAL_abs(phx*qhw - qhx*phw) * rhw*shw;
  const CGAL_Filtered_exact<CT,ET>  tbc2 = CGAL_abs(rhy*shw - shy*rhw) * phw*qhw;
  return (tbc2 < tbc1) ? CGAL_LARGER
                       : (tbc1 == tbc2) ? CGAL_EQUAL : CGAL_SMALLER;


template <class CGAL_Filtered_exact<CT,ET>>
// CGAL_KERNEL_INLINE
bool
CGAL_collinear_are_ordered_along_lineH2(
     const CGAL_Filtered_exact<CT,ET>& phx,
	const CGAL_Filtered_exact<CT,ET>& phy,
	const CGAL_Filtered_exact<CT,ET>& phw,
	     const CGAL_Filtered_exact<CT,ET>& qhx,
	const CGAL_Filtered_exact<CT,ET>& qhy,
	const CGAL_Filtered_exact<CT,ET>& qhw,
	     const CGAL_Filtered_exact<CT,ET>& rhx,
	const CGAL_Filtered_exact<CT,ET>& rhy,
	const CGAL_Filtered_exact<CT,ET>& rhw
                                       )
{
  CGAL_Comparison_result result;
  CGAL_FPU_set_rounding_to_infinity();
  try
  {
    result = CGAL_collinear_are_ordered_along_lineH2(
		phx.interval(),
		phw.interval(),
		qhx.interval(),
		qhw.interval(),
		rhy.interval(),
		rhw.interval(),
		shy.interval(),
		shw.interval(),
		phw.interval(),
		s.interval(),
		rhw.interval(),
		q.interval(),
		tbc1.interval(),
		CGAL_LARGER.interval(),
		tbc2.interval(),
		C.interval(),
		.interval());		CGAL_KERNEL_INLINE.interval(),
		bool.interval(),
		phx.interval(),
		phy.interval(),
		phw.interval(),
		qhx.interval(),
		qhy.interval(),
		qhw.interval(),
		rhx.interval(),
		rhy.interval(),
		rhw.interval(),
		.interval());
    CGAL_FPU_set_rounding_to_nearest();
  } 
  catch (CGAL_Interval_nt_advanced::unsafe_comparison)
  {
    CGAL_FPU_set_rounding_to_nearest();
    result = CGAL_collinear_are_ordered_along_lineH2(
		phx.exact(),
		phw.exact(),
		qhx.exact(),
		qhw.exact(),
		rhy.exact(),
		rhw.exact(),
		shy.exact(),
		shw.exact(),
		phw.exact(),
		s.exact(),
		rhw.exact(),
		q.exact(),
		tbc1.exact(),
		CGAL_LARGER.exact(),
		tbc2.exact(),
		C.exact(),
		.exact());		CGAL_KERNEL_INLINE.exact(),
		bool.exact(),
		phx.exact(),
		phy.exact(),
		phw.exact(),
		qhx.exact(),
		qhy.exact(),
		qhw.exact(),
		rhx.exact(),
		rhy.exact(),
		rhw.exact(),
		.exact());
  }
  return result;
}

template < class CT, class ET >
// CGAL_KERNEL_INLINE
bool
CGAL_collinear_are_strictly_ordered_along_lineH2(
     const CGAL_Filtered_exact<CT,ET>& phx,
	const CGAL_Filtered_exact<CT,ET>& phy,
	const CGAL_Filtered_exact<CT,ET>& phw,
	     const CGAL_Filtered_exact<CT,ET>& qhx,
	const CGAL_Filtered_exact<CT,ET>& qhy,
	const CGAL_Filtered_exact<CT,ET>& qhw,
	     const CGAL_Filtered_exact<CT,ET>& rhx,
	const CGAL_Filtered_exact<CT,ET>& rhy,
	const CGAL_Filtered_exact<CT,ET>& rhw)
{
  bool result;
  CGAL_FPU_set_rounding_to_infinity();
  try
  {
    result = CGAL_collinear_are_strictly_ordered_along_lineH2(
		phx.interval(),
		phy.interval(),
		phw.interval(),
		qhx.interval(),
		qhy.interval(),
		qhw.interval(),
		rhx.interval(),
		rhy.interval(),
		rhw.interval());
    CGAL_FPU_set_rounding_to_nearest();
  } 
  catch (CGAL_Interval_nt_advanced::unsafe_comparison)
  {
    CGAL_FPU_set_rounding_to_nearest();
    result = CGAL_collinear_are_strictly_ordered_along_lineH2(
		phx.exact(),
		phy.exact(),
		phw.exact(),
		qhx.exact(),
		qhy.exact(),
		qhw.exact(),
		rhx.exact(),
		rhy.exact(),
		rhw.exact());
  }
  return result;
}




#ifdef CGAL_ARITHMETIC_FILTER_H
#include <CGAL/Arithmetic_filter/predicates_on_rtH2.h>
#endif // CGAL_ARITHMETIC_FILTER_H

#endif // CGAL_ARITHMETIC_FILTER_PREDICATES_ON_RTH2_H
