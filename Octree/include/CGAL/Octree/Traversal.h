// Copyright (c) 2007-2020  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jackson Campolattaro, Cédric Portaneri, Tong Zhao

#ifndef CGAL_OCTREE_TRAVERSAL_CRITERION_H
#define CGAL_OCTREE_TRAVERSAL_CRITERION_H

#include <CGAL/license/Octree.h>

#include <iostream>
#include <boost/range/iterator_range.hpp>
#include "Node.h"
#include "Traversal_iterator.h"

namespace CGAL {

namespace Octree {

template<class Point_index, class Dimension>
const Node <Point_index, Dimension> *next_sibling(const Node <Point_index, Dimension> *n) {

  // Passing null returns the first node
  if (nullptr == n)
    return nullptr;

  // If this node has no parent, it has no siblings
  if (nullptr == n->parent())
    return nullptr;

  // Find out which child this is
  std::size_t index = n->index().to_ulong();

  constexpr static int degree = (2 << (Dimension::value-1));
  // Return null if this is the last child
  if (index == degree - 1)
    return nullptr;

  // Otherwise, return the next child
  return &((*n->parent())[index + 1]);
}

template<class Point_index, class Dimension>
const Node <Point_index, Dimension> *next_sibling_up(const Node <Point_index, Dimension> *n) {

  if (!n)
    return nullptr;

  auto up = n->parent();

  while (nullptr != up) {

    if (nullptr != next_sibling(up))
      return next_sibling(up);

    up = up->parent();
  }

  return nullptr;
}

template<class Point_index, class Dimension>
const Node <Point_index, Dimension> *deepest_first_child(const Node <Point_index, Dimension> *n) {

  if (!n)
    return nullptr;

  // Find the deepest child on the left
  auto first = n;
  while (!first->is_leaf())
    first = &(*first)[0];
  return first;
}

namespace Traversal {

/*!
 * \brief walker for preorder traversal
 */
struct Preorder {

  /*!
   * \brief retrieve the first node of a tree in a preorder traversal, given the root
   *
   * \tparam Point_index
   * \param root
   * \return
   */
  template<class Point_index, class Dimension>
  const Node <Point_index, Dimension> *first(const Node <Point_index, Dimension> *root) const {
    return root;
  }

  /*!
   * \brief retrieve the next node of a tree in a preorder traversal, given the current one
   *
   * \tparam Point_index
   * \param n
   * \return
   */
  template<class Point_index, class Dimension>
  const Node <Point_index, Dimension> *next(const Node <Point_index, Dimension> *n) const {

    if (n->is_leaf()) {

      auto next = next_sibling(n);

      if (nullptr == next) {

        return next_sibling_up(n);
      }

      return next;

    } else {

      // Return the first child of this node
      return &(*n)[0];
    }

  }
};

/*!
 * \todo
 */
struct Postorder {

  /*!
   * \brief retrieve the first node of a tree in a postorder traversal, given the root
   *
   * \tparam Point_index
   * \param root
   * \return
   */
  template<class Point_index, class Dimension>
  const Node <Point_index, Dimension> *first(const Node <Point_index, Dimension> *root) const {

    return deepest_first_child(root);
  }

  /*!
   * \brief retrieve the next node of a tree in a postorder traversal, given the current one
   *
   * \tparam Point_index
   * \param n
   * \return
   */
  template<class Point_index, class Dimension>
  const Node <Point_index, Dimension> *next(const Node <Point_index, Dimension> *n) const {

    auto next = deepest_first_child(next_sibling(n));

    if (!next)
      next = n->parent();

    return next;
  }
};

/*!
 * \brief walker for leaves-only traversal
 */
struct Leaves {

  /*!
   * \brief retrieve the first node of a tree in a leaves-only traversal, given the root
   *
   * \tparam Point_index
   * \param root
   * \return
   */
  template<class Point_index, class Dimension>
  const Node <Point_index, Dimension> *first(const Node <Point_index, Dimension> *root) const {

    return deepest_first_child(root);
  }

  /*!
   * \brief retrieve the next node of a tree in a leaves-only traversal, given the current one
   *
   * \tparam Point_index
   * \param n
   * \return
   */
  template<class Point_index, class Dimension>
  const Node <Point_index, Dimension> *next(const Node <Point_index, Dimension> *n) const {

    auto next = deepest_first_child(next_sibling(n));

    if (!next)
      next = deepest_first_child(next_sibling_up(n));

    return next;
  }
};

}

}

}

#endif //CGAL_OCTREE_TRAVERSAL_CRITERION_H
