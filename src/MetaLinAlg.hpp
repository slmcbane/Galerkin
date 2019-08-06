/*
 * Copyright (c) 2019, Sean McBane
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef METALINALG_HPP

/*!
 * @file MetaLinAlg.hpp
 * @brief A compile-time matrix facility intended to solve a linear system to
 * find coefficients of multinomial basis functions.
 */

#include "utils.hpp"

namespace Galerkin
{

namespace MetaLinAlg
{

/*!
 * @brief To store a row of a matrix we re-use the utility typeconst_list.
 * 
 * It is up to functions manipulating matrices to make sure that the types of
 * `Nums...` continue to make sense - they should all be rationals or integral
 * constants.
 */
template <class... Nums>
using MatrixRow = typeconst_list<Nums...>;

/*!
 * @brief A `Matrix` is just a `typeconst_list` of `MatrixRow`s.
 */
template <class... Rows>
using Matrix = typeconst_list<Rows...>;

/*!
 * @brief Make a row from a list of rationals.
 */
template <auto N, auto D, class... Nums>
constexpr auto make_row(Rationals::Rational<N, D>, Nums...) noexcept
{
    constexpr auto head = MatrixRow<Rationals::Rational<N, D>>();
    if constexpr (sizeof...(Nums) == 0)
    {
        return head;
    }
    else
    {
        return head.append(make_row(Nums()...));
    }
}

/*!
 * @brief Make a matrix from a list of rows.
 * 
 * This function requires and checks via `static_assert` that the length of each
 * row is the same.
 */
template <class Row, class... Rows>
constexpr auto make_matrix(Row, Rows...) noexcept
{
    constexpr auto head = Matrix<Row>();
    if constexpr (sizeof...(Rows) == 0)
    {
        return head;
    }
    else
    {
        constexpr auto tail = make_matrix(Rows()...);
        static_assert(Row().count() == tail.head().count(),
                      "All matrix rows must have the same length");
        return head.append(tail);
    }
}

/*!
 * @brief Return the specified row from `mat`, using a 0-based index.
 */
template <auto I, class... Rows>
constexpr auto get_row(Matrix<Rows...> mat) noexcept
{
    return get<I>(mat);
}

/*!
 * @brief Return the specified element from `mat`, using a 0-based index.
 */
template <auto I, auto J, class... Rows>
constexpr auto get_elt(Matrix<Rows...> mat) noexcept
{
    return get<J>(get_row<I>(mat));
}

/********************************************************************************
 * Test making matrices and accessing elements using this API.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

using namespace Galerkin::Rationals;

TEST_CASE("[Galerkin::MetaLinAlg] Construct a MetaMatrix, access elements and rows")
{
    constexpr auto mat = Matrix<
        MatrixRow<Rational<1, 2>, Rational<1, 3>>,
        MatrixRow<Rational<1, 3>, Rational<1, 2>>
    >();

    constexpr auto mat2 = make_matrix(make_row(rational<1, 2>, rational<2, 6>),
                                      make_row(rational<3, 9>, rational<1, 2>));

    REQUIRE(mat == mat2);
    REQUIRE(get_row<0>(mat) == make_row(rational<1, 2>, rational<2, 6>));
    REQUIRE(get_elt<0, 0>(mat) == rational<1, 2>);
    REQUIRE(get_elt<1, 0>(mat) == rational<1, 3>);

    // Should throw a static_assert; mismatched row lengths.
    /* constexpr auto discard = make_matrix(
        MatrixRow<Rational<1, 2>, Rational<1, 3>>(),
        MatrixRow<Rational<1, 1>, Rational<1, 1>, Rational<1, 1>>()
     ); */

    // Shouldn't compile; rows do not all consist of Rationals.
    // constexpr auto discard2 = make_matrix(make_row(std::integral_constant<int, 3>())):
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 * End of very basic matrix API tests.
 *******************************************************************************/

/// Replace row I (0-based) with given row in the matrix.
template <auto I, class Row, class... Rows>
constexpr auto replace_row(Matrix<Rows...> mat, Row) noexcept
{
    static_assert(I >= 0 && I < mat.count());
    if constexpr (I == 0)
    {
        return Matrix<Row>().append(mat.tail());
    }
    else
    {
        return make_matrix(mat.head()).append(replace_row<I-1>(mat.tail(), Row()));
    }
}

template <auto I, auto J, class... Rows>
constexpr auto swap_rows(Matrix<Rows...> mat)
{
    if constexpr (I > J)
    {
        return swap_rows<J, I>(Matrix<Rows...>());
    }
    else
    {
        if constexpr (I == 0)
        {
            return make_matrix(get_row<J>(mat)).append(replace_row<J-1>(mat.tail(), mat.head()));
        }
        else
        {
            return make_matrix(mat.head()).append(swap_rows<I-1, J-1>(mat.tail()));
        }
    }
    
}

/********************************************************************************
 * Test swapping 2 matrix rows; needed for LU pivoting when required.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::MetaLinAlg] Swap matrix rows")
{
    constexpr auto mat = make_matrix(
        make_row(rational<1, 1>, rational<1, 2>, rational<1, 3>),
        make_row(rational<2, 1>, rational<2, 2>, rational<2, 3>),
        make_row(rational<3, 1>, rational<3, 2>, rational<3, 3>)
    );

    constexpr auto swapped = swap_rows<1, 2>(mat);
    REQUIRE(get_row<2>(swapped) == get_row<1>(mat));
    REQUIRE(get_row<1>(swapped) == get_row<2>(mat));
    REQUIRE(get_row<0>(swapped) == get_row<0>(mat));
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 * End test of row swap.
 *******************************************************************************/

} /* namespace Galerkin */

} /* namespace MetaLinAlg */

#endif /* METALINALG_HPP */