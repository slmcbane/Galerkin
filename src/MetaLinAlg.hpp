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
#include "Rationals.hpp"

#include <tuple>

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

/// Swap the rows given by indices `I` and `J` in the given matrix.
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

/*!
 * @brief Find a row of the matrix with a non-zero on the given diagonal.
 * 
 * @tparam D The index of the column where we want a non-zero.
 * @tparam FIRST The index of the first row where we will look.
 */
template <auto COL, auto FIRST, class... Rows>
constexpr auto find_nonzero_entry(Matrix<Rows...>) noexcept
{
    static_assert(FIRST < sizeof...(Rows), "Non-zero pivot not found!");
    constexpr auto first_row = get_row<FIRST>(Matrix<Rows...>());
    if constexpr (get<COL>(first_row) != Rationals::rational<0>)
    {
        return FIRST;
    }
    else
    {
        return find_nonzero_entry<COL, FIRST+1>(Matrix<Rows...>());
    }
}

/********************************************************************************
 * Test the ability to find a row with non-zero diagonal element.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::MetaLinAlg] Find non-zero pivot element")
{
    constexpr auto mat = make_matrix(
        make_row(rational<1, 1>, rational<0, 2>, rational<1, 3>),
        make_row(rational<2, 1>, rational<0, 2>, rational<2, 3>),
        make_row(rational<3, 1>, rational<3, 2>, rational<3, 3>)
    );

    REQUIRE(find_nonzero_entry<0, 0>(mat) == 0);
    REQUIRE(find_nonzero_entry<1, 0>(mat) == 2);
    REQUIRE(find_nonzero_entry<0, 1>(mat) == 1);
    REQUIRE(find_nonzero_entry<1, 1>(mat) == 2);
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 *******************************************************************************/

template <auto I, auto N, auto D, class... Nums>
constexpr auto set_element(MatrixRow<Nums...> row,
                           Rationals::Rational<N, D>) noexcept
{
    static_assert(I >= 0 && I < sizeof...(Nums));
    if constexpr (I == 0)
    {
        return MatrixRow<Rationals::Rational<N, D>>().append(row.tail());
    }
    else
    {
        return make_row(row.head()).append(set_element<I-1>(row.tail(),
                                                            Rationals::Rational<N, D>()));
    }
}

template <auto I>
constexpr auto make_identity()
{
    [[maybe_unused]]
    constexpr auto zero_row = repeatedly<I>(Rationals::rational<0>);
    return static_reduce<0, I, 1>(
        [=](auto J) {
            return set_element<J()>(zero_row, Rationals::rational<1>);
        },
        MatrixRow<>(),
        [](auto rows, auto row) {
            return rows.append(make_list(row));
        }
    );
}

template <auto I>
constexpr auto eye = make_identity<I>();

template <auto I, class Num, class... Row1Nums, class... Row2Nums>
constexpr auto add_row_helper(Num, MatrixRow<Row1Nums...>,
                              MatrixRow<Row2Nums...>) noexcept
{
    if constexpr (I == sizeof...(Row1Nums))
    {
        return Matrix<Row2Nums...>();
    }
    else
    {
        constexpr auto y = Matrix<Row2Nums...>();
        constexpr auto elt = get<I>(MatrixRow<Row1Nums...>()) * Num();
        constexpr auto row = set_element<I>(y, get<I>(y) + elt);
        return add_row_helper<I+1>(Num(), MatrixRow<Row1Nums...>(), row);
    }
}

template <class Num, class... Row1Nums, class... Row2Nums>
constexpr auto add_row_with_mult(Num a, MatrixRow<Row1Nums...> x,
                                 MatrixRow<Row2Nums...> y) noexcept
{
    static_assert(sizeof...(Row1Nums) == sizeof...(Row2Nums));
    return add_row_helper<0>(a, x, y);
}

template <auto I, auto J, class... URows, class... LRows>
constexpr auto eliminate_row(Matrix<URows...> U, Matrix<LRows...> L)
{
    static_assert(J > I);
    auto mult = get_elt<J, I>(U) / get_elt<I, I>(U);
    auto newurow = add_row_with_mult(-mult, get_row<I>(U), get_row<J>(U));
    auto newu = replace_row<J>(U, newurow);

    auto newl = replace_row<J>(L, set_element<I>(get_row<J>(L), mult));
    return std::make_tuple(newu, newl);
}

/********************************************************************************
 * Test eliminating a row from a matrix (in the LU process). Rather than do this
 * in place I keep a separate L and U around for ease of implementation.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::MetaLinAlg] Test the identity function")
{
    REQUIRE(eye<3> == make_matrix(
        make_row(rational<1>, rational<0>, rational<0>),
        make_row(rational<0>, rational<1>, rational<0>),
        make_row(rational<0>, rational<0>, rational<1>)
    ));
}

TEST_CASE("[Galerkin::MetaLinAlg] Do row elimination on L and U factors")
{
    SUBCASE("A simple 2 x 2 matrix")
    {
        constexpr auto L = make_matrix(
            make_row(rational<2>, rational<3>),
            make_row(rational<1>, rational<4>)
        );

        auto [L1, U] = eliminate_row<0, 1>(L, eye<2>);
        REQUIRE(L1 == Matrix<
            MatrixRow<Rational<2, 1>, Rational<3, 1>>,
            MatrixRow<Rational<0, 1>, Rational<5, 2>>>());
        REQUIRE(U == Matrix<
            MatrixRow<Rational<1, 1>, Rational<0, 1>>,
            MatrixRow<Rational<1, 2>, Rational<1, 1>>>());
    }
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

} /* namespace Galerkin */

} /* namespace MetaLinAlg */

#endif /* METALINALG_HPP */