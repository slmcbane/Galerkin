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
        MatrixRow<Rational<1, 3>, Rational<1, 2>>>();

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
        return make_matrix(mat.head()).append(replace_row<I - 1>(mat.tail(), Row()));
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
            return make_matrix(get_row<J>(mat)).append(replace_row<J - 1>(mat.tail(), mat.head()));
        }
        else
        {
            return make_matrix(mat.head()).append(swap_rows<I - 1, J - 1>(mat.tail()));
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
        make_row(rational<3, 1>, rational<3, 2>, rational<3, 3>));

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
        return find_nonzero_entry<COL, FIRST + 1>(Matrix<Rows...>());
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
        make_row(rational<3, 1>, rational<3, 2>, rational<3, 3>));

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
        return make_row(row.head()).append(set_element<I - 1>(row.tail(), Rationals::Rational<N, D>()));
    }
}

template <auto I>
constexpr auto make_identity()
{
    [[maybe_unused]] constexpr auto zero_row = repeatedly<I>(Rationals::rational<0>);
    return static_reduce<0, I, 1>(
        [=](auto J) {
            return set_element<J()>(zero_row, Rationals::rational<1>);
        },
        MatrixRow<>(),
        [](auto rows, auto row) {
            return rows.append(make_list(row));
        });
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
        return add_row_helper<I + 1>(Num(), MatrixRow<Row1Nums...>(), row);
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

template <auto I, auto J, int K, class Row, class... Rows>
constexpr auto eliminate_row_helper(Matrix<Rows...>, Row, std::integral_constant<int, K>)
{
    [[maybe_unused]] constexpr auto A = Matrix<Rows...>();
    constexpr auto row = Row();
    [[maybe_unused]] constexpr auto mult = get_elt<J, I>(A) / get_elt<I, I>(A);

    if constexpr (K == I)
    {
        constexpr auto newrow = set_element<I>(row, mult);
        return eliminate_row_helper<I, J>(A, newrow, std::integral_constant<int, K+1>());
    }
    else if constexpr (K == sizeof...(Rows))
    {
        return row;
    }
    else
    {
        constexpr auto newrow = set_element<K>(row,
            get_elt<J, K>(A) - get_elt<I, K>(A) * mult);
        return eliminate_row_helper<I, J>(A, newrow, std::integral_constant<int, K+1>());
    }
}

template <auto I, auto J, class... Rows>
constexpr auto eliminate_row(Matrix<Rows...>)
{
    static_assert(J > I);
    constexpr auto A = Matrix<Rows...>();
    constexpr auto newrow = eliminate_row_helper<I, J>(A, get_row<J>(A),
                                                       std::integral_constant<int, I>());
    return replace_row<J>(A, newrow);
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
                          make_row(rational<0>, rational<0>, rational<1>)));
}

TEST_CASE("[Galerkin::MetaLinAlg] Do row elimination on L and U factors")
{
    SUBCASE("A simple 2 x 2 matrix")
    {
        constexpr auto A = make_matrix(
            make_row(rational<2>, rational<3>),
            make_row(rational<1>, rational<4>));

        auto LU = eliminate_row<0, 1>(A);
        REQUIRE(LU == Matrix<
                             MatrixRow<Rational<2, 1>, Rational<3, 1>>,
                             MatrixRow<Rational<1, 2>, Rational<5, 2>>>());
    }

    SUBCASE("A 3 x 3 matrix")
    {
        constexpr auto A = make_matrix(
            make_row(rational<2>, -rational<1>, rational<3>),
            make_row(rational<4>, rational<2>, rational<1>),
            make_row(-rational<6>, -rational<1>, rational<2>));

        auto LU = eliminate_row<0, 2>(eliminate_row<0, 1>(A));
        REQUIRE(LU == make_matrix(
            make_row(rational<2>, -rational<1>, rational<3>),
            make_row(rational<2>, rational<4>, -rational<5>),
            make_row(-rational<3>, -rational<4>, rational<11>)
        ));
        /* REQUIRE(L == make_matrix(
                         make_row(rational<2>, -rational<1>, rational<3>),
                         make_row(rational<0>, rational<4>, rational<-5>),
                         make_row(rational<0>, rational<-4>, rational<11>)));
        REQUIRE(U == make_matrix(
                         make_row(rational<1>, rational<0>, rational<0>),
                         make_row(rational<2>, rational<1>, rational<0>),
                         make_row(rational<-3>, rational<0>, rational<1>))); */
    }
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

/********************************************************************************
 *******************************************************************************/

template <int M, int N, class... URows, class... LRows>
constexpr auto do_row_elimination(Matrix<URows...>, Matrix<LRows...>)
{
    if constexpr (N == sizeof...(URows))
    {
        return std::make_tuple(Matrix<URows...>(), Matrix<LRows...>());
    }
    else
    {
        auto [U, L] = eliminate_row<M, N>(Matrix<URows...>(), Matrix<LRows...>());
        return do_row_elimination<M, N+1>(U, L);
    }
}

template <int M, int N, class... Rows>
constexpr auto do_row_elimination(Matrix<Rows...>) noexcept
{
    if constexpr (N == sizeof...(Rows))
    {
        return Matrix<Rows...>();
    }
    else
    {
        constexpr auto LU = eliminate_row<M, N>(Matrix<Rows...>());
        return do_row_elimination<M, N+1>(LU);
    }
}

template <class... Rows, int... swaps, int M>
constexpr auto factorize(Matrix<Rows...>,
                         typeconst_list<std::integral_constant<int, swaps>...>,
                         std::integral_constant<int, M>) noexcept
{
    if constexpr (M == sizeof...(Rows) - 1)
    {
        return std::tuple(Matrix<Rows...>(),
                          typeconst_list<std::integral_constant<int, swaps>...>());
    }
    else
    {
        constexpr auto swap = find_nonzero_entry<M, M>(Matrix<Rows...>());
        if constexpr (swap != M)
        {
            constexpr auto LU = swap_rows<M, swap>(Matrix<Rows...>());
            constexpr auto P = typeconst_list<std::integral_constant<int, swaps>...,
                                              std::integral_constant<int, swap>>();
            return factorize(LU, P, std::integral_constant<int, M>());
        }
        else
        {
            constexpr auto LU = do_row_elimination<M, M+1>(Matrix<Rows...>());
            // If number of swaps is the same as M we didn't do a swap for this
            // row yet; add the index of this row to swaps (indicates no swap).
            if constexpr (sizeof...(swaps) == M)
            {
                constexpr auto P = typeconst_list<std::integral_constant<int, swaps>...,
                                                  std::integral_constant<int, M>>();
                return factorize(LU, P, std::integral_constant<int, M+1>());
            }
            // otherwise, we did a swap and should not add anything to the swaps
            else
            {
                constexpr auto P = typeconst_list<std::integral_constant<int, swaps>...>();
                return factorize(LU, P, std::integral_constant<int, M+1>());
            }
        }
    }
}

template <class... URows, class... LRows, int... swaps, int M>
constexpr auto factorize(Matrix<URows...>, Matrix<LRows...>,
                         typeconst_list<std::integral_constant<int, swaps>...>,
                         std::integral_constant<int, M>) noexcept
{
    if constexpr (M == sizeof...(URows) - 1)
    {
        return std::make_tuple(Matrix<LRows...>(), Matrix<URows...>(),
                               typeconst_list<std::integral_constant<int, swaps>...>());
    }
    else
    {

        constexpr auto swap = find_nonzero_entry<M, M>(Matrix<URows...>());
        if constexpr (swap != M)
        {
            constexpr auto U = swap_rows<M, swap>(Matrix<URows...>());
            constexpr auto P = typeconst_list<std::integral_constant<int, swaps>...,
                                              std::integral_constant<int, swap>>();
            return factorize(U, Matrix<LRows...>(), P, std::integral_constant<int, M>());
        }
        else
        {
            auto [U, L] = do_row_elimination<M, M+1>(Matrix<URows...>(), Matrix<LRows...>());
            if constexpr (sizeof...(swaps) == M)
            {
                constexpr auto P = typeconst_list<std::integral_constant<int, swaps>...,
                                                  std::integral_constant<int, M>>();
                return factorize(U, L, P, std::integral_constant<int, M+1>());
            }
            else
            {
                constexpr auto P = typeconst_list<std::integral_constant<int, swaps>...>();
                return factorize(U, L, P, std::integral_constant<int, M+1>());
            }
        }
    }
}

template <class... Rows>
constexpr auto factorize(Matrix<Rows...>) noexcept
{
    constexpr auto P = typeconst_list<>();
    return factorize(Matrix<Rows...>(), P, std::integral_constant<int, 0>());
}

/********************************************************************************
 * Test LU factorization with pivoting to fix zero diagonals.
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[Galerkin::MetaLinAlg] Test LU factorization")
{
    SUBCASE("Test a matrix that doesn't require pivoting")
    {
        constexpr auto A = make_matrix(
            make_row(rational<2>, rational<-1>, rational<3>),
            make_row(rational<4>, rational<2>, rational<1>),
            make_row(rational<-6>, rational<-1>, rational<2>));

        // auto [L, U, P] = factorize(A);
        auto [LU, P] = factorize(A);

        REQUIRE(LU == make_matrix(
            make_row(rational<2>, -rational<1>, rational<3>),
            make_row(rational<2>, rational<4>, -rational<5>),
            make_row(-rational<3>, -rational<1>, rational<6>)
        ));

        REQUIRE(P == make_list(
            std::integral_constant<int, 0>(),
            std::integral_constant<int, 1>()));
    }

    SUBCASE("A more difficult case with pivoting")
    {
        
    }
}

#endif /* DOCTEST_LIBRARY_INCLUDED */

} // namespace MetaLinAlg

} // namespace Galerkin

#endif /* METALINALG_HPP */