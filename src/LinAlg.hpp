/*
 * Copyright (c) 2019, Sean McBane and The University of Texas at Austin.
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef GALERKIN_LINALG_HPP
#define GALERKIN_LINALG_HPP

/*!
 * @file MetaLinAlg.hpp
 * @brief A compile-time matrix facility intended to solve a linear system to
 * find coefficients of multinomial basis functions.
 */

#include "Rationals.hpp"

#include <array>
#include <cassert>
#include <tuple>

namespace Galerkin
{

namespace LinAlg
{

template <std::size_t N>
struct Matrix
{
    std::array<Rationals::Rational, N * N> data;

    constexpr const auto &operator()(int m, int n) const { return data[m * N + n]; }

    template <class T>
    constexpr Matrix set_entry(int m, int n, T val) const
    {
        std::array<Rationals::Rational, N * N> xs{data};
        xs[m * N + n] = Rationals::Rational(val);
        return Matrix(xs);
    }

    constexpr Matrix(const std::array<Rationals::Rational, N * N> &xs) : data{xs} {}

    constexpr Matrix(const Matrix &other) : data{other.data} {}

    constexpr Matrix() : data() {}

    constexpr Matrix &operator=(const Matrix &other) = default;
};

template <std::size_t N>
struct Vector
{
    std::array<Rationals::Rational, N> data;

    constexpr Vector(const std::array<Rationals::Rational, N> &xs) : data{xs} {}

    constexpr Vector(const Vector &other) : data{other.data} {}

    constexpr Vector() : data() {}

    constexpr const auto &operator()(int m) const { return data[m]; }

    template <class T>
    constexpr Vector set_entry(int m, T val) const
    {
        std::array<Rationals::Rational, N> xs{data};
        xs[m] = Rationals::Rational(val);
        return Vector(xs);
    }

    constexpr Vector &operator=(const Vector &other) = default;
};

template <auto N>
constexpr auto canonical(std::size_t i) noexcept
{
    if constexpr (N == 0)
    {
        return Vector<0>();
    }
    else
    {
        return Vector<N>().set_entry(i, 1);
    }
}

} // namespace LinAlg

/********************************************************************************
 * Test that canonical basis vector is correct.
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[LinAlg] Canonical basis vectors")
{
    constexpr auto c0 = LinAlg::canonical<3>(0);
    REQUIRE(c0(0) == Rationals::rational<1>);
    REQUIRE(c0(1) == Rationals::rational<0>);
    REQUIRE(c0(2) == c0(1));

    constexpr auto c1 = LinAlg::canonical<3>(1);
    REQUIRE(c1(0) == c0(1));
    REQUIRE(c1(1) == c0(0));
    REQUIRE(c1(2) == c0(1));

    constexpr auto c2 = LinAlg::canonical<3>(2);
    REQUIRE(c2(0) == c0(1));
    REQUIRE(c2(1) == c0(1));
    REQUIRE(c2(2) == c0(0));
}

#endif

/********************************************************************************
 *******************************************************************************/

namespace LinAlg
{

template <std::size_t N>
constexpr std::size_t find_pivot(const Matrix<N> &A, std::size_t col)
{
    for (std::size_t m = col; m < N; ++m)
    {
        if (A(m, col) != Rationals::rational<0>)
        {
            return m;
        }
    }
    return -1;
}

template <std::size_t N>
constexpr Matrix<N> eliminate_row(Matrix<N> A, std::size_t col, std::size_t row)
{
    assert(row > col);
    const auto mult = A(row, col) / A(col, col);

    A = A.set_entry(row, col, mult);
    for (std::size_t k = col + 1; k < N; ++k)
    {
        A = A.set_entry(row, k, A(row, k) - A(col, k) * mult);
    }
    return A;
}

template <std::size_t N>
constexpr Matrix<N> do_row_elimination(Matrix<N> A, std::size_t col)
{
    for (std::size_t m = col + 1; m < N; ++m)
    {
        A = eliminate_row(A, col, m);
    }
    return A;
}

template <std::size_t N>
constexpr Matrix<N> swap_rows(const Matrix<N> &A, std::size_t m, std::size_t n)
{
    Matrix<N> B = A;
    for (std::size_t j = 0; j < N; ++j)
    {
        B = B.set_entry(m, j, A(n, j))
            .set_entry(n, j, A(m, j));
    }
    return B;
}

template <std::size_t N>
constexpr std::tuple<Matrix<N>, std::array<std::size_t, N - 1>> factorize(Matrix<N> A)
{
    std::array<std::size_t, N - 1> swaps{0};
    for (std::size_t m = 0; m < N - 1; ++m)
    {
        auto swap = find_pivot(A, m);
        swaps[m] = swap;
        if (swap != m)
        {
            A = swap_rows(A, m, swap);
        }
        A = do_row_elimination(A, m);
    }
    return std::tuple(A, swaps);
}

} // namespace LinAlg

/********************************************************************************
 * Test row elimination
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[LinAlg] Do row elimination")
{

    SUBCASE("A simple 2 x 2 matrix")
    {
        constexpr auto A = LinAlg::Matrix<2>(std::array<Rationals::Rational, 4>{
            Rationals::rational<2>, Rationals::rational<3>, Rationals::rational<1>,
            Rationals::rational<4>});
        constexpr auto LU = LinAlg::eliminate_row(A, 0, 1);
        REQUIRE(LU(0, 0) == Rationals::rational<2>);
        REQUIRE(LU(0, 1) == Rationals::rational<3>);
        REQUIRE(LU(1, 0) == Rationals::rational<1, 2>);
        REQUIRE(LU(1, 1) == Rationals::rational<5, 2>);
    }

    SUBCASE("A 3 x 3 matrix")
    {
        constexpr auto A = LinAlg::Matrix<3>()
                               .set_entry(0, 0, 2)
                               .set_entry(0, 1, -1)
                               .set_entry(0, 2, 3)
                               .set_entry(1, 0, 4)
                               .set_entry(1, 1, 2)
                               .set_entry(1, 2, 1)
                               .set_entry(2, 0, -6)
                               .set_entry(2, 1, -1)
                               .set_entry(2, 2, 2);
        constexpr auto LU = eliminate_row(eliminate_row(A, 0, 1), 0, 2);
        REQUIRE(LU(0, 0) == 2);
        REQUIRE(LU(0, 1) == -1);
        REQUIRE(LU(0, 2) == 3);
        REQUIRE(LU(1, 0) == 2);
        REQUIRE(LU(1, 1) == 4);
        REQUIRE(LU(1, 2) == -5);
        REQUIRE(LU(2, 0) == -3);
        REQUIRE(LU(2, 1) == -4);
        REQUIRE(LU(2, 2) == 11);
    }

    SUBCASE("LU factorization with no pivoting")
    {
        constexpr auto A = LinAlg::Matrix<3>()
                               .set_entry(0, 0, 2)
                               .set_entry(0, 1, -1)
                               .set_entry(0, 2, 3)
                               .set_entry(1, 0, 4)
                               .set_entry(1, 1, 2)
                               .set_entry(1, 2, 1)
                               .set_entry(2, 0, -6)
                               .set_entry(2, 1, -1)
                               .set_entry(2, 2, 2);
        constexpr auto t = LinAlg::factorize(A);
        constexpr auto LU = std::get<0>(t);
        constexpr auto P = std::get<1>(t);
        REQUIRE(LU(0, 0) == 2);
        REQUIRE(LU(0, 1) == -1);
        REQUIRE(LU(0, 2) == 3);
        REQUIRE(LU(1, 0) == 2);
        REQUIRE(LU(1, 1) == 4);
        REQUIRE(LU(1, 2) == -5);
        REQUIRE(LU(2, 0) == -3);
        REQUIRE(LU(2, 1) == -1);
        REQUIRE(LU(2, 2) == 6);
        REQUIRE(P[0] == 0);
        REQUIRE(P[1] == 1);
    }

    SUBCASE("LU factorization with pivoting")
    {
        constexpr auto A = LinAlg::Matrix<4>()
            .set_entry(0, 0, 1)
            .set_entry(0, 1, 2)
            .set_entry(0, 2, 1)
            .set_entry(0, 3, 0)
            .set_entry(1, 0, 0)
            .set_entry(1, 1, 0)
            .set_entry(1, 2, 3)
            .set_entry(1, 3, 1)
            .set_entry(2, 0, 5)
            .set_entry(2, 1, 0)
            .set_entry(2, 2, 2)
            .set_entry(2, 3, 3)
            .set_entry(3, 0, 1)
            .set_entry(3, 1, 1)
            .set_entry(3, 2, 1)
            .set_entry(3, 3, 1);

        constexpr auto t = LinAlg::factorize(A);
        constexpr auto LU = std::get<0>(t);
        constexpr auto P = std::get<1>(t);
        REQUIRE(LU(0, 0) == 1);
        REQUIRE(LU(0, 1) == 2);
        REQUIRE(LU(0, 2) == 1);
        REQUIRE(LU(0, 3) == 0);
        REQUIRE(LU(1, 0) == 5);
        REQUIRE(LU(1, 1) == -10);
        REQUIRE(LU(1, 2) == -3);
        REQUIRE(LU(1, 3) == 3);
        REQUIRE(LU(2, 0) == 0);
        REQUIRE(LU(2, 1) == 0);
        REQUIRE(LU(2, 2) == 3);
        REQUIRE(LU(2, 3) == 1);
        REQUIRE(LU(3, 0) == 1);
        REQUIRE(LU(3, 1) == Rationals::rational<1, 10>);
        REQUIRE(LU(3, 2) == LU(3, 1));
        REQUIRE(LU(3, 3) == Rationals::rational<3, 5>);
        REQUIRE(P[0] == 0);
        REQUIRE(P[1] == 2);
        REQUIRE(P[2] == 2);
    }

} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED

/********************************************************************************
 *******************************************************************************/

namespace LinAlg
{

template <std::size_t N>
constexpr Vector<N> apply_permutation(const std::array<std::size_t, N-1> &P, Vector<N> b)
{
    for (std::size_t m = 0; m < N-1; ++m)
    {
        if (P[m] != m)
        {
            b = b.set_entry(P[m], b(m)).set_entry(m, b(P[m]));
        }
    }
    return b;
}

template <std::size_t N>
constexpr Vector<N> backsub_forward(const Matrix<N> &LU, const Vector<N> &b)
{
    auto soln = b;
    for (std::size_t m = 1; m < N; ++m)
    {
        for (std::size_t j = 0; j < m; ++j)
        {
            soln = soln.set_entry(m, soln(m) - soln(j) * LU(m, j));
        }
    }
    return soln;
}

template <std::size_t N>
constexpr Vector<N> backsub_back(const Matrix<N> &LU, const Vector<N> &b)
{
    auto soln = b.set_entry(N-1, b(N-1) / LU(N-1, N-1));
    if constexpr (N == 1)
    {
        return soln;
    }
    else
    {
        std::size_t m = N-2;
        while (true)
        {
            for (std::size_t j = m+1; j < N; ++j)
            {
                soln = soln.set_entry(m, soln(m) - soln(j) * LU(m, j));
            }
            soln = soln.set_entry(m, soln(m) / LU(m, m));
            if (m == 0)
            {
                return soln;
            }
            m -= 1;
        }
    }
}

template <std::size_t N>
constexpr Vector<N> backsub(const Matrix<N> &LU, const Vector<N> &b)
{
    return backsub_back(LU, backsub_forward(LU, b));
}

template <std::size_t N>
constexpr Vector<N> linear_solve(const Matrix<N> &A, const Vector<N> &b)
{
    const auto factorization = factorize(A);
    const auto LU = std::get<0>(factorization);
    const auto P  = std::get<1>(factorization);
    const Vector<N> rhs = apply_permutation(P, b);
    return backsub(LU, rhs);
}

} // namespace LinAlg

/********************************************************************************
 * Test a linear equation solve using LU factorization.
 *******************************************************************************/

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[LinAlg] Test full linear solve")
{
    constexpr auto A = LinAlg::Matrix<4>()
        .set_entry(0, 0, 1)
        .set_entry(0, 1, 2)
        .set_entry(0, 2, 1)
        .set_entry(0, 3, 0)
        .set_entry(1, 0, 0)
        .set_entry(1, 1, 0)
        .set_entry(1, 2, 3)
        .set_entry(1, 3, 1)
        .set_entry(2, 0, 5)
        .set_entry(2, 1, 0)
        .set_entry(2, 2, 2)
        .set_entry(2, 3, 3)
        .set_entry(3, 0, 1)
        .set_entry(3, 1, 1)
        .set_entry(3, 2, 1)
        .set_entry(3, 3, 1);
    constexpr auto b = LinAlg::Vector<4>()
        .set_entry(0, 1)
        .set_entry(1, 2)
        .set_entry(2, 3)
        .set_entry(3, 4);

    constexpr auto x = LinAlg::linear_solve(A, b);

    REQUIRE(x(0) == -2);
    REQUIRE(x(1) == 2);
    REQUIRE(x(2) == -1);
    REQUIRE(x(3) == 5);
} // TEST_CASE

#endif // DOCTEST_LIBRARY_INCLUDED

/********************************************************************************
 *******************************************************************************/

} // namespace Galerkin

#endif /* GALERKIN_LINALG_HPP */