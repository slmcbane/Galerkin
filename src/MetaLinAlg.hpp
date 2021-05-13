/*
 * Copyright (c) 2019, Sean McBane and The University of Texas at Austin.
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef METALINALG_HPP
#define METALINALG_HPP

/*!
 * @file MetaLinAlg.hpp
 * @brief A compile-time matrix facility intended to solve a linear system to
 * find coefficients of multinomial basis functions.
 */

// #include "Rationals.hpp"
// #include "utils.hpp"

#include <array>
#include <tuple>

namespace Galerkin
{

namespace LinAlg
{

template <std::size_t N>
struct Matrix
{
    std::array<double, N * N> data;

    constexpr double operator()(int m, int n) const { return data[m * N + n]; }

    template <class T>
    constexpr Matrix set_entry(int m, int n, T val) const
    {
        std::array<double, N * N> xs{data};
        xs[m * N + n] = val;
        return Matrix(xs);
    }

    constexpr Matrix(const std::array<double, N * N> &xs) : data{xs} {}

    constexpr Matrix(const Matrix &other) : data{other.data} {}

    constexpr Matrix() : data{0} {}
};

template <std::size_t N>
struct Vector
{
    std::array<double, N> data;

    constexpr Vector(const std::array<double, N> &xs) : data{xs} {}

    constexpr Vector(const Vector &other) : data{other.data} {}

    constexpr Vector() : data{0} {}

    constexpr double operator()(int m) const { return data[m]; }

    template <class T>
    constexpr Vector set_entry(int m, T val) const
    {
        std::array<double, N> xs{data};
        xs[m] = val;
        return Vector(xs);
    }
};

template <auto I, auto N>
constexpr auto canonical() noexcept
{
    static_assert(N >= 0 && I < N);
    if constexpr (N == 0)
    {
        return Vector<0>();
    }
    else
    {
        return Vector<N>().set_entry(I, 1);
    }
}

} // namespace LinAlg

/********************************************************************************
 * Test that canonical basis vector is correct.
 *******************************************************************************/
#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("[LinAlg] Canonical basis vectors")
{
    constexpr auto c0 = LinAlg::canonical<0, 3>();
    REQUIRE(c0(0) == 1.0);
    REQUIRE(c0(1) == 0.0);
    REQUIRE(c0(2) == 0.0);

    constexpr auto c1 = LinAlg::canonical<1, 3>();
    REQUIRE(c1(0) == 0.0);
    REQUIRE(c1(1) == 1.0);
    REQUIRE(c1(2) == 0.0);

    constexpr auto c2 = LinAlg::canonical<2, 3>();
    REQUIRE(c2(0) == 0.0);
    REQUIRE(c2(1) == 0.0);
    REQUIRE(c2(2) == 1.0);
}

#endif

/********************************************************************************
 *******************************************************************************/

namespace LinAlg
{

template <std::size_t N>
std::tuple<Matrix<N>, std::array<std::size_t, N-1>> factorize(Matrix<N> A)
{
    std::array<std::size_t, N-1> swaps{0};
    for (std::size_t m = 0; m < N-1; ++m)
    {
        auto swap = find_pivot(A, m);
        swaps[m] = swap;
        A = do_row_elimination(A, m, m+1);
    }
    return std::tuple(A, swaps);
}

} // namespace LinAlg

} // namespace Galerkin

#endif /* METALINALG_HPP */