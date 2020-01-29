/*
 * Copyright (c) 2019, The University of Texas at Austin & Sean McBane
 * This file is part of the Galerkin library; Galerkin is copyright under the
 * terms of the MIT license. Please see the top-level COPYRIGHT file for details.
 */

#ifndef ELEMENTBASE_HPP
#define ELEMENTBASE_HPP

/*!
 * @file ElementBase.hpp
 * @brief Defines the CRTP base class functionality for finite elements.
 */

#include "utils.hpp"

namespace Galerkin
{

/*!
 * @brief A tag type to specify the order of integration in ElementBase routines.
 */
template <int I>
struct IntegrationOrder
{
    constexpr static int order = I;
};

template <class T>
constexpr bool is_integration_order = false;

template <int I>
constexpr bool is_integration_order<IntegrationOrder<I>> = true;

/*!
 * @brief Specialize `DefaultIntegrationOrder` for derived element types to specify.
 * 
 * For a type `Derived` inheriting from `ElementBase<Derived>`, a specialization
 * of `DefaultIntegrationOrder` with a constexpr static data member `order` is
 * required to tell quadrature routines what order to use.
 */
template <class T>
struct DefaultIntegrationOrder
{};

/// Contains functionality relating to the abstract concept of a finite element.
namespace Elements
{

/// Specialize for your `Form` type if it is symmetric in its arguments.
template <class T>
struct IsSymmetric { static constexpr bool value = false; };

template <class T>
constexpr bool is_symmetric = IsSymmetric<T>::value;

/*!
 * @brief CRTP base class for elements.
 */
template <class Derived>
struct ElementBase
{
    /*!
     * @brief Take the partial derivative of `f` through the underlying coordinate transformation for the element.
     * @see `Transforms::TransformBase::partial`
     */
    template <auto I, class F>
    constexpr auto partial(const F &f) const
    {
        return derived().transform().template partial<I>(f);
    }

    /*!
     * @brief Integrate `f` in the coordinates given by the underlying coordinate transformation.
     * 
     * The optional second argument is an `IntegrationOrder` specifying the polynomial
     * degree up to which the quadrature scheme used should be accurate.
     * 
     * @see `Transforms::TransformBase::integrate`
     */
    template <class F, int Order = DefaultIntegrationOrder<Derived>::order>
    constexpr auto integrate(const F &f, 
                             IntegrationOrder<Order> order = IntegrationOrder<Order>{}) const noexcept
    {
        return derived().transform().template integrate<order.order>(f);
    }

    /*!
     * @brief Construct a matrix from a form applied to pairs of basis functions.
     * 
     * This is a utility in the sense that a user could create this functionality
     * for themselves using only `integrate` from above, but the metaprogramming
     * involved makes it a nuisance to do so. This function eases the creation of
     * common matrices like the `mass` and `stiffness` matrices. It accepts an
     * argument `form` which maps from a pair of scalar functions to a new scalar
     * function, and returns a "matrix" `K` with entries
     * `K(i, j) = integrate(form(basis[i], basis[j]), order)`.
     * 
     * The returned object is actually a closure that returns the appropriate
     * matrix element when called with indices. If `is_symmetric<Form>` is true,
     * the closure only uses storage for the upper triangular part of the matrix
     * and computation is only carried out for the upper triangular elements.
     * Elements are stored in row major order.
     */
    template <class Form, int Order = DefaultIntegrationOrder<Derived>::order>
    constexpr auto form_matrix(const Form &form,
                               IntegrationOrder<Order> order = IntegrationOrder<Order>{}) const noexcept
    {
        using eltype = decltype(integrate(
            form(get<0>(Derived::basis), get<0>(Derived::basis)), order));
        if constexpr (is_symmetric<Form>)
        {
            static_assert(false, "Symmetric specialization not implemented yet");
        }
        else
        {
            std::array<eltype, Derived::basis_size * Derived::basis_size> mat_storage {};
            static_for<0, Derived::basis_size, 1>(
                [&](auto i)
                {
                    // Written this way triggers an internal compiler error, instead of
                    // an (incorrect) error about constexpr
                    constexpr auto index1 = i();
                    static_for<0, Derived::basis_size, 1>(
                        [&](auto j)
                        {
                            mat_storage[index1*Derived::basis_size + j()] =
                                integrate(
                                    form(get<index1>(Derived::basis), get<j()>(Derived::basis)),
                                    order
                                );
                        }
                    );
                }
            );

            return [=](auto i, auto j)
            {
                return mat_storage[i * Derived::basis_size + j];
            };
        }
    }
private:
    ElementBase() = default;
    friend Derived;

    constexpr const auto &derived() const
    {
        return static_cast<const Derived&>(*this);
    }
};

} // namespace Elements

} // namespace Galerkin

#endif /* ELEMENTBASE_HPP */