#include "../Galerkin.hpp"
#include <iostream>

using namespace Galerkin;
using namespace Elements;
using namespace Rationals;
using namespace Polynomials;

#include <Eigen/Cholesky>
#include <Eigen/Core>

constexpr std::array<double, 5> scalings = {0.5, 1, 0.5, 0.25, 0.25};

constexpr std::array<double, 5> offsets = { 0.5, 2, 3.5, 4.25, 4.75 };

constexpr std::array<C1IntervalElement<double>, 5> elements = {
    C1IntervalElement(scalings[0], offsets[0]), C1IntervalElement(scalings[1], offsets[1]),
    C1IntervalElement(scalings[2], offsets[2]), C1IntervalElement(scalings[3], offsets[3]),
    C1IntervalElement(scalings[4], offsets[4])};

constexpr std::array<std::array<int, 4>, 5> dof_numbers = {
    std::array{0, 1, 2, 3}, std::array{3, 2, 4, 5}, std::array{5, 4, 6, 7}, std::array{7, 6, 8, 9},
    std::array{9, 8, 10, 11}};

constexpr std::size_t num_dofs = 12;

constexpr auto el_stiffness(const C1IntervalElement<double> &el)
{
    auto form = [&](auto f, auto g) { return el.template partial<0>(f) * el.template partial<0>(g); };

    return el.form_matrix(form);
}

Eigen::MatrixXd stiffness_matrix()
{
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(num_dofs, num_dofs);
    for (int i = 0; i < 5; ++i)
    {
        auto Ke = el_stiffness(elements[i]);
        for (int j = 0; j < 4; ++j)
        {
            for (int k = 0; k < 4; ++k)
            {
                K(dof_numbers[i][j], dof_numbers[i][k]) += Ke(j, k);
            }
        }
    }
    return K;
}

static_assert(elements[0].coordinate_map()(std::tuple(-1))[0] == 0);
static_assert(elements[0].coordinate_map()(std::tuple(1))[0] == 1);

constexpr auto u = make_poly(
    std::tuple(rational<3, 8>, rational<-47, 12>, rational<97, 8>, -rational<115, 12>),
    PowersList<Powers<4>, Powers<3>, Powers<2>, Powers<1>>{});

constexpr auto f = u.template partial<0>().template partial<0>();

Eigen::VectorXd forcing_vector()
{
    Eigen::VectorXd F = Eigen::VectorXd::Zero(num_dofs);
    for (int i = 0; i < 5; ++i)
    {
        const auto &el = elements[i];
        for (int j = 0; j < 4; ++j)
        {
            auto rule = el.coordinate_map().template quadrature_rule<6>();
            double detj = el.coordinate_map().detJ()();
            auto g = [&](auto xi) {
                return detj * f(el.coordinate_map()(std::tuple(xi))) * el.basis()[j](xi);
            };
            F[dof_numbers[i][j]] -= Quadrature::integrate(g, rule);
        }
    }

    return F;
}

void print_polynomial(const Polynomial<double, Powers<0>, Powers<1>, Powers<2>, Powers<3>> &p)
{
    printf("%f * x^3 + %f * x^2 + %f * x + %f", p.coeffs()[3], p.coeffs()[2], p.coeffs()[1], p.coeffs()[0]);
}

void eliminate_dirichlet_bc(Eigen::MatrixXd &K, Eigen::VectorXd &F, int which)
{
    for (int i = 0; i < K.rows(); ++i)
    {
        K(i, which) = 0;
        K(which, i) = 0;
    }
    K(which, which) = 1.0;
    F(which) = 0.0;
}

int main()
{
    auto K = stiffness_matrix();
    auto F = forcing_vector();
    eliminate_dirichlet_bc(K, F, 1);
    eliminate_dirichlet_bc(K, F, num_dofs - 2);
    std::cout << "K = " << K << '\n';
    std::cout << "F = " << F << '\n';
    Eigen::VectorXd U = K.ldlt().solve(F);

    for (int i = 0; i < 5; ++i)
    {
        const auto &el = elements[i];
        printf("Solution on element %d: ", i);
        auto p = 0.0 * el.basis()[0];
        for (int j = 0; j < 4; ++j)
        {
            p += el.basis()[j] * U[dof_numbers[i][j]];
        }
        print_polynomial(p);
        printf("\n");
    }

    return 0;
}