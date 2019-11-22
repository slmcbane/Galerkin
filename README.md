Galerkin
========
## Functional Function Spaces

Galerkin is a finite element toolbox tailored to re-use in research code that
may need to customize data structures and algorithms. In my experience,
fully-featured FEM toolboxes are very powerful __but__ suffer from one flaw.
That flaw is that, when the fundamental components of the finite element
method - elements, shape functions, quadrature, assembly - need to be re-used
in a context the library developers did not have in mind, the re-use of
library code may become extremely difficult. This may be an intrinsic fault -
a failure of software design - or just because the code is complex and poorly
documented. Galerkin exists to provide a toolbox for working with primitives
of the finite element method, and aims to stay out of the way of high-level
code.

I also believe that performance of existing finite element libraries can be
hampered by excessive use of OO programming, and the potential of modern
C++ for code generation is overlooked. The derivation of shape functions and
element mappings for common cases is mechanical and there is no need to
expend the effort to hard code particular cases; metaprogramming and
compile-time constructs can be used to specify reference finite elements
in an easily-readable domain specific language. Galerkin is a proof of
concept of this idea.

The library thus aims to:

- Provide FEM primitives that can be flexibly integrated into high-level code
- Generate highly efficient machine code; performance of templated code should
be on par with specialized hand-coded implementations.
- Eliminate the need to do mathematical derivations and hard-code common
constructs by hand; use modern C++ technique to make the compiler do it.

## Features
In general I will be adding features to Galerkin as I use them. This section is
also a roadmap for work that has yet to be completed.

As an initial goal, I will aim to solve problems in 1, 2, and 3 dimensions with
elements of Lagrange type (that is, the degrees of freedom are the values of
the shape functions at nodes of an element). These elements are appropriate for
solutions in $ C^0 $. To get to the point of solving a real problem, the
following will be implemented:

- Compile time generation of quadrature routines for polynomials of arbitrary
  order. Completed - up to polynomials of order 15. A modification of stopping
  tolerance for root finding of Legendre polynomials will allow to go higher.
- Compile time generation of multinomial shape functions with compile-time
  partial derivatives. Completed; implementation allows for deriving shape
  functions using arbitrary functions to specify the degrees of freedom. See
  example.
- Compile time generation of routines to compute reference -> physical mapping
  according to a specified form, the Jacobian and its determinant, etc.
- Reference finite elements may be defined using a DSL to specify the degrees
  of freedom, the vertices of the reference element and the map from reference
  to instantiated domains, and quadrature rules.

In order to implement a fully functional FEM, some extra tools will be needed
beyond what is laid out above. Primarily, I will implement simple mesh data
structures in the library as well, which may be used with or without the other
tools provided. Mesh tools will include I/O for VTK and perhaps other formats
for data exchange, or mesh generation. I have not decided what direction to
take as far as integration of any linear algebra tools in Galerkin.

## Examples
The below code is an example of using the primitives currently implemented. The
actual functionality shown will itself be abstracted into the library soon.

This example code implements a function computing the stiffness matrix
(discretized weak form of the Laplacian) on a quadrilateral element that is known
to be uniformly scaled from the reference square.

    #include "src/UniformScaling.hpp"
    #include "src/Elements.hpp"

    using Galerkin::Elements::make_form;
    using Galerkin::Metanomials::powers;

    template <auto i>
    constexpr auto intgr_constant = std::integral_constant<decltype(i), i>();

    /*
    * Only polynomial basis functions are supported right now. This constructs an
    * object representing a general polynomial of the form
    * p(x, y) = axy + bx + cy + d
    */
    constexpr auto form = make_form(powers(intgr_constant<1>, intgr_constant<1>),
                                    powers(intgr_constant<1>, intgr_constant<0>),
                                    powers(intgr_constant<0>, intgr_constant<1>),
                                    powers(intgr_constant<0>, intgr_constant<0>));

    using Galerkin::make_list;
    using Galerkin::Elements::evaluate_at;
    using Galerkin::Rationals::rational;

    /*
    * These are the degrees of freedom defining the shape function space. This
    * example derives the ordinary bilinear shape functions for a first-order
    * quadrilateral element; degrees of freedom are the value of the function at
    * each of the four nodes of the reference element.
    */
    constexpr auto DOF = make_list(evaluate_at(rational<-1>, rational<-1>),
                                evaluate_at(rational<-1>, rational<1>),
                                evaluate_at(rational<1>, rational<1>),
                                evaluate_at(rational<1>, rational<-1>));

    using Galerkin::Elements::derive_shape_functions;

    /*
    * Given the degrees of freedom and the form of shape functions to derive,
    * derive_shape_functions returns a list containing the canonical basis. All of
    * this computation happens at compile-time; shape functions are encoded in
    * types and contain no runtime data.
    */
    constexpr auto basis_fns = derive_shape_functions(form, DOF);

    // Make this function have C linkage for inspecting the assembly.
    extern "C" {
        void stiffness_matrix(double, double*);
    }

    using Galerkin::Transforms::UniformScaling;
    using Galerkin::get;

    /*
    * Compute the mass matrix for a quadrilateral element that is uniformly scaled
    * and translated from the reference element. Stores the upper triangular
    * elements in row order contiguously in the given array.
    */
    void stiffness_matrix(double scale, double* dst)
    {
        auto transform = UniformScaling(scale, std::array<double, 2>{0});
        
        // The index when taking a partial derivative is a template parameter,
        // so you cannot write an ordinary for loop. Galerkin::static_reduce is a
        // utility to do loop reductions with compile-time loop index, but I will
        // just write the code out here to avoid explaining it.
        
        // This defines a lambda that, given functions f1 and f2, computes the
        // integral over the transformed element of dot(grad(f1), grad(f2)).
        auto integral_form = [&](const auto &f1, const auto &f2)
        {
            // Taking partials returns another object that is itself a function.
            auto partial10 = transform.partial<0>(f1);
            auto partial11 = transform.partial<1>(f1);
            auto partial20 = transform.partial<0>(f2);
            auto partial21 = transform.partial<1>(f2);

            // Objects supporting the function protocol in Galerkin
            // (and inheriting from FunctionBase) can be multiplied, added,
            // subtracted, and divided to build expressions.
            return transform.integrate<2>(partial10 * partial20 + partial11 * partial21);
        };

        dst[0] = integral_form(get<0>(basis_fns), get<0>(basis_fns));
        dst[1] = integral_form(get<0>(basis_fns), get<1>(basis_fns));
        dst[2] = integral_form(get<0>(basis_fns), get<2>(basis_fns));
        dst[3] = integral_form(get<0>(basis_fns), get<3>(basis_fns));
        dst[4] = integral_form(get<1>(basis_fns), get<1>(basis_fns));
        dst[5] = integral_form(get<1>(basis_fns), get<2>(basis_fns));
        dst[6] = integral_form(get<1>(basis_fns), get<3>(basis_fns));
        dst[7] = integral_form(get<2>(basis_fns), get<2>(basis_fns));
        dst[8] = integral_form(get<2>(basis_fns), get<3>(basis_fns));
        dst[9] = integral_form(get<3>(basis_fns), get<3>(basis_fns));
    }

The above function compiles down to brief, vectorized floating point code as
shown below; this assembly was generated by clang-6.0.0 on Ubuntu for a Skylake
processor.

    stiffness_matrix:                       # @stiffness_matrix
        .cfi_startproc
    # %bb.0:
        vmovsd	.LCPI0_0(%rip), %xmm1   # xmm1 = mem[0],zero
        vdivsd	%xmm0, %xmm1, %xmm1
        vaddsd	%xmm1, %xmm1, %xmm2
        vmulsd	%xmm0, %xmm0, %xmm0
        vmulsd	%xmm0, %xmm2, %xmm3
        vmulsd	%xmm1, %xmm1, %xmm4
        vunpcklpd	%xmm0, %xmm3, %xmm0 # xmm0 = xmm3[0],xmm0[0]
        vmulpd	.LCPI0_1(%rip), %xmm0, %xmm0
        vunpcklpd	%xmm4, %xmm1, %xmm3 # xmm3 = xmm1[0],xmm4[0]
        vmulpd	%xmm3, %xmm0, %xmm3
        vmovupd	%xmm3, (%rdi)
        vpermilpd	$1, %xmm0, %xmm0 # xmm0 = xmm0[1,0]
        vmulsd	%xmm1, %xmm2, %xmm1
        vmulsd	%xmm1, %xmm0, %xmm0
        vmovsd	%xmm0, 16(%rdi)
        vmovhpd	%xmm3, 24(%rdi)
        vmovupd	%xmm3, 32(%rdi)
        vmovsd	%xmm0, 48(%rdi)
        vmovupd	%xmm3, 56(%rdi)
        vmovlpd	%xmm3, 72(%rdi)
        retq

### License
Copyright (c) 2019, Sean McBane and The University of Texas at Austin.

Galerkin is free software, licensed under the terms of the MIT license. Please
see the COPYRIGHT file for the formal terms and disclaimers.
