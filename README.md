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

### Generating shape functions

Generate the shape functions for a 1-dimensional element with quadratic shape
functions:

    constexpr auto form = make_form(
        powers(intgr_constant<2>),
        powers(intgr_constant<1>),
        powers(intgr_constant<0>));

    constexpr auto constraints = make_list(
        evaluate_at(rational<-1>),
        evaluate_at(rational<0>),
        evaluate_at(rational<1>));
        
    constexpr auto fns = derive_shape_functions(form, constraints);

    REQUIRE(get<0>(fns) ==
        multinomial(
            term(rational<1, 2>, powers(intgr_constant<2>)),
            term(-rational<1, 2>, powers(intgr_constant<1>))
        ));

    REQUIRE(get<1>(fns) ==
        multinomial(
            term(-rational<1>, powers(intgr_constant<2>)),
            term(rational<1>, powers(intgr_constant<0>))
        ));

    REQUIRE(get<2>(fns) ==
        multinomial(
            term(rational<1, 2>, powers(intgr_constant<2>)),
            term(rational<1, 2>, powers(intgr_constant<1>))
        ));

Generate shape functions for a 1-dimensional element with cubic shape
functions with the derivatives at the end nodes as DOF's; this element can be
used to construct a FEM that has solutions in $ C^1 $.

    constexpr auto form = make_form(
        powers(intgr_constant<3>),
        powers(intgr_constant<2>),
        powers(intgr_constant<1>),
        powers(intgr_constant<0>));

    constexpr auto constraints = make_list(
        evaluate_at(rational<-1>),
        evaluate_at(rational<1>),
        partial_at<0>(rational<-1>),
        partial_at<0>(rational<1>));

    constexpr auto fns = derive_shape_functions(form, constraints);

    REQUIRE(get<0>(fns) ==
        multinomial(
            term(rational<1, 4>, powers(intgr_constant<3>)),
            term(-rational<3, 4>, powers(intgr_constant<1>)),
            term(rational<1, 2>, powers(intgr_constant<0>))
        ));

    REQUIRE(get<1>(fns) ==
        multinomial(
            term(-rational<1, 4>, powers(intgr_constant<3>)),
            term(rational<3, 4>, powers(intgr_constant<1>)),
            term(rational<1, 2>, powers(intgr_constant<0>))
        ));

    REQUIRE(get<2>(fns) ==
        multinomial(
            term(rational<1, 4>, powers(intgr_constant<3>)),
            term(-rational<1, 4>, powers(intgr_constant<2>)),
            term(-rational<1, 4>, powers(intgr_constant<1>)),
            term(rational<1, 4>, powers(intgr_constant<0>))
        ));

    REQUIRE(get<3>(fns) ==
        multinomial(
            term(rational<1, 4>, powers(intgr_constant<3>)),
            term(rational<1, 4>, powers(intgr_constant<2>)),
            term(-rational<1, 4>, powers(intgr_constant<1>)),
            term(-rational<1, 4>, powers(intgr_constant<0>))
        ));

Generate bilinear shape functions for a quadrilateral:

    constexpr auto form = make_form(
        powers(intgr_constant<1>, intgr_constant<1>),
        powers(intgr_constant<1>, intgr_constant<0>),
        powers(intgr_constant<0>, intgr_constant<1>),
        powers(intgr_constant<0>, intgr_constant<0>)
    );

    constexpr auto constraints = make_list(
        evaluate_at(rational<-1>, rational<-1>),
        evaluate_at(rational<-1>, rational<1>),
        evaluate_at(rational<1>, rational<1>),
        evaluate_at(rational<1>, rational<-1>)
    );

    constexpr auto fns = derive_shape_functions(form, constraints);

    REQUIRE(get<0>(fns) ==
        multinomial(
            term(rational<1, 4>, powers(intgr_constant<1>, intgr_constant<1>)),
            term(-rational<1, 4>, powers(intgr_constant<1>, intgr_constant<0>)),
            term(-rational<1, 4>, powers(intgr_constant<0>, intgr_constant<1>)),
            term(rational<1, 4>, powers(intgr_constant<0>, intgr_constant<0>))
        ));

    REQUIRE(get<1>(fns) ==
        multinomial(
            term(-rational<1, 4>, powers(intgr_constant<1>, intgr_constant<1>)),
            term(-rational<1, 4>, powers(intgr_constant<1>, intgr_constant<0>)),
            term(rational<1, 4>, powers(intgr_constant<0>, intgr_constant<1>)),
            term(rational<1, 4>, powers(intgr_constant<0>, intgr_constant<0>))
        ));

    REQUIRE(get<2>(fns) ==
        multinomial(
            term(rational<1, 4>, powers(intgr_constant<1>, intgr_constant<1>)),
            term(rational<1, 4>, powers(intgr_constant<1>, intgr_constant<0>)),
            term(rational<1, 4>, powers(intgr_constant<0>, intgr_constant<1>)),
            term(rational<1, 4>, powers(intgr_constant<0>, intgr_constant<0>))
        ));

    REQUIRE(get<3>(fns) ==
        multinomial(
            term(-rational<1, 4>, powers(intgr_constant<1>, intgr_constant<1>)),
            term(rational<1, 4>, powers(intgr_constant<1>, intgr_constant<0>)),
            term(-rational<1, 4>, powers(intgr_constant<0>, intgr_constant<1>)),
            term(rational<1, 4>, powers(intgr_constant<0>, intgr_constant<0>))
        ));

### License
Copyright (c) 2019, Sean McBane and The University of Texas at Austin.

Galerkin is free software, licensed under the terms of the MIT license. Please
see the COPYRIGHT file for the formal terms and disclaimers.
