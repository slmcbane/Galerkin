Galerkin
========
## Functional Function Spaces

Galerkin is a finite element toolbox tailored to use in my own research code.
I often find existing, fully-featured finite element toolkits to be cumbersome to
work with precisely because they have **so many** features; just extracting the
few relevant bits of functionality can be very time consuming.

I also saw the potential to use the compiler to generate code for shape functions,
Jacobians of coordinate maps, and other mathematical tasks that are tedious to
derive and hard code. Galerkin is also a "because I can" proof of concept
demonstrating the use of metaprogramming to derive quadrature rules and shape
functions and define classes of finite elements at compile time, in a way that
results in high performance specialized code.

The library thus aims to:

- Provide FEM primitives that can be flexibly integrated into high-level code
- Generate highly efficient machine code using a high level of abstraction
- Provide an interface to define finite elements that removes as much work as possible from the implementer and moves it to the compiler.

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
- Compile time generation of polynomial shape functions with compile-time
  partial derivatives. Completed; implementation allows for deriving shape
  functions using arbitrary functionals to specify the degrees of freedom. See
  example.
- An abstract interface for reference -> instantiated coordinate maps and
  Jacobians; also compute integrals in reference coordinates. Completed, a few
  useful particular examples implemented (uniform scaling + translation, general
  quadrilateral, general triangle).
- Reference finite elements may be defined using a DSL. Completed - to create a
  "finite element", use `derive_shape_functions` to find the shape functions for
  your element at compile time given a form for the shape functions and a list
  of degrees of freedom, and define the coordinate map mapping reference to
  instantiated domain. See example code for a triangle element implementation.

In order to implement a fully functional FEM, some extra tools will be needed
beyond what is laid out above. Primarily, mesh data structures and linear algebra
tools are outside of Galerkin's scope. I will likely publish additional
repositories demonstrating a FEM implemented on top of Galerkin, as I use them in
research.

## Prerequisites
Galerkin has no external dependencies. It does make extensive use of C++17
features and so requires a compatible compiler. In addition, at time of this
writing (2/2020), GCC fails to compile the test suite with an internal compiler
error - some versions may have a cryptic error message instead, but the cause
is a compiler bug. Galerkin has been tested successfully with clang 6.0.0 and
clang 8.0.0.

### License
Copyright (c) 2020, Sean McBane and The University of Texas at Austin.

Galerkin is free software, licensed under the terms of the MIT license. Please
see the COPYRIGHT file for the formal terms and disclaimers.
