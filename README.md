Galerkin
========
### Functional Function Spaces

Galerkin is a toolbox inspired by my need to intrusively modify finite element
algorithms in research. Fully-featured finite element libraries like FeNICs or
MFEM are powerful, but if one needs to manipulate the underlying data structures,
you may quickly find yourself lost in undocumented source code. My research has
required intrusive modification of existing code so that it sometimes becomes
easier to just start from scratch and have full constrol of the data
structures underlying a problem. For this reason, I wanted a collection of my own
primitives for finite element methods that can composably be re-used in new code.

I also believe there are CPU cycles left on the table in existing codes; I would
prefer to avoid dynamic polymorphism more often than not, and modern features of
the C++ standard library as well as old tricks like the CRTP for static
polymorphism make this easier than ever. The derivations of shape functions and
basis functions are also very mechanical, and a prime candidate for attack via
metaprogramming.

The goal of this library is then threefold:
- Provide finite element primitives that can be flexibly integrated into
  specialized code, not the kitchen sink.
- Compile to machine code as economical as code hand-written for a particular
  task; and
- Don't hard-code quadrature rules, basis functions, etc.; use compile-time
  programming to generate them as needed.

### Features
In general I will be adding features to Galerkin as I use them. This section is
also a roadmap for work that has yet to be completed.

As an initial goal I will support Lagrange elements (appropriate for solutions in
H^1), in arbitrary dimensions. In practice, I will test on a small set of
practical problems in one, two, and three dimensions, with analytical solutions
for validation of the code. Auxiliary to these examples, the library will contain
a few simple mesh data structures appropriate for such simple problems.

To get to the point of writing code to solve a real problem, the following will
be implemented:
- Compile time generation of quadrature routines for polynomials of arbitrary
  order. Completed - up to polynomials of order 15. A modification of stopping
  tolerance for root finding of Legendre polynomials will allow to go higher.
- Compile time generation of multinomial shape functions with compile-time
  partial derivatives.
- Compile time generation of shape functions, Jacobians, and their inverses.
- Lagrange elements may be defined using a template DSL to specify nodes and the
  form of the basis and shape functions.
  
### License
Copyright (c) 2019, Sean McBane.

Galerkin is free software, licensed under the terms of the MIT license. Please
see the COPYRIGHT file for the formal terms and disclaimers.