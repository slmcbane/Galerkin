Galerkin
========
### Functional Function Spaces

Galerkin is a toolbox inspired by my experiences rolling my own finite element
code for miscellaneous problems when, for whatever reason, an off-the-shelf code
isn't satisfactory. I wanted a set of tools that made full use of modern C++
features and didn't force an object-oriented programming model. The tools should
be composable and not force the user into data structures they didn't choose -
you should easily be able to use only elements and integrators with your own
mesh data structure, etc.

Efficiency of the compiled code is a primary goal; however, this only addresses
performance on a single CPU. Parallelism is better addressed elsewhere, and
thanks to composability of the interface, you can flexibly hook your FEM code to
a distributed linear algebra backend, or implement your own.

The style of the code in Galerkin is heavily functional and influenced by the
[Julia language](https://julialang.org). Julia is probably my favorite language
to work in in any problem domain; however, I found that my Julia FEM code was not
as efficient as I knew it was capable of. Here, I hope I've leveraged to the
fullest the capabilities of C++ for code generation of specialized routines for
every kind of element.

### Features
This section is also a roadmap for work that has yet to be completed.

- Compile time generation of quadrature routines for polynomials of arbitrary
order. Almost completed.
- Compile time generation of multinomial shape functions with compile-time
partial derivatives. Requires compile time linear solve. Next priority.
- Compile time generation of a function evaluating the Jacobian of the reference
to instantiated element map, and the action of its inverse.
- Arbitrary multinomial shape functions may be combined with arithmetic
operators: `f1 * f2` returns a new multinomial object evaluating
`f1(xs...) * f2(xs...)`, and so forth. I do not anticipate adding multinomial
division.

Using these features many useful element types can be defined using a domain
specific language implemented entirely with templates. It will punish your
compiler, but the resulting code should be efficient. I will be examining the
generated assembly for common operations like computing an element stiffness
matrix and comparing to a more traditional approach using hardcoded functions for
each element type.

### License
Copyright (c) 2019, Sean McBane.

Galerkin is free software, licensed under the terms of the MIT license. Please
see the COPYRIGHT file for the formal terms and disclaimers.