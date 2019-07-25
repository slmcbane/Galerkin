Galerkin
========
### Functional Function Spaces

Galerkin is a toolbox inspired by my experiences rolling my own finite element
code for miscellaneous problems when, for whatever reason, an off-the-shelf code
isn't satisfactory. I wanted a set of tools that made full use of modern C++
features and didn't force an object-oriented programming model. The tools should
be composable and not force the user into data structures they didn't choose.
Efficiency of the compiled code is a primary goal; however, this only addresses
performance on a single CPU. Parallelism is better addressed elsewhere, and
thanks to composability of the interface, you can flexibly hook your FEM code to
a distributed linear algebra backend, or implement your own.

The style of the code in Galerkin is heavily functional and influenced by the
[Julia language](https://julialang.org). Julia is probably my favorite language
to work in in any problem domain; however, I found that when I wrote my FEM code
in Julia the way I wanted to write it, I couldn't quite get the performance that
it's capable of when hand-tuned. In addition, access to state-of-the-art
scientific libraries written in C and C++ is just simpler when working in the
host language (or its cousin). Galerkin is intended to be performant so that when
used in research code, I won't have to think twice to know the bottleneck isn't
in the low-level details of my finite element code.

### License
Copyright (c) 2019, Sean McBane.

Galerkin is free software, licensed under the terms of the MIT license. Please
see the COPYRIGHT file for the formal terms and disclaimers.