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

## Examples
The code below defines a tetrahedral element with linear shape functions using
Galerkin's interface. It demonstrates (C compatible) functions to compute the
upper triangular part of the mass and stiffness matrices.

    #include "src/ElementBase.hpp"
    #include "src/Elements.hpp"
    #include "src/TriangleTransform.hpp"

    using Galerkin::Metanomials::Powers;
    using Galerkin::Rationals::rational;
    using Galerkin::Elements::evaluate_at;

    class TetElement;

    // This template specialization must come *before* the definition of the
    // TetElement class because the instantiation of ElementBase depends on the
    // specialization of DefaultIntegrationOrder.
    template <>
    class Galerkin::DefaultIntegrationOrder<TetElement> : public Galerkin::IntegrationOrder<2>
    {};

    class TetElement : public Galerkin::Elements::ElementBase<TetElement>
    {
    public:
        template <class... Points>
        constexpr TetElement(const Points &...points) noexcept : m_trans(points...)
        {}

        static constexpr auto basis = Galerkin::Elements::derive_shape_functions(
            Galerkin::Elements::make_form(
                Powers<1, 0>{}, Powers<0, 1>{}, Powers<0, 0>{}
            ),
            Galerkin::make_list(
                evaluate_at(rational<-1>, rational<-1>),
                evaluate_at(rational<-1>, rational<1>),
                evaluate_at(rational<1>, rational<-1>)
            )
        );

        constexpr auto &coordinate_map() const noexcept { return m_trans; }

    private:
        Galerkin::Transforms::TriangleTransform<double> m_trans;
    };

    extern "C"
    {
        void mass_matrix(const double *verts, double *dst);
        void stiffness_matrix(const double *verts, double *dst);
    } // extern "C"

    void mass_matrix(const double *verts, double *dst)
    {
        const auto p1 = std::tuple(verts[0], verts[1]);
        const auto p2 = std::tuple(verts[2], verts[3]);
        const auto p3 = std::tuple(verts[4], verts[5]);
        const TetElement element(p1, p2, p3);

        constexpr auto form = [](auto f, auto g) { return f * g; };

        const auto matrix = element.form_matrix(form);
        int index = 0;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = i; j < 3; ++j)
            {
                dst[index++] = matrix(i, j);
            }
        }
    }

    void stiffness_matrix(const double *verts, double *dst)
    {
        const auto p1 = std::tuple(verts[0], verts[1]);
        const auto p2 = std::tuple(verts[2], verts[3]);
        const auto p3 = std::tuple(verts[4], verts[5]);
        const TetElement element(p1, p2, p3);

        auto form = [&](auto f, auto g)
        {
            return element.template partial<0>(f) * element.template partial<0>(g) +
                element.template partial<1>(f) * element.template partial<1>(g);
        };

        // Gradient of the first order shape functions is constant, so zero-th order
        // integration is fine.
        const auto matrix = element.form_matrix(form, Galerkin::IntegrationOrder<0>{});
        int index = 0;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = i; j < 3; ++j)
            {
                dst[index++] = matrix(i, j);
            }
        }
    }

The above compiles down to brief, vectorized floating point code as
shown below; this assembly was generated by clang 8.0.0 on Ubuntu for a Skylake
processor.

		.text
		.file	"test.cpp"
		.section	.rodata.cst16,"aM",@progbits,16
		.p2align	4               # -- Begin function mass_matrix
	.LCPI0_0:
		.quad	4590669220166325588     # double 0.083333333333333315
		.quad	4586165620538955094     # double 0.041666666666666671
		.section	.rodata.cst8,"aM",@progbits,8
		.p2align	3
	.LCPI0_1:
		.quad	4590669220166325589     # double 0.083333333333333329
		.text
		.globl	mass_matrix
		.p2align	4, 0x90
		.type	mass_matrix,@function
	mass_matrix:                            # @mass_matrix
		.cfi_startproc
	# %bb.0:
		vmovsd	(%rdi), %xmm0           # xmm0 = mem[0],zero
		vmovsd	8(%rdi), %xmm1          # xmm1 = mem[0],zero
		vmovsd	24(%rdi), %xmm2         # xmm2 = mem[0],zero
		vmovsd	16(%rdi), %xmm3         # xmm3 = mem[0],zero
		vmovsd	32(%rdi), %xmm4         # xmm4 = mem[0],zero
		vsubsd	%xmm0, %xmm4, %xmm4
		vsubsd	%xmm0, %xmm3, %xmm0
		vsubsd	40(%rdi), %xmm1, %xmm3
		vmulsd	%xmm0, %xmm3, %xmm0
		vsubsd	%xmm1, %xmm2, %xmm1
		vfmadd213sd	%xmm0, %xmm4, %xmm1 # xmm1 = (xmm4 * xmm1) + xmm0
		vmovddup	%xmm1, %xmm0    # xmm0 = xmm1[0,0]
		vmulpd	.LCPI0_0(%rip), %xmm0, %xmm0
		vpermpd	$84, %ymm0, %ymm2       # ymm2 = ymm0[0,1,1,1]
		vmulsd	.LCPI0_1(%rip), %xmm1, %xmm1
		vmovupd	%ymm2, -72(%rsp)
		vmovsd	%xmm1, -40(%rsp)
		vmovhpd	%xmm0, -32(%rsp)
		vpermilpd	$3, %xmm0, %xmm0 # xmm0 = xmm0[1,1]
		vmovupd	%xmm0, -24(%rsp)
		vmovsd	%xmm1, -8(%rsp)
		movq	-56(%rsp), %rax
		movq	%rax, 16(%rsi)
		vmovups	-72(%rsp), %xmm0
		vmovups	%xmm0, (%rsi)
		vmovups	-40(%rsp), %xmm0
		vmovups	%xmm0, 24(%rsi)
		movq	-8(%rsp), %rax
		movq	%rax, 40(%rsi)
		vzeroupper
		retq
	.Lfunc_end0:
		.size	mass_matrix, .Lfunc_end0-mass_matrix
		.cfi_endproc
											# -- End function
		.section	.rodata.cst16,"aM",@progbits,16
		.p2align	4               # -- Begin function stiffness_matrix
	.LCPI1_0:
		.quad	4602678819172646912     # double 0.5
		.quad	4602678819172646912     # double 0.5
	.LCPI1_2:
		.quad	-4620693217682128896    # double -0.5
		.quad	4602678819172646912     # double 0.5
	.LCPI1_3:
		.quad	4607182418800017408     # double 1
		.quad	4607182418800017408     # double 1
		.section	.rodata.cst8,"aM",@progbits,8
		.p2align	3
	.LCPI1_1:
		.quad	4607182418800017408     # double 1
		.text
		.globl	stiffness_matrix
		.p2align	4, 0x90
		.type	stiffness_matrix,@function
	stiffness_matrix:                       # @stiffness_matrix
		.cfi_startproc
	# %bb.0:
		vmovupd	(%rdi), %xmm0
		vmovupd	16(%rdi), %xmm1
		vmovupd	32(%rdi), %xmm2
		vsubpd	%xmm0, %xmm2, %xmm2
		vsubpd	%xmm0, %xmm1, %xmm0
		vmovapd	.LCPI1_0(%rip), %xmm1   # xmm1 = [5.0E-1,5.0E-1]
		vmulpd	%xmm1, %xmm2, %xmm2
		vmulpd	%xmm1, %xmm0, %xmm0
		vpermilpd	$1, %xmm0, %xmm1 # xmm1 = xmm0[1,0]
		vpermilpd	$1, %xmm2, %xmm3 # xmm3 = xmm2[1,0]
		vmulsd	%xmm0, %xmm3, %xmm4
		vfmsub231sd	%xmm1, %xmm2, %xmm4 # xmm4 = (xmm2 * xmm1) - xmm4
		vaddsd	%xmm3, %xmm2, %xmm3
		vmovsd	.LCPI1_1(%rip), %xmm5   # xmm5 = mem[0],zero
		vdivsd	%xmm4, %xmm5, %xmm5
		vmulsd	%xmm5, %xmm3, %xmm3
		vmulsd	%xmm3, %xmm3, %xmm3
		vaddsd	%xmm1, %xmm0, %xmm1
		vmulsd	%xmm5, %xmm1, %xmm1
		vfmadd213sd	%xmm3, %xmm1, %xmm1 # xmm1 = (xmm1 * xmm1) + xmm3
		vmovddup	%xmm4, %xmm3    # xmm3 = xmm4[0,0]
		vmulpd	.LCPI1_2(%rip), %xmm3, %xmm4
		vpermpd	$80, %ymm4, %ymm5       # ymm5 = ymm4[0,0,1,1]
		vpermilpd	$1, %xmm4, %xmm4 # xmm4 = xmm4[1,0]
		vmulsd	%xmm4, %xmm1, %xmm1
		vmovapd	.LCPI1_3(%rip), %xmm6   # xmm6 = [1.0E+0,1.0E+0]
		vdivpd	%xmm3, %xmm6, %xmm3
		vmulpd	%xmm3, %xmm2, %xmm2
		vmulpd	%xmm3, %xmm0, %xmm0
		vpermilpd	$1, %xmm2, %xmm3 # xmm3 = xmm2[1,0]
		vaddsd	%xmm2, %xmm3, %xmm3
		vpermilpd	$1, %xmm0, %xmm6 # xmm6 = xmm0[1,0]
		vaddsd	%xmm0, %xmm6, %xmm6
		vbroadcastsd	%xmm3, %ymm3
		vpermpd	$20, %ymm2, %ymm7       # ymm7 = ymm2[0,1,1,0]
		vblendpd	$12, %ymm7, %ymm3, %ymm3 # ymm3 = ymm3[0,1],ymm7[2,3]
		vpermpd	$84, %ymm2, %ymm7       # ymm7 = ymm2[0,1,1,1]
		vmulpd	%ymm7, %ymm3, %ymm3
		vbroadcastsd	%xmm6, %ymm6
		vpermpd	$20, %ymm0, %ymm7       # ymm7 = ymm0[0,1,1,0]
		vblendpd	$12, %ymm7, %ymm6, %ymm6 # ymm6 = ymm6[0,1],ymm7[2,3]
		vpermpd	$84, %ymm0, %ymm7       # ymm7 = ymm0[0,1,1,1]
		vfmadd213pd	%ymm3, %ymm6, %ymm7 # ymm7 = (ymm6 * ymm7) + ymm3
		vmulpd	%ymm5, %ymm7, %ymm3
		vmulsd	%xmm2, %xmm2, %xmm2
		vfmadd231sd	%xmm0, %xmm0, %xmm2 # xmm2 = (xmm0 * xmm0) + xmm2
		vmulsd	%xmm4, %xmm2, %xmm0
		vmovsd	%xmm1, -72(%rsp)
		vmovhpd	%xmm3, -64(%rsp)
		vmovupd	%ymm3, -56(%rsp)
		vextractf128	$1, %ymm3, %xmm1
		vblendpd	$1, %xmm3, %xmm1, %xmm1 # xmm1 = xmm3[0],xmm1[1]
		vmovupd	%xmm1, -24(%rsp)
		vmovsd	%xmm0, -8(%rsp)
		movq	-56(%rsp), %rax
		movq	%rax, 16(%rsi)
		vmovups	-72(%rsp), %xmm0
		vmovups	%xmm0, (%rsi)
		vmovups	-40(%rsp), %xmm0
		vmovups	%xmm0, 24(%rsi)
		movq	-8(%rsp), %rax
		movq	%rax, 40(%rsi)
		vzeroupper
		retq
	.Lfunc_end1:
		.size	stiffness_matrix, .Lfunc_end1-stiffness_matrix
		.cfi_endproc
											# -- End function

		.ident	"clang version 8.0.0-3~ubuntu18.04.2 (tags/RELEASE_800/final)"
		.section	".note.GNU-stack","",@progbits
		.addrsig
		.addrsig_sym __gxx_personality_v0


### License
Copyright (c) 2020, Sean McBane and The University of Texas at Austin.

Galerkin is free software, licensed under the terms of the MIT license. Please
see the COPYRIGHT file for the formal terms and disclaimers.
