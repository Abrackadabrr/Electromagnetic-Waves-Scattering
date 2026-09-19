# Electromagnetic-Waves-Scattering

## Purpose

This repository is the core C++ library for numerical solution of electromagnetic
surface integral equations (SIE), volume integral equations (VIE), and coupled
surface-volume integral equations (VSIE).

Treat this repository as production scientific code. Prefer mathematically clear,
testable, reusable implementations over experiment-specific shortcuts.

## Response language

Respond to the user primarily in Russian.

Mixed Russian-English technical prose is explicitly allowed and preferred when
translation could alter, narrow, or obscure meaning. Do not translate mathematical
symbols, operator names, function-space names, standard computational-
electromagnetics terminology, C++ identifiers, API names, library names, or source
notation merely to make the prose fully Russian.

## Shared Codex skills

Shared electromagnetic domain skills are expected under:

    .agents/skills/

The relevant skills include:

- `electromagnetics-notation`;
- `surface-integral-equations`;
- `volume-integral-equations`;
- `operator-discretization`;
- `uniform-grid-vie`.

Use them when the task matches their scope. The repository source tree remains
authoritative for the current implementation and API.

## Canonical electromagnetic convention

Use the repository convention

    exp(-i * omega * t)

throughout.

Therefore, in a homogeneous source-free isotropic medium,

    curl E = +i * omega * mu * H
    curl H = -i * omega * epsilon * E

and the outgoing scalar Helmholtz Green function is

    G_k(x,y) = exp(+i * k * |x-y|) / (4 * pi * |x-y|).

Do not import formulas written for the opposite harmonic convention by changing
only the sign in the exponential. Convert the complete convention bundle:

1. time factor;
2. Maxwell curl signs;
3. outgoing Green function;
4. `grad_x G` / `grad_y G` relations;
5. electric and magnetic current definitions;
6. normal orientation and trace-side definitions;
7. cross-product order;
8. field-representation prefactors;
9. jump terms;
10. resulting discretized matrix signs.

Project notation is authoritative. In particular, after complete convention
translation,

    K_project <-> k^2 * L_Gibson
    R_project <-> K_Gibson.

Never infer a mathematical sign only from a C++ type or function name.

## Repository structure

Public headers are under:

    include/EMW/

Implementations are under:

    src/

Tests are under:

    tests/

Benchmarks are under:

    benchmarks/

Research/prototype programs, when present, should remain separate from reusable
library functionality.

Before assuming a path, namespace, API, type alias, matrix layout, or ownership
rule, inspect the current working tree.

## Codebase-first rule

Before implementing a new operator, mesh abstraction, quadrature scheme, matrix
representation, fast algorithm, preconditioner, or solver:

1. inspect the existing implementation and nearby tests;
2. identify whether an equivalent abstraction already exists;
3. extend the existing design when it is technically appropriate;
4. avoid creating a parallel implementation only for convenience;
5. preserve public API and naming unless refactoring is explicitly requested.

Important existing concepts include, but are not limited to:

- `EMW::OperatorK`;
- `EMW::OperatorR`;
- `EMW::Helmholtz` helpers;
- `EMW::Operators::Volume`;
- `EMW::Mesh::VolumeMesh::CubeMesh`;
- Toeplitz / multi-level Toeplitz matrix abstractions;
- factored / low-rank structured blocks.

Local source is authoritative if it differs from skill reference notes.

## Scientific-code workflow

Keep the following layers conceptually separate:

1. continuous electromagnetic formulation;
2. definition and normalization of unknowns;
3. basis and testing spaces;
4. singular and near-singular integration;
5. discrete matrix element;
6. matrix representation and storage;
7. iterative solver and preconditioner;
8. acceleration / parallelization;
9. numerical validation.

Do not silently change several mathematical layers at once.

For every operator implementation or modification, explicitly verify:

- source variable versus observation variable;
- `grad_x G` versus `grad_y G`;
- cross-product ordering;
- normal orientation and trace side where relevant;
- basis and test normalization;
- self / near / far interaction classification;
- physical units and prefactors;
- degree-of-freedom and component ordering.

## Surface integral equations

Keep different discretization branches distinct.

Do not silently identify:

- project PWC/collocation surface discretization;
- RWG/Galerkin MoM;
- open-surface equations;
- closed-surface trace equations.

For PEC, dielectric, PMCHWT/Mueller, MFIE/CFIE, or junction formulas, derive or
translate jump terms using the active normal and time convention before coding.

## Volume integral equations

Before deriving or implementing a VIE, state the unknown explicitly: `E`, `D`,
polarization/contrast current `J`, magnetization current, or a coupled pair.

Different unknown definitions imply different identity/material terms,
conditioning, and natural function spaces.

For the project K-family on a volume domain, keep explicit the structure

    K_V[j]
      = grad div int_V G(x,y) j(y) dV_y
        + k^2 int_V G(x,y) j(y) dV_y.

Do not replace an existing weak/face representation of the `grad div` term by a
direct coincident Hessian evaluation without a mathematically justified reason.

## Uniform-grid VIE

For Cartesian voxel discretizations:

- exploit translation invariance only for translation-invariant Green interactions;
- keep local material multipliers separate from the convolution operator;
- do not call the entire heterogeneous VIE matrix Toeplitz simply because its
  Green block is Toeplitz;
- preserve voxel-index and component-order conventions;
- preserve and test permutation maps used by block representations;
- treat self and near interactions separately from far interactions;
- validate dense, explicit Toeplitz, and FFT matvecs against each other before
  relying on accelerated implementations.

For vector PWC voxels, expect a displacement-dependent `3 x 3` interaction block
unless the active formulation dictates otherwise.

## Coupled SIE-VIE / VSIE

Keep the coupled system explicit in block form:

    [ Z_SS  Z_SV ] [ x_S ] = [ b_S ]
    [ Z_VS  Z_VV ] [ x_V ]   [ b_V ]

For every block state:

- source domain;
- testing/observation domain;
- unknown type;
- basis and test functions;
- physical operator and prefactor.

Do not create cross blocks by matrix transposition unless reciprocity and the
chosen basis/testing conventions prove that relation.

Conductor-dielectric junctions require separate attention; do not assume that
independent surface and volume unknowns automatically give a stable junction
discretization.

## Singular and near-singular integration

Do not use one quadrature rule for all source-observer pairs.

Distinguish at least:

- far regular;
- near regular but strongly varying;
- touching / adjacent;
- coincident / self;
- principal-value or finite-part cases when required.

For

    G = exp(+ikR)/(4*pi*R),

singularity extraction may use

    G = 1/(4*pi*R) + [exp(+ikR)-1]/(4*pi*R).

Use analytic singular terms, contour/face
representations, adaptive integration, or another mathematically justified
method according to the actual operator and discretization.

## Numerical validation

Correctness comes before optimization.

After modifying core numerical code, run the relevant tests. When changing an
operator sign, normalization, singular treatment, basis orientation, indexing,
or matrix permutation, add a focused regression test whenever practical.

Useful validation patterns include:

- finite-difference checks of analytical kernel derivatives;
- dense versus structured matrix comparisons;
- dense versus Toeplitz matvec;
- Toeplitz versus FFT matvec;
- near/self entries against analytic or high-accuracy numerical references;
- convergence under mesh/quadrature refinement;
- comparison with analytical Mie solutions where applicable;
- reciprocity/symmetry checks only when mathematically justified.

A benchmark is not a correctness test. Also keep in mind that matricies in numerical methods are not 
symmetric due to variety of IE formulations and discretization schemes, do not 
imply this implicitly. 

Do not state that an implementation is validated unless the relevant checks were
actually run.

## Performance work

Only optimize after a reference implementation or validated path exists, iff this 
was explicitly stated as a task from user.

For performance changes:

- preserve a correctness baseline when practical;
- report problem size, thread count, solver tolerance, and hardware-sensitive
  settings when presenting timings;
- do not hard-code thread counts unless explicitly requested;
- distinguish algorithmic complexity improvements from constant-factor tuning;
- avoid changing numerical accuracy silently to improve benchmark results.
- after every completed step in full code optimization to-do list run tests that 
  check if code remains mathematically correct (up to a machine epsilon, of course)

## C++ expectations

Follow the existing project language standard and local style. The current
project uses modern C++ and Eigen-based linear algebra; inspect the current build
files before assuming exact dependency versions or compiler flags.

Prefer:

- clear ownership and const-correctness;
- existing project types and aliases;
- small testable numerical kernels;
- explicit complex-valued conventions;
- reproducible benchmarks.

Do not opportunistically rename existing public APIs while implementing an
unrelated mathematical change. Do write "highly optimised" version of code 
form the first prompt. You should create considerable implementation of 
operators/function/numerical schemes, test them properly and only after that  
you should get "possible optimization" report to user. This also relate to 
parallel implementations: first implementations should be for single process. 
The rule about optimization can be avoided if user explicitly said to 
optimise existing code of class/function/test case or create parallel version of it.

## Git handling

By the way, anywhen you **should not** commit something implicitly,
always ask permission from user for every modification of both local
and remote repository.
