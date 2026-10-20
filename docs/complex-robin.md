# Multiple potentials, complex problems and Robin conditions

## Build and run with deal.II 9.8.0

From the repository root, build both configurations in the supplied container:

```sh
docker run --rm -v "$PWD:/workspace" -w /workspace \
  dealii/dealii:v9.8.0-resolute-r1-arm64 bash -lc '
    cmake -S . -B build-docker/9.8.0 -GNinja \
      -DCMAKE_BUILD_TYPE=DebugRelease -DCMAKE_CXX_FLAGS=-Werror &&
    cmake --build build-docker/9.8.0 -j 4 &&
    cd build-docker/9.8.0 && ctest --output-on-failure
  '
```

Executables are `bem_fma_2d` and `bem_fma_3d` in the build directory;
Debug executables have the `.g` suffix. The first positional argument is the
number of threads, not a parameter filename. Each program reads
`parameters_bem_2.prm` or `parameters_bem_3.prm` in its working directory and
writes `used_parameters_bem_2.prm` or `used_parameters_bem_3.prm`.

The sphere examples use paths relative to `samples/`. Keep each run's output
separate, and copy both the parameter file and its grid directory:

```sh
mkdir -p build-docker/9.8.0/example-complex
cp -R samples/grids build-docker/9.8.0/example-complex/
cp samples/ref_params_sphere_complex_DR_direct.prm \
  build-docker/9.8.0/example-complex/parameters_bem_3.prm
docker run --rm -v "$PWD:/workspace" \
  -w /workspace/build-docker/9.8.0/example-complex \
  dealii/dealii:v9.8.0-resolute-r1-arm64 bash -lc '
    mpirun -np 2 ../bem_fma_3d 2
  '
```

The default container user can run MPI. If running as root, also pass
`-e OMPI_ALLOW_RUN_AS_ROOT=1 -e OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1` to Docker.
The samples cover real/complex potentials, DN (Dirichlet/Neumann), DR
(Dirichlet/Robin), NR (Neumann/Robin), and `direct`/`fma` methods.
`Direct` assembles integral matrices and uses an iterative solver; it does not
mean a dense direct factorization. FMA is supported in 3D only.

## Hosted CI environment

Local validation uses `dealii/dealii:v9.8.0-resolute-r1-arm64`, built from
the frozen recipes in [docker-files PR 90](https://github.com/dealii/docker-files/pull/90).
Until that image is published, hosted CI builds `docker/Dockerfile.ci` with the
exact deal.II 9.8.0 backports package on Ubuntu 26.04. This supports hosted AMD64
runners without depending on a tag that currently exists only locally.
The formatting job retains clang-format 16.0.6 from its existing image.

To reproduce the hosted environment locally, build it and substitute
`pi-bem-ci:9.8.0` for the image in the commands above:

```sh
docker build -f docker/Dockerfile.ci -t pi-bem-ci:9.8.0 .
```

The deal.II 9.8 update replaces the removed three-argument Trilinos vector
`reinit(source,false,true)` with assignment into the already initialized
ghost layout. The assignment imports both locally owned and ghost values;
it must participate collectively on the vector communicator.

## Components and boundary-function slots

`Potential components` counts real-valued storage components. Each complex
potential uses consecutive components: its real part at `i`, imaginary part at
`i+1`. `Complex components` lists the real indices, with no overlapping pairs:

```text
subsection Driver<3>
  set Potential components = 3
  set Complex components = 0
  set Set Global Refinement = true
end
```

This solves one complex potential in components 0 and 1, followed by a real
potential in component 2. Parameter declarations provide up to 14 storage
components. Every component has two boundary-function slots, numbered 0 and 1.

Under `BoundaryConditions<dim>`, old function subsection names must be replaced
by `Potential <component> <slot> <dim>d`, `Wind function <component> <slot> <dim>d`
and `Robin coefficients <component> <slot> <dim>d`. Potential expressions are
scalar; wind expressions have `dim` entries and supply the normal derivative
through their dot product with the outward normal. Functions used by the solve
are evaluated at time zero. Output uses the configured stem for component 0
and appends `_i` for later components.

## Robin conditions and geometry IDs

The 3D Robin function has three entries `a; b; g`, expressing
`a*phi + b*dphi_dn = g`. The solver stores `a/b` and `g/b`, so `b` must be
nonzero at every Robin node. For complex problems, the active component
contains the real parts of all three coefficients and the following component
contains their imaginary parts; the division is complex division.

For example, the condition `(2+3i)*phi + dphi_dn = 4+5i` in slot 0 is:

```text
subsection BoundaryConditions<3>
  subsection Robin coefficients 0 0 3d
    set Function expression = 2; 1; 4
  end
  subsection Robin coefficients 1 0 3d
    set Function expression = 3; 0; 5
  end
end
```

Boundary-type lists select mesh **manifold IDs**. Their names retain the word
"boundary" for compatibility. The three lists must be disjoint. The slot map
is a flat comma-separated list of `manifold, slot` pairs:

```text
subsection ComputationalDomain<3>
  subsection Boundary Conditions ID Numbers
    set Dirichlet boundary ids = 1
    set Neumann boundary ids =
    set Robin boundary ids = 2
    set Manifold to boundary slot map = 1,0, 2,1
  end
end
```

This uses slot 0 on manifold 1 and slot 1 on manifold 2, for every component.
Supply the corresponding potential/wind/Robin functions for each selected slot.
The bundled Robin samples are 3D: the present 2D parser allocates Robin
functions with two entries while Robin preparation requires three. Do not use
2D Robin conditions without resolving that implementation limitation.

## C++ API migration

All distributed operations must be invoked by every rank in the object's
communicator, with compatible vector ownership. The geometric dimension is
`dim`, and the boundary mesh dimension is `dim-1`.

| Interface | Contract and migration |
| --- | --- |
| `Driver(n_components)`, `BEMProblem(domain,comm,n_components)`, `BoundaryConditions(domain,bem,comm,n_components)` | Counts refer to real storage components; defaults preserve a single real potential. Keep referenced objects alive. |
| `n_phi_components()`, `set_n_phi_components(n)` | Query or grow component storage; shrinking is unsupported. Grow before initialization/parameter parsing. |
| `get_current_phi_component()`, `set_current_phi_component(i)` | Select a zero-based storage component. Selecting through `BoundaryConditions` also selects the bound BEM component. A complex solve requires `i+1<n`. |
| `get_phi([i])`, `get_dphi_dn([i])` | Mutable vector references for the current or explicit component replace direct access to single solution fields. Storage growth can invalidate references. |
| `get_wind(i,slot)`, `get_potential(i,slot)`, `get_robin_coeffs(i,slot)` | Parsed-function access requires initialized parameters, valid component and slot. Const public overloads expose the functions; mutable overloads are protected. |
| `get_phi_components([indices])` | Returns an owned block-vector copy collectively. Explicit indices retain their component block numbers; the no-argument overload selects all. |
| `solve_problem(reset_matrix)`, `solve_complex_problem(reset_matrix)` | Prepare data and solve current component/pair. Default `true` assembles the integral operators. Use `false` only while geometry, mapping, quadrature and discretization remain unchanged. |
| `BEMProblem::solve(...)` | Real overload gains `reset_matrix`; complex overload takes `phi, phi_imag, dphi_dn, dphi_dn_imag, rhs, rhs_imag, reset_matrix`. Boundary vectors and Robin data must already be prepared. |
| `prepare_robin_datastructs(...)` | Real overload fills diagonal/RHS; complex overload fills real and imaginary diagonal/RHS from consecutive component functions. |
| `vmult(...)`, `compute_rhs(...)`, `solve_system(...)` | Additional overloads operate on separate real/imaginary vectors. They require previously initialized operators and current boundary data. |
| `assemble_system_tbb()`, `assemble_system_singlethread()` | Explicit threaded/serial direct assembly. The high-level solve selects according to the thread limit. |
| `assemble_preconditioner_complex()`, `get_preconditioner(complex)` | Set up and select configured ILU/AMG storage. Selection returns a reference owned by the BEM problem. |
| `get_vector_gradients_solution([i])`, `get_vector_surface_gradients_solution([i])` | Mutable component-aware gradient access replaces the old single gradient fields. |
| `get_system_matrix()`, `get_system_matrix_complex()` | Collective dense diagnostic copies, size N or 2N, ordered real then imaginary. Require direct matrices, use quadratic memory per rank and omit constrained-row replacement. |
| `compute_errors(complex,current_is_real)` | Replace the old no-argument call: `(false,true)` for real, `(true,true)` for the real part of a pair, `(true,false)` for its imaginary part. |
| `output_results(stem)` | Writes the selected component after errors/gradients have been computed. |

For an initialized boundary object `bc`, component selection is explicit:

```cpp
bc.set_current_phi_component(0);
bc.solve_complex_problem(true);
bc.compute_errors(true, true);
bc.output_results("potential_real");
bc.set_current_phi_component(1);
bc.compute_errors(true, false);
bc.output_results("potential_imag");
```

The old `bem_problem_access.h` adapter and `time_integrator.h` are removed.
Use the component-aware public accessors; no time-integration replacement is
provided. `ComputationalDomain` now exposes `load_cad_objects()` (returns the
largest CAD tolerance), `refine_and_resize_by_aspect_ratio()`,
`refine_and_resize_by_cad_projections(max_tol)` and
`refine_and_resize_wrapup()` (partition and mesh output). Normal applications
should call `refine_and_resize(level)`, which runs these stages in order.
The option for attaching flat manifolds to IDs without CAD data preserves
import of grids whose flags also distinguish boundary conditions.

## Complex operators and FMA APIs

`ConstrainedComplexOperator<VEC,MATRIX>` accepts one flat vector ordered
`[real,imaginary]`. `BlockConstrainedComplexOperator<VEC,MATRIX>` accepts a
two-block vector in the same order. Their constructors take a matrix, real and
imaginary constraints, a constraint-relevant IndexSet and a communicator.
Both constraint objects must describe the same constrained rows with unshifted
indices. The matrix provides `this_cpu_set` and a four-vector `vmult`.
`distribute_rhs()` installs the corresponding inhomogeneities.
`BlockPreconditioner<P>` applies a referenced preconditioner independently to
each block through `vmult`/`Tvmult`. Referenced matrices, constraints and
preconditioners must outlive these adapters. The flat operator keeps mutable
scratch vectors and must not be used concurrently from multiple threads.

`BEMFMA` adds `_tbb` and, when compiled with `PBEM_OPENMP`, `_omp` backends for
direct integrals, multipole matrix-vector products and preconditioning.
Use the dispatching `direct_integrals()`, `multipole_matr_vect_products()` and
`FMA_preconditioner()` interfaces in application code.
`FMA_preconditioner_complex()` provides the complex-solve near-field ILU
preconditioner, with corresponding backend methods. Returned references belong
to the FMA object. TBB is the default; the optional OpenMP configuration is
selected by `PBEM_OPENMP=1` at CMake configure time.

`LocalExpansion` and `MultipoleExpansion` now own coefficient vectors and use
default value copy/assignment. `GetCoeffs()` returns a const vector reference,
not a writable pointer. `MultipoleExpansion::GetA_n_m()` and `GetCoeff()` also
return const references; use `SetCoeff()` or `AddToCoeff()` for coefficient
changes. Coefficients use triangular index `n*(n+1)/2+m`.
`Add()`, `Evaluate()` and `MultipoleExpansion::AddNormDer()` gain overloads with
a trailing `std::vector<std::complex<double>>&` scratch cache. An empty cache
is valid; calls resize/overwrite it. Keep one cache per concurrent operation.
The `AssLegFunction` pointer is non-owning and must outlive the expansion.

```cpp
AssLegFunction legendre;
MultipoleExpansion expansion(6, Point<3>(0, 0, 0), &legendre);
std::vector<std::complex<double>> scratch;
expansion.Add(1.0, Point<3>(0.1, 0, 0), scratch);
const double value = expansion.Evaluate(Point<3>(2, 0, 0), scratch);
const auto coefficients = expansion.GetCoeffs(); // independent snapshot
```

`LocalExpansion::equal(lhs,rhs,tolerance)` compares orders, centers and stored
coefficients with absolute tolerance (default `1e-10`).
`MultipoleExpansion::spherical_coords()` returns relative position, radius,
polar cosine and azimuth through output arguments; distinct points are required.

`OctreeBlock::small_set` is a sorted `boost::container::flat_set` of global IDs.
Node, quadrature and interaction-list getters return const references owned by
the block. Copy explicitly when a snapshot is required; mutation invalidates
iterators as dictated by the underlying container. `GetCenter()` returns the
box center. The custom copy/destructor and `CopyContent()` are removed in favor
of value copying; `NumIntListLevels()`/`NumNonIntListLevels()` are removed.
The complete interaction-list getter can be queried with `.size()`.

## Geometry and error-integration helpers

`QuasiSingularKernelIntegral<3>` searches a mapped quadrilateral for the closest
point to an external point. Its constructor initializes reference coordinates
`eta` and physical `min_distance`; the search considers the interior, edges and
vertices. The referenced iterator object, mapping and finite element must remain
alive. The 2D constructor does not initialize a projection or distance. This
helper currently provides closest-point geometry, not a public quadrature API.

`integrate_difference_based_on_material_id()` and its hp implementation
`do_integrate_difference_based_on_material_id()` select non-owning exact-function
pointers by **material ID**. They return cellwise differences and put zero on
unmapped/non-owned cells. Function component counts must match the finite
element. Compute the global norm separately with
`VectorTools::compute_global_error()`. These helpers depend on deal.II internal
integration interfaces and require checking when updating deal.II.
