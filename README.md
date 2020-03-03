# TriSym.jl

A [Nemo](http://nemocas.org/)/[Julia](https://julialang.org/) package to
compute trisymmetric decompositions.

## Installation

The package can be installed via the package manager `Pkg`, that is entered with
the key `]`.

```
(v1.2) pkg> add https://github.com/erou/TriSym.jl
```

## Usage

The package allows to work with Nemo's finite fields but also with *abstract*
commutative algebra over finite fields, *i.e.* algebras given by their
multiplication tensor.

### Finite fields

Here is an example with the finite field with 27 elements where we search for
decomposition of length 6 with two different *margins*.

```
julia> k, x = FiniteField(3, 3, "x")
(Finite field of degree 3 over F_3, x)

julia> multiplication_search(k, (0, 0, 0), 6)
1-element Array{Tuple{Vararg{Tuple{fq_nmod,Int64},N} where N},1}:
 ((x^2+1, 1), (x^2+2*x+1, 2), (2*x^2+x+1, 1), (x^2+x, 2), (2*x^2+x, 2), (x^2, 1))

julia> multiplication_search(k, (2, 1, 0), 6)
4-element Array{Tuple{Vararg{Tuple{fq_nmod,Int64},N} where N},1}:
 ((x+1, 1), (2*x+1, 2), (2*x^2+1, 1), (x, 1), (x^2+x, 2), (2*x^2+x, 1))                    
 ((x^2+1, 1), (x^2+2*x+1, 2), (2*x^2+x+1, 1), (x^2+x, 2), (2*x^2+x, 2), (x^2, 1))          
 ((x^2+1, 2), (x^2+x+1, 1), (2*x^2+1, 2), (2*x^2+2*x+1, 2), (2*x^2+x, 1), (x^2, 2))        
 ((x^2+x+1, 2), (x^2+2*x+1, 1), (2*x^2+1, 1), (2*x^2+x+1, 2), (2*x^2+2*x+1, 1), (x^2+x, 1))
```

### Abstract algebras

Here is the same example as before but in the case of *abstract algebras*, where
the user can specify the multiplication tensor. In this case we give the exact
same finite field `k`, but it also possible to call `AbstractCommutativeAlgebra`
with the matrix of the multiplication tensor.

```
julia> k
Finite field of degree 3 over F_3

julia> A = TriSym.AbstractCommutativeAlgebra(k)
TriSym.AbstractCommutativeAlgebra{3}([1  0  0  0  2  0]
[0  1  0  0  1  2]
[0  0  1  1  0  1], 3)

julia> L = MatrixSpace(ResidueRing(ZZ, 3), 1, 3)([ResidueRing(ZZ, 3)(coeff(tr(x^i), 0)) for i in 0:3-1])
[0  0  2]

julia> multiplication_search(A, L, (0, 0, 0), 6)
1-element Array{Tuple{Vararg{Tuple{nmod_mat,Int64},N} where N},1}:
 (([1]
[0]
[1], 1), ([1]
[2]
[1], 2), ([1]
[1]
[2], 1), ([0]
[1]
[1], 2), ([0]
[1]
[2], 2), ([0]
[0]
[1], 1))
```

## Tests

In the package manager, the package can be tested via

```
(v1.2) pkg> test TriSym
```

## Help

In Julia, *help* can be obtained by pressing the `?` key.

```
help?> multiplication_bilinear_map
search: multiplication_bilinear_map

  multiplication_bilinear_map(k::Nemo.FqNmodFiniteField)

  Compute the bilinear map associated with the multiplication of a finite field k.

  ─────────────────────────────────────────────────────────────────────────────────────────────────────

  multiplication_bilinear_map(A::AbstractCommutativeAlgebra{N}) where N

  Compute the bilinear map assciated with the multiplication in the algebra A.
```
