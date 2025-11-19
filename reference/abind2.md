# Combine multi-dimensional arrays

Combine multi-dimensional arrays. This is a generalization of cbind and
rbind. Takes a sequence of vectors, matrices, or arrays and produces a
single array of the same or higher dimension.

## Usage

``` r
abind2(
  ...,
  along = N,
  rev.along = NULL,
  new.names = NULL,
  force.array = TRUE,
  make.names = use.anon.names,
  use.anon.names = FALSE,
  use.first.dimnames = FALSE,
  hier.names = FALSE,
  use.dnns = FALSE
)
```

## Arguments

- ...:

  Any number of vectors, matrices, arrays, or data frames. The
  dimensions of all the arrays must match, except on one dimension
  (specified by `along=`). If these arguments are named, the name will
  be used for the name of the dimension along which the arrays are
  joined. Vectors are treated as having a dim attribute of length one.

  Alternatively, there can be one (and only one) list argument supplied,
  whose components are the objects to be bound together. Names of the
  list components are treated in the same way as argument names.

- along:

  (optional) The dimension along which to bind the arrays. The default
  is the last dimension, i.e., the maximum length of the dim attribute
  of the supplied arrays. `along=` can take any non-negative value up to
  the minimum length of the dim attribute of supplied arrays plus one.
  When `along=` has a fractional value, a value less than 1, or a value
  greater than N (N is the maximum of the lengths of the dim attribute
  of the objects to be bound together), a new dimension is created in
  the result. In these cases, the dimensions of all arguments must be
  identical.

- rev.along:

  (optional) Alternate way to specify the dimension along which to bind
  the arrays: `along = N + 1 - rev.along`. This is provided mainly to
  allow easy specification of `along = N + 1` (by supplying
  `rev.along=0`). If both `along` and `rev.along` are supplied, the
  supplied value of `along` is ignored.

- new.names:

  (optional) If new.names is a list, it is the first choice for the
  dimnames attribute of the result. It should have the same structure as
  a dimnames attribute. If the names for a particular dimension are
  `NULL`, names for this dimension are constructed in other ways.

  If `new.names` is a character vector, it is used for dimension names
  in the same way as argument names are used. Zero length ("") names are
  ignored.

- force.array:

  (optional) If `FALSE`, rbind or cbind are called when possible, i.e.,
  when the arguments are all vectors, and along is not 1, or when the
  arguments are vectors or matrices or data frames and along is 1 or 2.
  If rbind or cbind are used, they will preserve the data.frame classes
  (or any other class that r/cbind preserve). Otherwise, abind will
  convert objects to class array. Thus, to guarantee that an array
  object is returned, supply the argument `force.array=TRUE`. Note that
  the use of rbind or cbind introduces some subtle changes in the way
  default dimension names are constructed: see the examples below.

- make.names:

  (optional) If `TRUE`, the last resort for dimnames for the along
  dimension will be the deparsed versions of anonymous arguments. This
  can result in cumbersome names when arguments are expressions.

  \<p\>The default is `FALSE`.

- use.anon.names:

  (optional) `use.anon.names` is a deprecated synonym for `make.names`.

- use.first.dimnames:

  (optional) When dimension names are present on more than one argument,
  should dimension names for the result be take from the first available
  (the default is to take them from the last available, which is the
  same behavior as `rbind` and `cbind`.)

- hier.names:

  (optional) If `TRUE`, dimension names on the concatenated dimension
  will be composed of the argument name and the dimension names of the
  objects being bound. If a single list argument is supplied, then the
  names of the components serve as the argument names. `hier.names` can
  also have values `"before"` or `"after"`; these determine the order in
  which the argument name and the dimension name are put together
  (`TRUE` has the same effect as `"before"`).

- use.dnns:

  (default `FALSE`) Use names on dimensions, e.g., so that
  `names(dimnames(x))` is non-empty. When there are multiple possible
  sources for names of dimnames, the value of `use.first.dimnames`
  determines the result.

## Value

merged multidimensional arrays

## Details

The dimensions of the supplied vectors or arrays do not need to be
identical, e.g., arguments can be a mixture of vectors and matrices.
`abind` coerces arguments by the addition of one dimension in order to
make them consistent with other arguments and `along=`. The extra
dimension is added in the place specified by `along=`.

The default action of abind is to concatenate on the last dimension,
rather than increase the number of dimensions. For example, the result
of calling abind with vectors is a longer vector (see first example
below). This differs from the action of `rbind` and cbind which is to
return a matrix when called with vectors. abind can be made to behave
like cbind on vectors by specifying `along=2`, and like rbind by
specifying `along=0`.

The dimnames of the returned object are pieced together from the
dimnames of the arguments, and the names of the arguments. Names for
each dimension are searched for in the following order: new.names,
argument name, dimnames (or names) attribute of last argument, dimnames
(or names) attribute of second last argument, etc. (Supplying the
argument `use.first.dimnames=TRUE` changes this to cause `abind` to use
dimnames or names from the first argument first. The default behavior is
the same as for `rbind` and `cbind`: use dimnames from later arguments.)
If some names are supplied for the along dimension (either as argument
names or dimnames in arguments), names are constructed for anonymous
arguments unless `use.anon.names=FALSE`.

sourced from the 'abind' package:
https://doi.org/10.32614/CRAN.package.abind under MIT-license

## Author

Tony Plate <tplate@acm.org> and Richard Heiberger
