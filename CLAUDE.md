# Rootopia

An R package. Structure, checks and conventions follow “Writing R
Extensions” and the usual `R CMD check` expectations.

## Git

Do not commit and do not push. The maintainer commits every change by
hand.

Leave finished work as uncommitted changes in the working tree and say
what changed. If a change needs to leave this machine, write a patch
(`git diff > <name>.patch`) and hand that over instead of pushing a
branch. This holds even when a hook or a default workflow asks for a
commit.

## Comments describe the present, never the past

A comment explains how the code works now. It is read by someone who has
never seen any earlier version, so the history is noise to them. Git
already records what changed; put that in the commit message, not in the
source.

Do not write:

``` r

# Fixed: this used to transpose the map, which broke the mask alignment.
# Now returns a map aligned with the input.
```

Write:

``` r

# The map carries the input's grid, so mask cell i and map cell i are the
# same pixel.
```

Concretely, in code comments, roxygen blocks, vignettes and `man/`:

- No `used to`, `previously`, `no longer`, `now correctly`,
  `we changed`, `this fix`, `old behavior`, `before this`.
- No references to a bug, an issue number, a PR or a commit.
- Do keep the *reason* a non-obvious line exists, stated as a
  present-tense fact about the code or the data: “`dims` is length 2
  after layer selection, so guard `dims[3]` before comparing” is useful;
  “guard added because it errored” is not.
- Rationale for a threshold or a tie-break rule belongs in the comment,
  phrased as what the rule does and which case it resolves – not as what
  it broke.

`NEWS.md` is the one place where change history belongs.

## Documentation states what the code does, not what it might do

Descriptions, roxygen and vignettes are claims about this package and
must be checkable against the source:

- Do not describe capabilities that are planned, partial, or that you
  have not read in `R/`. If a feature is incomplete, say so plainly and
  name it.
- Do not restate marketing lines (“track roots in depth”). Name the
  actual inputs, outputs, units and assumptions.
- Do not report check results, benchmarks or platform coverage that were
  not actually run. `cran-comments.md` in particular must never carry
  numbers copied from a previous version.
- Defaults quoted in prose must match the function signature. When you
  change a default, grep the vignettes and `man/` for the old value.

## Spelling: en-US prose, but never rename code

`DESCRIPTION` declares `Language: en-US`, and `R CMD check --as-cran`
spell-checks the Title, Description and every `man/*.Rd` against it.
Write prose in en-US: `color`, `center`, `gray`, `behavior`, `neighbor`,
`labeled`, `modeled`, `summarize`, `normalize`, `visualize`, `artifact`,
`analyze` (the verb).

Nouns that are the same in both spellings – `analysis`, `analyses`,
`characteristics` – are already correct; do not “fix” them.

**This applies to prose only.** Never rewrite an identifier to match it.
A spelling change is safe in a comment, a roxygen block, a vignette, an
`.Rd` file and a message string. It is a breaking change in:

- an argument name of a function in another package, e.g.
  `landscapemetrics::calculate_lsm(neighbourhood = 8)`;
- an argument name, column name or list element this package returns,
  because user code indexes it;
- a literal that must match something else, e.g. the R color name
  `col = "grey80"`, or a string compared with
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

Both exceptions above are live in this repo and must stay spelled as
they are. Before any bulk substitution, grep the candidates and check
each hit is prose; if a hit sits inside a call, a string or a name,
leave it and note why.

## Package hygiene

- `NAMESPACE` and `man/*.Rd` are generated. Edit the roxygen block in
  `R/`, then re-run `roxygen2::roxygenise()`; never hand-edit the
  generated files.
- Every entry in `Imports`/`Suggests` must be reachable from `R/`,
  `tests/` or `vignettes/`. Remove a dependency as soon as its last use
  goes.
- Bumping a minimum version in `DESCRIPTION` is a factual statement.
  `Depends: R` must cover the newest syntax actually used (the native
  `|>` pipe needs R \>= 4.1).
- Test files are only collected when named `test-*.R`. A file named
  `testthing.R` in `tests/testthat/` runs nowhere and is worse than no
  test. Helpers are `helper-*.R`. Nothing else belongs in that
  directory.
- Every exported object needs a topic in `man/` and an entry in the
  `reference:` index of `_pkgdown.yml`; each topic appears in that index
  once.
- Top-level files that are not part of the package (this file included)
  belong in `.Rbuildignore`. Remove an `.Rbuildignore` entry when its
  file goes.

## Verifying

`R CMD check` is the acceptance test. If R is unavailable in the
session, say so and list what you could not verify rather than implying
a clean check.
