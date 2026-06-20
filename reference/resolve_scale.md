# Resolve the single \`scale\` argument, mapping the deprecated logical flags

\`normalize\`, \`binarize\`, and \`denormalize\` are kept only for
backward compatibility. Exactly one source of truth is allowed: either
\`scale\` or the old flags, never both.

## Usage

``` r
resolve_scale(scale, normalize, binarize, denormalize, scale_given)
```
