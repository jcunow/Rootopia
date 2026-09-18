# How It Works: Modules, Pipelines and Rules

## How It Works

The other vignettes show you *what to type*. This one explains *what
happens* — which function hands what to which, and the exact rule each
one applies when it has to make a decision.

### How to read this

Every module follows the same four-part shape:

- **What it is for** — one paragraph, no jargon.
- **The flow** — how the functions connect, and in what order.
- **The rules** — the decisions the code makes, stated exactly. If a
  number in your results surprises you, the explanation is here.
- **Collapsed sections** — the algorithm itself, and the edge cases.
  Open them when you need them; skip them otherwise.

Edge cases are deliberately kept out of the main text. They matter when
something goes wrong, and they get in the way when it does not.

Where Rootopia uses a standard method, this document names the method
and moves on. Detail is spent on the parts you will not find in a
textbook.

### The package pipeline

Most analyses run left to right through this chain. You can enter at any
point and stop at any point; nothing forces you through the whole thing.

     stitch  ->  clean  ->  (rotation censor)
                    |
                    v
               skeletonise  ---->  branching order
                    |
                    +-------->  size traits (length, diameter, pixels)
                    |
                    v
                depth map  ---->  root angle
                    |
                    v
            distribution indices

Turnover and soil & colour sit outside this chain: turnover compares two
images, soil & colour works on the colour channels rather than the root
shape.

[`root_depth_metrics()`](https://jcunow.github.io/Rootopia/reference/root_depth_metrics.md)
is not a module. It is a batch wrapper that runs branching order, size
traits, depth mapping, root angle and the distribution indices over a
folder of images. Read those sections to understand what it reports.

------------------------------------------------------------------------

### Preparing images

Optional steps that happen before any measurement.

#### Stitching

##### What it is for

A tube is often imaged as several overlapping frames. This module joins
them back into one long mosaic before analysis, so a root crossing a
frame boundary is measured once rather than twice.

##### The flow

``` r

stitch_root_scans("path/to/scans", pattern = ".tiff", tubes = "ask")
res <- stitch_root_scans("path/to/scans", pattern = ".tiff",
                         out_dir = "out", report = TRUE)
res$report          # per-join dx, dy, peak, overlap
```

##### The rules

**Files are grouped into tubes by a pattern in the filename.**
`group_regex` (default `"T0\\d{2}"`, matching labels like `T067`) pulls
a tube id from each path; files sharing an id form one sequence, sorted
by filename.

**Consecutive frames are aligned by FFT phase correlation** on a band
along the overlapping edge, then composited with a linear feather blend
across the overlap so the seam does not show.

**Frames are joined along the image width** by default. For sequences
acquired along the tube, set `direction = "vertical"`; frames are
transposed internally, stitched the same way, and transposed back.

**Check the `peak` column before trusting a mosaic.** It is the
normalised correlation height at each join — low values mean the
alignment was uncertain. Sorting the report by `peak` surfaces the weak
joins first.

**Edge cases**

| Situation | What happens |
|----|----|
| A tube has only one frame | Passed through unchanged |
| No file matches `group_regex` | Stops |
| `tubes` index out of range | Stops, naming how many tubes exist and what they are called |
| `tubes = "ask"` in a non-interactive session | Stops, suggesting indices or names instead |
| Poor alignment | Bring `edge_width` closer to the true overlap, raise `vertical_offset` past a header strip, or try `preprocess = "grad"` for uneven lighting |

#### Image input and cleaning

##### What it is for

Getting any image into a form the rest of the package can use, and
removing the speckles and pinholes that segmentation leaves behind. Both
matter more than they sound: a hole in a root becomes a fake loop in the
skeleton, and a speck becomes a fake root.

##### The flow

``` r

img     <- load_flexible_image(path, output_format = "spatrast",
                               scale = "binary", select_layer = 2)
cleaned <- clean_image(img, max_hole_size = 5, max_artifact_size = 5)
report_image_components(img)      # how big are the specks, before deciding
```

[`load_flexible_image()`](https://jcunow.github.io/Rootopia/reference/load_flexible_image.md)
is called internally by nearly every other function, so you rarely call
it yourself except when loading from a file path.

##### The rules

**One entry point, many input types.** File paths, `SpatRaster`,
`RasterBrick`, matrices, arrays, `cimg` and magick images all converge
to one internal representation. `scale` controls the value range:
`"binary"` maps to 0/1, `"to_01"` divides by 255, `"to_255"` multiplies,
`"none"` leaves values alone. Conversions use fixed factors, never a
per-image maximum, so two images stay comparable.

**A hole is not the same as a gap.**
[`clean_image()`](https://jcunow.github.io/Rootopia/reference/clean_image.md)
fills a background region only if it is fully enclosed — a background
region touching the image border is the outside world, not a hole, and
is never filled regardless of size.

**An artifact is a small connected foreground region.** By default
*every* foreground region is a candidate, including ones touching the
border. Set `protect_border = TRUE` to exempt border-touching regions,
on the assumption that they are roots leaving the frame rather than
specks.

**Size is in pixels**, counted as connected components. Use
[`report_image_components()`](https://jcunow.github.io/Rootopia/reference/report_image_components.md)
first to see the actual size distribution rather than guessing a
threshold.

**Under the hood**

Connected-component labelling and the optional edge smoothing are
standard morphology, done by `imager` (`label`, `dilate`, `erode`).
Rootopia adds the hole-versus-outside rule and the size thresholds; it
does not reimplement the morphology. See
[`?clean_image`](https://jcunow.github.io/Rootopia/reference/clean_image.md)
for the kernel options.

**Edge cases**

| Situation | What happens |
|----|----|
| Image has no cells | Stops |
| `max_hole_size` / `max_artifact_size` left `NULL` | That cleaning step is skipped entirely |
| Non-binary input to [`clean_image()`](https://jcunow.github.io/Rootopia/reference/clean_image.md) | Binarised on load |
| Unsupported file extension | Stops, listing the supported ones |

#### Rotation bias

##### What it is for

A minirhizotron tube can rotate slightly between sessions, and the
scanner does not cover the full 360°. Left alone, this means roots near
the edge of the scan are visible in one session and not the next — a
systematic bias that looks like root growth or death. This module
measures the rotation and crops it away.

The same slicing machinery also answers a different question: are roots
distributed evenly around the tube, or concentrated on one side?

##### The flow

``` r

r0    <- estimate_rotation_center(img)                 # once per tube
shift <- estimate_rotation_shift(img1, img2,           # per session pair
                                 cor_type = "phase")
kept  <- rotation_censor(img, center_offset = r0,
                         fixed_width = 800, fixed_rotation = TRUE)

slices <- slice_rotation(img, n = 48)                  # circumferential zones
```

##### The rules

**Finding the top of the tube.**
[`estimate_rotation_center()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_center.md)
looks for the white adhesive tape on the tube’s upper side. It builds a
bright reference band from a high quantile of the image, clusters the
result, and picks the cluster whose brightness exceeds `tape_brightness`
relative to the image. It returns a pixel **row** index. Tube geometry
is fixed, so this is done once per tube, not once per session.

**Measuring the shift.**
[`estimate_rotation_shift()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_shift.md)
compares two sessions by correlation over a shared depth window — either
cross-correlation (`"ccf"`) or phase correlation (`"phase"`, more robust
to brightness differences). It returns the depth-axis and rotation-axis
pixel offsets plus a peak height you can use as a confidence score.

**Cropping.**
[`rotation_censor()`](https://jcunow.github.io/Rootopia/reference/rotation_censor.md)
keeps a window of rows and discards the rest. Two modes:

- `fixed_rotation = FALSE` cuts proportionally to the measured offset,
  so the output width varies between images.
- `fixed_rotation = TRUE` centres a window of exactly `fixed_width` rows
  on a given row. Use this when comparing multiple sessions, because
  equal width is what makes counts comparable.

`center_offset` is read as a fraction of image height when between 0 and
1, and as an absolute row number when above 1.

**Slicing the circumference.** `slice_rotation(img, n)` cuts the image
into `n` equal bands of rows and returns them as a list, in rotation
order. Feeding a per-slice trait into
[`rhythmicity()`](https://jcunow.github.io/Rootopia/reference/rhythmicity.md)
(see Distribution indices) then tests whether roots favour one side of
the tube.

**Edge cases**

| Situation | What happens |
|----|----|
| `fixed_width` wider than the image allows | A message reports the maximum symmetric width and the crop is clamped to the image bounds — the output is then *not* the width you asked for |
| Cut would remove the whole image | Stops |
| Resulting window is empty | Warning, returns `NULL` |
| `n` greater than the number of rows | [`slice_rotation()`](https://jcunow.github.io/Rootopia/reference/slice_rotation.md) stops |
| Clustering fails in [`estimate_rotation_center()`](https://jcunow.github.io/Rootopia/reference/estimate_rotation_center.md) | Falls back to a direct brightness threshold |

> **Note on tube geometry.** Inner and outer tube diameters differ, so
> observed root length slightly underestimates true length in soil. No
> function here corrects for that; apply a resize coefficient yourself
> if you need it.

------------------------------------------------------------------------

### Core measurements

The shape of the root system, from a segmented image.

#### Skeletonisation

##### What it is for

A segmented root image is a solid shape. To measure architecture — where
roots branch, how they connect, which is a parent and which is a lateral
— you need the centre line rather than the shape. Skeletonisation thins
each root down to a line one pixel wide that runs along its middle,
keeping the connections intact.

##### The flow

``` r

skel   <- skeletonize_image(mask)            # solid mask -> 1-px centre line
points <- detect_skeleton_points(skel)       # find tips and branch points
clean  <- prune_skeleton(skel, mask,         # optional: remove short spurs
                         min_length = 10)
```

[`skeletonize_image()`](https://jcunow.github.io/Rootopia/reference/skeletonize_image.md)
is the entry point for everything downstream. The branching module will
call it for you if you pass a mask instead of a skeleton.

##### The rules

**Thinning.** Pixels are deleted from the outside inward, repeatedly,
until nothing more can go. Whether a pixel may be deleted depends only
on the pattern of its eight neighbours. Each of the 256 possible
patterns has a fixed verdict stored in a lookup table: keep, delete in
the first half-pass, delete in the second, or delete in both. Two
half-passes alternate so that the shape erodes evenly from both sides
rather than drifting; the neighbourhood is recomputed between them.

**Tips and branch points.** A skeleton pixel with exactly one neighbour
is a tip. One with three or more is a branch point. Pixels with exactly
two neighbours are ordinary line pixels and are neither.

**Spur pruning.**
[`prune_skeleton()`](https://jcunow.github.io/Rootopia/reference/prune_skeleton.md)
removes terminal branches shorter than `min_length` or thinner than
`min_diameter`. It works on the skeleton graph, because that is where
connectivity lives, but can return either the cleaned skeleton or the
cleaned solid mask.

**Under the hood — the thinning table**

The table is a Zhang–Suen *variant*, not the published table: it differs
at 25 of the 256 patterns. Twelve differences are more conservative —
the textbook deletes in both half-passes, this deletes in only one. The
other thirteen are staircase corners: a pixel whose only neighbours are
two orthogonal ones that already touch each other diagonally. The
textbook crossing-number test keeps those, which leaves a redundant
pixel on every 90° bend; deleting them gives a cleaner diagonal line and
cannot disconnect anything.

Two properties hold for every entry, and are checked in
`tests/testthat/test-skeleton.R` rather than taken on trust:

1.  No deletion can disconnect a pixel’s own neighbourhood.
2.  No pixel with a single neighbour — a root tip — is ever deleted, so
    thinning never shortens a root from its end.

For `output = "mask"`, each surviving skeleton pixel is regrown by its
local distance-transform radius and intersected with the original mask,
so real roots return to their true thickness while the spur’s body is
gone.

**Edge cases**

| Situation | What happens |
|----|----|
| Image has no foreground pixels | [`skeletonize_image()`](https://jcunow.github.io/Rootopia/reference/skeletonize_image.md) stops with “No foreground pixels” |
| Multi-layer raster | First layer used, after `select_layer` is applied |
| Input is not binary | `load_flexible_image(scale = "binary")` converts it first |
| Thinning does not converge | Stops after `max_iter` (default 200) passes |
| [`prune_skeleton()`](https://jcunow.github.io/Rootopia/reference/prune_skeleton.md) given no mask | Runs, but diameters collapse to ~1 px, so `min_diameter` is meaningless |

#### Branching order

##### What it is for

This module answers architectural questions: how many roots are there,
which are main axes and which are laterals, how often do they branch,
and how do length and diameter differ between generations. It turns a
skeleton into a network of segments and nodes, then labels every
segment.

This is the most involved module in the package, and the one where a
misunderstanding is most likely to produce a number you cannot explain.

##### The flow

One call does everything:

``` r

res <- branch_order_map(skel, mask = cleaned,
                        order = "branch_order",
                        unit = "cm", dpi = 300)
res$summary     # one row per order class
res$edges       # one row per segment
res$class_map   # raster, order value per root pixel
```

Underneath, in order:

    trace_segments()        skeleton -> segments and nodes
       |
    .splice_passthrough()   dissolve nodes that are not real branch points
       |
    resolve_crossings()     separate roots that merely cross
       |
    prune_terminal_segments()   optional, off by default
       |
    build_edge_table()      measure length and diameter per segment
       |
    compute_tip_order()     assign tip_order
       |
    assign_root_order()     group segments into roots; root_order, branch_order

Then
[`convert_root_units()`](https://jcunow.github.io/Rootopia/reference/convert_root_units.md)
rescales to cm,
[`summarize_orders()`](https://jcunow.github.io/Rootopia/reference/summarize_orders.md)
aggregates, and
[`order_classification_map()`](https://jcunow.github.io/Rootopia/reference/order_classification_map.md)
paints the result back onto the image.

Those seven steps are internal — they are named here so the order of
operations is visible, not because you call them. The two entry points
are
[`branch_order_map()`](https://jcunow.github.io/Rootopia/reference/branch_order_map.md)
and, if you want the segment table in pixels without the unit conversion
and summary,
[`root_graph_pipeline()`](https://jcunow.github.io/Rootopia/reference/root_graph_pipeline.md).

##### The rules

**Building the network.** Every skeleton pixel is classified by how many
neighbours it has: one means a tip, two means an ordinary line pixel,
three or more means a junction. *Adjacent junction pixels are merged
into a single node.* This matters — a real branch point is usually
several pixels wide in a thinned image, and treating each as its own
node would shatter the graph. Segments are the runs of two-neighbour
pixels between nodes.

**Segment length** is the distance along the pixel chain **plus** a
short stub at each end that lands on a merged junction, reaching from
the segment’s last pixel to the centre of the junction cluster. Without
that stub, every branch point would quietly eat about a pixel of root
from each arm that meets it.

**Three order schemes are always computed**, and they answer different
questions:

- `tip_order` — counted inward from the tips. Every terminal segment is
  1; peeling terminals away round by round, a segment’s order is one
  more than the highest of its children. This is centrifugal depth,
  *not* strict Strahler ordering: it increments at every junction, so a
  long chain of single laterals keeps climbing.
- `root_order` — segments are first grouped into continuous roots
  (below), then each root takes the highest `tip_order` along it. A
  thick main axis therefore keeps its high order all the way out to its
  tip.
- `branch_order` — generation counting. The thickest root in each
  disconnected piece of the image is order 1; its laterals are 2, theirs
  3, and so on. Independent of how deep any one subtree runs.

`order =` only chooses which of the three labels the summary and the
map; all three are always present in `res$edges`.

**Grouping segments into roots.** At a junction with three or more arms,
two of them are judged to be the same root continuing through, and the
rest are laterals. The pair chosen is the one maximising

> straightness + `diam_weight` × diameter similarity × relative
> thickness

where straightness compares the two arms’ directions, diameter
similarity is the thinner arm divided by the thicker, and relative
thickness is the thinner arm divided by the thickest arm at that node.
The last term is what stops two identical thin laterals from beating the
parent axis, which they otherwise would by scoring a perfect similarity
of 1.

**Crossings versus branches.** Where two roots cross in the image they
form a node with four arms, but they are not connected in reality. At
such nodes the three possible pairings are scored by straightness and
the straightest is chosen, but only if *both* of its pairs are
straighter than `crossing_straight` (default ≈ 120°). Otherwise the node
is left as a branch point.

**Under the hood — why crossings are genuinely ambiguous**

A thick axis crossed by a thin root, and a thick axis with two thin
laterals leaving in opposite directions, are the same shape in outline.
No threshold gets both right.

`crossing_diam_ratio` (off by default) lets you choose which error you
prefer: it requires the two candidate through-roots to be within that
thickness ratio of each other before accepting a crossing. Raising it
fixes the bilateral-branch case and breaks the crossing case. Set it
only if your material has more bilateral branching than fine-over-coarse
crossing.

Note also that only nodes with *exactly* four arms are examined. Nodes
with five or more are counted and reported under `verbose = TRUE`, but
left intact — they are then read as one continuation plus several
laterals.

**Under the hood — what “diameter” means here**

Diameters come from a distance transform of the **solid mask**, not the
skeleton: each pixel is labelled with its distance to the nearest
background pixel, so a skeleton pixel carries the radius of the largest
circle that fits inside the root there. Diameter is twice that.

This is why `mask =` matters. Pass only a skeleton and every diameter
collapses to about 1 px, and every diameter-dependent rule above — the
continuation choice, `branch_order`’s thickest-root seed — degenerates
with it.

**Edge cases**

| Situation | What happens |
|----|----|
| Image has no roots | Warning, empty result, no `class_map` |
| A ring with no tips (a genuine loop) | Those segments get `NA` order and a warning; they are excluded from summaries and counted in `attr(., "n_unordered")` |
| A solid mask passed as `skel` | Coverage check warns: nearly every pixel looks like a junction, so the graph comes out near-empty |
| Two separate root systems in one image | Each connected piece gets its own order-1 seed; they do not share numbering |
| A node with exactly two arms | Not a branch point — it is a bend or a thinning artefact, and the two arms are spliced back into one segment |
| Perfectly symmetric fork | The continuation score genuinely ties and the first pair enumerated wins. The `"fork"` validation phantom is excluded from per-root scoring for this reason |
| `class_map` pixel counts below the skeleton | Expected: merged junction interiors belong to no segment, leaving one unpainted pixel per branch point. Read lengths from `res$edges$length` |
| Pruning enabled | Thresholds use the same length definition the table reports, stub included |

#### Size traits

##### What it is for

The three quantities most analyses end up reporting: how much root there
is (pixels), how long it is (cm), and how thick it is (cm).

##### The flow

``` r

count_pixels(skel)                                  # raw foreground count
root_length(skel, unit = "cm", dpi = 300)           # length
root_diameter(cleaned, unit = "cm", dpi = 300)      # thickness
root_scape_metrics(img)                             # spatial pattern
```

Note which input each wants: **length from the skeleton, diameter from
the solid mask.** Swapping them silently gives wrong answers rather than
an error.

##### The rules

**Length is estimated, not counted.** Counting skeleton pixels
underestimates length, because a diagonal step covers √2 pixels of
distance but only one pixel of image. Rootopia counts orthogonal steps
(`No`) and diagonal steps (`Nd`) separately and offers four published
estimators:

| Method              | Formula                      |
|---------------------|------------------------------|
| `freeman_basic`     | √2·Nd + No                   |
| `freeman_corrected` | 0.948 × (√2·Nd + No)         |
| `kimura1`           | √(Nd² + (Nd + No)²)          |
| `kimura2` (default) | √(Nd² + (Nd + No/2)²) + No/2 |

Kimura2 is the default because it is the most stable across root
orientations. The Kimura estimators are **not additive**: the sum over
depth bins does not equal the whole-image value, which is why the
depth-profile wrapper uses a per-pixel Freeman length map instead.

**Diameter comes from a distance transform.** Every root pixel is
labelled with its distance to the nearest background pixel — the radius
of the largest circle that fits inside the root there. Diameter is twice
that. This is why
[`root_diameter()`](https://jcunow.github.io/Rootopia/reference/root_diameter.md)
needs the filled mask, not the skeleton.

**Unit conversion** is the same everywhere: `inch = px / dpi`,
`cm = px × 2.54 / dpi`.

**Under the hood**

Step counting uses four 3×3 focal kernels, one per direction pair. A
focal sum of exactly 2 means the centre pixel and one neighbour in that
direction are both root — that is, a connected step.

[`root_scape_metrics()`](https://jcunow.github.io/Rootopia/reference/root_scape_metrics.md)
is a thin wrapper around
[`landscapemetrics::calculate_lsm()`](https://r-spatialecology.github.io/landscapemetrics/reference/calculate_lsm.html),
which does the actual work; see that package for what each metric means.

**Edge cases**

| Situation | What happens |
|----|----|
| Image contains no roots | [`root_length()`](https://jcunow.github.io/Rootopia/reference/root_length.md) returns 0 without error |
| Non-binary image after preprocessing | Warning, computation continues |
| `unit = "cm"` or `"inch"` without a valid `dpi` | Stops |
| Skeleton passed to [`root_diameter()`](https://jcunow.github.io/Rootopia/reference/root_diameter.md) | Runs, but every diameter is ~1 px |
| Mask passed to [`root_length()`](https://jcunow.github.io/Rootopia/reference/root_length.md) | Set `skeletonize = TRUE` and it thins internally first |

------------------------------------------------------------------------

### Depth-resolved analysis

Placing those measurements on a real depth axis.

#### Depth mapping

##### What it is for

A minirhizotron tube goes into the soil at an angle, and the scanner
unrolls its curved surface into a flat rectangle. A pixel’s position in
that rectangle therefore does not equal its depth in the soil. This
module builds a map that gives every pixel its true soil depth in
centimetres, so traits can be reported per depth interval rather than
per pixel row.

##### The flow

``` r

dmap  <- create_depthmap(img, tilt = 45, dpi = 300,
                         tube_thicc = 7, start_soil = 0)
bins  <- binning(dmap, nn = 5)               # continuous cm -> 5 cm bins
slice <- depth_zoning(img, bins, depth = 20) # keep one bin
```

From there, two patterns cover everything downstream:

- **Whole profile at once** with
  [`terra::zonal()`](https://rspatial.github.io/terra/reference/zonal.html)
  — for traits that reduce to a per-zone sum or mean, such as pixel
  counts.
- **One slice at a time** with
  [`depth_zoning()`](https://jcunow.github.io/Rootopia/reference/depth_zoning.md)
  — for traits whose function needs a whole image, such as root length.
  Compute for one depth, then loop.

##### The rules

**Two axes, two corrections.** Depth increases along the image **width**
(columns). The tube’s curvature varies along the image **height**
(rows). The map is built in the image’s own orientation, so the result
lines up with the input cell for cell — no transposing needed
downstream.

**Depth along the tube.** Each column is converted from pixels to
centimetres by `2.54 / dpi`, then multiplied by `sin(tilt)`, because a
tube inserted at an angle covers less vertical depth than its length.
Finally `start_soil` is subtracted so that zero sits at the soil surface
rather than at the top of the image.

**Curvature across the tube.** With `sinoid = TRUE`, a cosine wave
spanning one full tube circumference is added across the rows: pixels
imaged from the upper side of the tube are nearer the surface than those
from the lower side. Its amplitude is the tilted tube diameter, and
`center_offset` rotates the wave to match where the top of your tube
actually sits. With `sinoid = FALSE` this term is zero — the right
choice for flat rhizotron windows.

**Masking.** Any pixel marked `1` in the optional `mask` becomes `NA` in
the depth map, which is how tape and other foreign objects are excluded
from every depth-resolved statistic.

**Edge cases**

| Situation | What happens |
|----|----|
| `tilt` outside 0–90°, exclusive | Stops. A tube at 0° or 90° has no defined depth gradient here |
| `center_offset` outside 0–1 | Stops |
| Image 1×1 or smaller | Stops |
| Mask dimensions differ from image | Stops |
| Image taller than one tube circumference | Stops with “Wrong Tube Diameter” — the sine wave cannot cover the image, which means `tube_thicc` or `dpi` is wrong |
| No mask supplied | An all-zero mask is used, so nothing is excluded |

#### Root angle

##### What it is for

[`deep_drive()`](https://jcunow.github.io/Rootopia/reference/deep_drive.md)
asks a single question: of all the root pixels in an image, what
fraction are growing in the direction that would take them deeper
fastest? It returns one number between 0 and 1 — a compact measure of
how directly a root system dives rather than spreading sideways.

##### The flow

``` r

dd <- deep_drive(DepthMap = dmap, RootMap = skel)          # just the number
dd <- deep_drive(DepthMap = dmap, RootMap = skel,
                 return = "all")                            # plus the maps
```

You can supply your own `AngleMap` instead of `RootMap` if you have
measured growth directions by another route.

##### The rules

**Two direction maps are compared, pixel by pixel.**

1.  The **actual** direction each root pixel is heading. When you supply
    a `RootMap`, this is derived by running a standard D8 flow-direction
    algorithm (`terra::terrain(v = "flowdir")`) over the depth map
    restricted to root pixels — water flowing downhill is the same
    problem as a root heading deeper. The eight D8 codes are relabelled
    as compass bearings, with 0° pointing up the image and angles
    increasing clockwise.

2.  The **optimal** direction: for each pixel, which of its eight
    neighbours lies deepest. This is computed directly from the depth
    map, and diagonal neighbours are divided by √2 so that a diagonal
    step is not unfairly favoured over an orthogonal one.

**The score** is simply the count of pixels where the two agree exactly,
divided by the number of root pixels with a defined direction. There is
no partial credit: a pixel heading 45° away from optimal counts the same
as one heading 180° away.

**Edge cases**

| Situation | What happens |
|----|----|
| Neither `AngleMap` nor `RootMap` supplied | Stops |
| Depth map is entirely `NA` | Stops |
| No root pixels with a defined direction | Warning, returns `NA` |
| D8 code 0 (a pixel with no downhill neighbour) | Becomes `NA` and is excluded from both counts |
| Depth values negative | The absolute value is used, so sign conventions do not matter |

#### Distribution indices

##### What it is for

Summarising *where* roots are rather than how many there are: how deep
the centre of mass sits, whether two profiles differ, whether there is a
repeating pattern around the tube.

##### The flow

These are ordinary numeric functions — they take vectors, not images, so
they come after a depth profile has been built.

``` r

MRD(w = depth, roots = rootlength.density)          # mean rooting depth
compare_depth_distribution(P, Q, metric = "wasserstein")
root_accumulation(df, group = "Tube", depth = "depth",
                  variable = "rootlength", stdrz = "relative")
rhythmicity(x = slice_index, y = trait, fix_period = 48)
circular_mean(angles)
modal_peaks(diameters)
```

##### The rules

**Mean rooting depth** is the depth average weighted by root amount:
`sum(depth × roots) / sum(roots)`. One number, in the unit of `depth`.

**Comparing two profiles.**
[`compare_depth_distribution()`](https://jcunow.github.io/Rootopia/reference/compare_depth_distribution.md)
measures the distance between profiles `P` and `Q` by one of three
metrics: Wasserstein (earth-mover), Jensen–Shannon, or Kullback–Leibler.
Identical profiles give zero. Optional tail weighting lets deep layers
count for more than shallow ones, which matters when the interesting
differences are at depth but most root mass is near the surface.

**Accumulation** converts a per-bin variable into a cumulative profile
down the depth axis, per group. `stdrz` chooses the scale: raw counts,
additive, or relative to each group’s total.

**Rhythmicity** fits
`y = amplitude × sin(2π/period × (x + phase)) + offset` and tests
whether the amplitude differs from zero — an F-test comparing the fitted
model against a flat line, or a likelihood-ratio variant. It reports
amplitude, phase, period, peak position, R² and a p-value.

The function is agnostic about what `x` means. Applied to slice index
with the period fixed to the number of slices, it tests whether roots
are distributed evenly around the tube circumference. That is a
*spatial* use of a tool normally introduced for time series.

**Circular mean** averages angles correctly. The ordinary mean of 350°
and 10° is 180°, which is exactly wrong; the circular mean gives 0°.

**Modal peaks** finds local maxima in a kernel density estimate above a
prominence threshold — used mainly to ask whether a diameter
distribution has one mode or several. Optionally refines the result with
model-based clustering via `mclust`.

**Edge cases**

| Situation | What happens |
|----|----|
| [`MRD()`](https://jcunow.github.io/Rootopia/reference/MRD.md) given unequal-length vectors, any `NA`, negative roots, or all-zero roots | Stops in each case |
| [`circular_mean()`](https://jcunow.github.io/Rootopia/reference/circular_mean.md) given `NA` or infinite values | Warns and drops them; stops if nothing remains |
| [`root_accumulation()`](https://jcunow.github.io/Rootopia/reference/root_accumulation.md) given a missing column name | Stops, naming the column |
| [`rhythmicity()`](https://jcunow.github.io/Rootopia/reference/rhythmicity.md) with `fix_period = NULL` | Estimates the period too, costing one degree of freedom |
| Tail weighting off | The weighting argument is ignored and set to constant |

------------------------------------------------------------------------

### Specialised

Questions that need a second image or the colour channels.

#### Turnover

##### What it is for

How much root was produced, how much died, and how much stayed, between
two points in time.

##### The flow

Two methods, for two kinds of input:

``` r

# two separate sessions
root_turnover(img1, img2, method = "tc", tc_method = "kimura",
              unit = "cm", dpi = 300)

# one multi-layer RootDetector image
root_turnover(dpc_img, method = "dpc")
```

##### The rules

**Temporal comparison (`"tc"`)** measures root amount in each image and
takes the difference. `tc_method` chooses the measure: `"kimura"` (root
length) or `"rootpx"` (pixel count). It reports standing root at each
time point, production as the difference, and new root as a percentage
of each time point. This method sees only the *net* change — root that
grew and died between the two scans is invisible.

**Decay–production–constant (`"dpc"`)** decomposes a single three-layer
image in which a segmentation tool has already colour-coded each pixel
as produced, decayed, or unchanged. This resolves gross change rather
than net.

The three layers are separated by thresholding each against the
unchanged layer at `blur_capture` (default 0.95) of its maximum, which
tolerates the soft edges that image compression leaves. A tape layer is
identified and removed first.

**Two ratio conventions**, and the choice changes the numbers:

- `include_virtualroots = FALSE` (default) — new growth is production
  divided by production plus constant; decay is decay divided by decay
  plus constant.
- `include_virtualroots = TRUE` — both are divided by the total of all
  three, so they share one denominator and count roots present at *any*
  time point.

**Edge cases**

| Situation | What happens |
|----|----|
| Fewer than 3 layers for `"dpc"` | Stops |
| `product_layer` equals `decay_layer` | Stops |
| Any layer entirely `NA` | Warning, returns `NULL` |
| All three classes empty | Warning, ratios returned as `NA` |
| `im_return = TRUE` | Returns the four classified rasters instead of the numbers, for visual checking |

#### Soil and colour

##### What it is for

Everything that uses colour rather than shape: classifying what each
pixel is made of, summarising tube colour, quantifying surface texture,
and building a rhizosphere zone around roots.

##### The flow

``` r

result <- classify_soil_rgb(rgb_img)                 # per-pixel class
plot_soil_classification(result)
result$metrics                                       # per-class statistics

cents <- build_soil_centroids(picks, max_dist)       # calibrate to your scanner
tube_coloration(rgb_img)                             # whole-image colour summary
analyze_soil_texture(rgb_img)                        # GLCM texture
halo <- create_root_buffer(seg, width = 3, halo_only = TRUE)
```

##### The rules

**Classification is nearest-centroid in CIE LAB.** Each pixel is
converted from RGB to LAB — a colour space where Euclidean distance
roughly matches perceived difference — and assigned to the closest class
centroid.

**Each class has its own distance limit.** A pixel further than
`MAX_DIST` from *every* centroid is left `"unclassified"` rather than
forced into the nearest class. The limits are per class because some
materials are tighter in colour than others.

**Every material in your scan needs its own class.** Anything without
one is snapped into whichever class happens to sit closest. The shipped
centroids were calibrated on one scanner at one site;
[`build_soil_centroids()`](https://jcunow.github.io/Rootopia/reference/build_soil_centroids.md)
derives your own from colour picks, and can blend them with the existing
ones via `alpha` if you are refining progressively.

**The halo is dilation.**
[`create_root_buffer()`](https://jcunow.github.io/Rootopia/reference/create_root_buffer.md)
grows the root mask outward by `width` iterations using either an
8-neighbour (`"circle"`) or 4-neighbour (`"diamond"`) kernel. With
`halo_only = TRUE` the roots themselves are subtracted, leaving only the
ring.

**Under the hood**

RGB→LAB uses the standard D65 white point. Assignment is vectorised over
all pixels at once via the expansion ‖p − c‖² = ‖p‖² − 2p·c + ‖c‖²,
which is why a full-resolution scan classifies in seconds;
`downsample_fact` speeds it up further for previews.

Texture metrics come from
[`glcm::glcm()`](https://rdrr.io/pkg/glcm/man/glcm.html) (grey-level
co-occurrence matrix) — a standard method, documented in that package.

[`tube_coloration()`](https://jcunow.github.io/Rootopia/reference/tube_coloration.md)
uses the exact Rec. 709 luma coefficients (0.2126/0.7152/0.0722). Note
that
[`rgb2gray()`](https://jcunow.github.io/Rootopia/reference/rgb2gray.md)
uses rounded ones (0.21/0.72/0.07), so the two do not produce identical
greyscale values; this is deliberate and the difference is documented at
both sites in the source.

**Edge cases**

| Situation | What happens |
|----|----|
| Pixel beyond every `MAX_DIST` | Class 0, reported as unclassified. Watch this fraction — a high value means your centroids do not match your images |
| Two centroids very close together | Assignment between them becomes unstable. The shipped defaults note that `dark_soil` and `red_soil` sit about 10 LAB units apart |
| Input is not 3-band | [`tube_coloration()`](https://jcunow.github.io/Rootopia/reference/tube_coloration.md) stops |
| RGB weights not summing to 1 | [`tube_coloration()`](https://jcunow.github.io/Rootopia/reference/tube_coloration.md) stops; [`rgb2gray()`](https://jcunow.github.io/Rootopia/reference/rgb2gray.md) only warns |

------------------------------------------------------------------------

### How the branching numbers are checked

On a real scan, nothing is known in advance, so a wrong length looks
exactly like a right one. The test suite therefore draws synthetic root
images whose length, width and topology are prescribed, runs the whole
branching pipeline over them, and scores the output against what was
drawn.

Five designs are used, each stressing something different: laterals on
one axis, 45° laterals, three generations of nested branching, two roots
crossing without touching, and a symmetric fork. Each is scored twice.
Feeding the exact one-pixel centre line isolates the graph — tracing,
junction handling, crossing resolution, ordering, length integration.
Feeding the filled image instead also carries the thinning error,
chiefly the erosion of about one root radius at every tip, so it gets
more room: 6% on length against 3% for the centre-line route. Tip,
branch-point and root counts must be exact. Diameters allow 10%, scored
against the package’s own `2 × EDT` convention rather than a formula.

The fork is scored on totals only. Its two arms are equally straight and
equally thick, so which one continues the parent is genuinely undefined,
and a per-root answer would be arbitrary rather than wrong.

There is nothing here for you to call: `devtools::test()` re-runs all of
it, and a change that shifts a branching number fails a test instead of
quietly landing in your results.

------------------------------------------------------------------------

### Where the numbers can surprise you

A short index of the rules most likely to explain an unexpected result.

| Symptom | Rule responsible | Module |
|----|----|----|
| Diameters are all about 1 px | A skeleton was passed where a mask was needed | Branching order, Size traits |
| Depth bins sum to less than the whole image | Kimura length estimators are not additive | Size traits |
| `class_map` has fewer pixels than the skeleton | Merged junction interiors belong to no segment | Branching order |
| A symmetric fork is ordered arbitrarily | The continuation score genuinely ties | Branching order |
| Tip and branch-point counts disagree with the graph | A four-way crossing is one branch point, not two | Skeletonisation |
| Rotation crop is not the width requested | `fixed_width` did not fit and was clamped | Rotation bias |
| Many pixels come back unclassified | Every centroid was beyond its `MAX_DIST` | Soil and colour |
| Turnover ratios differ between runs | `include_virtualroots` changes the denominator | Turnover |
| [`detect_skeleton_points()`](https://jcunow.github.io/Rootopia/reference/detect_skeleton_points.md) output will not overlay | Its rasters are rebuilt without the input’s extent | Skeletonisation |
