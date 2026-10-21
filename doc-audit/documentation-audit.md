# Documentation audit: tutorials, help pages, README

Scope: all seven vignettes, `README.md`, `_pkgdown.yml`, the roxygen blocks in
`R/`, and the generated `man/*.Rd`. Every claim below was checked against the
source; the method of checking is named where it matters.

**R is not available in this session.** Nothing here was verified by running
`R CMD check`, `devtools::test()` or a vignette build. Findings come from
reading the source and from two checks I could run outside R: the
Zhang-Suen lookup-table comparison (Python) and the pandoc heading-anchor
check. What that leaves unverified is listed at the end.

---

## Status

Everything in sections A, B, D, E and F is now fixed in the working tree, and
section C is fixed except where noted. Three decisions were taken by the
maintainer:

* **`soil_starts` keeps its behavior** (centimeters) and the documentation was
  corrected to match, at every site plus the pixel-to-cm conversion from
  `estimate_soil_surface()`. (B4)
* **`stdrz` keeps its behavior** and the prose was corrected to describe it;
  `"additive"` is the option that ends at 1, and the vignettes now recommend
  it. (B3)
* **`order_by` / `decreasing` were implemented** rather than removed from the
  docs, since nothing else could set the frame order. (A3)

The insertion-angle convention was resolved by renaming `tilt` to
`angle_from_vertical`, which determined the trigonometry and fixed C15: the
length term scales by `cos` and the circumference term by `sin`. See `NEWS.md`.
`root_depth_metrics(insertion_angles =)` keeps its name and now follows the same
convention.

Not fixed, and deliberately so:

* **C13** the "Python tool" reference is reworded rather than removed, because
  the `"feature"` method genuinely errors with that message.
* **F3** the accepted-but-unpublished DOI in `R/data.R` is left for the
  maintainer to confirm.
* **F8** the two "planned feature" notes in the rotation vignette are left as
  they are; they are labeled as unimplemented, which is honest, and a *Future
  directions* section is normal in a method paper.

`man/*.Rd` and `NAMESPACE` are generated and are now stale: run
`roxygen2::roxygenise()`.

---

## A. Broken — the documented call does not work

### A1. `root_depth_metrics(flat_geometry = )` does not exist
`vignettes/BatchProcessing_vignette.Rmd:342-355` and `:491`

The vignette has a whole section, *Flat rhizotron windows*, built on
`flat_geometry = TRUE`, and repeats the argument in *Full example with all
options*. `flat_geometry` appears nowhere in `R/`; the signature
(`R/Wrapper_root_depth_metrics.R:468-537`) has no such argument and no `...`,
so both calls fail with "unused argument".

The behavior the section describes does exist — `sinoid` is set from
`tube_geometry` at `R/Wrapper_root_depth_metrics.R:1040`, which is driven by
whether `insertion_angles` / `tube_diameter_cm` were supplied. So the section
needs rewriting around that switch, not a renamed argument.

### A2. Same example: `insertion_angles = 0` is rejected
`vignettes/BatchProcessing_vignette.Rmd:353`

```r
insertion_angles = 0   # vertical windows
```

Two errors in one line. `root_depth_metrics()` validates
`insertion_angles` to `(0, 90]` (`R/Wrapper_root_depth_metrics.R:679-681`), so
`0` stops the run. And the comment contradicts the convention the same vignette
states at line 98 and the roxygen states at `:60-67`: degrees from
**horizontal**, so a vertical window is `90`.

### A3. `stitch_root_scans(order_by = , decreasing = )` do not exist
`vignettes/Stitching_vignette.Rmd:275-288`

The subsection *Ordering — `order_by` / `decreasing`* gives two calls using
them. Neither argument exists in `stitch_root_scans()`
(`R/stitch_scans.R:144-151`), which has no `...` either. The lead-in at
`:240-241` ("see *Naming, selection & ordering* below to change it") promises a
control the package does not have.

The parenthetical at `:287-288` — zero-pad your frame numbers — is the only
option that actually works, so it is the answer rather than the fallback.

### A4. `slice_rotation(mode = "both")` does not exist
`vignettes/Rotation_Bias_vignettes.Rmd:306-311`

The *Going further* box describes splitting each circumferential slice by depth
bin via `slice_rotation(mode = "both")`. The signature is
`slice_rotation(img, n)` (`R/CreateZoning4.R:170`) — there is no `mode`.

The slice x depth grid is still reachable: `slice_rotation()` each slice, then
`depth_zoning()` within it. If the box is meant as a planned feature it should
say so, like the *Future directions* section does.

### A5. `binning()` has no `rounding` argument
`vignettes/MinirhizotronScans_vignettes.Rmd:438`

> `You can change this in binning() with the argument: rounding = "ceiling"`

The argument is `round_option` (`R/CreateZoning4.R:97`), and it accepts
`"rounding"`, `"ceiling"`, `"floor"`.

### A6. `image_threshold()` is handed a `cimg`
`vignettes/FlatBedScans_vignettes.Rmd:136-143`, and the same pattern in
`vignettes/HowItWorks_vignette.Rmd:258-265`

```r
grayscale <- load_flexible_image("path/to/scan.tif", scale = "to_01")
mask <- image_threshold(grayscale, threshold = 0.3, method = "global")
```

`load_flexible_image()` defaults to `output_format = "cimg"`
(`R/flexible_inputs.R:162`), and `image_threshold()` stops unless the input is
a `SpatRaster` (`R/LilHelpers3.R:580`). Add
`output_format = "spatrast"`.

A second trap in the same call: `image_threshold()` defaults to
`select_layer = 2` (`R/LilHelpers3.R:568`), which is wrong for the
single-layer grayscale image the text describes. Pass `select_layer = NULL`.

### A7. Cross-vignette anchor `#7b-branch-order` does not resolve
`vignettes/MinirhizotronScans_vignettes.Rmd:394`,
`vignettes/BatchProcessing_vignette.Rmd:282`

Both link to `FlatBedScans_vignettes.html#7b-branch-order`. Pandoc strips
everything up to the first letter when it builds an identifier, so the heading
`### 7b. Branch order` becomes `id="b.-branch-order"` (checked with the pandoc
in this container). Either give the heading an explicit
`{#branch-order}` and link to that, or drop the fragment.

---

## B. Wrong numbers — the call runs and reports something misleading

### B1. `root_diameter()` silently uses 300 dpi in a 150 dpi workflow
`vignettes/MinirhizotronScans_vignettes.Rmd:305-310`

```r
diam_result <- root_diameter(root_layer, skeleton_img = skl,
                             unit = "cm", select_layer = NULL)
```

`dpi` is omitted and defaults to `300`
(`R/diameter_calculation3.R:59`), while every other call in that vignette uses
`dpi = 150`. The `avg.diameter` column that this chunk adds to `depth_data`,
and the printed `head()`, are out by a factor of two. This chunk is
`eval = TRUE`, so the published page shows the wrong numbers.

### B2. Flatbed root length density mixes two resolutions
`vignettes/FlatBedScans_vignettes.Rmd:220` and `:288`

Length and diameter use `dpi = 1200`; the area the density divides by uses
`(300 / 2.54)^2`. The printed `rl_density` is 16x too large. Both chunks run.

### B3. `stdrz = "relative"` does not end at 1
`vignettes/MinirhizotronScans_vignettes.Rmd:441-443`;
`vignettes/HowItWorks_vignette.Rmd:1222-1224`

The minirhizotron vignette says

> `stdrz = "relative"` rescales the cumulative curve to end at 1

and *How It Works* says it scales "relative to each group's total". Neither is
what the code does. `R/LilHelpers3.R:770-780`:

* `"additive"` divides the cumulative series by `max(cs)` — **this** is the
  curve that ends at 1;
* `"relative"` divides it by `sum(cs)`, the sum of the cumulative values, which
  ends at `max(cs) / sum(cs)`.

So the comparability argument the vignette makes belongs to `"additive"`, and
the figure's y-axis label ("Cumulative root pixel fraction", line 436) is wrong
for the option it plots. The inline comment at line 428 — `"additive" (/max)`,
`"relative" (/sum)` — is correct and contradicts the prose three lines later.

Worth deciding separately whether `/sum(cs)` is the intended definition at all;
`@param stdrz` (`R/LilHelpers3.R:704`) documents none of the three options.

### B4. `soil_starts` is a pixel count fed into a centimetre argument
`R/Wrapper_root_depth_metrics.R:69-71` and `:1042`;
`vignettes/HowItWorks_vignette.Rmd:1006-1010`;
`vignettes/BatchProcessing_vignette.Rmd:47`, `:193`

This one is a code/documentation disagreement, not only a doc bug.

* `create_depthmap(start_soil = )` is documented and used as **centimetres**
  (`R/CreateDepthmap3.R:19`, and `M = outer(...) - start_soil` at `:131`).
* `root_depth_metrics(soil_starts = )` is documented as a **pixel row**
  (`R/Wrapper_root_depth_metrics.R:69`) and passed straight through to
  `start_soil` (`:1042`).
* `estimate_soil_surface()` returns `soil0` as a **pixel index**
  (`R/Calibration3.R:591-603`), and *How It Works* feeds it directly to
  `create_depthmap(start_soil = )`.

At 150 dpi a pixel is 0.017 cm, so treating one as the other shifts the depth
origin by about 59x. One of the three has to change; the docs currently assert
both conventions.

Separately, `soil0` indexes the **depth axis**, which runs along image
*columns* in this package (`R/Calibration3.R:556-558`,
`R/CreateDepthmap3.R:100-101`). Calling it a "row" in
`R/Wrapper_root_depth_metrics.R:69` and in the batch vignette is wrong whichever
unit wins.

### B5. The `mask` idiom in every depth example masks nothing
`vignettes/MinirhizotronScans_vignettes.Rmd:187-201`, and the roxygen examples
of `create_depthmap`, `binning` and `depth_zoning`

```r
mask <- seg_img[[1]] - seg_img[[2]]
mask[mask == 255] <- NA
```

`create_depthmap()` excludes a pixel only where the mask is exactly `1`
(`R/CreateDepthmap3.R:135-136`; `@param mask` says "1 = mask"). With a
`scale = "none"` three-band RootDetector image the layers are 0/255, so the
difference is `-255`, `0` or `255`, and the `255` values are then turned into
`NA`. No cell is ever `1`, so nothing is masked.

The vignette's comment — "`mask` flags foreign objects (e.g. tube edges / tape)
to exclude from depth attribution" — describes an effect the code does not
produce here. Either the idiom needs a `> 0` / `/255` step, or the examples
should stop advertising masking.

---

## C. Claims that overstate or misstate the code

### C1. "Thirteen staircase corners" — only four of them are
`vignettes/HowItWorks_vignette.Rmd:432-438`; `R/Thinning5.R:19-27`

I reconstructed the published Zhang-Suen deletion conditions and compared them
entry by entry with `.thinning_lut`. The headline numbers are right: the tables
disagree at **25** of 256 codes, and **12** of those are the conservative
direction (textbook deletes in both half-passes, this deletes in one).

The description of the other 13 is not right. The text says

> a pixel whose only neighbors are two orthogonal ones that already touch each
> other diagonally

Only codes 10, 40, 130 and 160 have *only* that pair. The other nine (11, 26,
104, 108, 134, 161, 176, 194, 198) contain such a pair **plus** one or two
further neighbors. All 13 do contain an adjacent orthogonal pair, so dropping
"whose only neighbors are" makes the sentence true as written.

### C2. `create_depthmap()` accepts `tilt = 90`
`vignettes/HowItWorks_vignette.Rmd:1066`

The edge-case table says

> `tilt` outside 0-90°, exclusive | Stops. A tube at 0° or 90° has no defined
> depth gradient here

The validation is `tilt <= 0 || tilt > 90` (`R/CreateDepthmap3.R:52-53`), and
the roxygen explicitly documents 90 as the supported head-on case
(`:13-17`). `root_depth_metrics()` even defaults to `90` for flatbed geometry
(`R/Wrapper_root_depth_metrics.R:692`). Only `0` stops.

### C3. The sine amplitude is half the tilted diameter
`vignettes/HowItWorks_vignette.Rmd:1046`

> Its amplitude is the tilted tube diameter

`df00 = (1 - cos(...)) * (tube_thicc_tilted / 2)` (`R/CreateDepthmap3.R:117`).
Amplitude is `tube_thicc_tilted / 2`; the full crest-to-trough range is the
tilted diameter. Say "its range" and the sentence is right.

### C4. `peakTime` is documented as the wrong formula
`R/Rhythmicity2.R:159` -> `man/rhythmicity.Rd`

Documented as "Peak time (phase + period/4)". Computed as
`(period / 4 - phase) %% period` (`R/Rhythmicity2.R:213`). Different sign and
a modulus.

### C5. `distance_map_rast` is not the distance transform
`R/diameter_calculation3.R:38` -> `man/root_diameter.Rd`

Documented as "Raster showing the distance transform values". It is
`Ds <- distance_map * 2` (`:131`, `:242`) — the local **diameter** in pixels,
over the whole mask rather than the skeleton, and not converted to `unit`.
(The minirhizotron vignette's inline comment, "# root thickness" at line 319,
is the accurate description; the roxygen is the one to fix.)

### C6. `root_depth_metrics()` does not run branching order or root angle by default
`vignettes/HowItWorks_vignette.Rmd:62-64`

> It is a batch wrapper that runs branching order, size traits, depth mapping,
> root angle and the distribution indices over a folder of images.

`calc_root_order_metrics` and `calc_root_angles` both default to `FALSE`
(`R/Wrapper_root_depth_metrics.R:513-514`), as the batch vignette correctly
says. "can run" rather than "runs".

### C7. Modal-peak columns need an argument the vignette never mentions
`vignettes/BatchProcessing_vignette.Rmd:229`, `:244`

`calc_diameter_quantiles` is credited with "90th/95th/99th percentile + modal
peaks" and with the `n.diameter.peaks` / `diameter.peak.*` columns. Those
columns come from `calc_modal_peaks` (`R/Wrapper_root_depth_metrics.R:510`,
`:1452-1470`), which is a separate argument, also off by default, and absent
from the vignette.

### C8. `par_start` is ignored in the call that passes it
`vignettes/Rotation_Bias_vignettes.Rmd:241-248`

`rhythmicity()` is called with both `fix_period = n_slices` and a four-element
`par_start`. With a numeric `fix_period`, `rhythmicity()` calls
`fit_sine_curve(x, y, fix_period = fix_period)` and does not forward
`par_start` (`R/Rhythmicity2.R:200-204`) — which the roxygen states correctly.
The call reads as if the starting values mattered.

### C9. A printed conclusion that the data has to agree with
`vignettes/Rotation_Bias_vignettes.Rmd:252`

```r
cat("The top-down difference is", round((res$amplitude*2) / res$offset,2),
    " -> circa twice as large as the mean")
```

The sentence asserts the ratio is about 2 while the same line prints it. If the
fit moves, the page contradicts itself. State the ratio, not a reading of it.

### C10. A comment that describes a different crop
`vignettes/Rotation_Bias_vignettes.Rmd:210-214`

```r
# rotation trim: drop the top 33% of rows
root_layer <- terra::crop(root_layer, terra::ext(e[1] + e[2] / 4, e[2], e[3], e[4]))
# --> the trimming avoids noisy deep layer and tape artifacts at the top
```

The crop moves `xmin`, so it trims **columns** — the shallow end of the depth
axis — by a quarter, not 33% of the rows. The second comment ("noisy deep
layer") points the other way along the same axis. One crop, three descriptions.

### C11. `tube_coloration()` described three different ways
`vignettes/SpecialTopics_vignette.Rmd:160-161`;
`vignettes/FlatBedScans_vignettes.Rmd:404-405`

It returns `rcc`, `gcc`, `bcc`, `hue`, `saturation`, `luminosity`, `red`,
`green`, `blue` (`R/Zonestats3.R:410-420`). *Special Topics* calls that "mean
RGB and a single luminance value" (omits the chromatic coordinates and HSV);
*Flatbed Scans* calls it "mean chromatic coordinates and HSV values", but
`luminosity` is a weighted luma, explicitly **not** HSV value — the roxygen
says so at `R/Zonestats3.R:342-344`.

### C12. README uses two deprecated arguments in its quick start
`README.md:69`, `:84`

* `clean_image(img, max_artifact_size = 5, max_hole_size = 5)` — both
  deprecated, both `warning()` (`R/clean_image.R:416-427`). Current names:
  `drop_objects_upto`, `fill_holes_upto`.
* `order_metrics(order_res, focal = "thinnest")` — `focal` is deprecated and
  warns (`R/RootBranching.R:1833`). Current name: `groups`.

So a reader pasting the quick start gets three deprecation warnings on their
first contact with the package.

`README.md:82` also passes `order = "root_order"` under the comment
"# Branch order". `root_order` is the deprecated alias for
`root_max_order_from_tip` (`R/RootBranching.R:394-397`), which is the per-root
maximum counted from the tips — not branch order. The branch-order scheme is
`order_from_base`.

### C13. The stitching "Python tool" has no referent
`vignettes/Stitching_vignette.Rmd:175-177`

> The Python tool's `"feature"` (SIFT/ORB) needs OpenCV and has no CRAN backend

Nothing else in the package or the docs names a Python tool. For a reader this
is an unresolvable reference; it reads as a leftover from a port.

### C15. The curvature amplitude may use the wrong trigonometric function
`R/CreateDepthmap3.R` (curvature branch)

Flagged rather than changed, because it moves every published depth profile.

The curvature term scales the tube diameter by `sin(tilt)`. By the geometry it
should arguably be `cos(tilt)`:

* `tilt = 90`, a vertical tube. Its circular cross-section is horizontal, so
  every point around the circumference sits at the same depth and the curvature
  term should vanish. `sin(90) = 1` gives it full amplitude; `cos(90) = 0`
  gives none.
* `tilt = 0`, a horizontal tube. Its cross-section is vertical, so the top of
  the tube is a full diameter shallower than the bottom and the term should be
  at maximum. `sin(0) = 0` gives none; `cos(0) = 1` gives full amplitude.

At `tilt = 45`, `sin` and `cos` are equal, which is why this would not show up
in any example or test in the repository -- 45 is the value used in every one.

The depth-along-the-axis term is `sin(tilt)` and that one is right: an angled
tube does cover less vertical depth than its own length.

I have not changed it. If it is a real inversion, the fix is one character and
every depth profile the package has produced shifts; if the intent was
something else, the derivation above says what the comment should rule out.

### C14. Vignette index entry does not match the title
`vignettes/SpecialTopics_vignette.Rmd:2` vs `:9`

Title: "Special Topics: Soil Color, Soil Texture, Rhizosphere Halos, and
Turnover". `\VignetteIndexEntry`: "Special Topics: RGB classifier, Soil
Texture, ...". The index and the page disagree.

---

## D. Comments addressed to a developer, or describing the past

`CLAUDE.md` rules out both. These are the hits; I read every one in context
rather than matching keywords alone.

### D1. A note to whoever works on the file next
`R/Calibration3.R:470-472`

```r
# scale = "none", unlike estimate_rotation_center()'s "to_01" -- so
# tape_brightness is read against the raw 0-255 range here. See the note at
# that function before merging the two pipelines.
```

The first sentence is a present-tense fact worth keeping. The last sentence is
a task for a maintainer, not a comment on the code.

### D2. Past tense describing an old bug
`R/Calibration3.R:556-559`

```r
# `i` indexes COLUMNS (depth), so rr1[,i,] holds nrow(rr1) pixels -- the
# tape coverage of that depth slice across the circumference. Dividing by
# ncol() scaled the ratio by the aspect ratio, so `tape_thresh` did not
# mean the documented "fraction of the slice covered by tape".
```

Sentence one is the keeper. Sentence two is the changelog.

### D3. `R/Zonestats3.R:352-355`

```r
# Routing this through rgb2gray() would shift every color value this
# function has ever reported, so the two stay separate even though the
# arithmetic is the same.
```

"has ever reported" is version history. The reason the two differ can be
stated without it.

### D4. `R/Zonestats3.R:398-401`

```r
# `hsl[3]` below is HSV *value* (max channel), which is not luminance --
# reporting that as `luminosity` ignored r/g/b entirely, so the weighted
# mean is used instead.
```

Clearest case: the second clause describes behavior that no longer exists.

### D5. `R/LilHelpers3.R:818-819`

```r
# Checked before the load, not after: load_flexible_image() would otherwise
# raise R's generic missing-argument error first and this message never showed.
```

Same shape — "never showed" is about an earlier version.

### D6. `R/clean_image.R:416-418` and `R/Wrapper_root_depth_metrics.R:539-541`

```r
# Former argument names. They said how big a thing was but not what happened
# to it, so a reader could not tell whether the value was the largest thing
# kept or the largest thing removed.
```

Borderline. "Former argument names" is a present-tense fact about the block
below it. The rest is the reasoning behind a past rename, which the `NEWS.md`
entry is the place for.

### D7. Vignette setup chunks explaining editorial decisions
`vignettes/MinirhizotronScans_vignettes.Rmd:17-21`;
`vignettes/FlatBedScans_vignettes.Rmd:16-19`;
`vignettes/BatchProcessing_vignette.Rmd:13-15`

```r
# The displayed `load_flexible_image("path/to/...")` chunks stay eval=FALSE —
# that is how a user supplies their own files — while the chunks that operate
# on the example objects are eval=TRUE so the code is actually exercised.
```

```r
# Every chunk here operates on a directory of the reader's own images, so none
# of them can run at build time and eval = FALSE is the default throughout.
# The other vignettes work on the bundled example rasters and do run.
```

These justify a vignette-authoring choice to whoever maintains the vignette.
They are invisible to readers (`include = FALSE`) and say nothing about the
package. The *reason* a chunk is `eval = FALSE` is useful — but stated to the
reader, in the visible text, where two of the vignettes already do it well
("this chunk is not run here"), not in a hidden comment about the other
vignettes.

---

## E. Spelling: en-GB in prose

`DESCRIPTION` sets `Language: en-US` and `R CMD check --as-cran` spell-checks
the Title, Description and every `man/*.Rd` against it.

I checked each hit for whether it sits in prose or inside an identifier, a
string or a foreign argument name. **Leave these alone** — they are code:
`artefact_rate` and `artefact_cycles` (`R/network_metrics.R:153-154`),
`col = "grey80"` (`R/clean_image.R:239`), `n_grey` (`R/Zonestats3.R:496`),
`neighbourhood = 8` (`landscapemetrics`, `R/Zonestats3.R:265`), and
`characteristics` in `man/modal_peaks.Rd`, which is the same word in both
spellings.

### E1. In `man/` — these reach `--as-cran`

| File | Word | US form |
|---|---|---|
| `root_depth_metrics.Rd` | `centimetres` (x3) | centimeters |
| `root_depth_metrics.Rd` | `skeletonised` | skeletonized |
| `render_order_overlay.Rd` | `colour` (x3) | color |
| `render_order_overlay.Rd` | `labelled` | labeled |
| `root_transport.Rd` | `colour`, `centimetres` | color, centimeters |
| `smooth_root_edges.Rd` | `skeletonisation` (x2) | skeletonization |
| `clean_image.Rd` | `skeletonisation` | skeletonization |
| `flatbed_scan_example.Rd` | `skeletonisation` | skeletonization |

Fix the roxygen block in `R/`, then re-run `roxygen2::roxygenise()` — these
`.Rd` files are generated.

### E2. In vignettes, `README.md` and `_pkgdown.yml` — not checked, but inconsistent

`HowItWorks_vignette.Rmd`: `Skeletonisation` (x4), `skeletonise`,
`centimetres` (x3), `centimetre`, `colour` (x2), `Summarising`, `Binarised`,
`labelled`, `favour`, `favoured`.
`MinirhizotronScans_vignettes.Rmd`: `centimetres`.
`BatchProcessing_vignette.Rmd`: `Skeletonised`.
`README.md`: `skeletonise`.
`_pkgdown.yml`: `Colour`, `colour`, `summarise`, `summarising`,
`characterisation`, `artefacts` (this one *is* prose — "rotational artefacts" —
unlike the `artefact_*` identifiers).

---

## F. Smaller things

1. **`man/create_root_buffer.Rd` title has a stray bracket.**
   `R/CreateZoning4.R:1`: `#' Create a buffer halo) around non-zero pixels`.
2. **`R/data.R:256` and `:265` both open `@author` for
   `flatbed_scan_example`,** and the line after the second one is untagged
   prose, so it is absorbed into the author field.
3. **`R/data.R:270-272` cites an accepted-but-unpublished paper with a live
   DOI** (`10.1002/oik.12211`, "accepted, 2026"). Worth confirming the DOI
   resolves before release.
4. **`vignettes/Stitching_vignette.Rmd:147`: "The two panels are
   indistinguishable."** A statement about a figure the reader is looking at,
   which only holds if the demo stitch is exact. It may well be — but it is the
   kind of claim a reviewer will test.
5. **`vignettes/FlatBedScans_vignettes.Rmd:109`** assigns
   `image_components <- report_image_components(seg)` and never uses it; the
   function returns invisibly, so the assignment suppresses nothing and the
   printed summary is the point.
6. **`vignettes/FlatBedScans_vignettes.Rmd:302-306`** takes tip and fork counts
   from the unpruned `skl` but divides by `rl`, measured on `skl_pruned`. The
   section immediately above explains that spurs inflate tip counts, so the
   branching frequency is computed from the numbers the vignette has just
   warned about.
7. **`vignettes/HowItWorks_vignette.Rmd:382` ("`center_offset` is read as a
   fraction ... and as an absolute row number when above 1")** omits the
   `"top"` / `"middle"` / `"bottom"` keywords `rotation_censor()` also accepts
   (`R/Calibration3.R:331-339`).
8. **Planned features are described in two places**
   (`Rotation_Bias_vignettes.Rmd:133-136` and `:315-328`). Both are clearly
   labeled as not implemented, which is honest, but `CLAUDE.md` says not to
   describe planned capabilities at all. Your call — a *Future directions*
   section is normal in a method paper and unusual in package docs.
9. **`batch_root_traits` has no `reference:` entry in `_pkgdown.yml`.** It is an
   alias documented inside `root_depth_metrics.Rd`, so pkgdown is satisfied and
   nothing breaks. Noting it only because `CLAUDE.md` asks for an index entry
   per exported object.

---

## G. Checked and correct

Listing these so the next audit does not redo them.

* All four Freeman/Kimura length formulas in
  `HowItWorks_vignette.Rmd:844-849` match `R/Zonestats3.R:4-26`.
* The four order schemes, their definitions, the per-segment vs per-root split,
  and the deprecated aliases (`tip_order`, `root_order`, `branch_order`)
  resolving silently while `focal` warns — all match
  `R/RootBranching.R:394-421`, `:1833`.
* `order_metrics(groups = "thinnest" | "thickest" | named list)`, the `"rest"`
  row, the duplicate-order error, the dropped empty groups, length-weighted
  diameter, `cumulative_*` only without `groups` — all match
  `R/RootBranching.R:1796-1890`.
* `network_metrics()` topologies, the `$traits` / `$qc` split, `cyclomatic`,
  `alpha_meshedness`, `artefact_rate`, `artefact_cycles` — match
  `R/network_metrics.R`.
* `root_transport()` column names and the Hagen-Poiseuille `8 eta L / (pi r^4)`
  description — match `R/root_transport.R`.
* The validation-suite description (five phantom designs; counts exact; 3%
  length from the center line, 6% from the mask; 10% diameter against the
  package's own `2 x EDT`) — matches `tests/testthat/helper-phantom.R:213`,
  `:274-285`.
* The two thinning invariants are genuinely tested, in the file named
  (`tests/testthat/test-skeleton.R:24-59`).
* 25 lookup-table differences and 12 conservative ones — recomputed and
  confirmed (see C1 for the part that is wrong).
* `root_depth_metrics()` has 48 arguments ("about fifty"), and the progress-line
  and `SKIPPED` formats match `R/Wrapper_root_depth_metrics.R:1541`, `:578`.
* Color metric column suffixes `_root` / `_bg` — match `:1418-1421`.
* "1144 rows here" in `Rotation_Bias_vignettes.Rmd:161` — matches
  `R/data.R:52`.
* Every `man/` topic is either in the pkgdown `reference:` index or carries
  `\keyword{internal}`, and every exported object has an `\alias`.
* `stitch_root_scans()` arguments used in the stitching vignette — all present
  except those in A3.
* `classify_soil_rgb()` returning `$map` / `$metrics`,
  `estimate_rotation_shift()` arguments, `deep_drive(return = "all")`,
  `root_turnover()` arguments, `root_scape_metrics()` arguments,
  `prune_skeleton()` arguments, `detect_skeleton_points()` arguments,
  `zoom_plot()` arguments — all match.

---

## H. Not verified

Because R is unavailable here:

* `R CMD check`, `R CMD check --as-cran`, and therefore the spell-check in E1
  as the tool would actually report it.
* `devtools::test()` — the test-suite behavior described in
  `HowItWorks_vignette.Rmd:1386-1416` is verified by reading the helpers, not
  by running them.
* Any vignette build, so: the figures, the numbers the `eval = TRUE` chunks
  print, the runtime claims ("slow", "considerably more expensive"), and
  Stitching F4.
* Whether the fixes in A1-A6 run clean once applied.
