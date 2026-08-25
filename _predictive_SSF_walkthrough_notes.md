# Predictive SSF Walkthrough — design notes

Working notes for `predictive_SSF_walkthrough.qmd`. This file records *why* the walkthrough is
built the way it is — the decisions, the equations, the provenance of each section, and the
things that were checked rather than assumed. The `.qmd` itself is the teaching artefact; this
is the record behind it.

---

## 1. What the walkthrough is for

A simplified, heavily-commented version of the pipeline from
[`swforrest/dynamic_SSF_sims`](https://github.com/swforrest/dynamic_SSF_sims), the code
accompanying Forrest et al. (2025, *Ecography*), "Predicting fine-scale distributions and emergent
spatiotemporal patterns from temporally dynamic step selection simulations".

That repo is six scripts. This is one, and it goes:

```
fit  ->  simulate  ->  aggregate  ->  hourly summaries  ->  assess predictions
```

The teaching payload is a **contrast between two models**:

| | |
|---|---|
| **Model A — static** | An ordinary iSSF. One selection coefficient per covariate, constant across the day. |
| **Model B — temporally dynamic** | A GAM where each coefficient is a cyclic smooth function of hour, `s(hour, by = x, bs = "cc")`. |

Both are simulated from, both are aggregated into utilisation distributions, and both are
validated against the observed data. The point a reader should leave with: a model whose
*coefficients* look more interesting also has to *predict* better, and you can check whether it
does.

---

## 2. Data and covariates

| | |
|---|---|
| Individual | Buffalo **2158** |
| Period | **2018-07-25 → 2018-10-31** (2018 dry season) |
| Covariates | **NDVI**, **canopy cover**, **slope** |
| Rasters | `mapping/ndvi_aug_2018.tif`, `canopy_cover.tif`, `slope_raster.tif` — all 25 m, EPSG:3112, identical extents |

**Why a single individual.** Keeps the script readable and fast, and makes validation
unambiguous: predicted use is compared against *that animal's* observed locations. Extending to
multiple individuals is a modelling discussion (pooling, random effects, two-step estimation)
that would swamp the pipeline this script is about.

**Why the dry-season subset.** Buffalo 2158 is not range-stationary across the full 15 months of
data — it shifts roughly 20 km east during 2019. A long-term utilisation distribution validated
against the whole track would score badly for reasons that have nothing to do with the model.
Restricting to Jul–Oct 2018 keeps the animal in one ~9.5 × 2.5 km band. It also matches the NDVI
layer, which is August 2018 only.

**Canopy cover is rescaled to 0–1** (`canopy_cover / 100`), as in the source repo, so its
coefficient is on a comparable scale to NDVI's.

### Measured baselines (sanity checks — these should reproduce)

- 2,235 observed locations, ~1 h fix interval
- Tentative gamma: **shape 0.4266, scale 786.8**
- Tentative von Mises: **κ 0.398**
- `n_control = 10` → 24,563 rows / 2,233 strata; `n_control = 25` → 58,058 rows / 2,233 strata
- Model A corrected movement parameters: **shape 0.4308141, scale 790.3658, κ 0.3924898**
  (hand-computed and `amt::update_sl_distr()` agree to 1.7e-9)
- Model B hourly parameters stay valid at every hour: shape 0.34–0.81, scale 219–1413,
  κ −0.39–1.30, mean step length **93 m to 1127 m** across the day

### Measured timings

| Step | Time |
|---|---|
| Model A (`fit_issf`) | 0.4 s |
| Model B (`cox.ph` GAM, 6 smooths, `n_control = 25`) | ~213 s |
| Manual simulation, 5 trajectories × 720 steps | 1.5 s (so ~2.5 min per 500-trajectory set) |
| `amt::simulate_path`, 10 × 500 steps | ~68 s |

`n_control = 25` roughly doubles the Model B fit time relative to 10 (94 s → 213 s) but is worth
it: the canopy smooth is only marginally significant at 10 (p = 0.17) and clearly significant at
25 (p = 0.017). All six smooths are significant at `n_control = 25`, with edf 5.4–8.6.

---

## 3. The parameter block

Everything scalable lives in one chunk at the top of the `.qmd`, so a reader can re-run the whole
thing bigger or smaller without hunting through it:

| Parameter | Default | What it does |
|---|---|---|
| `which_buffalo` | `2158` | Individual to fit |
| `season_start` / `season_end` | 2018-07-25 / 2018-11-01 | Fitting window |
| `n_control` | `25` | Available steps per used step |
| `n_traj` | `500` | Simulated trajectories per model |
| `n_steps` | `720` | Steps per trajectory (720 h = 30 days) |
| `n_ch` | `25` | Candidate steps evaluated at each simulated step |
| `burn_in` | `24` | Leading steps discarded from each trajectory |
| `start_mode` | `"random"` | Seeding strategy — see §6 |
| `boundary` | `"wrapped"` | Landscape boundary — see §6 |
| `ud_res` | `100` m | Resolution of the all-hours UD |
| `ud_res_hourly` | `250` m | Resolution of the 24 hourly UDs (coarser: each holds ~1/24 of the locations) |

---

## 4. Model formulations

### Model A — static iSSF

```r
amt::fit_issf(case_ ~ ndvi + canopy + slope + sl_ + log_sl_ + cos_ta_ + strata(step_id_), ...)
```

Also fitted as the equivalent `mgcv::gam(cbind(times, step_id_) ~ ..., family = cox.ph)` so that
both models share one interface for prediction, and so the reader can see the two are the same
model.

Movement terms are named `sl_`, `log_sl_`, `cos_ta_` — **amt's convention, with the trailing
underscore**. This matters: `amt::update_sl_distr()` and `update_ta_distr()` look for exactly
these names, and they're used in §5 to cross-check the hand-computed movement parameters.

### Model B — temporally dynamic GAM

```r
mgcv::gam(cbind(times, step_id_) ~
            s(hour, by = ndvi,    bs = "cc") +
            s(hour, by = canopy,  bs = "cc") +
            s(hour, by = slope,   bs = "cc") +
            s(hour, by = sl_,     bs = "cc") +
            s(hour, by = log_sl_, bs = "cc") +
            s(hour, by = cos_ta_, bs = "cc"),
          knots = list(hour = c(0, 24)),
          family = cox.ph, weight = case_, ...)
```

This is a **varying-coefficient model**: `s(hour, by = ndvi)` fits β_NDVI(hour), a coefficient
that is a smooth function of time of day. `bs = "cc"` (cyclic cubic regression spline) makes the
curve join up at midnight. It's the same idea as the harmonic terms in the paper —
β(τ) = α₀ + Σ αⱼ sin(2jπτ/T) + Σ αⱼ₊ₚ cos(2jπτ/T) — with a spline basis instead of a Fourier one,
following Klappstein et al. (2024).

Movement terms get the same treatment, which is what lets step length and directional persistence
vary across the day in the simulations.

### Two things that were checked, not assumed

**(a) Numeric `by=` smooths are NOT centred — no parametric main effect is needed.**

The obvious worry with `s(hour, by = ndvi)` and no `ndvi` term alongside it is that mgcv's
sum-to-zero identifiability constraint would force the average of β_NDVI(hour) across the day to
be zero, throwing away the overall selection strength. That worry is unfounded: mgcv skips the
centring constraint when the `by` variable is numeric. Checked two ways:

- Basis dimension. `smoothCon(s(hour, by = z, bs = "cc", k = 6))` keeps **5** columns after
  constraints; plain `s(hour, bs = "cc", k = 6)` keeps **4**. The constraint isn't applied.
- Full fit on this data, with and without a parametric `ndvi` term: β_NDVI(hour) spans
  −6.635 → 5.606 **either way**, mean −0.754 either way, edf 8.064 either way. The main effect is
  redundant and costs ~2 AIC (10553.45 → 10555.44) for the extra parameter.

So β_x(τ) is the smooth alone, with nothing added back. The `.qmd` keeps a short box on this,
because it's a natural thing for a reader to be uneasy about, and the basis-dimension check is a
satisfying one-line answer.

In the end the coefficients are extracted with `predict()` rather than `gratia::smooth_estimates()`,
via a small `beta_by_hour()` helper. Since neither model has an intercept, predicting the linear
predictor at a design where one covariate is 1 and all others are 0 returns that covariate's
coefficient directly. The same call works unchanged on the static model (returning a constant) and
the dynamic one (returning the smooth), and `se.fit = TRUE` gives confidence intervals — so one
function serves both models and the downstream code never needs to know which it is holding.

**(b) `knots = list(hour = c(0, 24))` is genuinely needed.**

Without it, mgcv places the cyclic endpoints at `range(hour)` — which is 0 and 23, not 0 and 24:

```
no knots argument   ->  knots at 0, 4.6, 9.2, 13.8, 18.4, 23
knots = c(0, 24)    ->  knots at 0, 4.8, 9.6, 14.4, 19.2, 24
```

The first wraps hour 23 onto hour 0 and squeezes an hour out of the daily cycle. Small, but real.
The same two-line fix has been applied to the temporal GAMs in `ssf_responses.qmd`.

---

## 5. From coefficients to a simulator

Both models are reduced to the **same object**: a 24-row table with columns
`hour, ndvi, canopy, slope, shape, scale, kappa`. Model A's is flat; Model B's cycles. Everything
downstream is model-agnostic, which is the whole reason the source repo used this format
(`TwoStep_2pDaily_coefs_dry_*.csv`).

Habitat coefficients come from the fitted model. Movement parameters use the iSSF update rules
(Avgar et al. 2016; Fieberg et al. 2021) — the tentative distributions fitted to the observed
steps are *corrected* by the fitted movement coefficients:

```
shape(τ) = shape_tentative + β_log_sl_(τ)
scale(τ) = 1 / ( 1/scale_tentative − β_sl_(τ) )
kappa(τ) = kappa_tentative + β_cos_ta_(τ)
```

For Model A these are checked against `amt::update_sl_distr()` / `update_ta_distr()`, which
implement exactly this. That cross-check is the highest-value assertion in the script: if the
update rules are wrong, every simulated step length is wrong and nothing downstream means
anything.

**Habitat selection surfaces** are then precomputed once per hour, on the log scale:

```
log w(x, τ) = β_ndvi(τ)·NDVI(x) + β_canopy(τ)·canopy(x) + β_slope(τ)·slope(x)
```

giving a 24-layer raster stack per model. The simulator indexes into it by hour rather than
re-evaluating the linear predictor at every candidate step.

---

## 6. Simulation

The sampler is the redistribution kernel written out longhand: **movement kernel × habitat
kernel**. At each step, draw `n_ch` candidate step lengths from `rgamma(shape(τ), scale(τ))` and
turning angles from a von Mises with κ(τ), project the candidate endpoints, look up `log w` in
the hour-τ layer, and choose one with `sample(n_ch, 1, prob = exp(log_w))`.

Two performance decisions, both of which matter because the loop runs ~360,000 times per
500-trajectory set:

- **All random numbers are drawn up front**, not inside the loop.
- **Raster lookup avoids `terra::extract()`.** The stack is digested once into a plain matrix via
  `terra::values()` (one row per cell, row-major, one column per layer), and lookups are done with
  cell arithmetic: `col = floor((x - xmin)/resx) + 1`, `row = floor((ymax - y)/resy) + 1`,
  `cell = (row - 1) * ncol + col`. Verified against `terra::extract()` on 1,000 random points
  (max absolute difference 0, NA handling identical) and against `cellFromXY()`. This is what
  brings a 500 × 720-step run down to ~2.5 minutes.

The source repo calls `terra::extract()` inside the loop; on the scale it was run at (an HPC
cluster) that was fine, but it is the main reason a naive port feels slow on a laptop.

### Two bugs found during the first render

**1. Turning angles were misaligned with hours (real correctness bug).**

The first version built the von Mises inputs per *candidate*:

```r
kappa_vec <- rep(coefs$kappa[step_hours], each = n_ch)   # length n_steps * n_ch
ta <- as.vector(mapply(Rfast::rvonmises, n = n_ch, m = mu_vec, k = abs(kappa_vec))) - pi
```

`mapply` then made `n_steps * n_ch` calls, each returning `n_ch` draws — 450,000 values instead
of 18,000, and block *i* of the result came from `kappa[step_hours[ceiling(i / n_ch)]]` rather
than `kappa[step_hours[i]]`. Step 2 was getting hour 1's concentration, step 30 hour 2's, and so
on. The fix is to pass a per-*step* vector, so `mapply` makes one call per step:

```r
kappa_step <- coefs$kappa[step_hours]                    # length n_steps
ta <- as.vector(mapply(Rfast::rvonmises, n = n_ch, m = mu_step, k = abs(kappa_step))) - pi
```

This is silent — no error, just wrong turning angles — so it is the more dangerous of the two.

**2. `lookup()` crashed on a non-finite coordinate — root cause: `Rfast::rvonmises`.**

`inside <- col >= 1 & col <= ncol & ...` evaluates to `NA` (not `FALSE`) when the coordinate is
`NaN`, and `NA` is not permitted as a subscript in `out[inside] <- ...`. This aborted the first
render at trajectory 53. Guarded with an explicit `!is.na(col) & !is.na(row) & ...`, which is the
correct semantics anyway: a non-finite coordinate is not inside the raster, so it gets the same
negligible weight as one that falls outside.

Tracked the `NaN` to its source by running trajectories until one failed and inspecting the state:
the failing step had exactly one non-finite **turning angle**, with all step lengths finite. Two
problems with `Rfast::rvonmises()`, both confirmed directly:

- **It occasionally returns `NaN`** — 6 non-finite in 4.8 × 10⁷ draws (~1 in 8 million), at small
  positive concentrations (κ ≈ 0.16, 0.25). Negligible-sounding, but a production run makes
  `n_traj × n_steps × n_ch` ≈ 9 million draws *per model*, so roughly one per run. That is exactly
  why it killed one render and then resisted reproduction across 6,000 trajectories.
- **It ignores `set.seed()`.** `Rfast` uses its own RNG: two calls under the same seed return
  different values. The simulations were therefore never reproducible, seed at the top or not.

The second problem is the more serious one for a teaching script, and it is inherited straight from
the source repo. **Switched to `circular::rvonmises()`**, which is reproducible under `set.seed()`,
produced zero non-finite values in 4.8 × 10⁷ draws with the same parameters, and is marginally
faster. Drawing is now done once per *hour* (24 calls) rather than once per step, scattering into
position with `outer()` — verified for length, hour alignment (correlation between mean |turning
angle| and κ across hours = −0.994) and reproducibility.

`circular` is called with `::` rather than attached, because attaching it masks `stats::sd` and
`stats::var`.

A related, harmless edge case in the wrapped boundary: `y %% y_extent` can return exactly 0, which
maps to `row = nrow + 1` and is therefore treated as outside. Measure-zero and non-fatal, so left
as is.

### Starting locations (`start_mode`)

Three strategies, because they answer different questions:

- **`"observed"`** — every trajectory starts at the animal's first observed location. Asks *where
  would this animal have gone from where it actually was?* The spread that emerges is
  home-range-like. The right choice for short simulations, and the one to use if the question is
  about space use rather than landscape-scale distribution.
- **`"observed_sample"`** — starts drawn at random from the observed locations. Seeded in habitat
  the animal actually used, but not all from one point, so the resulting UD isn't dominated by a
  single origin.
- **`"random"`** — uniform across the raster extent. Asks *what is the stationary distribution of
  this movement process over this landscape?* This is what the long-term UD and the validation
  need, and what the source repo used.

### Boundary

- **`"wrapped"`** (toroidal, via `%%`) — no edge effects, so the process has a well-defined
  stationary distribution. Used for everything quantitative. Requires the raster origin at (0,0),
  restored afterwards.
- **`"reflective"`** — candidate steps outside the extent get near-zero weight. Used only for the
  illustrative trajectory plot, because paths don't jump across the frame.

Burn-in interacts with both: seeding from a single observed point needs a real burn-in before the
trajectory forgets where it started; seeding uniformly at random needs much less.

---

## 7. Aggregation, summaries, validation

- **Aggregation** — `terra::rasterize(coords, template, fun = "sum")`, all-hours at 100 m and per
  hour at 250 m, normalised to sum to 1 so the two models are comparable. Plus a convergence
  check built from subsets of the already-simulated trajectories (25/50/100/250/500), correlated
  against the full result — this is how a reader decides whether `n_traj` was enough for *their*
  system.
- **Hourly summaries** — hourly means of step length and each covariate, simulated vs observed.
  With one individual there's no between-animal spread to plot, so the observed track is split
  into weekly blocks to give a comparable distribution of observed hourly curves.
- **Validation** — continuous Boyce index (`ecospat::ecospat.boyce`), Spearman and Pearson,
  overall and separately for each of the 24 hours. The hourly panel is the miniature version of
  the key result in Forrest et al. (2025).

---

## 8. The `amt` route

Near the end, the same thing done with the package: `amt::redistribution_kernel()` +
`simulate_path()` (Signer et al. 2023), shown for Model A.

**Two things that are needed to make this work**, neither of them obvious from the resulting error
messages:

1. The default extraction function uses `where = "both"`, producing columns named `ndvi_start` and
   `ndvi_end`. A model referring to plain `ndvi` then fails with
   `object 'ndvi' not found`. A custom `fun` with `where = "end"` fixes it — but must *also*
   compute `log_sl_` and `cos_ta_`, since the kernel supplies only `sl_` and `ta_` (otherwise the
   next error is `object 'log_sl_' not found`).
2. The kernel has no toroidal option and terminates a path once too many proposed steps fall
   outside the map. On the 14 × 14 km simulation crop, paths stop after ~150 of 500 steps. The
   walkthrough therefore hands `amt` the *full* covariate rasters, and uses the opportunity to
   point out that this is exactly the edge effect the wrapped boundary in our own sampler avoids.

The deeper limitation, and why the manual sampler exists here: the kernel is built from one fitted
model against a fixed `map`, so coefficients that change with hour need either a custom extraction
function that indexes the hour-τ layer by the point's timestamp — at the cost of a per-row
`terra::extract()`, precisely the overhead the matrix lookup was written to avoid — or a kernel
rebuilt per hour.

---

## 9. Provenance

| Walkthrough section | Source |
|---|---|
| Steps, random steps, model fitting | `DynamicSSF_1_Step_generation.qmd`, `DynamicSSF_2a_Model_fit_dry_season.qmd` |
| Hourly coefficients, RSF surfaces, `simulate_ssf()` | `DynamicSSF_3_Simulating_trajectories.qmd` |
| Rasterising simulated locations | `DynamicSSF_4a_Aggregating_simulations.R` |
| Convergence check | `DynamicSSF_4b_Simulation_convergence.R` |
| Hourly summary statistics and plots | `DynamicSSF_5_Hourly_summaries.qmd` |
| Boyce index, hourly validation | `DynamicSSF_6_Assessing_predictions.qmd` |

---

## 9a. Results, and how they should be read

The first complete render produced a more interesting outcome than expected, and the walkthrough
text was rewritten to report it honestly rather than to tell the tidy story.

| Comparison | Static | Dynamic |
|---|---|---|
| AIC | 14523.35 (edf 6) | **13311.24** (edf 43.1), ΔAIC 1212 |
| Boyce index, whole-day UD | **0.729** | 0.692 |
| Spearman vs observed KDE | 0.025 | **0.223** |
| Boyce index, mean across 24 hours | 0.275 | **0.529** (better in 21 of 24 hours) |

So the dynamic model wins decisively on fit, on the KDE comparison and on the hourly validation —
but is *marginally worse* on the single whole-day Boyce index. That is not a contradiction: the
whole-day UD pools every hour, which averages away exactly the advantage the dynamic model has.
The two models put the animal in much the same places over 30 days and disagree about the
schedule; a pooled metric cannot see a schedule.

This is now the walkthrough's closing point — the validation has to match the claim being made,
and a better-fitting model has not automatically earned the claim that it predicts better.

**The convergence check fails, and that is left visible on purpose.** Spearman correlation against
the full 500-trajectory UD runs 0.26 / 0.36 / 0.46 / 0.62 / 0.86 / 1.00 at 10 / 25 / 50 / 100 /
250 / 500 trajectories — still climbing steeply, nowhere near flat. An earlier draft of the prose
claimed the curve "flattens well before the full set"; that was wrong and was corrected. The
diagnostic now runs at two aggregation resolutions to show the second lever (coarser cells
converge much faster on the same simulations), and a callout states plainly that the results are
an illustration of the method rather than a settled statement about this animal.

## 10. Status

- [x] Design decisions settled and verified against the data
- [x] `by=` centring question resolved empirically
- [x] `knots` issue confirmed and fixed in `ssf_responses.qmd`
- [x] `predictive_SSF_walkthrough.qmd` written
- [x] Model fitting timings confirmed acceptable for rendering (~3.5 min for Model B)
- [x] Simulations run at full scale (500 x 720, three sets, ~2.5 min each)
- [x] Gifs generated (hourly selection surfaces; hourly UDs)
- [x] Turning-angle alignment bug found and fixed
- [x] `Rfast::rvonmises` NaN + non-reproducibility diagnosed; switched to `circular::rvonmises`
- [x] Validation results checked, and prose corrected where it contradicted them
- [ ] Final render with the `circular::rvonmises` switch (in progress)
- [ ] Whole-site render with the new navbar entry

### Known limitations, stated in the walkthrough rather than hidden

- The utilisation distributions have **not converged** at 500 trajectories.
- Fitted to a single individual over a single season, so nothing here generalises to the
  population; there is no random-effects or two-step estimation.
- NDVI is a single August 2018 layer held fixed, while the fitting window spans late July to
  October.
- The Boyce index validates against the same animal's observed locations — in-sample in the sense
  that the model was fitted to those same locations. A genuine out-of-sample test would hold out a
  period or an individual.

## 11. References

- Forrest, S.W. et al. (2025) *Ecography* 2025(2). doi:10.1111/ecog.07421
- Signer, J. et al. (2023) *Methods Ecol Evol*. Simulating animal space use from fitted iSSFs.
- Klappstein, N.J. et al. (2024) *Methods Ecol Evol*. Step selection functions with non-linear and
  random effects.
- Avgar, T. et al. (2016) *Methods Ecol Evol*. Integrated step selection analysis.
- Fieberg, J. et al. (2021) *J Anim Ecol*. A 'How to' guide for interpreting parameters in
  habitat-selection analyses.
