# Stanza handling in the Atlantis to EwE conversion

`R-scripts/atlantis_to_ewe_funcs.R` converts Atlantis output into the EwE style csv files read by
`EBFM_indicator_calculation_template.qmd`. By default each Atlantis functional group becomes a single EwE
group: the age classes (cohorts) are only kept separate inside the calculations and are then summed.

Setting `stanzas = TRUE` instead splits the age structured groups into stanzas, like EwE multi-stanza groups.
Each stanza is written as its own group in all the output files (biomass, catch, Q/B, mortality, TL, the
consumption and diet matrices and the Species_ID template).

## Using it

In `R-scripts/run_atlantis_to_ewe.R` (or the `atlantis.conversion` chunk of the qmd, where the settings are
called `atl_stanzas`, `atl_stanza_groups`, and so on):

```r
stanzas <- TRUE          # Split the age structured groups
stanza_groups <- NULL    # Codes of the groups to split - NULL for all age structured groups, e.g. c("FTS", "FLG")
stanza_breaks <- NULL    # Optional - first cohort of each stanza after the first, e.g. list(FTS = c(3, 6))
stanza_labels <- NULL    # Optional - labels of the stanzas, e.g. list(FTS = c("larval", "juv", "adult"))
```

Or call the function directly:

```r
convert_atlantis_to_ewe(..., stanzas = TRUE, stanza_groups = NULL, stanza_breaks = NULL, stanza_labels = NULL)
```

The runner script and the qmd pass the same stanza settings to both the scenario run and the unfished
run. If you call `convert_atlantis_to_ewe()` yourself, do the same, so that both runs have the same groups.

With `stanzas = FALSE` (the default) the output is identical to the output before stanza handling was added.

## Which groups are split, and how

- Only groups with per-cohort tracers (`<Name><cohort>_Nums`, `_StructN` and `_ResN`) in the main output
  `.nc` can be split. These are the vertebrates, including `FISH_INVERT` groups such as squid and octopus.
  Requesting any other group in `stanza_groups` gives a warning and the group is left whole.
- **Default split:** the Atlantis juvenile / adult split. A cohort is an adult once its age class
  (cohort / `NumGeneTypes`) reaches `<Code>_age_mat` in the biology prm file, as set up in `atbiolsetup.c`.
  The stanzas are labelled `juv` and `adult`.
- **Custom splits:** `stanza_breaks` gives the first cohort of each stanza after the first, counting
  cohorts from 1 as in the Atlantis variable names. For example, `list(FTS = c(3, 6))` gives cohorts 1-2,
  3-5 and 6-10. Custom stanzas are labelled `s1`, `s2`, ... unless `stanza_labels` is given.
- Empty stanzas are dropped. A group left with only one stanza (e.g. `age_mat` of 0) is not split.
- Stanza groups are named `<Code>_<label>` (e.g. `FTS_juv`) and `<Name>_<label>` (e.g. `Striped_tuna_juv`).
  They take the place of the Atlantis group in the group order, so the detritus groups stay at the end as
  `funcs.R` expects.

## How each quantity is derived for a stanza

| Quantity | Derivation |
|---|---|
| Biomass | Biomass per cohort = `Nums * (StructN + ResN) * X_CN * k_wetdry / 1e9` (tonnes) from the main output `.nc`, excluding boundary boxes, taken as the annual mean of the snapshots in each model year. These cohort biomasses are used to share the group's `BiomIndx.txt` biomass across its stanzas, so group totals match the unsplit conversion. |
| Consumption (predator side) | Each cohort's consumption (`<Name><cohort>_Eat` in `PROD.nc`) and its own diet (`DietCheck.txt`) go to the stanza the cohort is in. |
| Consumption (prey side) | See [Sharing consumption across prey stanzas](#sharing-consumption-across-prey-stanzas) below. |
| Q/B, M2, TL | Calculated from the stanza-level consumption matrix and stanza biomass, as for unsplit groups. |
| M0 | Biomass-weighted mean, across the cohorts in the stanza, of `-log(1 - M1)` from the `<Code>-<cohort>-S<stock>-M1` columns of `SpecificMort.txt`. Unsplit groups keep the simple (unweighted) mean across cohorts. |
| Catch, landings, discards | Each group's catch per fleet is shared across its stanzas using the catch and discards per cohort in `CATCH.nc` (`<Name><cohort>_Catch` and `<Name><cohort>_Discards` numbers, times the mean individual weight). If these are not in the output, the catch is shared in proportion to the stanza biomass, the same share for every fleet. |
| Species_ID: Linf | Length of the oldest cohort in the stanza (initial conditions weight with `li_a`, `li_b`). |
| Species_ID: MaxAge | The species longevity (`NumCohorts * NumAgeClassSize`), the same for every stanza. |
| Species_ID: other columns | `PelDemID` and `Classification` come from the parent group. `PISCZOOPL` is recalculated from each stanza's own diet. `pisczoopl_groups` accepts stanza codes or names, or the parent Atlantis code or name (which applies to all of its stanzas). |

## Sharing consumption across prey stanzas

`DietCheck.txt` reports only the prey group eaten, not the prey cohort. So when the prey group is split, each
predator cohort's consumption of that prey group is shared across the prey cohorts (and so the prey stanzas).
This follows `Avail_Fish()` in `atecology/atvertprocesses.c`. The availability of prey cohort *c* to a
predator cohort is:

```
pPREY[prey stage of c][predator][predator stage]  x  size window scalar  x  biomass of prey cohort c
```

- **`pPREY`:** read from `pPREY<prey stage><Pred><pred stage>` in the biology prm file (`pPREY<Pred>` for
  invertebrate predators). Where a group has only stage 1 values, they are used for stage 2 too.
- **Size window:** the `sizeScalar` of `Avail_Fish()`, applied only when the prey is a vertebrate. The form
  depends on the flags in the biology prm file:
  - `UseHardFeedingWindow = 1`: available if `KLP * predator StructN <= prey StructN <= KUP * predator StructN`.
  - Otherwise, with `UseBiLogisticFeedingWindow = 0`: the humped curve on StructN using `Kmax_coefft_`.
  - Otherwise, with `UseBiLogisticFeedingWindow = 1`: the bi-logistic curve on length (cm), from
    `(StructN + ResN)` and `li_a_`, `li_b_`.

  The scalar is 1 when the predator has no cohort sizes (e.g. invertebrate predators) or no window parameters.
- **Sizes:** the annual mean, weighted by numbers, of each cohort's individual `StructN` and `ResN`.
- **Spatial overlap:** the sharing is done box by box, with the boxes weighted by how much the predator
  cohort ate in each (its `Eat` per box). Consumption in boxes where none of the prey was present is shared
  using the whole-model-domain availability.
- **Fallback:** if no prey cohort is available at the annual mean sizes (although `DietCheck.txt` shows the
  prey was eaten), the consumption is shared by prey cohort biomass. The run reports how often this happened.

This only redistributes consumption among a prey group's stanzas. Summed back over stanzas, the consumption
matrix is identical to the unsplit conversion.

## Limitations

- The refuge, stock, turbidity and pH scalars in `Avail_Fish()` are not included.
- Sizes and biomasses are annual means, so cohorts that grow into or out of a feeding window during the year
  are only approximated. Feeding is matched at the box level, not by layer, because `Eat` is written per box.
- `flag_olddiet = 1` is not supported for the prey split, and a warning is given.
- Juvenile Linf is the length at the end of the juvenile stanza. For late-maturing groups it is close to the
  adult value. Check and edit the Species_ID template as needed.

## Checks and messages during a run

- **Cohort biomass vs `BiomIndx.txt`:** the range of the ratio is reported. They agree at matching snapshot
  times. They differ because `BiomIndx.txt` is often written only yearly (`tsumout`), while the `.nc` is
  written every `toutinc`, and because `BiomIndx.txt` includes migrants outside the model domain. For
  example, in the SEAP NESP run humpbacks are 264,900 t in `BiomIndx.txt` but 1,373 t inside the domain.
- The number of predator cohort - prey pairs that used the biomass fallback for the prey split.
- Any groups whose catch had to be shared by biomass because `CATCH.nc` had no per-cohort catch.

## Testing (SEAP NESP wind run, 20 years, no fisheries output)

- With `stanzas = FALSE`, all 12 output files are byte-identical to the previous version.
- With `stanzas = TRUE`, 73 groups were split into 146 stanzas (174 output groups). The conversion took
  about 1.3 minutes.
- Summed over stanzas, the consumption matrix matched the unsplit conversion to within 1.5e-8 t (of 1.9e13 t
  in total), and biomass matched to rounding error.
- The catch split was not exercised, as that run has no fisheries output.
