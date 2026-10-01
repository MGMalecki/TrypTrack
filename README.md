<p align="center"><img src="www/logo.svg" alt="TrypTrack" width="420"></p>

**TrypTrack** is a Shiny app for growth curves and doubling times of
*Trypanosoma brucei* cultures that are counted and diluted repeatedly, for
example after RNAi induction with tetracycline.

**Live app:** https://michamaleckis.shinyapps.io/tryptrack/

## What it does

1. **Experiment design.** Choose cell lines, conditions (e.g. tet- / tet+),
   number of biological replicates or clones, the counting schedule (every day,
   every 2 days or custom days), the usual counting time and default dilution
   factors. The app builds a form to fill in.
2. **Data.** Enter the average count for each culture and day. Only the
   dilution factor is needed; the density after dilution is calculated.
   Each count has a day and a clock time, so counting at 11:00 one day and
   16:00 the next is handled correctly. Save the table as CSV and upload it
   again the next day to continue, or add more days at any time.
3. **Plots.**
   - cumulative growth corrected for dilution (log scale)
   - measured densities with dilutions (the sawtooth)
   - doubling time or doublings per 24 h for each interval
   - growth relative to a reference condition (e.g. tet+ as % of tet-)

   Replicates can be shown as mean +/- SD or as separate panels (useful for
   picking the clone with the strongest phenotype). Plots download as PDF/PNG.
4. **Doubling time.** Log-linear fit over a chosen time window, per replicate
   and as mean +/- SD.

## How the numbers are calculated

| quantity | formula |
|---|---|
| density after dilution | `count / dilution_factor` (or `count_after_dilution` if measured) |
| growth in an interval | `count / density after the previous dilution` |
| cumulative cells | `first count x product of all growth factors` |
| doubling time of an interval | `interval (h) x ln 2 / ln(growth)` |
| elapsed time | `24 x (day - first day) + (clock - clock of first count)` |

Densities are shown as cells/ml = count x 10^4 (adjustable).

## Data format

One row per count. `example_data.csv` is a complete example.

| column | meaning |
|---|---|
| `cell_line` | RNAi line / genotype |
| `condition` | e.g. `tet-`, `tet+` |
| `replicate` | biological replicate or clone |
| `day` | day of the experiment (0, 1, 2, ...) |
| `clock` | time of counting, e.g. `11:00` (empty = same as the first count) |
| `count` | average count from the counting device |
| `dilution_factor` | dilution made right after this count (empty on the last day) |
| `count_after_dilution` | optional, only if counted again after diluting |

Files from the first version of the app (`cells_diluted` / `cells_24h`) are
converted automatically.

## Run it on your own computer

```r
install.packages(c("shiny", "bslib", "dplyr", "tidyr", "ggplot2",
                   "rhandsontable", "zip"))
shiny::runGitHub("TrypTrack", "YOUR-GITHUB-USERNAME")
```

or download the repository and open `app.R` in RStudio -> **Run App**.
