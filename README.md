# sweep_regions
command line tool to define outliers and merge them into regions

This is an experimental project. 
Use with an abundance of caution.

## Motivation 
In population genetics, we often want to define outliers for signals of adaptation
based on some summary statistic -- like FST -- that is above and beyond what we can expect from neutral demographic processes.
However, outliers maybe defined at single SNP positions or a handful of them with considerable local autocorrelation.
We don't really want to treat those autocorrelated regions as independent
as they are likely being driven by the same process plus linkage.
Some methods like Baypass deal with this formally and nicely, but not all do.
This is my pretty simple duct tape style attempt to remidy this using R's `smooth.spline()` function.
The basic idea is to use a smoothing spline to identify outlier regions with a statistic that is consitently elevated.
The outliers can be found based on quantiles from neutral simulations.
We can then merge those regions that are within a certain distance 
and exclude those that are too large or small. 
Needless to say this all takes considerable thinking and experimetation ahead of time to make reasonable choices, 
ideally with lots of simulations where the truth is known. 

Requirements

- This repo or at least the `sweep_regions.R`

- up to date version  of R

- packages `argparse`, `dplyr`, `tidyr`, and `vroom`



with RAiSD-like files as input

```bash
./sweep_regions.R neutral_cutoff --quantile 0.99 -f example_data/raisd2_neutral_example.txt --input_raisd2 TRUE --delimeter "\t"

./sweep_regions.R sweep_regions --data_frame example_data/raisd8_example.txt -R TRUE --delimeter "\t" --merge_size 1e5 --cutoff 100 --min_size 1e3 --max_size 1e8 --out_file example_data/raisd_out.csv
```


using generic input files with headers to identify the `positions` and statistic `values` columns

```bash
./sweep_regions.R neutral_cutoff --quantile 0.99 -f example_data/generic_neutral_example.txt --delimeter "\t" --positions position --values stat


./sweep_regions.R sweep_regions --data_frame example_data/generic_example.txt --delimeter "\t" --merge_size 1e5 --cutoff 100 --min_size 1e3 --max_size 1e8 --out_file example_data/raisd_generic_out.csv --positions position --values stat --chromsome chr1
```


If true is known about specific regions, in the case of simulated inout, this table can be provided that lists a position as a sweep or neutral. The output regions will then be labeled with a `truth class` column that lists each regions as `true_positive` or `false_positive`. `true_negatives` and `false_negative`  will also be added to table if known neutral positions are in the `truth_data` input.

run unit tests 
```bash
sweep_regions.R --do_tests
```