# Finite-strata weak-IV reference implementation

`weak_iv_finite_strata.py` implements the conservative projection described in Section E.3 of the Technical Supplement. It is deliberately narrow: one pre-event period, one post-event period, a flat continuation, fixed finite strata, and fixed external reference weights. It requires Python 3.10 or later and uses only the Python standard library.

The package contains the runnable script, a worked CSV input, its expected JSON result, and unit tests. The script addresses only the narrow inferential case described here.

Run the included example from the project directory:

```bash
python3 implementation/weak_iv_finite_strata.py \
  implementation/example_finite_strata.csv \
  --output implementation/example_finite_strata_result.json
```

With the bundled inputs, the point estimate is 1.00 and the conservative 95 percent set is approximately `[-1.37, 3.50]`. Every simultaneous first-stage interval is consistent with the declared positive orientation, so the script does not trigger its causal stop. This check confirms the arithmetic; it is not a calibration for another application. The set is a statistical statement about the specified contrast. Calling that contrast an event-induced change in LATE still requires well-defined intervention and observation processes, valid period-specific IVs, support and a sampling bridge, no anticipation and factual/no-event exclusion, complier-population links, and the declared continuation.

The input CSV must contain:

- a stratum label and strictly positive fixed reference weight (omit zero-weight rows, and make the retained weights sum to one);
- the pre-event reduced form and its standard error;
- the pre-event first stage and its standard error;
- the post-event reduced form and its standard error; and
- the post-event first stage and its standard error.

The standard errors are inputs. They must already reflect the application's sampling weights, clustering, and dependence, and the primitive estimators must have credible asymptotic normal approximations. The script's normal critical value is not a few-cluster correction. It uses Bonferroni-simultaneous intervals for the primitive moments, projects them through each IV ratio, and combines the cell intervals with the fixed weights. If a required first-stage interval contains zero, it returns the whole real line.

The JSON output also reports `first_stage_orientation`. Point-estimated signs are a warning, not a test of a population assumption. The script therefore distinguishes them from the Bonferroni-simultaneous first-stage intervals. An interval wholly at or below zero contradicts the declared positive orientation at the reported confidence level and triggers `causal_interpretation_stop: true`. An interval that spans zero with a positive upper endpoint is labeled inconclusive and yields weak or unbounded ratio evidence, but sampling noise alone is not called a causal-assumption failure. The script never silently reverses the instrument or assigns a causal label.

Run the included checks with:

```bash
python3 -m unittest discover -s implementation -p 'test_*.py'
```

Do not use this script unchanged for estimated reference weights, data-dependent support or strata, flexible covariates, nonlinear continuations, overlapping cohort stacks, or estimated upper-level/cohort weights. It also does not choose between an average of upper-level-unit-by-stratum conditional ratios and a pooled moment ratio. Those cases require an inference procedure for the exact target and dependence structure.
