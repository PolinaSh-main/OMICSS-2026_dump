# 03b_fst_equal_n — pairwise FST with the group sizes held equal

```
cd pipelines/03b_fst_equal_n
nextflow run .
```

Ten independent draws of 12 samples from each of the seven K = 7 groups,
all 21 comparisons rerun on each draw with the same `vcftools` call as
`03_fst`, plus within-group diversity measured on both the full groups
and the equal-n draws.

`03_fst` is untouched. It stays the reference run behind the group's
figures; this is a check on how its numbers should be read.

## Why it was run

In `03_fst` the 21 FST values correlate with group size at Spearman
−0.79. The five highest all involve a group of 12 or 13; the five lowest
all involve one of 22 to 49. The obvious suspicion was that this is a
sampling artefact: FST is between-group variance over total variance, so
a small sample would give a noisy, inflated estimate.

**That suspicion was wrong**, and the run says so clearly.

## What came out

| | full groups | equal n = 12 |
|---|---|---|
| Spearman(mean group size, FST) | −0.770 | **−0.756** |
| Spearman(smaller of the two, FST) | −0.573 | −0.575 |

Equalising the sample sizes changed nothing. Nor is it one odd group
carrying the effect — dropping K3 leaves −0.696, dropping K6 leaves
−0.582, dropping K4 leaves −0.775.

The ranking is also perfectly stable, which was the other worry:

```
median SD between the ten draws : 0.0019
median gap between neighbouring pairs : 0.0036
Kendall tau, old ranking vs new : +0.952
```

So the 21 numbers *can* be ranked. They reproduce.

## What the ranking actually measures

Nucleotide diversity within each group (mean per-site π over the
segregating sites, so "how variable are the members of this group",
not genome-wide π):

| group | n | π |
|---|---|---|
| K4 | 12 | 0.2403 |
| **K3** | 13 | **0.2214** |
| K1 | 15 | 0.2395 |
| K7 | 17 | 0.2489 |
| K2 | 22 | 0.2440 |
| K5 | 32 | 0.2417 |
| K6 | 49 | 0.2482 |

π barely moves with sample size — it is a per-site average, so cutting a
group from 49 to 12 leaves it almost unchanged. That is exactly why
equalising n did nothing.

Now put Hs, the mean π of the two groups in a pair, against their FST:

```
Spearman(Hs of the pair, FST)   : -0.840  (p = 1.9e-06)
Spearman(size of the pair, FST) : -0.756
implied Ht = Hs / (1 - FST)     : 0.2647 +/- 0.0035   (0.2553 to 0.2703)
```

**Ht — the diversity of the two groups pooled — is the same for every
one of the 21 pairs**, to within 1.3%. Hold it at its mean and move only
Hs and the whole observed spread falls out:

| | Hs | predicted FST |
|---|---|---|
| least diverse pair | 0.2305 | 0.1294 |
| most diverse pair | 0.2486 | 0.0610 |
| observed range | | 0.0557 – 0.1408 |

A 12% spread in Hs produces a threefold spread in FST because FST is one
minus a ratio close to one, so small moves in the numerator of that
ratio are amplified.

## How to say it

The 21 FST values are close to a **restatement of within-group
diversity**. Because the pooled diversity of any two groups is constant
across the dataset, FST here is measuring how internally uniform the two
groups are, not how far apart they sit.

The constancy of Ht is itself the finding worth showing: no pair of
Caucasian groups is more diverged than any other in absolute terms. It
is one gene pool.

K3 is the group to talk about — π = 0.221 against ~0.245 everywhere
else, and it appears in the top five pairs. Thirteen Armenian and
Azerbaijani cultivated accessions that are unusually similar to each
other. The reading is close kinship among a few cultivars, not an
anciently divergent lineage. That is a claim about K3 that a
relatedness check (KING, as in the reference paper) would settle.

## What this does not fix

Equal n does not equalise diversity, and diversity is what drives the
number. If the goal were groups that are comparable as populations, the
lever is the Q ≥ 0.75 threshold that defines them, not the sample size.
