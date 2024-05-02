# MagTetris-GA
ezyMRI Hackathon 2024

Implementation of a genetic algorithm for optimizing sparse Halbach array
lowfield MRI designs. Uses
[MagTetris](https://github.com/BioMed-EM-Lab/MagTetris) as the core field estimator.

## Layout
`code` - Main code files.

Main script at `halbach_ga.m`.

### Remanence-Robust Array Optimization
`robust_halbach_ga.m`, `robust_halbach_ga_obj.m`

Small variations in the remanence of permanent magnets may lead to the actual
inhomogeneity being larger than the naive optimization would suggest. The
remanence-robust version of this algorithm samples `N` random remanence
distributions over the entire array and takes either the worst-performing one
(`robust_mode="worst"`) or the average performances (`robust_mode="mean"`) as
its field objective.


`deviation.m`

Assess the expected performance and robustness of an array design by Monte-Carlo
sampling from a remanence distribution (assumed Gaussian).

#### Results
The non-robust design outperforms the robust design in B0 field strength.
The robust design outperforms the non-robust design on average in B0 field
inhomogeneity, for large standard deviation, but the improvement is relatively
small. For the 6 layer design, the final numbers were:

##### zy-plane
Field Strength
- Non-robust: 161 mT
- Robust: 136 mT

B0 inhomogeneity (max deviation / mean)
- @std. dev. = 1e-2 (realistic)
  - Non-robust: 9836ppm
  - Robust: 6047ppm
  - Robust non-best (population): 5607ppm
- @std. dev. = 1e-6 (randomness basically negligible)
  - Non-robust: 165 ppm
  - Robust: 333 ppm
  - Robust non-best (population): 5432 ppm
     - Basically the same as the high-noise case. But at least it's the same!

I also checked what happens in the average case (e.g. for the robust GA, maybe
the "best" solution is actually just a lucky solution). In fact, the random
robust population member that I picked DID have a robust B0 inhomogeneity
curve. The curve was flat! However, the actual inhomogeneity was about the same
as the non-robust solution's inhomogeneity in the 1e-2 case, and was strictly
worse for all other levels of randomness.



## Implementation notes
### Force calculation
Force calculation is very expensive. By default, the force objective is disabled
and the genetic algorithm is allowed to run until convergence without it. If the
force constraint is violated, the final population from the first stage can be
used to initialize the GA solve with the force objective active.




## References
If you use this code, please cite the original MagTetris Paper:

```
Ting-Ou Liang, Yan Hao Koh, Tie Qiu, Erping Li, Wenwei Yu, Shao Ying Huang,
MagTetris: A simulator for fast magnetic field and force calculation for
permanent magnet array designs, Journal of Magnetic Resonance, Volume 352, 2023,
107463, ISSN 1090-7807, https://doi.org/10.1016/j.jmr.2023.107463.
```
