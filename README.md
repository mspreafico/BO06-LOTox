# BO06-LOTox

Code for implementing **longitudinal Latent Overall Toxicity (LOTox) profiles** analysis applied to MRC BO06 trial in Osteosarcoma.

Data are not publicly available due to confidentiality and privacy restrictions.

## Reference

Spreafico M., Ieva F. & Fiocco M. (2021). Longitudinal Latent Overall Toxicity (LOTox) profiles in osteosarcoma: a new taxonomy based on latent Markov models. https://arxiv.org/abs/2107.12863

## Description

- Files:
  - **01_LMmodel_selection.R**: Selection of the best latent Markov model for longitudinal toxicity data.
  - **02_LOTox_profiles.R**: Computation of both longitudinal Probability profiles of LOTox (P-LOTox) and longitudinal Relative Risk profiles of LOTox (RR-LOTox).
  - **utils.R**: Utils function for longitudinal LOTox profiles computation.
  - **legend_covariates**: Covariates legend.
- Sub-folder **./data/** contains pre-processed data.
- Sub-folder **./results/** contains results related to the best LM model (M5), LOTox sequences, P-LOTox and RR-LOTox profiles.


## Software
- R software.
- Packages: compositions, data.table, Formula, LMest, xtable.

(Last update: January 25th, 2021)
