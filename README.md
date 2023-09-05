# BO06-LOTox

Code for implementing **longitudinal Latent Overall Toxicity (LOTox) profiles** analysis applied to MRC BO06 trial in Osteosarcoma.

### Reference
Spreafico M., Ieva F. & Fiocco M. (2023). Longitudinal Latent Overall Toxicity (LOTox) profiles in osteosarcoma: a new taxonomy based on latent Markov models. https://arxiv.org/abs/2107.12863

### Data Availability
Data are not publicly available due to confidentiality and privacy restrictions.
Access to the full dataset of MRC BO06 trial can be requested to MRC Clinical Trials Unit at UCL, Institute of Clinical Trials and Methodology, UCL, London.

Along with the code a fake dataset is provided in order to allow researchers who want to replicate the same analysis to properly get how the code has to be run and how results are displayed.

## Description

- Files:
  - **01_LMmodel_selection.R**: Selection of the best latent Markov model for longitudinal toxicity data.
  - **02_LOTox_profiles.R**: Computation of both longitudinal Probability profiles of LOTox (P-LOTox) and longitudinal Relative Risk profiles of LOTox (RR-LOTox).
  - **utils.R**: Utils function for longitudinal LOTox profiles computation.
- Sub-folder **./data/** contains pre-processed data. [NOT PUBLICY ACCESSIBLE due to privacy restrictions]
- Sub-folder **./results/** contains results related to the best LM model (M5), LOTox sequences, P-LOTox and RR-LOTox profiles. [NOT PUBLICY ACCESSIBLE]
- Sub-folder **./fake_data_results/** contains the fake dataset along with its legend and results:
	- **fake_dataset.Rdata**: fake dataset related to 75 subjects to emulate pre-processed information used within the study;
	- **legend_covariates.txt**: variables legend;
  - **fake_LM_final_M5.Rdata**: LM model M5 fitted on the fake dataset;
  - **fake_LOTox_profiles.Rdata**: LOTox profiles obtained by LM model M5 fitted on the fake dataset.


## Software
- R software.
- Packages: compositions, data.table, Formula, LMest, xtable.

(Last update: September 5th, 2023)
