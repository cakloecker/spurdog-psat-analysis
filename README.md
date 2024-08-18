# Seasonal habitat use and diel vertical migration in female spurdogs
This repository contains R code used for analysing PSAT data and producing figures linked to the journal article entitled 'Seasonal habitat use and diel vertical migration in female spurdogs (Squalus acanthias Linnaeus, 1758) in Nordic waters' published in Movement Ecology. It processes archival PSAT data, incl. temperature, depth, light level and accelleration, of 19 female spurdogs tagged around the fjord system around Bergen, Norway for 82 - 366 days between 2019 and 2023. Amongst others it inspects the depth-temperature niche as well as periodicity in depth use across sub-daily to seasonal cycles using continous wavelet analysis.

### Contained files
- The `00_functions.R` script defines functions, i.e. used in `04_analysis_paper-supp.Rmd`, for calculating weights for niche plot, plotting time series, wavelet scalogramme etc.
- The 5s resolution archival data is read in and summarised on a level of minutes, hours and days in `01_archive-preprocessing.R`. 
- Data is inspected and cleaned in `02_data-cleaning.R` to remove potential tagging or capture effects. 
- A wavelet analysis is performed in `03_wavelet-analysis.R` for the hourly depth timeseries and adds a 'DVM' covariate indicating the presence of a significant 24h period in this timeseries.
- The Rmd-file `04_analysis_paper-supp.Rmd` is used to inspect and perform descriptive analysis on the data as well as to generate all figures contained in the main paper and the supplement. We also added the knitted `04_analysis_paper-supp.html` to include respective output.
- We also added the associated poster `BLS8-Poster_4PI-21_20240221.pdf` presented at the 8th International Bio-logging Science Symposium in Tokyo, Japan 2024.

### Visual Abstract
![VisualAbstract](https://github.com/user-attachments/assets/964c144b-ac52-4608-9404-700e47111240)

### Citation

Klöcker, C. A., Albert, O., Albert, O., Ferter, K. P., Bjelland, O., Lennox, R. J., Albretsen, J., Pohl, L., Dahlmo, L. S., Queiroz, N., & Junge, C. (2024). Seasonal habitat use and diel vertical migration in female spurdogs in Nordic waters.
