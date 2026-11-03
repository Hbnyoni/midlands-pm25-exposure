# Midlands PM2.5 exposure and its sources: code and derived data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23342800.svg)](https://doi.org/10.5281/zenodo.23342800)

Code and derived data for *Fires, industry and regional background: a source-attributed model of PM2.5 exposure in an unmonitored African province* (Nyoni, Chinwadzimba, Sibanda, Tigere and Moyo; manuscript under review at GeoHealth). Version 1.0.0 accompanied the first submission, *Population exposure and health risk profiles for fine particulate matter in Midlands Province, Zimbabwe*; version 1.1.0 adds the source-attributed exposure model described below.

The analysis combines modelled daily and annual PM2.5 (GHAP, with ACAG as a comparison product), year-matched WorldPop residential population and geoBoundaries district boundaries for 2017 to 2022. It reports population-weighted exposure, WHO guideline comparisons, cumulative and persistent exposure across eight district reporting areas, and conditional health contrasts transferred from published concentration–response functions. Version 1.1.0 adds a source-attributed exposure model that splits every 1 km cell-day of GHAP PM2.5 into a fire component (MODIS active fires within 50 km, estimated within district and year-month with ERA5-Land weather), a local industrial and urban component (sources identified from TROPOMI nitrogen dioxide) and a regional background, and attributes population-days above the WHO 24 hour guideline to each source with the Shapley value. No survey, clinical or individual-level data are used.

## Contents

| Folder | What it holds |
|---|---|
| `analysis/` | Python and R scripts for data acquisition, population allocation, exposure metrics, health profiles, independent R verification and all figures |
| `data/metadata/` | Source collection identifiers, download manifest with SHA-256 hashes, boundary provenance and the GHAP unit audit |
| `data/public/` | Boundary subsets (GeoJSON) used in the analysis |
| `data/sources/` | Source-model inputs: daily fire-pixel grids (0.05°, 2017 to 2022), TROPOMI NO2 and SO2 period means, monthly site series, daily fire and ERA5-Land series, ACAG V5.GL monthly district means and the site coordinates |
| `evidence/` | Published health functions, their source verification and the analysis plans |
| `results/` | Derived numerical tables (CSV and JSON) behind every table and figure |
| `qa/` | Independent R verification results and supporting check tables |
| `figures/` | Final figures as PNG, numbered as in the manuscript (`figure_number_mapping.json` links them to the scripts' file names) |

Large raster extracts (about 90 public GeoTIFFs) are not stored here. They are reacquired from the public sources with `analysis/01_acquire_province_gee.py`, and every file is checked against the hashes in `data/metadata/gee_download_manifest.json`.

## Reproducing the analysis

Requirements: Python 3.12 with `analysis/requirements_scientific.txt`, and base R 4.4.2 (plus ggplot2, patchwork, ragg and svglite for the statistical figures). Clear any conflicting `PROJ_LIB` and `PROJ_DATA` variables; the scripts do this for geospatial steps. Earth Engine downloads need your own authenticated account and registered Cloud project; no credentials are distributed.

From the repository root:

1. `python analysis/00_prepare_province.py`: boundary setup (the supplied boundaries can be kept)
2. `python analysis/01_acquire_province_gee.py --project YOUR_EARTH_ENGINE_PROJECT --workers 3`
3. `python analysis/02_build_population_overlays.py`
4. `python analysis/03_analyse_province.py`
5. `Rscript analysis/04_verify_in_R.R`
6. `python analysis/09_register_health_evidence.py`
7. `python analysis/10_build_health_profiles.py`
8. `Rscript analysis/11_verify_health_profiles.R .`
9. `python analysis/15_strengthen_inference.py`
10. `Rscript analysis/16_verify_inference.R .`
11. `python analysis/18_exposure_distributions.py`
12. `Rscript analysis/19_verify_distributions.R .`
13. `python analysis/24_fetch_osm_locator_tiles.py`: caches the OpenStreetMap tiles for the locator map
14. `python analysis/20_publication_maps.py`: map figures
15. `Rscript analysis/21_publication_graphs.R .`: statistical figures
16. `python analysis/22_submission_workflow.py`: workflow figure (Figure 1)

Source-attributed exposure model (version 1.1.0), after steps 1 to 4:

17. `python analysis/28_source_signatures.py --project YOUR_EARTH_ENGINE_PROJECT`: TROPOMI, fire and weather series
18. `python analysis/29_source_statistics.py --project YOUR_EARTH_ENGINE_PROJECT`: site interrupted time series and fire statistics
19. `python analysis/33_model_inputs.py --project YOUR_EARTH_ENGINE_PROJECT`: daily fire grids and ACAG monthly district means
20. `python analysis/34_source_model.py`: fire model with checks, local excess, decomposition and Shapley attribution (`results/model_*`)
21. `python analysis/35_model_figures.py`: Figures 1 to 11 of the revised article

`30_source_figures.py` draws the intermediate source figures of the first revision and is kept for transparency.

`05_make_figures.py`, `12_make_health_figures.py` and `17_refine_publication_figures.py` are earlier figure versions kept for transparency; steps 14 to 16 produce the published figures. Unit tests for the allocation and scientific calculations are in `analysis/test_*.py`.

## Key method choices

- GHAP V1 concentrations are used as supplied (multiplier 1.0); the unit audit is in `data/metadata/ghap_units_audit.json`.
- Population counts are allocated to the native exposure cells by area fraction in EPSG:6933, conserving totals. Municipal units are grouped with their district, and Redcliff is included in Kwekwe.
- Annual means require 90% of days and 75% in every quarter. Missing dates stay missing; nothing is interpolated.
- Fire effect: ln PM2.5 regressed on ln(1 + fire pixels within 50 km over three days) with district by year-month fixed effects, precipitation and wind; standard errors clustered by year-month; placebo with fires 3 to 5 days later.
- Local sources: TROPOMI NO2 excess (10 km disc minus 20 to 50 km ring) above zero at the 95% level; PM2.5 local excess from GHAP and ACAG annual means.
- Health contrasts keep each published endpoint, exposure window, increment and uncertainty. They are conditional scenarios, not local health estimates, and overlapping endpoints are not added.

## Data sources and licences

Third-party data keep their own licences; please cite the original sources:

- **GHAP daily and annual PM2.5**: Wei et al. (2023), *Nature Communications* 14:8349, https://doi.org/10.1038/s41467-023-43862-3; accessed through the Earth Engine Community Catalog.
- **ACAG V6.GL.02 annual PM2.5**: Shen et al. (2024), *ACS ES&T Air* 1:332–345, https://doi.org/10.1021/acsestair.3c00054.
- **WorldPop Global2 R2025A v1** constrained population, Zimbabwe, 2017 to 2022 (CC BY 4.0).
- **geoBoundaries** gbOpen Zimbabwe ADM1 and ADM2, commit 9469f09 (CC BY 3.0 IGO; attribution ZIMSTAT/OCHA). Runfola et al. (2020), *PLoS ONE* 15:e0231866.
- **MODIS active fires (FIRMS)**: Giglio et al. (2016), *Remote Sensing of Environment* 178:31–41; Earth Engine `FIRMS`.
- **ERA5-Land daily aggregates**: Muñoz-Sabater et al. (2021), *Earth System Science Data* 13:4349–4383.
- **Sentinel-5P TROPOMI NO2 and SO2, offline L3**: Veefkind et al. (2012), *Remote Sensing of Environment* 120:70–83.
- **ACAG V5.GL monthly PM2.5**: van Donkelaar et al. (2021), *Environmental Science & Technology* 55:15287–15300.
- **OpenStreetMap** tiles in the locator map and site coordinates: © OpenStreetMap contributors (ODbL).

## Licence

Code is released under the MIT licence (`LICENSE`). Derived data and figures created by the authors are released under CC BY 4.0 (`LICENSE-DATA`). Third-party inputs keep their original licences listed above.

## Citation

Please cite the archived release and the associated article:

Nyoni, H. B., Moyo, P. L., Chinwadzimba, Z., Sibanda, R., & Tigere, D. (2026). *Midlands PM2.5 exposure and its sources: code and derived data* (Version 1.1.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.23342800 (all versions; version 1.0.0: https://doi.org/10.5281/zenodo.23342801)

## Contact

Handsome Bongani Nyoni (nyonih@staff.msu.ac.zw) and Perez Livias Moyo (perezlivias1997@gmail.com).
