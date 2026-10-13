# Midlands PM2.5 population exposure and health risk profiles

Code and derived data for *Population exposure and health risk profiles for fine particulate matter in Midlands Province, Zimbabwe* (Nyoni, Moyo, Chinwadzimba, Sibanda and Tigere; manuscript under review).

The analysis combines modelled daily and annual PM2.5 (GHAP, with ACAG as a comparison product), year-matched WorldPop residential population and geoBoundaries district boundaries for 2017 to 2022. It reports population-weighted exposure, WHO guideline comparisons, cumulative and persistent exposure across eight district reporting areas, and conditional health contrasts transferred from published concentration–response functions. No survey, clinical or individual-level data are used.

## Contents

| Folder | What it holds |
|---|---|
| `analysis/` | Python and R scripts for data acquisition, population allocation, exposure metrics, health profiles, independent R verification and all figures |
| `data/metadata/` | Source collection identifiers, download manifest with SHA-256 hashes, boundary provenance and the GHAP unit audit |
| `data/public/` | Boundary subsets (GeoJSON) used in the analysis |
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

`05_make_figures.py`, `12_make_health_figures.py` and `17_refine_publication_figures.py` are earlier figure versions kept for transparency; steps 14 to 16 produce the published figures. Unit tests for the allocation and scientific calculations are in `analysis/test_*.py`.

## Key method choices

- GHAP V1 concentrations are used as supplied (multiplier 1.0); the unit audit is in `data/metadata/ghap_units_audit.json`.
- Population counts are allocated to the native exposure cells by area fraction in EPSG:6933, conserving totals. Municipal units are grouped with their district, and Redcliff is included in Kwekwe.
- Annual means require 90% of days and 75% in every quarter. Missing dates stay missing; nothing is interpolated.
- Health contrasts keep each published endpoint, exposure window, increment and uncertainty. They are conditional scenarios, not local health estimates, and overlapping endpoints are not added.

## Data sources and licences

Third-party data keep their own licences; please cite the original sources:

- **GHAP daily and annual PM2.5**: Wei et al. (2023), *Nature Communications* 14:8349, https://doi.org/10.1038/s41467-023-43862-3; accessed through the Earth Engine Community Catalog.
- **ACAG V6.GL.02 annual PM2.5**: Shen et al. (2024), *ACS ES&T Air* 1:332–345, https://doi.org/10.1021/acsestair.3c00054.
- **WorldPop Global2 R2025A v1** constrained population, Zimbabwe, 2017 to 2022 (CC BY 4.0).
- **geoBoundaries** gbOpen Zimbabwe ADM1 and ADM2, commit 9469f09 (CC BY 3.0 IGO; attribution ZIMSTAT/OCHA). Runfola et al. (2020), *PLoS ONE* 15:e0231866.
- **OpenStreetMap** tiles in the locator map: © OpenStreetMap contributors (ODbL).

## Licence

Code is released under the MIT licence (`LICENSE`). Derived data and figures created by the authors are released under CC BY 4.0 (`LICENSE-DATA`). Third-party inputs keep their original licences listed above.

## Citation

Please cite this repository using the archived release DOI (see `CITATION.cff`) and the associated article.

## Contact

Handsome Bongani Nyoni (nyonih@staff.msu.ac.zw) and Perez Livias Moyo (perezlivias1997@gmail.com).
