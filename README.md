# Gulf of Lion Bird Migration Data Processing

## Overview

This project processes GPS tracking data from multiple bird migration studies to analyze flight patterns over the Gulf of Lion region in the Mediterranean Sea. The workflow identifies individuals that cross the Gulf of Lion, processes their tracking data, classifies their behavior, and generates summary visualizations.

**Author:** Victor Cat  
**Date:** December 2024

---

## Project Structure

```
project/
├── Analysis.R                          # Main analysis script
├── R/
│   ├── 01_formatting_functions.R       # Data cleaning and formatting
│   ├── 02_spatial_functions.R          # Spatial analysis functions
│   ├── 03_read_data.R                  # Data reading functions
│   └── 04_visualization_functions.R    # Plotting and visualization
├── data/
│   ├── MB_1552820885.csv              # Movebank data
│   ├── none_MB/                        # Non-Movebank data sources
│   │   ├── Tourterelle&Oedicnemes/
│   │   ├── Bondrée_apivore_patrick_byholm/
│   │   ├── BIRDMAP/
│   │   └── schlaich_et_al.2016/
│   └── supp/
│       ├── europe.rds                  # European coastline shapefile
│       ├── GOL_extended.rds            # Extended study area
│       ├── GOL_restricted.rds          # Restricted Gulf of Lion area
│       └── EC.xlsx                     # Ecological groups mapping
├── data/raw_individuals/               # Temporary: individual CSV files
└── output/
    ├── gps_migration.csv               # Final processed dataset
    └── plots/                          # Generated visualizations
```

---

## Workflow Overview

The analysis follows a **four-phase workflow**:

### **PHASE 1: Data Reading and Splitting**
- Reads data from multiple sources (Movebank + 5 other datasets)
- Splits large datasets by individual bird
- Saves one CSV file per individual in `data/raw_individuals/`

### **PHASE 2: Individual Processing**
- Applies formatting pipeline to each individual
- Filters, downsamples, classifies, and enriches data
- Appends processed data to final output CSV
- Removes temporary individual files

### **PHASE 3: Summary and Visualization**
- Loads final processed dataset
- Adds ecological group classifications
- Generates summary statistics and plots

### **PHASE 4: Bayesian Flight Height Modeling (Optional)**
- Hierarchical Bayesian analysis of flight altitudes
- Models post-breeding migration patterns
- Calculates proportions by altitude range
- Generates ecological group-specific distribution plots
- **Note:** Computationally intensive (10-30 minutes)

---

## Phase 1: Data Reading

### Data Sources

1. **Movebank** (`MB_1552820885.csv`)
   - Large dataset processed in chunks (100,000 rows at a time)
   - Memory-efficient streaming approach

2. **Turtle Dove & Stone Curlew** (Frédéric Jiguet)
   - 6 CSV files with different individuals
   - Species: *Streptopelia turtur*, *Burhinus oedicnemus*

3. **Honey Buzzard** (Patrick Byholm)
   - 3 files with different formats
   - Species: *Pernis apivorus*

4. **BIRDMAP** (Adojaan K et al.)
   - Occurrence database with multiple species
   - Tab-delimited format

5. **Montagu's Harrier** (Schlaich et al. 2016)
   - Published study data
   - Species: *Circus pygargus*

### Output
- One CSV file per individual bird in `data/raw_individuals/`
- Standardized column names across all sources

---

## Phase 2: Processing Pipeline

Each individual file goes through the following **11-step pipeline**:

### 1. **Initial Filtering**
```r
clean_row()
```
- Removes rows with missing coordinates or timestamps
- Removes duplicate entries
- Sorts by timestamp

### 2. **Equipment Location**
```r
first_loc()
```
- Records first GPS location (assumed equipment deployment site)
- Stores in `equip_loc` column

### 3. **Spatial Filtering**
```r
select_pts(area_res, area_ext)
```
- **Requirement:** Individual must cross the restricted Gulf of Lion area
- **Action:** Keeps only points within the extended study area
- **Result:** Excludes individuals that never crossed Gulf of Lion
- **Result:** Removes points outside extended area for qualifying individuals

### 4. **Skip Conditions (Applied Here)**
- Skip if no `height_above_msl` data
- Skip if excluded species (marine residents, etc.)
- Skip if no points after spatial filter

### 5. **Temporal Downsampling**
```r
downsample_tracks(time_interval = 10)
```
- Reduces temporal resolution to 10-minute intervals
- Keeps middle point per interval
- Reduces data volume while maintaining trajectory structure

### 6. **Trajectory Assignment**
```r
assign_trajectories(max_gap_days = 5)
```
- Splits tracks into discrete trajectories
- **Splitting criteria:**
  - July 1st boundary (pre/post breeding season)
  - Gaps > 5 days between consecutive points
- **Naming:** `{individual_id}_traj_{n}`
- Example: `bird_123_traj_1`, `bird_123_traj_2`

### 7. **Trajectory Filtering**
```r
remove_short_trajectories(min_points = 5)
```
- Removes trajectories with fewer than 5 points
- Ensures sufficient data for analysis

### 8. **Skip Conditions (Applied Here)**
- Skip if no valid trajectories remain

### 9. **Season Classification**
```r
add_season()
```
- Classifies each point by month:
  - **January-June:** "pre-breeding"
  - **July-December:** "post-breeding"

### 10. **Restricted Area Classification**
```r
check_restricted_area(area_res)
```
- Checks if each point falls within Gulf of Lion
- Adds `in_GOL` column (TRUE/FALSE)

### 11. **Behavior Classification**
```r
classify_staging_migration(eps = 0.1, minPts = 10)
```
- Uses DBSCAN clustering on coordinates
- **Parameters:**
  - `eps = 0.1`: Maximum distance between points (degrees)
  - `minPts = 10`: Minimum points to form cluster
- **Classification:**
  - Points in clusters → **"Staging"** (stopover sites)
  - Points not in clusters → **"Migration"** (active movement)
- Adds columns: `db_cluster`, `state`

### 12. **Distance to Shore**
```r
calculate_distance_to_shore(coastline = eu)
```
- Calculates distance to nearest shoreline
- Uses projected CRS (EPSG:3035) for accuracy
- Adds `distance_to_shore` column (meters)

### 13. **Type Conversion**
```r
class_col()
```
- Converts all columns to character for CSV storage
- Ensures compatibility when appending data

### Final Output
- Appends processed data to `output/gps_migration.csv`
- Deletes temporary individual file

---

## Phase 3: Summary and Visualization

### Data Loading
- Reads final `output/gps_migration.csv`
- Converts numeric columns back to proper types
- Adds ecological group classifications from `EC.xlsx`

### Ecological Groups

Five groups based on taxonomy and behavior:

1. **Ardeid** (Herons)
   - *Ardea purpurea*, *Ardeola ralloides*, *Ixobrychus minutus*

2. **Large wader** (Large wading birds)
   - *Grus grus*, *Phoenicopterus roseus*, *Plegadis falcinellus*, *Platalea leucorodia*

3. **Larid** (Gulls and terns)
   - *Larus fuscus*, *Chroicocephalus* spp., *Chlidonias hybrida*, etc.

4. **Shorebird** (Shorebirds and waders)
   - *Limosa limosa*, *Numenius arquata*, *Burhinus oedicnemus*, etc.

5. **Terrestrial** (Terrestrial migrants)
   - *Circus pygargus*, *Falco* spp., *Pernis apivorus*, *Streptopelia turtur*, etc.

### Generated Plots

Six high-resolution plots (300 dpi PNG) saved to `output/plots/`:

1. **`01_gol_migration_both_seasons.png`**
   - Spatial map of migration through Gulf of Lion
   - Faceted by season (pre-breeding | post-breeding)
   - Colored by individual, no legend

2. **`02_gol_migration_vs_staging.png`**
   - Spatial map showing migration vs staging behavior
   - Colored by state (orange = Migration, blue = Staging)
   - All data shown (not filtered to Gulf of Lion only)
   - All seasons combined

3. **`03_gol_ecological_groups_pre_breeding.png`**
   - Pre-breeding migration faceted by ecological group
   - 3-column layout

4. **`04_gol_ecological_groups_post_breeding.png`**
   - Post-breeding migration faceted by ecological group

5. **`05_flight_height_distribution.png`**
   - Histogram of flight heights (meters above MSL)
   - Only includes migration points within Gulf of Lion

6. **`06_flight_height_by_ecological_group.png`**
   - Flight height histograms faceted by ecological group
   - Free y-axis scales for different sample sizes

---

## Phase 4: Bayesian Flight Height Modeling (Optional)

This optional phase performs hierarchical Bayesian analysis of bird flight altitudes during post-breeding migration.

### Model Description

**Hierarchical lognormal model** with random effects by ecological group:

```
Observation level:
  obs_alt[i] ~ Normal(true_alt[i], sigma_obs)
  
Latent process:
  true_alt[i] ~ Lognormal(mu[i], sigma_lognorm[i])
  mu[i] = alpha_0_re[group] + alpha_1_re[group] * distance_to_shore
  sigma_lognorm[i] = beta_0_re[group]

Random effects:
  alpha_0_re[j] ~ Normal(alpha_0_mu, alpha_0_sigma)
  alpha_1_re[j] ~ Normal(alpha_1_mu, alpha_1_sigma)
  beta_0_re[j] ~ Normal(beta_0_sigma, beta_0_sigma_sigma)
```

### Data Filtering

- Post-breeding season only
- Migration state (not staging)
- Points in Gulf of Lion (`in_GOL == TRUE`)
- Height < 5000 m
- Removes bottom 2% of observations

### MCMC Settings

- **Iterations:** 50,000 per chain
- **Burn-in:** 10,000 (20%)
- **Thinning:** 5
- **Chains:** 2
- **Convergence:** Monitored via Gelman-Rubin diagnostics

### Outputs

**Model files** (`output/model/`):
- `model_lognorm_pstB.rds` - Full MCMC samples and WAIC
- `convergence_diagnostics.txt` - Summary statistics

**Proportion tables** (`output/proportions/`):
- `proportions_[group_name]_pstB.rds` - By ecological group (ardeid, large_wader, larid, shorebird, terrestrial, population)
- `proportion_data_pstB.rds` - Combined dataset

**Visualizations** (`output/plots/`):
- `07_model_validation_pstB.png` - Observed vs simulated distributions
- `08-13_altitude_distribution_[group_name]_pstB.png` - Group-specific distributions (ardeid, large_wader, larid, shorebird, terrestrial, population)
- `proportions_table_pstB.html/png` - Summary table

### Altitude Range Proportions

The model calculates the proportion of birds flying in each altitude range:
- **0-20 m:** Very low altitude flight
- **20-200 m:** Low altitude flight (collision risk zone)
- **20-300 m:** Extended low altitude zone
- **Above 300 m:** High altitude flight

Results include median estimates with 95% credible intervals.

### Running the Model

By default, Phase 4 is **disabled** (computationally intensive).

To enable:
```r
# In Analysis.R, line ~395:
run_modeling <- TRUE  # Change from FALSE to TRUE
```

Expected runtime: **10-30 minutes** depending on hardware.

### Required Packages

Additional packages for Phase 4:
- `nimble` - MCMC sampling
- `coda` - MCMC diagnostics
- `MCMCvis` - Visualization
- `gt` - Table formatting
- `cowplot` - Plot composition

---

## Data Columns in Final Output

| Column | Type | Description |
|--------|------|-------------|
| `individual_id` | character | Unique bird identifier |
| `study_id` | character | Source study/project |
| `individual_taxon_canonical_name` | character | Scientific species name |
| `location_long` | numeric | Longitude (WGS84) |
| `location_lat` | numeric | Latitude (WGS84) |
| `height_above_msl` | numeric | Flight height (meters above mean sea level) |
| `timestamp` | POSIXct | Date and time (UTC) |
| `equip_loc` | character | First location (equipment site) |
| `trajectory` | character | Trajectory ID (format: `id_traj_n`) |
| `season` | character | "pre-breeding" or "post-breeding" |
| `in_GOL` | logical | TRUE if point is in Gulf of Lion restricted area |
| `state` | character | "Migration" or "Staging" |
| `db_cluster` | integer | DBSCAN cluster ID (0 = noise/migration) |
| `distance_to_shore` | numeric | Distance to nearest shoreline (meters) |
| `ecological_group` | character | Ecological group classification |

---

## Excluded Species

The following species are excluded from analysis (marine residents, non-migratory):

- *Puffinus yelkouan* (Yelkouan Shearwater)
- *Calonectris diomedea* (Scopoli's Shearwater)
- *Thalasseus sandvicensis* (Sandwich Tern)
- *Larus michahellis* (Yellow-legged Gull)
- *Aquila clanga* (Greater Spotted Eagle)
- *Aquila chrysaetos* (Golden Eagle)
- *Ciconia ciconia* (White Stork)

---

## Configuration Parameters

### Spatial
- **Restricted area:** Gulf of Lion core zone (`GOL_restricted.rds`)
- **Extended area:** Broader study region (`GOL_extended.rds`)
- **Coastline:** European shoreline for distance calculations (`europe.rds`)

### Temporal
- **Downsampling interval:** 10 minutes
- **Trajectory gap threshold:** 5 days
- **Season boundary:** July 1st

### Trajectory Filtering
- **Minimum trajectory length:** 5 points

### Behavior Classification (DBSCAN)
- **eps:** 0.1 degrees (~11 km at equator)
- **minPts:** 10 points minimum for cluster

---

## Required R Packages

### Core Packages (Phases 1-3)
```r
tidyverse    # Data manipulation and visualization
sf           # Spatial data handling
readr        # Fast CSV reading
readxl       # Excel file reading
lubridate    # Date/time manipulation
dbscan       # Clustering for behavior classification
DT           # (optional) Interactive tables
parallel     # (optional) Parallel processing
pbapply      # (optional) Progress bars
```

### Additional Packages for Phase 4 (Bayesian Modeling)
```r
nimble       # MCMC sampling and Bayesian inference
coda         # MCMC diagnostics and convergence
MCMCvis      # MCMC visualization and summary
gt           # Table formatting and styling
cowplot      # Plot composition
patchwork    # Plot arrangement
transport    # Optimal transport distances
grid         # Graphics devices
png          # PNG image handling
```

---

## Memory Efficiency

The workflow is designed for **memory-efficient processing** of large datasets:

1. **Chunked reading:** Movebank data processed in 100k row chunks
2. **Individual file strategy:** Never loads entire dataset into memory
3. **Streaming append:** Writes processed data incrementally to CSV
4. **Temporary file cleanup:** Deletes individual files after processing

This approach allows processing datasets with **millions of GPS points** on standard computers.

---

## Usage

### Basic Execution

```r
# Run the entire workflow
source("Analysis.R")
```

### Step-by-Step Execution

```r
# Load packages and functions
source("R/01_formatting_functions.R")
source("R/02_spatial_functions.R")
source("R/03_read_data.R")
source("R/04_visualization_functions.R")

# Run Phase 1: Data reading
# [Phase 1 code from Analysis.R]

# Run Phase 2: Individual processing
# [Phase 2 code from Analysis.R]

# Run Phase 3: Visualization
# [Phase 3 code from Analysis.R]
```

### Custom Processing

You can modify parameters in the main script:

```r
# Change downsampling interval
downsample_tracks(time_interval = 30)  # 30 minutes instead of 10

# Change trajectory gap threshold
assign_trajectories(max_gap_days = 7)  # 7 days instead of 5

# Change DBSCAN parameters
classify_staging_migration(eps = 0.2, minPts = 15)

# Change minimum trajectory length
remove_short_trajectories(min_points = 10)
```

---

## Output Summary Statistics

The script reports:
- Total records processed
- Number of unique individuals
- Number of unique species
- Records per species (top 20)
- Processing status for each individual
- Plot generation status

---

## Troubleshooting

### Common Issues

**Problem:** "File not found" errors  
**Solution:** Ensure all data files are in correct directory structure

**Problem:** Memory errors with large datasets  
**Solution:** Already handled via chunked processing; check available RAM

**Problem:** Missing height_above_msl data  
**Solution:** Automatically skipped; individuals without altitude data are excluded

**Problem:** No points after spatial filtering  
**Solution:** Normal - not all individuals cross Gulf of Lion; they are excluded

**Problem:** DBSCAN clustering errors  
**Solution:** Automatically handled; trajectories with < minPts get "Migration" state

---

## Citation

If you use this workflow, please cite the relevant data sources:

- Movebank data: [Include Movebank study citation]
- BIRDMAP: Adojaan K et al. [Include full citation]
- Montagu's Harrier: Schlaich et al. 2016, Journal of Applied Ecology
- [Other data sources as appropriate]

---

## Contact

**Victor Cat**  
[Contact information]

---

## License

[Include appropriate license information]

---

## Version History

- **v1.0** (December 2024): Initial workflow implementation
  - 6 data sources integrated
  - 11-step processing pipeline (Phases 1-2)
  - 6 summary visualizations (Phase 3)
  - Hierarchical Bayesian modeling (Phase 4, optional)
  - Ecological group classification