Pacific Ocean Oxygen Subduction Analysis - Project Repository

Project Overview:
This repository contains the complete code and data for the study "Wind–Stress Curl and Eddies Control the Seasonal Oxygen Subduction and Obduction in the Pacific". The project analyzes oxygen subduction (Sox) processes in the Pacific Ocean using ECCO model data and GOBAI-O2 reanalysis data.

Last updated: 2025/11/8 
Contact: zhan_xyue00@sjtu.edu.cn / ytzhou@sjtu.edu.cn


==== Repository Structure ====

### data/ - Research Data ###
Comprehensive dataset for oxygen subduction analysis. See 'README_data.txt' for detailed variable descriptions, units, and processing workflows.
上传在zenodo

- ** grid/ ** - Spatial reference data
  - basin_mask_pacific.mat - Pacific Ocean basin definitions
  - grid_area_pacific.mat - Grid cell surface areas  
  - grid_coordinates.mat - Longitude/Latitude coordinates

- ** input/ ** - Primary input data
  - oxygen_concentration.mat - Oxygen concentration fields (preprocessed from GOBAI-O2)
  - physical_variables.mat - Comprehensive physical forcing fields (preprocessed from ECCOv4r4)
  - velocity_eastward_350m.mat - Eastward velocity for Equatorial Undercurrent calculations (preprocessed from ECCOv4r4)

- ** processed/ ** - Derived physical quantities
  - transport_euc.mat - Equatorial Undercurrent transport
  - wind_stress_curl.mat - Wind stress curl fields
  - wind_stress_curl_zonal_climatology.mat - Zonal mean wind stress curl

- ** results/ ** - Analysis outputs
  - oxygen_subduction_results.mat - Sox component calculations
  - analysis_eddy_regions.mat - Eddy-dominated region analyses

### figures/ - Plotting Scripts for Manuscript Figures ###
MATLAB scripts generating all figures for the publication.

- ** FIG1_cde.m ** - Figure 1c,d,e: Spatial distributions
  - Climatological mean Sox
  - Seasonal standard deviation  
  - Seasonal correlation (r²)

- ** FIG1_fgh.m ** - Figure 1f,g,h: Latitudinal distributions
  - Subduction area fraction (DJF vs JJA)
  - Sox and Oox rates by latitude
  - Zonally integrated Sox

- ** FIG2.m ** - Figure 2: Equatorial region analysis
  - Seasonal Sox components (170E-80W, ±2.5°)
  - Hovmöller diagrams of Sox and vertical velocity
  - ENSO-Sox relationship analysis

- ** FIG3.m ** - Figure 3: Zonal mean analysis
  - Sox and vertical velocity climatology/anomalies
  - Wind stress curl relationships
  - Seasonal region identification

- ** FIG4.m ** - Figure 4: Eddy-dominated regions
  - Relative percentage of eddy contribution map
  - Seasonal cycles of Sox components
  - Three case study regions analysis

### scripts/ - Core Analysis Scripts ###
Fundamental calculation scripts used across the analysis.

- ** analyze_eddy_dominated_regions.m ** - Identifies and analyzes eddy-dominated regions
- ** calculate_EUC_transport.m ** - Computes Equatorial Undercurrent transport
- ** calculate_oxygen_subduction.m ** - Main Sox calculation with 4 components
- ** calculate_wind_stress_curl.m ** - Wind stress curl calculations

Numerical Core Functions:
- ** curldffs.m ** - 2D curl calculations using gradient methods
- ** distcoordinate.m ** - Converts lat/lon to distance coordinates
- ** divdffs.m ** - 2D divergence calculations
- ** sw_dist.m ** - Seawater toolbox distance calculations

### toolbox/ - External Dependencies ### 
Third-party toolboxes required for analysis and visualization.

- ** m_map1.4.tar/ ** - Mapping toolbox for MATLAB
- ** Mapcolors/ ** - Custom colormap functions

