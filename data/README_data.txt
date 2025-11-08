Oxygen Subduction Project - Data Description
Last updated: 2025/11/8 
Contact: zhan_xyue00@sjtu.edu.cn / ytzhou@sjtu.edu.cn

This document describes the data files used in the oxygen subduction analysis project.
All data are organized by category for better reproducibility.

===GRID DATA===
Files containing spatial grid information and domain definitions.

[1] basin_mask_pacific.mat
- Variables: BASIN、BASIN1、basin_str
- Description: Pacific Ocean basin mask defining regional boundaries
- Usage: Used for regional analysis and masking calculations to Pacific Ocean


[2] grid_area_pacific.mat  
- Variables: Areaall
- Description: Grid cell surface areas for Pacific Ocean
- Units: km²
- Usage: Used in flux calculations and regional integration


[3] grid_coordinates.mat
- Variables: oceanlon, oceanlat
- Description: Longitude and latitude coordinates for ECCO model grid
- Units: degrees
- Spatial resolution: 0.5 degree
- Domain: Global
- Usage: Spatial referencing for all calculations


===INPUT DATA===
Raw or primary input data used in calculations.

[1] oxygen_concentration.mat
- Variables: oxy, oxy10
- Description: Oxygen concentrations at the base of the permanent-maximum mixed layer and 10 m below it [μmol/kg]
- Source: GOBAI-O₂ reanalysis data
- Usage: Primary tracer fields for oxygen subduction calculations

[2] physical_variables.mat
- Variables: Hml, Hst, Hmax, evel1, nvel1, evel2, nvel2, wb, estar1, estar2, nstar1, nstar2, taue, taun, rhob10, rho
- Description: A preprocessed physical forcing fields
- Variables breakdown:
  - Hml: Thickness of the mixed layer [m]
  - Hst: Thickness of the seasonal thermocline [m] 
  - Hmax: Thickness of the permanent-maximum mixed layer during 2004-2017 [m]
  - evel1, nvel1: Velocity components (eastward, northward) in the mixed layer [m/s]
  - evel2, nvel2: Velocity components (eastward, northward) in the seasonal thermocline [m/s]
  - wb: Vertical velocity at the base of the permanent-maximum mixed layer [m/s]
  - estar1, nstar1: Eddy-induced velocity components (eastward, northward) in the mixed layer [m/s]
  - estar2, nstar2: Additional eddy components (eastward, northward) in the seasonal thermocline [m/s]
  - taue, taun: Wind stress components [N/m²]
  - rho, rhob10: Potential density at the base of the permanent-maximum mixed layer and 10 m below it  [kg/m³]
- Source: ECCOv4r4 model output
- Usage: Physical drivers for subduction calculations

[3] velocity_eastward_350m.mat
- Variables: EVEL, dep
- Description: Eastward velocity in upper 350m and depth levels
- Units: EVEL - m/s; dep - m
- Source: ECCOv4r4 model output
- Usage: Specifically for Equatorial Undercurrent (EUC) transport calculations

===PROCESSED DATA===
Intermediate calculated fields and derived quantities.

[1] transport_euc.mat
- Variables: Trspt, long
- Description: Equatorial Undercurrent transport time series
- Units: Trspt - Sverdrups (Sv, 10⁶ m³/s); long - degrees
- Calculated from: velocity_eastward_350m.mat using calculate_euc_transport.m
- Usage: Analysis of EUC variability and its role in oxygen transport
- Method: Integrated eastward flow between ±2.5° latitude and 0-350m depth

[2] wind_stress_curl.mat
- Variables: WC
- Description: Wind stress curl field
- Units: N/m³
- Calculated from: taue, taun in physical_variables.mat using calculate_wind_stress_curl.m
- Usage: Analysis of wind-driven ocean circulation
- Method: ∇ × τ = ∂τ_n/∂x - ∂τ_e/∂y

[3] wind_stress_curl_zonal_climatology.mat  
- Variables: WCm
- Description: Zonally-averaged wind stress curl climatology
- Units: N/m³
- Calculated from: wind_stress_curl.mat using calculate_wind_stress_curl.m
- Usage: Analysis of seasonal wind patterns
- Method: Zonal mean of monthly climatology over Pacific basin

===RESULTS DATA===
Final analysis outputs and scientific results.

[1] oxygen_subduction_results.mat
- Variables: OSLavoX, OSLavoY, OSLavo, OSwb, OSeddy
- Description: Oxygen subduction rate components
- Units: mol O₂ m⁻² month⁻¹
- Components:
  - OSLavoX: Lateral induction (east-west direction)
  - OSLavoY: Lateral induction (north-south direction) 
  - OSLavo: Total lateral induction
  - OSwb: Vertical velocity contribution
  - OSeddy: Eddy-induced contribution
- Calculated from: oxygen_concentration.mat + physical_variables.mat using calculate_oxygen_subduction.m
- Usage: Quantification of oxygen subduction mechanisms
- Method: Upwind scheme applied to water mass subduction rates

[2] analysis_eddy_regions.mat
- Variables: BASIN_ed, EddyRegion, EddyLine, OSbar_ed, CpnPrp_ed
- Description: Eddy-dominated region analysis results
- Variables breakdown:
  - BASIN_ed: Mask of eddy-dominated grid cells
  - EddyRegion: Region names ('Kuroshio Extension', 'North Equatorial Pacific', 'Coastal Peru')
  - EddyLine: Polygon boundaries for each region
  - OSbar_ed: Monthly mean oxygen subduction by term and region [Tmol/month]
  - CpnPrp_ed: Relative percentage contribution of each term [%]
- Calculated from: oxygen_subduction_results.mat using analyze_eddy_regions.m
- Usage: Regional analysis of eddy contributions to oxygen subduction


===Citation===
When using this data, please cite: ***
