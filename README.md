# ALBClimateShift
## Overview
**ALBClimateShift** is a research-based project that explores the climate-driven shifts in the suitable areas of Alternaria leaf blotch (ALB) on apples in China. Leveraging multiple species distribution models (SDMs) and climate change scenarios, this repository provides the data and code used in our analysis to enhance the understanding of ALB dynamics under future climatic conditions.

This project is based on the study:
> "Climate-Driven Shifts in the Suitable Areas of Alternaria Leaf Blotch (Alternaria mali Roberts) on Apples: Projections and Uncertainty Analysis in China."

## Features
- **Species Distribution Modeling (SDM):** Implementation of five SDMs to analyze the relationship between environmental variables and ALB distribution.
- **Climate Data Projections:** Utilization of five Global Climate Models (GCMs) from CMIP6 under four Shared Socioeconomic Pathways (SSPs).
- **Uncertainty Analysis:** Quantification of prediction uncertainties across GCMs, SDMs, and scenarios.
- **Interactive Maps:** Visualization of current and future ALB-suitable regions under different climate scenarios.
- **Reproducible Analysis:** Well-documented scripts and datasets for replicating results.

---

## Getting Started

### Prerequisites
1. R version 4.4.2 
2. Main R Packages:
   - `sdm(1.2.52)`
   - `usdm(2.1.7)`
   - `raster(3.6.30)`
   - `tidyverse(2.0.0)`

### Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/chenbin98/ALBClimateShift.git
2.	Install the required R packages.

### Data Description

The dataset includes:
1. **ALB presence Records**: Collected from orchard surveys, public databases, and literature.
2. **Environmental Variables**: Climatic data (temperature, precipitation, etc.) from CMIP6. The topographic variables were derived from a digital elevation model (DEM)  from the Data Center of Resources and Environmental Sciences (http://www.resdc.cn/).

### Project Structure
```
ALBClimateShift/
├── data/                  # Raw and processed data
├── src/                   # Analysis scripts
│   ├── 01_data_preprocessing.R
│   ├── 02_model_training.R
│   ├── 03_projection_analysis.R
│   └── 04_visualization.py
├── figs/                  # Output figures
├── README.md              # Project documentation
└── LICENSE                # License file
```
### Results

Key findings:
1. ALB suitable areas are projected to shift northwestward (137–263 km) and to higher elevations (288–680 m) by 2090s under high-emission scenarios.
![Change in suitable areas of ALB](figs/climate_changing_compare_to_baseline.jpg)

2. The uncertainty in predictions is primarily driven by differences among GCMs (42.2%).
![](figs/uncertainty_analysis.jpg)

### License

This project is licensed under the MIT License. See the LICENSE file for details.

### Acknowledgments

This work was supported by:
	•	College of Soil and Water Conservation Science and Engineering, Northwest A&F University
	•	State Key Laboratory of Soil Erosion and Dryland Farming on the Loess Plateau
	•	China Agricultural University

For inquiries, contact:
	•	Prof. Qiang Yu: yuq@nwafu.edu.cn
	•	Prof. Gang Zhao: gang.zhao@nwafu.edu.cn

