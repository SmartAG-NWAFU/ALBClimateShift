![fig5_climate_changing_for_baseline](https://github.com/user-attachments/assets/1db522c1-404f-46eb-9232-72baaff76690)# ALBClimateShift-Predicting-Alternaria-Leaf-Blotch-Distribution-Under-Climate-Change-in-China

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
2. R (>=4.1)
4. R Packages:
   - `biomod2`
   - `dismo`
   - `ggplot2`

### Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/ALBClimateShift.git
2.	Install the required R packages.

### Data Description

The dataset includes:
	1.	ALB Occurrence Records: Collected from orchard surveys, public databases, and literature.
	2.	Environmental Variables: Climatic data (temperature, precipitation, etc.) from CMIP6.
	3.	Baseline and Future Climate Scenarios: GCMs data for the baseline (1970–2000) and projections (2030s, 2050s, 2070s, 2090s) under SSP126, SSP245, SSP370, and SSP585.


### Project Structure

ALBClimateShift/
├── data/                      # Raw and processed data
├── scripts/                   # Analysis scripts
│   ├── 01_data_preprocessing.py
│   ├── 02_model_training.R
│   ├── 03_projection_analysis.R
│   └── 04_visualization.py
├── results/                   # Output files and figures
├── notebooks/                 # Jupyter/RMarkdown notebooks for exploratory analysis
├── README.md                  # Project documentation
└── LICENSE                    # License file


### Results

Key findings:
	•	ALB suitable areas are projected to shift northwestward (137–263 km) and to higher elevations (288–680 m) by 2090s under high-emission scenarios.

 
	•	The uncertainty in predictions is primarily driven by differences among GCMs (42.2%).

 ### Acknowledgments

This work was supported by:
	•	College of Soil and Water Conservation Science and Engineering, Northwest A&F University
	•	State Key Laboratory of Soil Erosion and Dryland Farming on the Loess Plateau
	•	China Agricultural University

For inquiries, contact:
	•	Prof. Qiang Yu: yuq@nwafu.edu.cn
	•	Prof. Gang Zhao: gang.zhao@nwafu.edu.cn

