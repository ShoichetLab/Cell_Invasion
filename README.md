# Cell_Invasion
This project provides a script to determine depth of invasion of cells into a 3D hydrogel based on a surface parabola

# Setup and Installation
Download all .m files and save in the same file directory

# How to Use
1. Imaging
Objective: 10x
Stains: - Nuclei stain
        - Fluorescent Beads for the surface 
Z-stack, recommeneded imaging increments: 20 um

2. Imaris pre-processing (v.
- Import images to Imaris
- Open representative image
- Create spots for cells
    - "Add new spots"
    - Source Channel: Nuclei Stain (e.g. Channel 1 - Hoeschst33258)
    - Adjust filter to capture entire cell population
    - Execute
- Repeat to create spots for surface beads
    - Source Channel: Beads (e.g. Channel 2 - Alexa Fluor 568)
- Select the spots
- Select the "Creation" Tab
- "Store Parameters for Batch"
- Back in the Arena Tab, Right click the spot creation tool and click "Run Batch"
- Right click generated statistics and click "Export Statistics". Save as .xslx
- Save Beads and Cells in their own individual folders for use in MATLAB analysis (example below)

3. Matlab Analysis
- Open InvasionScript.m 
- Run the script
    - A module should open
    - Select file folder containing beads
    - Select file folder containing cells
    - Select desired threshold (typically 30-75 um, selected based on non-invading control)
    - Optional: Select folder to save generated MATLAB images
    - Hit 'Run'
- An excel file containing average invasion depth and the percent of invading cells / sample is generated and can be used for plotting


# Sample Folder Structure
|> Home
|  > User
|    > Experiment 1
|      > Cells
|        - CellStatistics.xlsx
|      > Beads
|        - BeadsStatistics.xlsx
