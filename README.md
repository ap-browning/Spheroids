# Spheroids

Repository for the preprint "Quantitative analysis of tumour spheroid structure" available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.08.05.455334v1).

## Figures
 - [Figure 1](https://htmlpreview.github.io/?https://github.com/ap-browning/Spheroids/blob/main/Figures/Fig1/fig1.html)
 - [Figure 3](https://htmlpreview.github.io/?https://github.com/ap-browning/Spheroids/blob/main/Figures/Fig3/fig3.html)
 - [Figure 4](https://htmlpreview.github.io/?https://github.com/ap-browning/Spheroids/blob/main/Figures/Fig4/fig4.html)
 - [Figure 5](https://htmlpreview.github.io/?https://github.com/ap-browning/Spheroids/blob/main/Figures/Fig5/fig5.html)
 - [Figure 6](https://htmlpreview.github.io/?https://github.com/ap-browning/Spheroids/blob/main/Figures/Fig6/fig6.html)
 - [Figure S7](https://htmlpreview.github.io/?https://github.com/ap-browning/Spheroids/blob/main/Figures/FigS7/figS7.html)
 - [Figure S8](https://htmlpreview.github.io/?https://github.com/ap-browning/Spheroids/blob/main/Figures/FigS8/figS8.html)
 - [Figure S9](https://htmlpreview.github.io/?https://github.com/ap-browning/Spheroids/blob/main/Figures/FigS9/figS9.html)

## Data

All data used in the analysis is available in `Data`. Data collected from the Incucyte S3 (Satorius) (i.e., outer radius measurements only) is available in `IncucyteData.csv`. Data summarised from confocal microscopy images (find the repository containing image processing codes here) is available in `AllConfocalData.csv`.

The script `SubsetData.jl` randomly subsets the complete confocal data set in `Data/AllConfocalData.csv` so that the sample size for each condition is approximately the same (a random seed is used to ensure reproducability). 

The datasets contain the following variables
| Variable          | Values                    | Description  |
| ----------------- |---------------------------| -------------|
| CellLine          | `983b` or `793b`          | Cell line
| InitialCondition  | `2500`, `5000`, `10000`   | Initial number of cells in each spheroid.
| Day               | `Number`                  | Days elapsed between spheroid seeding and measurement.
| R                 | `Number`                  | Spheroid outer radius.
| φ                 | `Number` ∈ [0,1]          | Proportion of outer radius containing inhibited or dead cells.
| η                 | `Number` ∈ [0,1]          | Proportion of outer radius void of living cells.


## Running the code
 
### Installation

Code to produce all results is contained within this repository, including the Julia modules `Greenspan` and `Inference`. To download, first clone this Github repo. Next, add the two module folders to your `LOAD_PATH`:
```
push!(LOAD_PATH,"/path/to/repo/Module/Greenspan")
push!(LOAD_PATH,"/path/to/repo/Module/Inference")
```
Next, run the following code (press `]` to enter the `Pkg` REPL) to install all required packages and activate the project
```
(v1.6) pkg> activate .
(Spheroids) pkg> instantiate
```

### Results

Code used to produce each figure in the main document and supplementary material is available in the `Figures` folder. For example, to reproduce Figure 4, run `Figures/Figure 4/fig4.jl`. This will create a plot stored in the variable `fig4` that contains Figure 4. 

Note that most figures contain a version saved as `html`. While these figures may differ stylistically from those in the main document, opening these in a web browser allows for interaction with the figure.

## Citation
If you use the data or software in this repository in your own work, please cite the following:

AP Browning, JA Sharp, RJ Murphy, G Gunasingh, B Lawson, K Burrage, NK Haass, MJ Simpson. 2021 Quantitative analysis of tumour spheroid structure. _bioRxiv_ [doi:test](https://dx.doi.org/test)

```
@article{Browning.2021gh,
	title        = {Quantitative analysis of tumour spheroid structure},
	author       = {Browning, Alexander P and Sharp, Jesse A and Murphy, Ryan J and Gunasingh, Gency and Lawson, Brodie and Burrage, Kevin and Haass, Nikolas K and Simpson, Matthew J},
	year         = 2021,
	journal      = {bioRxiv},
	doi          = {10.1101/2021.08.05.455334}
}
```
