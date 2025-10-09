This contains code and data for manuscript titled "Epigenetic variation can promote adaptation by smoothing rugged fitness landscapes" - https://doi.org/10.1101/2025.06.06.658353

The scripts to generate fitness landscapes are - 

CreatingLandscapes_Main.R (to create a set of correlated fitness landscapes which form the main part of the paper), and CreatingLandscapes_HoC.R (following the house of cards model to create a set of uncorrelated landscapes)

We focus on landscapes with 5-loci in this work, but the number of loci can be easily changed by changing the parameter L in either script. The degree of ruggedness of the landscape (another critical feature in the work) can also be tuned by changing the strength of epistatic interactions in either script. 

We use the Wright-Fisher model to simulate adaptation on these landscapes. The code for our simulations is in WFmodel.R

All specific code to generate each figure are included in MainText_Code and Supplementary_Code, and the corresponding data files are included in MainText_Data and Supplementary_Data. 
