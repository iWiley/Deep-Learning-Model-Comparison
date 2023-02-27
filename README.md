## Introduction

This repository contains the source code and data of all programs used in the manuscript "Exploring the use of deep learning to assess immune infiltration in patients with hepatocellular carcinoma".

## How to get the original data？
The raw data of R scripts are stored in the /R/Data directory, and the deep learning training data used by C# programs are available on the release page. In addition, the tile data of the Xijing Hospital cohort used for validation is stored in the RELEASE page.

## Folder Structure

    DL-in-HCC
    ├── C#                                            C# related projects
    │   ├── DL-in-HCC.sln                               Solution file
    │   ├── Tile Classifier                             Function as folder name
    │   └── Wistu.Lib.ClassifyModel                     Public library, 'wistu' is the author's blog name
    ├── LICENSE
    ├── R                                             R related scripts
    │   ├── 00.Functions.R                              Basic functions of the script
    │   ├── 01.CalculateCutoff.R                        The name is its function
    │   ├── 02.DrawOS.R                                 The name is its function
    │   ├── 03.CoxRegression.R                          The name is its function
    │   ├── 04.Nomogram.R                               The name is its function
    │   ├── 05.CalibrationCurve.R                       The name is its function
    │   ├── 06.Time-dependentROCCurve.R                 The name is its function
    │   ├── 07.CIBERSORT.R                              The name is its function
    │   ├── 08.EnrichmentAnalysis.R                     The name is its function
    │   ├── 09.ImmuneCheckpointExpression.R             The name is its function
    │   ├── Data                                        Contains the data needed for the script
    │   └── R.Rproj                                     R project file
    └── README.md

## C# Projects

### Tile Classifier

The main role of this tool is to assist users in fast and efficient tile classification, while it can use already trained models to assist in classification.

#### Usage

Use keyboard shortcuts to sort the tiles. For more details, please refer to the source code.

## R Projects
The R scripts have been split into separate units according to their functions, and each R script runs interdependently but independently of each other. Note: If you cannot output the image by executing the corresponding R script directly, please press Ctrl+Enter in R studio to run the code line by line and get the image in the Plot window.

## Python Projects

Python项目主要用于模型的训练以及评估。
