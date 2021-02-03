# Estimating-under-reporting-of-leprosy-in-Brazil

## Table of contents
* [General info](#general-info)
* [Compilation](#compilation)
* [Input data](#data)
* [References](#references)

## General info
In this directory, we present descriptive analysis and the code to perform a Bayesian hierarchical modeling to estimate under reported cases of leprosy in Brazil.

The code originated from the results described in [1]. This was an interdisciplinary work aimed to apply statistical thecniches to better understand the dissemination of leprosy in Brazil. 

## Compilation
Our model was implemented using NIMBLE in conjunction with the software R [2], following the supplementary materials provided in [3]. The code is located at the "scripts" folder.

## Input data
On the folder "Data_and_results" the reader can find a CSV file containing all necessary variables to perform the analysis as well as the obtained results for each microregion. This file also contains the raw data used in the figures presented in the main manuscript [1].

## References 
[1] Oliveira GL, Oliveira JF, Pescarini JM, Andrade RF, Nery JS, Ichihara MY, Smeeth L, Brickley EB, Barreto ML, Penna GO, Penna ML. Estimating Under Reporting of Leprosy in Brazil using a Bayesian Approach. medRxiv. 2020 Jan 1.

[2] R Core Team (2015). R: A Language and Environment for Statistical Computing, R Foundation for Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, (2015). URL:https://www.R-project.org/

[3] Stoner, O; Economou, T; Drummond, G. (2019). A Hierarchical Framework for Correcting Under-Reporting in Count Data. Journal of the American Statistical Association , DOI 10.1080/01621459.2019.1573732.
