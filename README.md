# This Fork...
## Feature
- Support parallelized `ForestIV`  
## Usage

- Local dev dir  
    ```R
    local_folder <- "/Users/rey/Developer/EnsembleIV/"
    devtools::install(local_folder, dependencies = TRUE)
    ```

- Install remotely  
    ```R
    remote_repo <- "https://github.com/reycn/EnsembleIV-Parallel.git"
    devtools::install_github(remote_repo, dependencies = TRUE)
    library(EnsembleIV)
    ```
## Speed Comparison
TBD...
# EnsembleIV and ForestIV
An R package that implements [EnsembleIV](https://arxiv.org/abs/2303.02820) and [ForestIV](https://pubsonline.informs.org/doi/abs/10.1287/ijds.2022.0019).

To install, run 
```r
devtools::install_github("mochenyang/EnsembleIV")
```
