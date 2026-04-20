# HNSCclassifier: An R Package to Predict Molecular Subtypes of Head and Neck Squamous Cell Carcinoma

**Install package**

```{r}
  # Required packages: run if not already installed
  if(!requireNamespace('BiocManager')){
    install.packages('BiocManager')
  }
  if(!requireNamespace('devtools')){
    install.packages('devtools')
  }

  ## Users need to install these packages before using GCclassifier
  BiocManager::install(c(
    'impute', 'dplyr', 'magrittr', 'randomForest', 
   'shiny', 'DT', 'shinyjs', 'BiocStyle', 'survminer'), force = T)
  
  ### install: latest version (R version >= 4.1.0 required)
  ### build_vignettes = T, if installing using RStudio
  devtools::install_github("Ronlee12355/HNSCclassifier", build_vignettes = T)
  
  ### if not installing from RStudio IDE, no vignette creating is recommended since it requires Pandoc and other dependancies
  devtools::install_github("Ronlee12355/HNSCclassifier", build_vignettes = F)
```
