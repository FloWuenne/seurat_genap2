# seurat_genap2
Implementation of Seurat pipeline for Galaxy in GenAP2. 

Code is built upon the seurat package in the Galaxy tool shed: https://github.com/galaxyproject/tools-iuc/blob/master/tools/seurat/

All credit to the basic implementation of the Seurat code goes to the authors of the original package.
The following features were added to this version:
  - Support of Salmon - Alevin as input via internal conversion to cell x gene matrix (ongoing)
  - Dimensional reduction using UMAP instead of tSNE (ongoing)
  - Export of data as single RDS object which is compatible with downstream shiny visualisation
