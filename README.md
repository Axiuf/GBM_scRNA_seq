# GBM_scRNA_seq

This is my single cell RNA-seq analysis pipeline for the glioblastoma (GBM) project.

The file structure shows as below :

```
.
└── GBM_scRNA_seq
    ├── src # All source code files are stored here
    ├── plots # Figures out put from the project
    ├── data 
    │   └── Cellranger # Raw data from Cell ranger (10X)
    ├── .gitignore
    ├── GBM_scRNA_seq.Rproj
    ├── README.md
    └── requirements.R # The environment required by the project
```

- Subdirs in `Plots` and intermediate objects in `data` are accordant with the Rscripts in `src` from 
    which they are derived.
- Figures in `Plots` are named as `contents_to_plot_function_to_use.tiff`.
- When running Rscript in `src`, the workdir should be `/your/path/to/GBM_scRNA_seq`, 
    not `/your/path/to/GBM_scRNA_seq/src`.
- Seurat objects in all source files are named `GBM` for convenient of code reuse. Check which file is loaded
    to figure out the exact object being processed.