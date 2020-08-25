# GBM_scRNA_seq

This is my single cell RNA-seq analysis pipline for the glioblastoma (GBM) project.

The file structure shows as below :

```
.
└── GBM_scRNA_seq
    ├── src # All source code files are stored here.
    ├── plots # Figures out put from the project.
    ├── data 
    │   └── Cellranger # Raw data from Cell ranger (10X).
    ├── .gitignore
    ├── GBM_scRNA_seq.Rproj
    ├── README.md
    └── requirements.R # The environment required by the project.
```

- Subdirs in `Plots` are accordant with the Rscripts in `src` from which they are derived.
- Figures in `Plots` are named as `contents_to_plot_function_to_use.tiff`.
- When running Rscript in `src`, the work dir should be `/your/path/to/GBM_scRNA_seq`, 
    not `/your/path/to/GBM_scRNA_seq/src`.