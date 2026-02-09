# Nextflow Single-Cell RNA-seq Analysis Pipeline

**Version:** 1.0  
**Last Updated:** February 8, 2026  
**Author:** Pia Lagler

---

## 📋 Overview

This Nextflow pipeline provides a comprehensive workflow for single-cell RNA-seq analysis, from raw FASTQ files to biological insights. The pipeline is optimized for high-performance computing (HPC) environments and uses conda environments for reproducible analysis.

### Pipeline Steps

| Step | Script | Description | Output |
|------|--------|-------------|--------|
| 00 | STARsolo | Read alignment and quantification | Count matrices |
| 01 | QC & Integration | Quality control, SoupX, doublet removal, STACAS integration | Integrated Seurat object |
| 02 | Clustering & Annotation | PCA, UMAP, clustering, SingleR annotation | Annotated clusters |
| 03 | Differential Analysis | Differential abundance (MiloR), differential state (Seurat + Muscat) | DE genes, DA regions |
| 04a | CellChat Analysis | Cell-cell communication inference, comparison mode | Interaction networks, merged CellChat objects |
| 04b | NicheNet Analysis | Ligand-receptor prioritization, target prediction | Ligand activity rankings, target networks |
| 05 | Gene Programs | GeneNMF meta-programs, GSEA enrichment | Gene program activities |
| 06 | Trajectory Analysis | Slingshot trajectory inference, tradeSeq DE | Pseudotime trajectories |
| 07 | Benchmarking | Performance metrics, memory profiling | QC metrics |
| 08 | Dashboard | Interactive HTML reports | Visualization dashboard |

---

## Start

### 1. Clone Repository
```bash
git clone <repository-url>
cd nextflow_test
```

### 2. Install Conda Environments
```bash
cd envs

# Create each environment
conda env create -f seurat_env_3_soupx.yml -n seurat_env_3_soupx
conda env create -f annotation_singler2.yml -n annotation_singler2
conda env create -f cellcomm_env.yml -n cellcomm_env
conda env create -f cellcomm_genenmf.yml -n cellcomm_genenmf
conda env create -f trajectory_analysis.yml -n trajectory_analysis
conda env create -f benchmarking_env.yml -n benchmarking_env
conda env create -f starsolo_env.yml -n starsolo_env

# Install GitHub packages (see Environment Setup section)
```

### 3. Configure Pipeline
Edit `nextflow.config` and update paths:
- Input data directory
- Output directory  
- Conda environment paths
- Reference genome paths

### 4. Run Pipeline
```bash
# Full pipeline from FASTQ
nextflow run main.nf

# Alternative workflows - start from specific steps:
nextflow run main.nf -entry from_qc          # Skip STARsolo, start from QC
nextflow run main.nf -entry from_clustering  # Start from clustering
nextflow run main.nf -entry from_differential # Start from differential analysis
nextflow run main.nf -entry from_cellcomm    # Start from cell communication
nextflow run main.nf -entry from_nichenet    # Only run NicheNet (needs CellChat results)
nextflow run main.nf -entry from_geneprog    # Start from gene programs
nextflow run main.nf -entry from_trajectory  # Start from trajectory
nextflow run main.nf -entry from_benchmark   # Only run benchmarking
```

---

## 📦 Environment Setup

### Required Conda Environments

The pipeline uses **8 conda environments**, each optimized for specific analysis steps:

| Environment | Purpose | R Version | Key Packages |
|-------------|---------|-----------|------------|
| `starsolo_env` | Read alignment | Python 3.9 | STAR 2.7.10b |
| `seurat_env_3_soupx` | QC, integration | 4.1.3 | Seurat 4.3.0, SoupX, STACAS |
| `annotation_singler2` | Cell annotation | 4.3.3 | Seurat 5.0.3, SingleR |
| `cellcomm2_env` | CellChat analysis | 4.3.3 | CellChat 1.6.1, patchwork |
| `cellcomm_env` | NicheNet analysis | 4.3.3 | nichenetr, multinichenetr |
| `cellcomm_genenmf` | Gene programs | 4.3.3 | GeneNMF, msigdbr |
| `trajectory_analysis` | Trajectory inference | 4.4.3 | Slingshot, tradeSeq |
| `benchmarking_env` | Performance metrics | 4.1.3 | microbenchmark |

### Installation

#### Step 1: Create Conda Environments

```bash
cd envs

# Install all environments
for env_file in *.yml; do
  env_name=$(basename $env_file .yml)
  echo "Creating environment: $env_name"
  mamba env create -f $env_file -n $env_name
done
```

Or create each environment individually:

```bash
mamba env create -f starsolo_env.yml -n starsolo_env
mamba env create -f seurat_env_3_soupx.yml -n seurat_env_3_soupx
mamba env create -f annotation_singler2.yml -n annotation_singler2
mamba env create -f cellcomm2_env.yml -n cellcomm2_env
mamba env create -f cellcomm_env.yml -n cellcomm_env
mamba env create -f cellcomm_genenmf.yml -n cellcomm_genenmf
mamba env create -f trajectory_analysis.yml -n trajectory_analysis
mamba env create -f benchmarking_env.yml -n benchmarking_env
```

#### Step 2: Install GitHub Packages

Some packages must be installed manually from GitHub:

```bash
# For seurat_env_3_soupx (QC & integration)
conda activate seurat_env_3_soupx
R -e 'devtools::install_github("carmonalab/STACAS", upgrade = "never")'
R -e 'BiocManager::install("miloR", update = FALSE, ask = FALSE)'

# For cellcomm2_env (CellChat)
conda activate cellcomm2_env
R -e 'devtools::install_github("sqjin/CellChat", ref = "v1.6.1", upgrade = "never")'

# For cellcomm_env (NicheNet)
conda activate cellcomm_env
R -e 'devtools::install_github("saeyslab/nichenetr", ref = "v2.2.1.1", upgrade = "never")'
R -e 'devtools::install_github("saeyslab/multinichenetr", upgrade = "never")'
R -e 'install.packages("umap", repos = "https://cloud.r-project.org/")'

# For cellcomm_genenmf 
conda activate cellcomm_genenmf
R -e 'install.packages("GeneNMF", repos = "https://cloud.r-project.org/")'
```


#### Step 3: Update nextflow.config

Find your conda environment paths:
```bash
conda info --envs
```

Edit `nextflow.config` and update:
```groovy
params {
  conda_starsolo_env = "/path/to/.conda/envs/starsolo_env"
  conda_seurat_env = "/path/to/.conda/envs/seurat_env_3_soupx"
  conda_annotation_env = "/path/to/.conda/envs/annotation_singler2"
  conda_cellchat_v2 = "/path/to/.conda/envs/cellcomm2_env"
  conda_cellcomm_multi = "/path/to/.conda/envs/cellcomm_env"
  conda_nmf_env = "/path/to/.conda/envs/cellcomm_genenmf"
  conda_trajectory_env = "/path/to/.conda/envs/trajectory_analysis"
  conda_benchmarking_env = "/path/to/.conda/envs/benchmarking_env"
}
```

---

## 📂 Directory Structure

```
nextflow_test/
├── main.nf
├── nextflow.config
├── README.md
├── envs/
│   ├── *.yml (7 environments)
│   └── README.md
├── Scripts/
│   ├── 01_qc_integration.R
│   ├── 02_clustering_annotation.R
│   ├── 03_differential_analysis.R
│   ├── 04a_cellchat.R
│   ├── 04b_nichenet.R
│   ├── 05_gene_program_discovery.R
│   ├── 06_trajectory.R
│   ├── 07_benchmarking.R
│   ├── 08_dashboard.R
│   └── utils/
├── ref_mouse/
└── outputs/
    ├── 00_starsolo/
    ├── 01_qc_integration/
    ├── 02_clustering/
    ├── 03_differential/
    ├── 04_cellcomm/
    ├── 05_gene_programs/
    ├── 06_trajectory/
    ├── 07_benchmarking/
    └── 08_reports/
```

---

## ⚙️ Pipeline Configuration

### Key Parameters in `nextflow.config`

#### Input/Output
```groovy
params {
  fastq_dir = "data/"                 // Input FASTQ directory
  output_dir = "outputs/"             // Output directory
}
```

#### Reference Genome
```groovy
params {
  genome_dir = "ref_mouse/STAR_index"
  whitelist = "ref_mouse/3M-february-2018.txt"
  gtf_file = "ref_mouse/genes.gtf"
}
```

#### Conda Environments
```groovy
params {
  conda_seurat_env = "/path/to/seurat_env_3_soupx"
  conda_annotation_env = "/path/to/annotation_singler2"
  // ... (see Environment Setup section)
}
```

#### Analysis Parameters
```groovy
params {
  qc {
    min_features_per_cell = 200
    n_variable_features = 3000
    integration_method = "STACAS"
    random_seed = 42
  }
  
  clustering {
    resolution = [0.1, 0.25, 0.5, 0.75, 1.0]
    n_neighbors = 30
    n_pcs = 30
  }
  
  differential {
    min_pct = 0.1
    logfc_threshold = 0.25
  }
}
```

---

## 📊 Interactive Dashboard

The pipeline generates a comprehensive HTML dashboard (`08_reports/comprehensive_analysis_dashboard.html`) that aggregates all results in an interactive, tab-based interface.

### Features
- **Tab-based navigation** through all analysis modules
- **Before/After comparisons** for QC steps (side-by-side)
- **All plot types supported:** HTML (interactive), PNG, PDF
- **Self-contained:** No internet connection required after generation
- **Portable:** Download the `08_reports/` folder and open anywhere

### Dashboard Sections
1. **Key Results** - Highlights from your analysis
2. **QC & Integration** - Quality metrics, ambient RNA, doublet removal
3. **Clustering & Annotation** - PCA, UMAP, cell type assignments
4. **Differential Analysis** - DA (MiloR), DS (Seurat/Muscat), volcano plots
5. **Cell Communication** - CellChat networks, NicheNet predictions
6. **Gene Programs** - Metaprograms, enrichment analysis
7. **Trajectory Analysis** - Slingshot pseudotime, GAM plots
8. **Benchmarking** - Performance metrics, memory usage

### Usage
```bash
# Open dashboard after pipeline completes
firefox outputs/08_reports/comprehensive_analysis_dashboard.html

# The entire 08_reports/ folder is portable
zip -r dashboard.zip outputs/08_reports/
# Share dashboard.zip with collaborators
```

---

## 🔗 References

### Software
- **Nextflow:** https://www.nextflow.io/
- **Seurat:** https://satijalab.org/seurat/
- **CellChat:** https://github.com/sqjin/CellChat
- **NicheNet:** https://github.com/saeyslab/nichenetr
- **STACAS:** https://github.com/carmonalab/STACAS
- **MiloR:** https://bioconductor.org/packages/miloR
- **Slingshot:** https://bioconductor.org/packages/slingshot



## 👤 Author

**Pia Lagler**  

---

**Last Updated:** February 9, 2026
