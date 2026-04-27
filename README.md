# 🧬 Spatial Transcriptomics Analysis Tutorials

A comprehensive collection of four spatial transcriptomics tutorials using **Scanpy** and **Squidpy**, covering Visium (H&E and Fluorescence), MERFISH, and Xenium platforms.

---

## 📋 Table of Contents

1. [Tutorial 1 — Scanpy: Basic Spatial Analysis (Visium + MERFISH)](#tutorial-1--scanpy-basic-spatial-analysis-visium--merfish)
2. [Tutorial 2 — Squidpy: Visium Fluorescence Analysis](#tutorial-2--squidpy-visium-fluorescence-analysis)
3. [Tutorial 3 — Squidpy: Visium H&E Analysis](#tutorial-3--squidpy-visium-he-analysis)
4. [Tutorial 4 — Squidpy: Xenium Analysis](#tutorial-4--squidpy-xenium-analysis)
5. [Environment Setup](#environment-setup)
6. [Quick Reference: Key Functions](#quick-reference-key-functions)

---

## Tutorial 1 — Scanpy: Basic Spatial Analysis (Visium + MERFISH)

**Source:** [scanpy-tutorials.readthedocs.io](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html)  
**Author:** Giovanni Palla  
**Platform:** 10x Genomics Visium + MERFISH

### Overview

Introduces the fundamentals of working with spatial transcriptomics data in Scanpy. Covers loading, quality control, normalization, clustering, and spatial visualization — all on a real human lymph node Visium dataset. Ends with a MERFISH demonstration.

---

### Section 1.1 — Reading the Data

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sc.datasets.visium_sge()` | scanpy | Download and load 10x Visium dataset |
| `adata.var_names_make_unique()` | anndata | Ensure unique variable (gene) names |
| `sc.pp.calculate_qc_metrics()` | scanpy | Compute standard QC statistics per spot |

#### Inputs
- `sample_id = "V1_Human_Lymph_Node"` — Sample identifier string for the 10x public dataset
- Downloaded automatically from the 10x Genomics portal

#### Outputs
- `adata` — AnnData object with:
  - **4035 spots × 36601 genes**
  - `obs`: `in_tissue`, `array_row`, `array_col`, `total_counts`, `n_genes_by_counts`, `pct_counts_mt`, etc.
  - `var`: `gene_ids`, `feature_types`, `mt` (boolean for mitochondrial genes)
  - `uns['spatial']` — tissue image and scale factors
  - `obsm['spatial']` — 2D spatial coordinates for each spot

---

### Section 1.2 — QC and Preprocessing

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sns.histplot()` | seaborn | Visualize count and gene distributions |
| `sc.pp.filter_cells()` | scanpy | Remove spots with too few/many counts |
| `sc.pp.filter_genes()` | scanpy | Remove genes detected in very few spots |
| `sc.pp.normalize_total()` | scanpy | Normalize total counts per spot |
| `sc.pp.log1p()` | scanpy | Log-transform normalized counts |
| `sc.pp.highly_variable_genes()` | scanpy | Select top 2000 variable genes |

#### Inputs
- Raw AnnData object from Section 1.1
- Filter thresholds: `min_counts=5000`, `max_counts=35000`, `pct_counts_mt < 20`, `min_cells=10`
- `n_top_genes=2000`, `flavor="seurat"` for HVG selection

#### Outputs
- Filtered dataset: **3861 spots × ~19685 genes** (after filtering)
- `adata.var['highly_variable']` — boolean column marking top 2000 HVGs
- `adata.var['means']`, `adata.var['dispersions']` — gene-level statistics
- 4-panel QC histogram plot of count/gene distributions

#### Key Result
> Spots with fewer than 5000 or more than 35000 total counts and those with >20% mitochondrial reads are removed. Genes detected in fewer than 10 spots are also excluded.

---

### Section 1.3 — Manifold Embedding and Clustering

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sc.pp.pca()` | scanpy | Principal Component Analysis (dimensionality reduction) |
| `sc.pp.neighbors()` | scanpy | Build k-nearest neighbors graph on PCA space |
| `sc.tl.umap()` | scanpy | Non-linear 2D embedding of neighbor graph |
| `sc.tl.leiden()` | scanpy | Community detection clustering algorithm |
| `sc.pl.umap()` | scanpy | Plot UMAP colored by cluster/covariates |

#### Inputs
- Preprocessed, normalized, log-transformed AnnData
- `n_comps=50` for PCA
- Leiden parameters: `flavor="igraph"`, `directed=False`, `n_iterations=2`

#### Outputs
- `adata.obsm['X_pca']` — 50-dimensional PCA embedding
- `adata.uns['neighbors']` — neighbor graph metadata
- `adata.obsm['X_umap']` — 2D UMAP coordinates
- `adata.obs['clusters']` — **10 Leiden clusters** per spot
- UMAP plot colored by `total_counts`, `n_genes_by_counts`, and `clusters`

---

### Section 1.4 — Visualization in Spatial Coordinates

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sc.pl.spatial()` | scanpy | Overlay spot data onto H&E tissue image |

#### Key Parameters of `sc.pl.spatial()`
| Parameter | Description |
|-----------|-------------|
| `img_key` | Key for image in `adata.uns` (e.g., `"hires"`) |
| `color` | Variable(s) to color spots by |
| `crop_coord` | `[left, right, top, bottom]` for zoomed view |
| `alpha_img` | Transparency of background image |
| `alpha` | Transparency of spot overlay |
| `size` | Scaling factor for spot sizes |
| `groups` | Subset clusters to highlight |

#### Inputs
- AnnData with spatial coordinates and cluster labels
- H&E high-resolution image stored in `adata.uns['spatial']`

#### Outputs
- Spatial scatter plots with spots colored by:
  - `total_counts` and `n_genes_by_counts` (continuous)
  - `clusters` (10 Leiden clusters overlaid on tissue)
  - Zoomed inset showing clusters 5 and 9 with alpha blending

#### Key Result
> Spots in the same transcriptional cluster frequently co-localize spatially, revealing tissue organization (e.g., cluster 5 surrounded by cluster 0).

---

### Section 1.5 — Cluster Marker Genes

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sc.tl.rank_genes_groups()` | scanpy | Differential expression — find cluster markers |
| `sc.pl.rank_genes_groups_heatmap()` | scanpy | Heatmap of top marker gene expression |
| `sc.pl.spatial()` | scanpy | Spatial expression plot for marker genes |

#### Inputs
- AnnData with cluster labels
- `method="t-test"`, `groupby="clusters"`
- `n_genes=10` — top 10 markers per cluster

#### Outputs
- `adata.uns['rank_genes_groups']` — ranked marker gene lists with scores, log fold-changes, adjusted p-values
- Heatmap of top 10 markers for cluster 9 across all clusters
- Spatial plots of `CR2`, `COL1A2`, `SYPL1` expression

#### Key Result
> Gene *CR2* recapitulates the spatial structure of cluster 9. Genes like *COL1A2* and *SYPL1* show distinct spatial expression patterns corresponding to tissue regions.

---

### Section 1.6 — MERFISH Example

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `pd.read_excel()` | pandas | Load cell coordinate table |
| `sc.read_csv()` | scanpy | Load gene expression count matrix |
| `sc.pp.normalize_per_cell()` | scanpy | Normalize to 1e6 counts per cell |
| `sc.pl.embedding()` | scanpy | Plot cells in spatial coordinate system |

#### Inputs
- `pnas.1912459116.sd15.xlsx` — cell spatial coordinates
- `pnas.1912459116.sd12.csv` — gene expression count matrix
- From publication: Xia et al. 2019 (cultured U2-OS cells, 12,903 genes)

#### Outputs
- `adata_merfish` — AnnData with **645 cells × 12903 genes**
- `obsm['spatial']` — manually set XY coordinates
- **6 Leiden clusters** (corresponding to cell-cycle stages)
- UMAP and spatial embedding side-by-side plots

#### Key Result
> MERFISH clusters correspond to cell-cycle states. No spatial structure is expected (cultured cells), confirming method correctness.

---
---

## Tutorial 2 — Squidpy: Visium Fluorescence Analysis

**Source:** [squidpy.readthedocs.io/tutorial_visium_fluo](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_fluo.html)  
**Platform:** 10x Genomics Visium (Fluorescence)  
**Dataset:** Mouse brain coronal section — pre-processed, pre-annotated

### Overview

Demonstrates Squidpy's image analysis capabilities on a **fluorescence Visium** dataset (3-channel: DAPI, anti-NEUN, anti-GFAP). Shows how to segment nuclei, extract multi-scale image features, and compute image-based clusters to complement gene expression clusters.

---

### Section 2.1 — Import Packages & Data

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.datasets.visium_fluo_image_crop()` | squidpy | Load pre-processed fluorescence image crop |
| `sq.datasets.visium_fluo_adata_crop()` | squidpy | Load pre-annotated AnnData object |
| `sq.pl.spatial_scatter()` | squidpy | Visualize spot clusters on tissue |
| `img.show(channelwise=True)` | squidpy | Display each fluorescence channel separately |

#### Inputs
- Pre-processed dataset downloaded automatically (~303MB image + ~65MB AnnData)
- Mouse brain coronal section — smaller crop for efficiency

#### Outputs
- `img` — `squidpy.im.ImageContainer` with 3-channel fluorescence image:
  - Channel 0: **DAPI** (DNA marker)
  - Channel 1: **anti-NEUN** (neuronal marker)
  - Channel 2: **anti-GFAP** (glial cell marker)
- `adata` — pre-annotated AnnData with cluster labels in `adata.obs['cluster']`
- Spatial scatter plot of gene-space clusters

---

### Section 2.2 — Image Segmentation

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.im.process()` | squidpy | Apply image processing (e.g., Gaussian smoothing) |
| `sq.im.segment()` | squidpy | Segment nuclei using watershed algorithm |
| `img.crop_corner()` | squidpy | Crop a region of the image for visual inspection |
| `img_crop.show()` | squidpy | Visualize original vs. segmented image |

#### Inputs
- `img` — raw fluorescence ImageContainer
- `layer="image"`, `method="smooth"` — Gaussian smoothing on original image
- `layer="image_smooth"`, `method="watershed"`, `channel=0` (DAPI), `chunks=1000`

#### Outputs
- `img['image_smooth']` — Gaussian-smoothed version of the image
- `img['segmented_watershed']` — Label image where each nucleus is assigned a unique integer ID
- Side-by-side comparison plot: raw DAPI channel vs. segmented label image

#### Key Result
> Each nucleus in the DAPI channel is isolated as a distinct integer-labeled region. This enables downstream counting and morphological analysis per spot.

---

### Section 2.3 — Segmentation Features

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.im.calculate_image_features()` | squidpy | Compute per-spot image features from segmentation |
| `sq.pl.extract()` | squidpy | Temporarily move obsm features to obs for plotting |
| `sq.pl.spatial_scatter()` | squidpy | Spatial scatter with feature overlay |

#### Inputs
- `features="segmentation"`
- `label_layer="segmented_watershed"` — use the watershed label image
- `key_added="features_segmentation"` — where to store results in `adata.obsm`
- `n_jobs=1`

#### Outputs
- `adata.obsm['features_segmentation']` — per-spot matrix containing:
  - `segmentation_label` — number of segmented nuclei per spot (cell count estimate)
  - `segmentation_ch-0_mean_intensity_mean` — mean DAPI intensity in segmented regions
  - `segmentation_ch-1_mean_intensity_mean` — mean anti-NEUN intensity (neuron density)
  - `segmentation_ch-2_mean_intensity_mean` — mean anti-GFAP intensity (glial density)
- 4-panel spatial scatter plot comparing cell count, gene cluster, and per-channel intensities

#### Key Result
> The pyramidal layer of the Hippocampus has noticeably more cells per spot. Cortex_1 and Cortex_3 clusters show high anti-NEUN intensity (more neurons). Fiber_tracts and lateral ventricles show elevated anti-GFAP (more glial cells).

---

### Section 2.4 — Extract and Cluster Image Features

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.im.calculate_image_features()` | squidpy | Compute summary/histogram/texture features |
| `sc.pp.scale()` | scanpy | Scale features before PCA |
| `sc.pp.pca()` | scanpy | PCA on extracted features |
| `sc.tl.leiden()` | scanpy | Leiden clustering on image feature space |
| `sq.pl.spatial_scatter()` | squidpy | Visualize image-based clusters spatially |

#### Feature Extraction Configurations
| Config Name | Features | Scale | Context |
|-------------|----------|-------|---------|
| `features_orig` | summary, texture, histogram | 1.0 | Spot only (masked circle) |
| `features_context` | summary, histogram | 1.0 | With surrounding context |
| `features_lowres` | summary, histogram | 0.25 | Low resolution, more context |

#### Inputs
- `adata` and `img`
- Three different scale/context parameter sets per above table
- Combined into `adata.obsm['features']` via `pd.concat()`

#### Outputs
- `adata.obsm['features']` — concatenated multi-scale feature matrix
- `adata.obs['features_summary_cluster']` — Leiden clusters from summary features
- `adata.obs['features_histogram_cluster']` — Leiden clusters from histogram features
- `adata.obs['features_texture_cluster']` — Leiden clusters from texture features
- 4-panel spatial plot comparing all three image-based clusterings vs. gene-space clusters

#### Key Result
> Image-based clusters are spatially coherent and reveal finer-grained structure than gene-space clusters — especially in the Hippocampus and cortex layers. Texture, summary, and histogram features each capture different aspects of tissue morphology.

---
---

## Tutorial 3 — Squidpy: Visium H&E Analysis

**Source:** [squidpy.readthedocs.io/tutorial_visium_hne](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html)  
**Platform:** 10x Genomics Visium (H&E stain)  
**Dataset:** Mouse brain coronal section — pre-processed, pre-annotated

### Overview

Full pipeline for Visium H&E data using Squidpy. Combines image feature extraction with spatial graph analysis tools including **neighborhood enrichment**, **co-occurrence scoring**, **ligand-receptor interaction analysis**, and **Moran's I** spatially variable gene detection.

---

### Section 3.1 — Import Packages & Data

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.datasets.visium_hne_image()` | squidpy | Load pre-processed H&E tissue image |
| `sq.datasets.visium_hne_adata()` | squidpy | Load pre-annotated AnnData |
| `sq.pl.spatial_scatter()` | squidpy | Visualize cluster annotation spatially |

#### Inputs
- Pre-processed mouse brain coronal section Visium dataset (~314MB)
- Cluster annotations derived using Allen Brain Atlas + Linnarsson lab resources

#### Outputs
- `img` — ImageContainer with H&E image
- `adata` — AnnData with **2688 spots** and pre-annotated `adata.obs['cluster']`
- Spatial scatter of cluster annotation (Cortex, Hippocampus, Fiber_tracts, etc.)

---

### Section 3.2 — Image Features

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.im.calculate_image_features()` | squidpy | Extract summary features at multiple scales |
| `sc.pp.scale()` | scanpy | Feature standardization before PCA |
| `sc.pp.pca()`, `sc.pp.neighbors()`, `sc.tl.leiden()` | scanpy | Feature-space clustering |
| `sq.pl.spatial_scatter()` | squidpy | Compare image vs. gene cluster annotations |

#### Inputs
- Summary features computed at two scales: `scale=1.0` (local) and `scale=2.0` (more context)
- `key_added`: `"features_summary_scale1.0"` and `"features_summary_scale2.0"`
- `n_jobs=4`

#### Outputs
- `adata.obsm['features']` — combined multi-scale summary feature matrix
- `adata.obs['features_cluster']` — image-space Leiden clusters
- Side-by-side spatial scatter of `features_cluster` vs. `cluster` (gene-space)

#### Key Result
> Image clusters partially recapitulate gene clusters (e.g., *Fiber_tract*, Hippocampus region). However, cortex gene clusters show layered structure while image clusters reveal different cortical regions — complementary information.

---

### Section 3.3 — Spatial Statistics and Graph Analysis

#### Subsection: Neighborhood Enrichment

##### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.gr.spatial_neighbors()` | squidpy | Build spatial connectivity matrix (adjacency graph) |
| `sq.gr.nhood_enrichment()` | squidpy | Compute permutation-based neighborhood enrichment score |
| `sq.pl.nhood_enrichment()` | squidpy | Heatmap of enrichment z-scores between cluster pairs |

##### Inputs
- `adata` with spatial coordinates
- `cluster_key="cluster"` — annotation to use for enrichment
- Default `n_perms=1000` permutations

##### Outputs
- `adata.uns['spatial_neighbors']` — spatial connectivity and distance matrices
- `adata.obsp['spatial_connectivities']` — adjacency matrix
- `adata.uns['nhood_enrichment']` — enrichment z-score matrix (clusters × clusters)
- Heatmap visualization of enrichment scores

##### Key Result
> *Pyramidal_layer_dentate_gyrus* and *Pyramidal_layer* clusters are strongly enriched as neighbors of the *Hippocampus* cluster — consistent with their anatomical co-localization.

---

#### Subsection: Co-occurrence Across Spatial Dimensions

##### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.gr.co_occurrence()` | squidpy | Compute co-occurrence probability across radii |
| `sq.pl.co_occurrence()` | squidpy | Line plot of co-occurrence score vs. distance |

##### Score Formula
$$\text{score} = \frac{p(\text{exp} | \text{cond})}{p(\text{exp})}$$

Where:
- $p(\text{exp}|\text{cond})$ = conditional probability of seeing cluster `exp` given cluster `cond` is nearby
- $p(\text{exp})$ = baseline probability of observing cluster `exp`

##### Inputs
- `cluster_key="cluster"` — cluster annotation column
- `clusters="Hippocampus"` — condition cluster to analyze

##### Outputs
- `adata.uns['co_occurrence']` — co-occurrence score matrix across distance bins
- Line plot showing co-occurrence of all clusters relative to *Hippocampus* across increasing radii

##### Key Result
> *Pyramidal_layer* cluster co-occurs at short distances with *Hippocampus*, confirming tight spatial proximity. Distance units are in pixels of the source Visium image.

---

#### Subsection: Ligand-Receptor Interaction Analysis

##### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.gr.ligrec()` | squidpy | CellPhoneDB-inspired LR interaction analysis |
| `sq.pl.ligrec()` | squidpy | Dot plot of significant LR pairs |

##### Inputs
- `n_perms=100` — permutations for significance testing
- `cluster_key="cluster"`
- Filter parameters for visualization:
  - `source_groups="Hippocampus"`
  - `target_groups=["Pyramidal_layer", "Pyramidal_layer_dentate_gyrus"]`
  - `means_range=(3, np.inf)` — only highly expressed pairs
  - `alpha=1e-4` — stringent significance threshold

##### Outputs
- `adata.uns['ligrec']` — dictionary with mean expression and p-values for all LR pairs
- Dot plot where dot size = significance, dot color = mean expression
- Filtered to show LR pairs from Hippocampus → Pyramidal layers

##### Key Result
> Several candidate ligand-receptor pairs with potential roles in Hippocampal cellular communication are identified. Results could be further refined by integrating cell-type deconvolution.

---

#### Subsection: Spatially Variable Genes with Moran's I

##### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.gr.spatial_autocorr()` | squidpy | Compute Moran's I spatial autocorrelation |
| `sq.pl.spatial_scatter()` | squidpy | Visualize top spatially variable genes on tissue |

##### Inputs
- `mode="moran"` — Moran's I statistic (vs. `"geary"` for Geary's C)
- `genes` — subset of 1000 highly variable genes
- `n_perms=100`, `n_jobs=1`

##### Outputs
- `adata.uns['moranI']` — DataFrame with columns:
  - `I` — Moran's I score (higher = more spatially clustered)
  - `pval_norm`, `pval_sim` — p-values (normal approx. and simulation-based)
  - `pval_norm_fdr_bh` — Benjamini-Hochberg adjusted p-values
- Top genes ranked by I score: *Olfm1* (0.76), *Plp1* (0.75), *Itpka* (0.73), *Snap25* (0.72)
- Spatial scatter of top spatially variable genes

##### Key Result
> Top spatially variable genes are associated with pyramidal layers and fiber tracts, confirming their role as region-specific markers within the mouse brain.

---
---

## Tutorial 4 — Squidpy: Xenium Analysis

**Source:** [squidpy.readthedocs.io/tutorial_xenium](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_xenium.html)  
**Platform:** 10x Genomics Xenium (in-situ transcriptomics)  
**Dataset:** Human Lung Cancer FFPE tissue — 161,000 cells × 480 genes

### Overview

End-to-end pipeline for **Xenium** data — a high-resolution, single-cell in-situ platform. Uses the `spatialdata` framework for multi-modal data management (images, shapes, labels, transcripts). Covers QC, clustering, spatial graph analysis, co-occurrence, neighborhood enrichment, Moran's I, and interactive visualization with napari.

---

### Section 4.1 — Loading Xenium Data

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `xenium()` (from `spatialdata_io`) | spatialdata-io | Parse raw Xenium output directory |
| `sdata.write()` | spatialdata | Convert and save to Zarr format |
| `sd.read_zarr()` | spatialdata | Read from Zarr store for efficient access |

#### Inputs
- Raw Xenium output directory (`./Xenium/`) containing:
  - Transcript detection files
  - Cell segmentation files
  - Morphology focus image
- `zarr_path = "./Xenium.zarr"` — output Zarr store path

#### SpatialData Object Structure
```
SpatialData object:
├── Images:    'morphology_focus'  → MultiscaleSpatialImage [cyx]  (5 resolution levels)
├── Labels:    'cell_labels'       → MultiscaleSpatialImage [yx]   (cell segmentation)
│              'nucleus_labels'    → MultiscaleSpatialImage [yx]   (nucleus segmentation)
├── Points:    'transcripts'       → DataFrame (40,257,199 × 11)   (3D transcript coords)
├── Shapes:    'cell_boundaries'   → GeoDataFrame (161,000 cells)
│              'cell_circles'      → GeoDataFrame (161,000 cells)
│              'nucleus_boundaries'→ GeoDataFrame (161,000 cells)
└── Tables:    'table'             → AnnData (161,000 × 480)
```

#### Outputs
- `sdata` — unified SpatialData object with all Xenium modalities
- `adata = sdata.tables["table"]` — AnnData with:
  - **161,000 cells × 480 genes**
  - `obs`: `cell_id`, `transcript_counts`, `control_probe_counts`, `total_counts`, `cell_area`, `nucleus_area`, `z_level`
  - `obsm['spatial']` — 2D XY cell centroid coordinates

---

### Section 4.2 — Quality Control Metrics

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sc.pp.calculate_qc_metrics()` | scanpy | Compute per-cell QC statistics |
| `sns.histplot()` | seaborn | Distribution plots of QC metrics |
| `sc.pp.filter_cells()` | scanpy | Filter low-quality cells |
| `sc.pp.filter_genes()` | scanpy | Filter lowly detected genes |

#### Inputs
- `percent_top=(10, 20, 50, 150)` — top-gene percentage thresholds
- Control probe and codeword percentages calculated from `adata.obs`
- Filter thresholds: `min_counts=10` (cells), `min_cells=5` (genes)

#### Outputs
- Control probe % : **~0.005%** (very low → good data quality)
- Control codeword % : **~0.0025%** (very low → good data quality)
- 4-panel QC histogram:
  - Total transcripts per cell distribution
  - Unique transcripts per cell distribution
  - Segmented cell area distribution
  - Nucleus-to-cell area ratio distribution

#### Key Result
> Xenium data shows extremely low control probe rates, indicating high-quality transcript detection. Cell area and nucleus ratio distributions help identify atypical cells for exclusion.

---

### Section 4.3 — Preprocessing, Dimensionality Reduction, and Clustering

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sc.pp.normalize_total()` | scanpy | Normalize total counts per cell |
| `sc.pp.log1p()` | scanpy | Log-transform normalized counts |
| `sc.pp.pca()` | scanpy | PCA for dimensionality reduction |
| `sc.pp.neighbors()` | scanpy | Build k-NN graph |
| `sc.tl.umap()` | scanpy | UMAP embedding |
| `sc.tl.leiden()` | scanpy | Leiden community detection |

#### Inputs
- `adata.layers["counts"] = adata.X.copy()` — raw counts saved before normalization
- Standard Scanpy normalization → log1p → PCA → neighbors → UMAP → Leiden

#### Outputs
- `adata.obsm['X_umap']` — 2D UMAP coordinates
- `adata.obs['leiden']` — Leiden cluster labels per cell
- UMAP plots colored by `total_counts`, `n_genes_by_counts`, and `leiden`
- Spatial scatter plot of Leiden clusters over tissue

---

### Section 4.4 — Computation of Spatial Statistics

#### Subsection: Building a Spatial Neighborhood Graph

##### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.gr.spatial_neighbors()` | squidpy | Compute cell connectivity graph via Delaunay triangulation |

##### Inputs
- `coord_type="generic"` — for non-Visium (non-grid) spatial data
- `delaunay=True` — use Delaunay triangulation for neighbor detection

##### Outputs
- `adata.obsp['spatial_connectivities']` — connectivity matrix
- `adata.obsp['spatial_distances']` — distance matrix

---

#### Subsection: Centrality Scores

##### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.gr.centrality_scores()` | squidpy | Compute graph centrality measures per cluster |
| `sq.pl.centrality_scores()` | squidpy | Bar plots of centrality scores per cluster |

##### Centrality Metrics Computed
| Metric | Description |
|--------|-------------|
| **Closeness centrality** | How close a cluster's nodes are to all other nodes |
| **Degree centrality** | Fraction of non-cluster nodes connected to cluster members |
| **Clustering coefficient** | How densely cluster nodes interconnect with each other |

##### Inputs
- `cluster_key="leiden"` — Leiden cluster labels

##### Outputs
- `adata.uns['leiden_centrality_scores']` — DataFrame of per-cluster centrality scores
- Bar plot visualization of all three centrality metrics across Leiden clusters

---

#### Subsection: Co-occurrence Probability

##### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sc.pp.subsample()` | scanpy | Randomly subsample 50% of cells for efficiency |
| `sq.gr.co_occurrence()` | squidpy | Compute pairwise co-occurrence across radii |
| `sq.pl.co_occurrence()` | squidpy | Plot co-occurrence score curves |
| `sq.pl.spatial_scatter()` | squidpy | Spatial scatter of subsampled dataset |

##### Inputs
- `adata_subsample` — 50% random subsample (~80,000 cells)
- `cluster_key="leiden"`
- `clusters="12"` — condition cluster to analyze

##### Outputs
- `adata_subsample.uns['co_occurrence']` — co-occurrence score across radii
- Line plot of co-occurrence probability ratio for cluster "12" vs. all others
- Spatial scatter plot of Leiden clusters on subsampled data

---

#### Subsection: Neighborhood Enrichment Analysis

##### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.gr.nhood_enrichment()` | squidpy | Permutation-based enrichment score |
| `sq.pl.nhood_enrichment()` | squidpy | Heatmap of z-scores |

##### Inputs
- `cluster_key="leiden"`
- Default `n_perms=1000`

##### Outputs
- `adata.uns['nhood_enrichment']` — enrichment z-score matrix
- Heatmap of neighborhood enrichment
- Side-by-side comparison with spatial scatter of subsampled data

---

#### Subsection: Moran's I Spatial Autocorrelation

##### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `sq.gr.spatial_autocorr()` | squidpy | Global spatial autocorrelation per gene |
| `sq.pl.spatial_scatter()` | squidpy | Spatial expression plot for top SVGs |

##### Inputs
- `mode="moran"` — Moran's I statistic
- `n_perms=100`, `n_jobs=1`
- All 480 Xenium genes evaluated

##### Outputs
- `adata_subsample.uns['moranI']` — ranked DataFrame:
  - Top SVGs: *AREG* (I=0.696), *MET* (I=0.683), *ANXA1* (I=0.667), *EPCAM* (I=0.633)
- Spatial scatter of top spatially variable genes
- `spatialdata-plot` visualization overlaying gene expression on morphology image

---

### Section 4.5 — Interactive Visualization with napari-spatialdata

#### Tools Used
| Tool | Package | Purpose |
|------|---------|---------|
| `napari_spatialdata.Interactive` | napari-spatialdata | Interactive multi-layer viewer |

#### Inputs
- Complete `sdata` SpatialData object (images, labels, shapes, points, tables)

#### Outputs
- Interactive napari GUI with:
  - Morphology focus image layer
  - Cell/nucleus boundary overlays
  - Gene expression coloring (e.g., AREG expression per cell)
  - Leiden cluster annotation overlay

---
---

## Environment Setup

### Required Packages

```bash
# Core packages
pip install scanpy squidpy anndata

# For Xenium tutorial
pip install spatialdata spatialdata-io spatialdata-plot

# For interactive visualization
pip install napari napari-spatialdata

# Additional utilities
pip install pandas seaborn matplotlib openpyxl
```

### Conda Environment

```bash
conda env create -f environment.yml
conda activate squidpy-env
```

---

## Quick Reference: Key Functions

| Function | Package | Used In | Purpose |
|----------|---------|---------|---------|
| `sc.datasets.visium_sge()` | scanpy | T1 | Load 10x Visium dataset |
| `sc.pp.calculate_qc_metrics()` | scanpy | T1, T4 | QC statistics |
| `sc.pp.filter_cells/genes()` | scanpy | T1, T4 | Basic filtering |
| `sc.pp.normalize_total()` | scanpy | T1, T4 | Count normalization |
| `sc.pp.log1p()` | scanpy | T1, T4 | Log transformation |
| `sc.pp.highly_variable_genes()` | scanpy | T1 | Select HVGs |
| `sc.pp.pca()` | scanpy | T1–T4 | Dimensionality reduction |
| `sc.tl.leiden()` | scanpy | T1–T4 | Graph clustering |
| `sc.tl.umap()` | scanpy | T1, T4 | 2D embedding |
| `sc.pl.spatial()` | scanpy | T1 | Spatial visualization |
| `sc.tl.rank_genes_groups()` | scanpy | T1 | Differential expression |
| `sq.im.process()` | squidpy | T2 | Image smoothing |
| `sq.im.segment()` | squidpy | T2 | Nucleus segmentation |
| `sq.im.calculate_image_features()` | squidpy | T2, T3 | Image feature extraction |
| `sq.gr.spatial_neighbors()` | squidpy | T3, T4 | Spatial graph construction |
| `sq.gr.nhood_enrichment()` | squidpy | T3, T4 | Neighborhood enrichment |
| `sq.gr.co_occurrence()` | squidpy | T3, T4 | Co-occurrence probability |
| `sq.gr.ligrec()` | squidpy | T3 | Ligand-receptor analysis |
| `sq.gr.spatial_autocorr()` | squidpy | T3, T4 | Moran's I / Geary's C |
| `sq.gr.centrality_scores()` | squidpy | T4 | Graph centrality metrics |
| `sq.pl.spatial_scatter()` | squidpy | T2–T4 | Spatial scatter plots |
| `xenium()` | spatialdata-io | T4 | Load Xenium data |
| `sd.read_zarr()` | spatialdata | T4 | Read Zarr store |

---

## Data Summary

| Tutorial | Platform | Tissue | Cells/Spots | Genes | Key Analysis |
|----------|----------|--------|-------------|-------|-------------|
| T1 (Scanpy) | Visium + MERFISH | Human Lymph Node + U2-OS cells | 3861 spots + 645 cells | 36601 + 12903 | Clustering, marker genes, spatial visualization |
| T2 (Squidpy Fluo) | Visium Fluorescence | Mouse Brain | 704 spots | Pre-processed | Image segmentation, feature clustering |
| T3 (Squidpy H&E) | Visium H&E | Mouse Brain | 2688 spots | Pre-processed | Nhood enrichment, co-occurrence, LR, Moran's I |
| T4 (Squidpy Xenium) | Xenium | Human Lung Cancer | 161,000 cells | 480 | Full spatial stats, centrality, napari viewer |

---

*Tutorials sourced from [Scanpy Tutorials](https://scanpy-tutorials.readthedocs.io) and [Squidpy Documentation](https://squidpy.readthedocs.io). Part of the scverse ecosystem.*
