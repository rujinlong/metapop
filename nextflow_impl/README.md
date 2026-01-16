# MetaPop Nextflow Implementation

## Overview
This is a Nextflow-based re-implementation of the [MetaPop](https://github.com/rujinlong/metapop) pipeline. It is designed for the macro- and micro-diversity analysis of metagenomic-derived populations (MAGs/genomes).

**Key Improvements over Original:**
*   **Scalability**: Uses Nextflow to enable distributed execution on HPC clusters (Slurm, etc.) and cloud environments.
*   **Parallelism**: Replaces single-node `multiprocessing` with efficient process-level parallelization.
*   **Reproducibility**: Containerized execution (Docker/Singularity) and strict version control.
*   **Performance**: Optimized resource usage with configurable profiles (memory/CPU).

---

## 🏗 Architecture

The pipeline uses **Nextflow DSL2** for modularity and reusability. It is organized into:

```text
nextflow_impl/
├── main.nf              # Orchestrator workflow
├── nextflow.config      # Configuration & Profiles (HPC/Docker)
├── modules/             # Process definitions
│   ├── preprocessing.nf # QC, Coverage, Genome Filtering
│   ├── snp_calling.nf   # Variant detection & annotation
│   └── microdiversity.nf# Population genetics stats (pN/pS, Fst)
└── bin/                 # Helper scripts (Python)
    ├── filter_reads.py  # Custom BAM filtering
    ├── calc_coverage.py # Breadth/TAD calculation
    └── ...
```

---

## 🧬 Detailed Workflow Description

### 1. Reference Preparation (`STAGE 0`)

Before processing samples, the reference genome is prepared:

*   **indexing**: `samtools faidx` creates an index for random access.
*   **gene_prediction**: 
    *   Uses **Pyrodigal** (a highly parallel implementation of Prodigal).
    *   **Mode Logic**: Defaults to `single` mode (optimized for finished genomes). 
    *   **Fallback Strategy**: If sequences are too short (<20kb) for `single` mode training, the pipeline automatically retries with `meta` mode (metagenomic anonymous mode).
    *   Generates `.gff` (annotations), `.fna` (genes), and `.faa` (proteins).

### 2. Preprocessing (`STAGE 1`)

Ensures only high-quality reads mapping to targeted genomes are used.

*   **FILL_MD** (`samtools calmd`): 
    *   Recomputes MD tags to accurately track mismatches.
    *   *Why?* Many aligners do not generate perfect MD tags, which are crucial for calculating percent identity.

*   **FILTER_READS** (`bin/filter_reads.py`):
    *   Filters reads based on **Percent Identity** (defaults to ≥95%) and **Alignment Length**.
    *   *Bio Logic*: Removes reads that may originate from related but distinct organisms (e.g., different species/strains), reducing noise in SNP calling.
    *   Supports both Global Identity (matches / read_len) and Local Identity (matches / align_len).

*   **CALCULATE_COVERAGE** (`bin/calc_coverage.py`):
    *   Computes two key metrics per contig:
        1.  **Breadth**: % of bases covered by at least 1 read.
        2.  **TAD (Truncated Average Depth)**: Average depth after removing the top/bottom X% (default 10%) of positions.
    *   *Bio Logic*: TAD is more robust than mean depth as it removes outlier peaks (e.g., repeat regions, transposons) and deep troughs.

*   **FILTER_GENOMES** (`bin/filter_genomes.py`):
    *   Filters the BAM files to retain *only* reads mapping to contigs that pass coverage thresholds (default: ≥20% breadth, ≥10x depth).
    *   *Bio Logic*: Ensures that we only call SNPs on genomes that are sufficiently abundant in the sample to give reliable population statistics.

### 3. SNP Calling (`STAGE 2`)

Identifies single nucleotide polymorphisms within the population.

*   **CALL_VARIANTS** (`bcftools`):
    *   Uses `bcftools mpileup` + `bcftools call` to identify potential variant sites.
    *   Filters: `QUAL > 20`.
    *   *Note*: Variants are called per-sample but can be aggregated.

*   **REFINE_SNPS** (`refine_snps.py`):
    *   Counts A, T, C, G alleles at every identified variant position.
    *   *Bio Logic*: Essential for calculating allele frequencies, which drive microdiversity metrics (pi, theta).

*   **CALL_CONSENSUS** (`call_consensus.py`):
    *   Determines the major allele at every position across the population.

*   **ASSOCIATE_GENES** (`associate_genes.py`):
    *   Maps SNPs to predicted genes.
    *   Determines **Codon Position** (1st, 2nd, 3rd base).
    *   *Bio Logic*: Crucial for distinguishing **Synonymous (pS)** vs. **Non-Synonymous (pN)** mutations, a key indicator of selection pressure.

### 4. Microdiversity (`STAGE 3`)

(Currently partially implemented placeholders structure)
*   **Goal**: Calculate population genetic metrics:
    *   **$\pi$ (Nucleotide Diversity)**: Average difference between any two sequences.
    *   **$\theta$ (Watterson's Estimator)**: Based on number of segregating sites.
    *   **Tajima's D**: $\pi - \theta$ (Detects selection/demography).
    *   **Fixation Index ($F_{ST}$)**: Population differentiation.
    *   **Linked SNPs**: Analysis of SNPs on the same read to phase haplotypes.

---

## 💻 HPC & Resource Management

This workflow is optimized for HPC execution.

### Profiles (`nextflow.config`)
*   `standard`: Local execution.
*   `slurm`: Uses Slurm scheduler + Singularity.
*   `docker`: Uses local Docker engine.

### Resource Labels
Each process is tagged with a label (`process_low`, `process_medium`, `process_high`) which maps to specific resources in `nextflow.config`.
*   **Example**: `FILTER_READS` is `process_medium` (4 CPUs, 8GB RAM).
*   **Dynamic Resources**: Can be adjusted in custom configs without changing code.

### Containers
*   Uses **Docker/Singularity** to encapsulate all dependencies (`samtools`, `bcftools`, `python`, `pyrodigal`).
*   Configured to auto-mount paths for Singularity on HPCs.

---

## 🛠 For Developers: How to Extend

1.  **Adding a new module**:
    *   Create `modules/new_module.nf`.
    *   Define `process NEW_PROCESS`.
    *   Import it in `main.nf`.

2.  **Updating dependencies**:
    *   Modify `env.yml` (Conda env file).
    *   Rebuild Docker image: `docker build -t jinlongru/metapop-latest .`

3.  **Testing**:
    *   Use the `test` profile: `nextflow run main.nf -profile test`.
    *   It uses bundled test data in `test_data/`.
