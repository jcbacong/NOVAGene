# NOVAGene

Nonlinear Overlapping Variability Analysis of Gene Networks or NOVAGene is a gene coexpression analysis framework that extends WGCNA by incorporating nonlinear association measures. This repository contains all scripts and notebooks for reproducing the analyses presented in our study.

Bacong et al (2026) NOVAGene: a nonlinear approach to improve biological resolution in weighted gene co-expression networks for low-variability datasets. Front. Bioinform. 6:1813626. doi: 10.3389/fbinf.2026.1813626

---
## Requirements

### 💻 Computational Environment

The majority of heavy computations were performed in a **SLURM-based HPC (High-Performance Computing) Linux environment**. I personally recommend using an HPC environment for full-scale analyses; however, a desktop machine with a capable GPU may be used as an alternative.

R and Python notebooks include step-by-step computations for subsequent analyses and can be run on personal computers.

### Python
Install the required Python packages:

```sh
# Create and activate conda environment
conda create -n novagene python=3.10
conda activate novagene

# Install dependencies
pip install -r requirements.txt
```

### R
Install the required R packages by running the following in your R console:

```r
install.packages(c("WGCNA", "DGCA", "DiffCoEx", "energy"))
```

---

## 📁 Repository Files & Structure

### `dataset/`
Contains CSV files where:
- **Rows** represent genes of interest from microarray data (but can be applied to RNA-seq data as well)
- **Columns** represent patient identifiers

---

### `wgcna_R/`
Contains R Markdown (`.Rmd`) notebooks for:
- Standard **WGCNA** analysis
- **dcor_WGCNA** analysis
- **NOVAGene** (our proposed method)

Please refer to the individual notebooks for step-by-step procedures.

---

### `calc_dcor.py`
Computes the **N×N distance correlation matrix** from gene expression data.

```python
python calc_dcor.py --input <input.txt> --output <output.csv>
```

```bash
sbatch --array=<no of compute nodes> run_dcor.sh
```

---

### `correxp-v2.py`
Calculates the **correlation matrix** between two phenotypes (or groups).

```python
python correxp-v2.py \
    --moduleA <moduleA.csv> \
    --moduleB <moduleB.csv> \
    --exprA <exprA.csv> \
    --exprB <exprB.csv> \
    --adjA <adjA.csv> \
    --adjB <adjB.csv> \
    --output <output.csv>
```

```bash
sbatch --array=<no of compute nodes> run_correxpv2.sh
```


---

### `benchmarking/`
Contains benchmarking scripts for comparing against:
- **DGCA**
- **DiffCoEx**
---

## 🚀 Getting Started

Clone the repository:

```sh
git clone https://github.com/jcbacong/NOVAGene.git
cd NOVAGene
```

---

## 📄 Citation

If you use NOVAGene in your research, please cite our paper accordingly.

Bacong JRC, Nevado JJB, Alejandria MM, Salamat MSS and Buensalido JAL (2026) NOVAGene: a nonlinear approach to improve biological resolution in weighted gene co-expression networks for low-variability datasets. Front. Bioinform. 6:1813626. doi: 10.3389/fbinf.2026.1813626

---

## 📬 Contact

For questions or issues, please open a GitHub Issue or contact me directly at bacong.junelle@gmail.com
