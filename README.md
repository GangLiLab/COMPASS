
# COMPASS

**COMPASS** (*Comprehensive Multi-omics Polycomb Analysis Scoring System*) is a context-aware framework for prioritising candidate Polycomb target genes by integrating PcG-perturbation RNA-seq and H3K27me3 ChIP-seq evidence.

COMPASS converts heterogeneous transcriptomic and epigenomic evidence into a continuous, interpretable **Polycomb Target Score (PTS)**, enabling gene prioritisation across tissues, cell types and disease contexts.

## Overview

Polycomb target identification is often based on either chromatin occupancy or perturbation-responsive expression. However, these assays capture complementary but incomplete aspects of Polycomb regulation.

COMPASS addresses this by:

- scoring RNA-seq perturbation evidence and H3K27me3 ChIP-seq evidence independently;
- integrating evidence within biological contexts;
- applying biologically interpretable pattern weights;
- retaining evidence provenance and incomplete-evidence information;
- generating a continuous, context-specific Polycomb Target Score.

## Key Features

- Context-aware Polycomb target prioritisation
- Integration of perturbation transcriptomics and H3K27me3 evidence
- Interpretable scoring framework
- Support for incomplete single-modality evidence
- Evidence-traceable gene ranking
- Robustness, ablation and benchmark analysis modules

## Main Output

The main output is a gene-context level **Polycomb Target Score (PTS)**.

A high PTS indicates strong available evidence for Polycomb-associated regulation in a defined context. It should be interpreted as prioritisation evidence, not proof of direct PcG binding or causal repression.

## Repository Contents

```text
COMPASS/
├── data/          # input metadata and processed evidence tables
├── scripts/       # COMPASS scoring and analysis scripts
├── results/       # output scores, benchmarks and figures
└── README.md
```

## Citation

Wait to update


## License

License information will be added before public release.
