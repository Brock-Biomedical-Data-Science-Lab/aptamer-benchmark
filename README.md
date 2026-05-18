# Protein-Aptamer Complex Benchmark
This repository contains the datasets and selected scripts used in the paper:

“Comprehensive Evaluation of Artificial Intelligence-Empowered Approaches for Protein-Aptamer Complex Prediction” (Zhao et al.)

## Directory Structure
```text
data/
├── processed_no_ions/        # All processed protein–aptamer complex PDB files without ions
├── processed_with_ions/      # All processed protein–aptamer complex PDB files with ions
└── raw/                      # PDB files directly downloaded from the PDB (Biological Assembly 1)

scripts/
├── aptamer-only/             # Scripts for aptamer-only
└── protein-aptamer-complex/  # Scripts for protein–aptamer complex
```
## Citing This Work
@Article{Zhao2026,
author={Zhao, Jiani
and Tram, Kha
and Yan, Hongbin
and Li, Yifeng},
title={Comprehensive evaluation of artificial intelligence-empowered approaches for protein--aptamer complex prediction},
journal={Briefings in Bioinformatics},
year={2026},
month={May},
day={01},
volume={27},
number={3},
pages={bbag206},
abstract={Drug discovery is a time-consuming, expensive, and high-risk process. Recent advances in artificial intelligence (AI) have enabled major breakthroughs in small-molecule and protein therapeutics. However, AI-driven design of aptamer drugs remains largely unexplored. Aptamers are short (15--100 nt) single-stranded DNAs or RNAs that exhibit high binding affinity, high specificity, and low immunogenicity, making them promising candidates for disease (such as cancer) therapeutics. Compared with protein--ligand or protein--protein systems, protein--aptamer complexes are under-represented in public structural databases, and aptamers themselves are highly flexible and relatively large molecules. These characteristics present distinct challenges for AI-based structural modeling. Here, we systematically evaluate recent AI frameworks, including AlphaFold3, Chai-1, Boltz-2, and RoseTTAFold2NA, along with a template-based approach, in predicting protein--aptamer complex structures and estimating binding free energies. We establish an independent benchmark to assess their performance in structural accuracy, stability, and energetic consistency. This study provides a foundation for the application of AI in aptamer drug design and offers a reference framework for future research in nucleic-acid therapeutics and biomolecular modeling.},
issn={1477-4054},
doi={10.1093/bib/bbag206},
url={https://doi.org/10.1093/bib/bbag206}
}
