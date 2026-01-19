## Computational Design of Ductile Additively Manufactured WTaNb Alloys

This repository contains the computational data and code supporting the paper on the design of ductile W-Ta-Nb alloys. The core of the work involves accelerating the computational screening of compositions by correlating the predicted Pugh ratio ($B/G$) with experimental ductility (crack fraction).

### Key Concepts

- **Ductility Predictor:** The Pugh ratio ($K/G$), the ratio of bulk modulus to shear modulus) is used as a computationally efficient predictor for material ductility.
- **Computational Acceleration:** Machine Learning Interatomic Potentials (MLIPs), specifically MACE, are used to accelerate Density Functional Theory (DFT) calculations of elastic constants, enabling efficient screening along the melting point-Pugh ratio Pareto front.
- **W-Ta-Nb System:** Focuses on designing tungsten-based alloys with improved printability for additive manufacturing.
- **Electronic Basis:** The trends in the Pugh ratio are explained via the electronic density of states at the Fermi level.
- **Experimental Validation:** The computed Pugh ratio shows a strong correlation with experimentally observed crack fractions in additively manufactured alloys, successfully identifying the most ductile compositions (e.g., W20Ta70Nb10 and W30Ta60Nb10).

### Repository Contents

| Directory/File               | Purpose                                       | Key Content                                                                                                                                                                                                                                                |
| :--------------------------- | :-------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `miscellaneous_code/`        | Computational scripts and analysis notebooks. | Structure generation (`generate_sqs.ipynb`), structural relaxation (`relaxing_structures.py`), MLIP/MACE calculations (`mace_ordered_sqs.ipynb`), DFT VASP setup (`run_vasp_internal_elastic_calc.ipynb`), and final analysis (`sqs_elastic_trend.ipynb`). |
| `paper_figures/`             | Final figures for the publication.            | Plots showing the Pareto front, electronic analysis (DOS), and correlation between computed Pugh ratio and experimental crack fraction.                                                                                                                    |
| `sqs_structures*.json`       | Initial configuration data.                   | Special Quasirandom Structures (SQS) used for simulating disordered alloys.                                                                                                                                                                                |
| `enumerated_structures.json` | Additional alloy configurations.              | Enumerated alloy structures for potential MLIP training or testing.                                                                                                                                                                                        |
| `relaxed_structures*.json`   | Final structural data.                        | Atomically relaxed structures, including calculated total energies, ready for property extraction.                                                                                                                                                         |
