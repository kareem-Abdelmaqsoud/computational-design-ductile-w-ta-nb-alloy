"""Compute UMA (uma-s-1p2p1) elastic properties for the W-Ta-Nb SQS structures.

Outputs (in this directory):
  - sqs_structures_full_triangle_with_elastic_uma_s_1p2p1.json
      63 full-triangle SQS structures + elastic_properties (pugh_ratio = K/G)
  - uma_s_1p2p1_dft_sqs_elastic.json
      UMA on the exact 10% Nb SQS cells used for DFT, plus pure W, Ta, Nb

uma-s-1p2p1 needs fairchem-core >= 2.22; run with the fairchem223 environment:
  /home/jovyan/shared-scratch-kabdelma-pvc/kabdelma/envs/fairchem223/bin/python generate_uma_elastic_data.py
"""
import json
import os
import re

from ase.filters import FrechetCellFilter
from ase.io import read
from ase.lattice.cubic import BodyCenteredCubic
from fairchem.core import FAIRChemCalculator, pretrained_mlip
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from tqdm import tqdm

from elastic_utils import calculate_elasticity

MODEL = "uma-s-1p2p1"
EV_TO_GPA = 160.2176621
DFT_SQS_DIR = "../vasp_internal_elastic/sqs/elastic"
# Max 1% shear strain, comparable to the small DFT finite-difference distortions
NORM_STRAINS = (-0.01, -0.005, 0.005, 0.01)
SHEAR_STRAINS = (-0.01, -0.005, 0.005, 0.01)

predictor = pretrained_mlip.get_predict_unit(MODEL, device="cuda")
calc = FAIRChemCalculator(predictor, task_name="omat")


def elastic_properties(atoms):
    results = calculate_elasticity(
        atoms,
        norm_strains=NORM_STRAINS,
        shear_strains=SHEAR_STRAINS,
        calculator=calc,
        cell_filter_cls=FrechetCellFilter,
        fix_symmetry=False,
    )
    bulk_gpa = results["bulk_modulus_vrh"] * EV_TO_GPA
    shear_gpa = results["shear_modulus_vrh"] * EV_TO_GPA
    voigt = results["elastic_tensor"].voigt * EV_TO_GPA
    return {
        "model": MODEL,
        "max_shear_strain": max(SHEAR_STRAINS),
        "bulk_modulus_vrh": bulk_gpa,
        "shear_modulus_vrh": shear_gpa,
        "pugh_ratio": bulk_gpa / shear_gpa,  # K/G
        "c11": voigt[0][0],
        "c12": voigt[0][1],
        "c44": voigt[3][3],
    }


def w_fraction_from_dft_name(name):
    """DFT folder names like 'ta5nbw4', 'tanbw8', 'ta8nbw', 'nbw9' (10 atoms)."""
    match = re.search(r"w(\d+)", name)
    return (int(match.group(1)) if match else 1) / 10.0


# --- 1. Full-triangle SQS structures ---
with open("sqs_structures_full_triangle.json") as f:
    sqs_dicts = json.load(f)

full_triangle = {}
for key, struct_dict in tqdm(sqs_dicts.items(), desc="full triangle"):
    atoms = AseAtomsAdaptor.get_atoms(Structure.from_dict(struct_dict))
    entry = dict(struct_dict)
    entry["elastic_properties"] = elastic_properties(atoms)
    full_triangle[key] = entry

with open("sqs_structures_full_triangle_with_elastic_uma_s_1p2p1.json", "w") as f:
    json.dump(full_triangle, f, indent=2)

# --- 2. DFT-matched 10% Nb SQS cells and pure elements ---
dft_matched = {"model": MODEL, "sqs": {}, "pure": {}}
for name in tqdm(sorted(os.listdir(DFT_SQS_DIR)), desc="DFT SQS cells"):
    poscar = os.path.join(DFT_SQS_DIR, name, "POSCAR")
    if not os.path.exists(poscar):
        continue
    props = elastic_properties(read(poscar))
    props["W_fraction"] = w_fraction_from_dft_name(name)
    dft_matched["sqs"][name] = props

for element in tqdm(["W", "Ta", "Nb"], desc="pure elements"):
    dft_matched["pure"][element] = elastic_properties(BodyCenteredCubic(element).repeat(2))

with open("uma_s_1p2p1_dft_sqs_elastic.json", "w") as f:
    json.dump(dft_matched, f, indent=2)

print("Done.")
for element, p in dft_matched["pure"].items():
    print(f"{element}: K={p['bulk_modulus_vrh']:.1f} G={p['shear_modulus_vrh']:.1f} "
          f"C44={p['c44']:.1f} K/G={p['pugh_ratio']:.3f}")
