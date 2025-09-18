import json
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

# load from JSON
with open("enumerated_structures.json", "r") as f:
    loaded = json.load(f)

# convert JSON dicts → pymatgen Structures
pmg_structures = [Structure.from_dict(d) for d in loaded]

# convert pymatgen Structures → ASE Atoms
structures = [AseAtomsAdaptor.get_atoms(pmg) for pmg in pmg_structures]


from fairchem.core import pretrained_mlip, FAIRChemCalculator

# load predictor
predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cpu")
calc = FAIRChemCalculator(predictor, task_name="omat")


from ase.optimize import LBFGS
from ase.constraints import ExpCellFilter  # alias for FrechetCellFilter in new ASE
from tqdm import tqdm

relaxed_structures = []
for atoms in tqdm(structures):
    atoms.calc = calc
    opt = LBFGS(ExpCellFilter(atoms))
    opt.run(fmax=0.05, steps=200)
    relaxed_structures.append(atoms.copy())



import json
from pymatgen.io.ase import AseAtomsAdaptor

# convert ASE Atoms → pymatgen Structures
relaxed_pmg = [AseAtomsAdaptor.get_structure(atoms) for atoms in relaxed_structures]

# convert to dicts
relaxed_dicts = [s.as_dict() for s in relaxed_pmg]

# save as JSON
with open("relaxed_structures.json", "w") as f:
    json.dump(relaxed_dicts, f, indent=2)