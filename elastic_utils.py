from ase.constraints import FixSymmetry
from ase.optimize import FIRE

from ase import Atoms
from ase.filters import Filter

def relax_atoms(
    atoms: Atoms,
    steps: int = 500,
    fmax: float = 0.02,
    optimizer_cls: None = None,
    fix_symmetry: bool = False,
    cell_filter_cls: type[Filter] | None = None,
) -> Atoms:
    """Simple helper function to run relaxations and return the relaxed Atoms

    Args:
        atoms: ASE atoms with a calculator
        steps: max number of relaxation steps
        fmax: force convergence threshold
        optimizer_cls: ASE optimizer. Default FIRE
        fix_symmetry: fix structure symmetry in relaxation: Default False
        cell_filter_cls: An instance of an ASE filter.

    Returns:
        Atoms: relaxed atoms
    """

    if fix_symmetry:
        atoms.set_constraint(FixSymmetry(atoms))

    if cell_filter_cls is not None:
        _atoms = cell_filter_cls(atoms)
    else:
        _atoms = atoms

    optimizer_cls = FIRE if optimizer_cls is None else optimizer_cls
    opt = optimizer_cls(_atoms, logfile=None)
    opt.run(fmax=fmax, steps=steps)

    atoms.info |= {"opt_nsteps": opt.nsteps}
    return atoms


from pymatgen.analysis.elasticity import DeformedStructureSet, ElasticTensor, Strain
from pymatgen.io.ase import AseAtomsAdaptor
from collections.abc import Sequence
from ase import Atoms
from ase.calculators.calculator import Calculator


def calculate_elasticity(
    atoms: Atoms,
    calculator: Calculator,
    norm_strains: Sequence[float] | float = (-0.01, -0.005, 0.005, 0.01),
    shear_strains: Sequence[float] | float = (-0.06, -0.03, 0.03, 0.06),
    relax_initial: bool = True,
    relax_strained: bool = True,
    use_equilibrium_stress: bool = True,
    **relax_kwargs,
):
    """Calculate elastic tensor, bulk, shear moduli following MP workflow

    Will not run a relaxation. We do that outside in order to be able to have more control.

    Args:
        atoms: ASE atoms object
        calculator: an ASE Calculator
        norm_strains: sequence of normal strains
        shear_strains: sequence of shear strains
        relax_initial: relax the initial structure. Default is True
        relax_strained: relax the atomic positions of strained structure. Default True.
        use_equilibrium_stress: use equilibrium stress in calculation. For relaxed structures this
            should be very small

    Returns:
        dict of elasticity results
    """
    atoms.calc = calculator
    if relax_initial:
        atoms = relax_atoms(atoms, **relax_kwargs)

    eq_stress = atoms.get_stress(voigt=False) if use_equilibrium_stress else None

    deformed_structure_set = DeformedStructureSet(
        AseAtomsAdaptor.get_structure(atoms),
        norm_strains,
        shear_strains,
    )

    strains = [
        Strain.from_deformation(deformation)
        for deformation in deformed_structure_set.deformations
    ]

    stresses = []
    for deformed_structure in deformed_structure_set:
        atoms = deformed_structure.to_ase_atoms()
        atoms.calc = calculator

        if relax_strained is True:
            relax_kwargs.update({"cell_filter_cls": None})
            atoms = relax_atoms(atoms, **relax_kwargs)

        stresses.append(atoms.get_stress(voigt=False))

    elastic_tensor = ElasticTensor.from_independent_strains(
        strains=strains, stresses=stresses, eq_stress=eq_stress
    )

    # results are in eV/A^3
    results = {
        "elastic_tensor": elastic_tensor,
        "shear_modulus_vrh": elastic_tensor.g_vrh,
        "bulk_modulus_vrh": elastic_tensor.k_vrh,
    }

    return results