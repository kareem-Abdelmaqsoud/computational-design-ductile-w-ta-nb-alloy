import json
import numpy as np
import subprocess

# Base template for W–Ta SQS
template = {
    "structure": {
        "lattice": [
            [3.20, 0.0, 0.0],
            [0.0, 3.20, 0.0],
            [0.0, 0.0, 3.20]
        ],
        "coords": [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5]
        ],
        "species": ["W", "Ta"],
        "supercell": [3, 3, 3]
    },
    "iterations": 5000000,
    "shell_weights": {
        "1": 1.0,
        "2": 0.5
    },
    "composition": {
        "W": 27,
        "Ta": 27
    },
    "max_results_per_objective": 10
}

total_atoms = 54

# Loop over W fraction (skip pure W or Ta)
for w_frac in np.linspace(0.1, 0.9, 9):
    w_atoms = int(round(total_atoms * w_frac))
    ta_atoms = total_atoms - w_atoms

    # Update composition in template
    template["composition"]["W"] = w_atoms
    template["composition"]["Ta"] = ta_atoms

    # Save JSON for this composition
    json_file = f"sqs_w{w_atoms}_ta{ta_atoms}.json"
    with open(json_file, "w") as f:
        json.dump(template, f, indent=4)

    # Run SQS generation
    subprocess.run(["sqsgen", "run", "-i", json_file], check=True)