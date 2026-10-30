from pathlib import Path
import json
import numpy as np
json_dir = Path(__file__).parent / "scaled_designs"
output_fp = Path(__file__).parent / "pcc_designs-comsol_scale-map.json"
wavelengths = np.arange(-9, 10, dtype=int) * 5 + 2923
if __name__ == "__main__":
    temp = dict()
    for i in range(5):
        design_name = f"design-{i}"
        temp[design_name] = {}
        json_file = json_dir / f"{design_name}.json"
        with open(json_file, "r") as f:
            data = json.load(f)["designs"]
        for d in data:
            w = d['name'][-4:]
            s = d['search']['scale']
            temp[design_name][w] = s
    with open(output_fp, "w") as f:
        json.dump(temp, f, indent=4)


