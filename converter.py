import json
import os
import sys

INPUT_JSON = 'tree.json'
OUTPUT_TSV = 'targets.tsv'
MAX_Y_HG38 = 57227415

if __name__ == "__main__":
    if not os.path.exists(INPUT_JSON):
        sys.exit(1)

    with open(INPUT_JSON, 'r', encoding='utf-8') as f:
        try:
            data = json.load(f)
        except json.JSONDecodeError:
            sys.exit(1)

    unique_positions = set()
    all_nodes = data.get('allNodes', {})

    for node_id in all_nodes:
        node = all_nodes[node_id]
        variants = node.get('variants', [])

        for v in variants:
            pos_val = v.get('position')
            if pos_val is not None:
                try:
                    pos = int(pos_val)
                    if 0 < pos <= MAX_Y_HG38:
                        unique_positions.add(pos)
                except (ValueError, TypeError):
                    continue

    sorted_positions = sorted(list(unique_positions))

    if not sorted_positions:
        sys.exit(1)

    with open(OUTPUT_TSV, 'w', encoding='utf-8') as f_out:
        for pos in sorted_positions:
            f_out.write(f"chrY\t{pos}\n")
