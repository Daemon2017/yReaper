import os


def prepare_data(json_data):
    nodes = json_data.get('allNodes', {})
    pos_map, node_to_snps, anc_map = {}, {}, {}

    for node_id, node in nodes.items():
        name = node['name']
        node_to_snps[name] = []

        for v in node.get('variants', []):
            pos_val = v.get('position')
            if pos_val is not None:
                pos = str(pos_val)
                pos_map[pos] = {
                    'der': v.get('derived', ''),
                    'anc': v.get('ancestral', ''),
                    'var': v.get('variant', 'Unk'),
                    'node': name
                }
                node_to_snps[name].append(pos)

        path, curr = [], node
        while curr is not None:
            path.append(curr['name'])
            parent_id = curr.get('parentId')
            curr = nodes.get(str(parent_id)) if parent_id else None
        anc_map[name] = path[::-1]

    return pos_map, node_to_snps, anc_map


def get_vcf_calls(vcf_path):
    calls = {}
    if not os.path.exists(vcf_path):
        return calls

    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            cols = line.strip().split('\t')
            if len(cols) < 10:
                continue

            pos = cols[1]
            alts = cols[4].split(',')
            fmt_keys = cols[8].split(':')
            sample_vals = cols[9].split(':')
            fmt_dict = dict(zip(fmt_keys, sample_vals))

            if 'AD' in fmt_dict and fmt_dict['AD'] != '.':
                ad_values = fmt_dict['AD'].split(',')
                try:
                    ref_count = int(ad_values[0])
                    for i, alt_base in enumerate(alts):
                        alt_idx = i + 1
                        if alt_idx < len(ad_values):
                            alt_count = int(ad_values[alt_idx])
                            if pos not in calls:
                                calls[pos] = {}
                            calls[pos][alt_base] = {
                                'cnt': alt_count,
                                'tot': ref_count + alt_count
                            }
                except (ValueError, IndexError):
                    continue
    return calls


def analyze_sample(vcf_path, pos_map, node_to_snps, anc_map, out_file):
    vcf_calls = get_vcf_calls(vcf_path)
    found_nodes = set()

    for pos, alts_data in vcf_calls.items():
        if pos in pos_map:
            target = pos_map[pos]['der']
            if target in alts_data and alts_data[target]['cnt'] > 0:
                found_nodes.add(pos_map[pos]['node'])

    if not found_nodes:
        out_file.write(f"Образец {os.path.basename(vcf_path)}: Релевантных SNP не найдено.\n")
        return

    node_scores = {}
    for node in found_nodes:
        path = anc_map.get(node, [])
        score = 0
        for ancestor in path:
            for p in node_to_snps.get(ancestor, []):
                if p in vcf_calls:
                    target = pos_map[p]['der']
                    if target in vcf_calls[p] and vcf_calls[p][target]['cnt'] > 0:
                        score += 1
                        break
        node_scores[node] = score

    term_node = max(node_scores, key=node_scores.get)

    out_file.write(f"{'=' * 95}\n")
    out_file.write(f"Образец: {os.path.basename(vcf_path)} | Терминальный субклад: {term_node}\n")
    out_file.write(f"{'=' * 95}\n")
    out_file.write(f"{'Узел/уровень':<20} | {'SNP Статус (+/+!/-/?) [Мутация/Всего]'}\n")
    out_file.write(f"{'-' * 20}-|-{'-' * 72}\n")

    for node in anc_map[term_node]:
        res_list = []
        for pos in node_to_snps.get(node, []):
            info = pos_map[pos]
            if pos in vcf_calls:
                if info['der'] in vcf_calls[pos]:
                    data = vcf_calls[pos][info['der']]
                    dmg = ""
                    if data['cnt'] > 0:
                        status = "+"
                        if (info['anc'] == 'C' and info['der'] == 'T') or (info['anc'] == 'G' and info['der'] == 'A'):
                            dmg = "!"
                    else:
                        status = "-"
                    res_str = f"{info['var']}{status}{dmg} [{data['cnt']}/{data['tot']}]"
                else:
                    any_v = list(vcf_calls[pos].values())[0]
                    res_str = f"{info['var']}- [0/{any_v['tot']}]"
            else:
                res_str = f"{info['var']}?"
            res_list.append(res_str)

        sorted_res = sorted(res_list, key=lambda x: ('+' not in x, '-' not in x))
        out_file.write(f"{node:<20} | {' '.join(sorted_res)}\n")
