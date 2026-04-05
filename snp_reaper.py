import json
import os
import sys

from utils import prepare_data, analyze_sample

JSON_PATH = 'tree.json'
VCF_DIR = './output'
RESULTS_DIR = './results'

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(1)

    vcf_path = sys.argv[1]

    if not os.path.exists(JSON_PATH):
        print(f"Ошибка: {JSON_PATH} не найден!")
        sys.exit(1)

    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)

    with open(JSON_PATH, 'r', encoding='utf-8') as f:
        try:
            tree_data = json.load(f)
        except json.JSONDecodeError:
            print("Ошибка: Некорректный формат JSON в tree.json!")
            sys.exit(1)

    pos_m, node_s, anc_m = prepare_data(tree_data)

    if not os.path.exists(vcf_path):
        print(f"Ошибка: Файл {vcf_path} не найден!")
        sys.exit(1)

    sample_name = os.path.basename(vcf_path).replace("_Y.vcf", "").replace(".vcf", "")
    report_path = os.path.join(RESULTS_DIR, f"{sample_name}_report.txt")

    print(f"Анализ: {sample_name}...")

    try:
        with open(report_path, 'w', encoding='utf-8') as f_report:
            analyze_sample(vcf_path, pos_m, node_s, anc_m, f_report)
        print(f"Готово. Отчет: {report_path}")
    except Exception as e:
        print(f"Ошибка при обработке {sample_name}: {e}")
