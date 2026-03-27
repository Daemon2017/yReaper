import json
import os
import sys

from utils import prepare_data, analyze_sample

JSON_PATH = 'tree.json'
VCF_DIR = './output'
RESULTS_DIR = './results'

if __name__ == "__main__":
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

    if not os.path.exists(VCF_DIR):
        print(f"Ошибка: Директория {VCF_DIR} не найдена!")
        sys.exit(1)

    vcf_files = sorted([f for f in os.listdir(VCF_DIR) if f.endswith("_Y.vcf")])

    if not vcf_files:
        print("VCF файлы для анализа не найдены.")
        sys.exit(1)

    print(f"Найдено файлов для анализа: {len(vcf_files)}")

    for file in vcf_files:
        sample_name = file.replace("_Y.vcf", "").replace(".vcf", "")
        report_path = os.path.join(RESULTS_DIR, f"{sample_name}_report.txt")
        vcf_path = os.path.join(VCF_DIR, file)

        print(f"Анализ: {sample_name}...")

        try:
            with open(report_path, 'w', encoding='utf-8') as f_report:
                analyze_sample(vcf_path, pos_m, node_s, anc_m, f_report)
        except Exception as e:
            print(f"Ошибка при обработке {sample_name}: {e}")

    print(f"Готово. Отчеты сохранены в папке: {RESULTS_DIR}")
