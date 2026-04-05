#!/bin/bash

THREADS=$(nproc)
TOTAL_RAM=$(free -g | awk '/^Mem:/ {print $2}')
if [ "$TOTAL_RAM" -lt 16 ]; then
  MEM_SORT="1G"
else
  MEM_SORT="2G"
fi

REF="hg38.fa"
INPUT_DIR="./input"
OUTPUT_DIR="./output"
RESULTS_DIR="./results"
TARGETS="targets.tsv"
REF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz"
TREE_URL="https://www.familytreedna.com/public/y-dna-haplotree/get"
Y_MAX=57227415

if [ ! -d "$INPUT_DIR" ]; then
  mkdir -p "$INPUT_DIR"
fi
if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi
if [ ! -d "$RESULTS_DIR" ]; then
  mkdir -p "$RESULTS_DIR"
fi

for pkg in python3 bwa samtools bcftools wget gunzip fastp; do
  if ! command -v "$pkg" &>/dev/null; then
    sudo apt update && sudo apt install -y "$pkg"
  fi
done

if [ ! -f "$REF" ]; then
  if [ ! -f "${REF}.gz" ]; then
    wget -c -L "$REF_URL" -O "${REF}.gz"
  fi
  gunzip "${REF}.gz"
fi

if [ ! -f "${REF}.fai" ]; then
  samtools faidx "$REF"
fi
if [ ! -f "${REF}.bwt" ]; then
  bwa index "$REF"
fi
if [ ! -f "tree.json" ]; then
  wget -q -L "$TREE_URL" -U "Mozilla/5.0" -O "tree.json"
fi

if [ ! -f "$TARGETS" ]; then
  if [ -f "converter.py" ]; then
    python3 converter.py
  else
    exit 1
  fi
fi

Y_NAME=$(grep ">" "$REF" | grep -m 1 -oE "chrY|Y" | head -n 1)
if [ -z "$Y_NAME" ]; then
  Y_NAME="chrY"
fi
sed -i "s/^[^[:space:]]*/$Y_NAME/" "$TARGETS"

shopt -s nullglob
for FASTQ in "$INPUT_DIR"/*.fastq.gz; do
  SAMPLE=$(basename "${FASTQ%.fastq.gz}")
  VCF_OUT="$OUTPUT_DIR/${SAMPLE}_Y.vcf"
  REPORT_FILE="./results/${SAMPLE}_report.txt"
  if [ -f "$REPORT_FILE" ]; then
    echo "Пропуск: Отчет для $SAMPLE уже существует."
    continue
  fi

  echo "Обработка образца: $SAMPLE"

  if ! fastp -i "$FASTQ" -o "${OUTPUT_DIR}/${SAMPLE}_tmp.fq.gz" \
    --thread "$THREADS" \
    --length_required 25 \
    --trim_poly_g \
    --json /dev/null --html /dev/null 2>"${OUTPUT_DIR}/${SAMPLE}_f.log"; then
    continue
  fi

  if [ ! -s "${OUTPUT_DIR}/${SAMPLE}_tmp.fq.gz" ]; then
    rm -f "${OUTPUT_DIR}/${SAMPLE}_tmp.fq.gz" "${OUTPUT_DIR}/${SAMPLE}_f.log"
    continue
  fi

  bwa mem -t "$THREADS" -k 17 -T 20 -M "$REF" "${OUTPUT_DIR}/${SAMPLE}_tmp.fq.gz" |
    samtools view -@ "$THREADS" -u -h -L <(echo -e "$Y_NAME\t0\t$Y_MAX") - |
    samtools sort -@ "$THREADS" -m "$MEM_SORT" - |
    samtools markdup -@ "$THREADS" -r - "${OUTPUT_DIR}/${SAMPLE}_tmp.bam"

  samtools index "${OUTPUT_DIR}/${SAMPLE}_tmp.bam"

  bcftools mpileup -Ou -A -B -Q 13 -a FORMAT/AD,FORMAT/DP -f "$REF" -T "$TARGETS" -r "$Y_NAME" "${OUTPUT_DIR}/${SAMPLE}_tmp.bam" |
    bcftools call -m -f GQ -A -Ov -o "$VCF_OUT"

  if [ -f "snp_reaper.py" ] && [ -f "$VCF_OUT" ]; then
    python3 snp_reaper.py "$VCF_OUT"
  fi

  rm -f "${OUTPUT_DIR}/${SAMPLE}_tmp.fq.gz" "${OUTPUT_DIR}/${SAMPLE}_tmp.bam" "${OUTPUT_DIR}/${SAMPLE}_tmp.bam.bai" "${OUTPUT_DIR}/${SAMPLE}_f.log" "$VCF_OUT"
done
