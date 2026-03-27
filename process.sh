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
  if [ -f "$VCF_OUT" ]; then
    continue
  fi

  echo "Обработка образца: $SAMPLE"

  if ! fastp -i "$FASTQ" -o "${SAMPLE}_tmp.fq.gz" \
    --thread "$THREADS" \
    --length_required 25 \
    --trim_poly_g \
    --json /dev/null --html /dev/null 2>"${SAMPLE}_f.log"; then
    continue
  fi

  if [ ! -s "${SAMPLE}_tmp.fq.gz" ]; then
    rm -f "${SAMPLE}_tmp.fq.gz" "${SAMPLE}_f.log"
    continue
  fi

  bwa mem -t "$THREADS" -k 17 -T 20 -M "$REF" "${SAMPLE}_tmp.fq.gz" |
    samtools view -@ "$THREADS" -u -h -L <(echo -e "$Y_NAME\t0\t$Y_MAX") - |
    samtools sort -@ "$THREADS" -m "$MEM_SORT" - |
    samtools markdup -@ "$THREADS" -r - "${SAMPLE}_tmp.bam"

  samtools index "${SAMPLE}_tmp.bam"

  bcftools mpileup -Ou -A -B -Q 13 -a FORMAT/AD,FORMAT/DP -f "$REF" -T "$TARGETS" -r "$Y_NAME" "${SAMPLE}_tmp.bam" |
    bcftools call -m -f GQ -A -Ov -o "$VCF_OUT"

  rm -f "${SAMPLE}_tmp.fq.gz" "${SAMPLE}_tmp.bam" "${SAMPLE}_tmp.bam.bai" "${SAMPLE}_f.log"
done

if [ -f "reaper.py" ]; then
  python3 reaper.py
fi
