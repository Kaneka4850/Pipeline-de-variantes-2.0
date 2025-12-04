#!/bin/bash
# --- SCRIPT DE ANOTAÇÃO DE VARIANTES (VEP) ---
# Autor: Cleber - Versão Final Corrigida
# Data: 26/11/2025

set -e

# --- 1. CONFIGURAÇÃO DE CAMINHOS ---
VEP_EXEC="$HOME/ensembl-vep-release-115/vep"
CACHE_DIR="$HOME/.vep"

# CAMINHO EXATO DO FASTA (Baseado no seu 'find')
FASTA_FILE="/root/.vep/homo_sapiens/115_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"

# --- 2. ARGUMENTOS ---
if [ "$#" -ne 2 ]; then
    echo "Uso: $0 <NomeDaAmostra> <Arquivo_VCF_Entrada>"
    exit 1
fi

SAMPLE_NAME="$1"
INPUT_VCF="$2"

OUTPUT_DIR="07_annotated_variants"
mkdir -p $OUTPUT_DIR

OUTPUT_VCF="${OUTPUT_DIR}/${SAMPLE_NAME}.annotated.vcf"
STATS_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_vep_summary.html"

echo ">>> [VEP] Iniciando anotação para: $SAMPLE_NAME"
echo ">>> FASTA utilizado: $FASTA_FILE"

# --- 3. EXECUÇÃO ---
# Aqui está a correção: adicionamos a linha --fasta $FASTA_FILE explicitamente

$VEP_EXEC \
    --input_file $INPUT_VCF \
    --output_file $OUTPUT_VCF \
    --format vcf \
    --vcf \
    --compress_output bgzip \
    --force_overwrite \
    --assembly GRCh37 \
    --cache \
    --dir_cache $CACHE_DIR \
    --merged \
    --offline \
    --fasta $FASTA_FILE \
    --everything \
    --fork 4 \
    --stats_file $STATS_FILE \
    --verbose

echo ">>> [VEP] Sucesso! Relatório: $STATS_FILE"
