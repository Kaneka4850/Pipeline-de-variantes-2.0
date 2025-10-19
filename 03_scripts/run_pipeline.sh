# --- CONFIGURAÇÃO DE SEGURANÇA ---
# (Para o script imediatamente se qualquer comando falhar)
set -e

# --- 1. DEFINIÇÃO DE VARIÁVEIS PRINCIPAIS ---

# Caminho para o genoma de referência
GENOME_REF="/root/pipeline_variantes/02_genome_reference/hg19/hg19.fa"

# Ferramentas a serem utilizadas
GATK="gatk"
SAMTOOLS="samtools"
BWA="bwa"
FASTQC="fastqc"

# --- 2. DEFINIÇÃO DAS PASTAS DO PROJETO ---

RAW_DATA_DIR="01_raw_data"
QC_DIR="04_quality_control"
ALIGN_DIR="05_alignment_results"
VARIANT_DIR="06_variant_calling"
ANNOT_DIR="07_annotated_variants"
REPORT_DIR="08_final_reports"
LOG_DIR="logs"

echo ">>> [SETUP] Verificando e criando diretórios de saída..."
mkdir -p $QC_DIR
mkdir -p $ALIGN_DIR
mkdir -p $VARIANT_DIR
mkdir -p $ANNOT_DIR
mkdir -p $REPORT_DIR
mkdir -p $LOG_DIR
echo ">>> [SETUP] Estrutura de diretórios garantida."
echo ""

# --- 3. INFORMAÇÕES DA AMOSTRA (Recebidas via Argumentos) ---
echo ">>> Lendo informações da amostra da linha de comando..."

# Verificação de segurança: Checa se os 3 argumentos foram passados
# $# é a contagem de argumentos. -ne 3 significa "não é igual a 3"
if [ "$#" -ne 3 ]; then
    echo "---------------------------------------------------------------------"
    echo "Erro: Uso incorreto do script."
    echo "Este pipeline precisa de 3 argumentos para rodar:"
    echo ""
    echo "Uso: $0 <NomeDaAmostra> <Caminho_FASTQ_R1> <Caminho_FASTQ_R2>"
    echo ""
    echo "Exemplo:"
    echo "$0 Amostra_ABC 01_raw_data/amostra_R1.fq.gz 01_raw_data/amostra_R2.fq.gz"
    echo "---------------------------------------------------------------------"
    exit 1 # Para o script com um código de erro
fi

# $0 é o nome do script. $1, $2, $3 são os argumentos.
SAMPLE_NAME="$1"
FASTQ_R1="$2"
FASTQ_R2="$3"

# --- Informações de Read Group (RG) - ESSENCIAL para o GATK ---

RG_ID="${SAMPLE_NAME}_RG"  # ID único do grupo (ex: Amostra_ABC_RG)
RG_SM="$SAMPLE_NAME"       # Nome da Amostra (ex: Amostra_ABC)
RG_PL="ILLUMINA"           # Plataforma de sequenciamento
RG_LB="lib1"               # Nome da Biblioteca (pode ser um placeholder)

# Construindo a string completa do Read Group para o BWA
RG_HEADER="@RG\\tID:${RG_ID}\\tSM:${RG_SM}\\tPL:${RG_PL}\\tLB:${RG_LB}"

echo ">>> [CONFIG] Amostra: $SAMPLE_NAME"
echo ">>> [CONFIG] FASTQ R1: $FASTQ_R1"
echo ">>> [CONFIG] FASTQ R2: $FASTQ_R2"
echo ">>> [CONFIG] Read Group: $RG_HEADER"

# --- FIM DA SEÇÃO 3 ---

echo ">>> Variáveis do pipeline definidas."
echo ">>> Genoma de Referência: $GENOME_REF"

# --- 4. PASSO 1: CONTROLE DE QUALIDADE (FASTQC) ---
echo ""
echo ">>> [PASSO 1] Iniciando FastQC nos dados brutos..."

$FASTQC -o $QC_DIR $FASTQ_R1
$FASTQC -o $QC_DIR $FASTQ_R2

echo ">>> [PASSO 1] FastQC concluído. Relatórios em $QC_DIR"


# --- 5. PASSO 2: ALINHAMENTO E PÓS-PROCESSAMENTO ---
echo ""
echo ">>> [PASSO 2] Iniciando Alinhamento (BWA) e Processamento (Samtools)..."

# Definindo nomes dos arquivos de saída
BAM_SAIDA="${ALIGN_DIR}/${SAMPLE_NAME}.bam"
BAM_ORDENADO="${ALIGN_DIR}/${SAMPLE_NAME}.sorted.bam"

# O 'pipe' (|) é a mágica aqui.
# 1. BWA-MEM alinha e gera um SAM (que vai para o 'pipe')
# 2. SAMTOOLS view converte o SAM (que vem do 'pipe') para BAM (que vai para o 'pipe')
# 3. SAMTOOLS sort ordena o BAM (que vem do 'pipe') e salva no arquivo final.
#
# O arquivo .sam gigante NUNCA é salvo no disco.
# O -t 4 usa 4 threads/núcleos. Ajuste se necessário.

echo ">>> Alinhando $SAMPLE_NAME com BWA-MEM..."
$BWA mem -t 4 -R "$RG_HEADER" $GENOME_REF $FASTQ_R1 $FASTQ_R2 | \
    $SAMTOOLS view -@ 4 -bS - | \
    $SAMTOOLS sort -@ 4 -o $BAM_ORDENADO -

echo ">>> Alinhamento e ordenação concluídos: $BAM_ORDENADO"

# 6. Indexando o BAM ordenado (necessário para os próximos passos)
echo ">>> Indexando o arquivo BAM ordenado..."
$SAMTOOLS index $BAM_ORDENADO

echo ">>> [PASSO 2] Alinhamento e pós-processamento concluídos."

# --- 6. PASSO 3: MARCAR DUPLICATAS (GATK) ---
echo ""
echo ">>> [PASSO 3] Marcando duplicatas com GATK MarkDuplicates..."

# O BAM ordenado da etapa anterior
BAM_ORDENADO="${ALIGN_DIR}/${SAMPLE_NAME}.sorted.bam"

# Definindo arquivos de saída para esta etapa
BAM_MARCADO="${ALIGN_DIR}/${SAMPLE_NAME}.sorted.marked.bam"
METRICS_FILE="${ALIGN_DIR}/${SAMPLE_NAME}.marked_dup_metrics.txt"

$GATK MarkDuplicates \
    -I $BAM_ORDENADO \
    -O $BAM_MARCADO \
    -M $METRICS_FILE

echo ">>> Duplicatas marcadas. Saída: $BAM_MARCADO"

# 7. Indexar o BAM com duplicatas marcadas
# (O HaplotypeCaller precisa que o BAM final esteja indexado)
echo ">>> Indexando o BAM com duplicatas marcadas..."
$SAMTOOLS index $BAM_MARCADO

echo ">>> [PASSO 3] Marcação de duplicatas e indexação concluídas."

# --- FIM DO CABEÇALHO ---
# --- 8. PASSO 4: CHAMADA DE VARIANTES (GATK HaplotypeCaller) ---
echo ""
echo ">>> [PASSO 4] Iniciando chamada de variantes com GATK HaplotypeCaller..."

# BAM final da etapa anterior
BAM_MARCADO="${ALIGN_DIR}/${SAMPLE_NAME}.sorted.marked.bam"

# Arquivo de saída VCF (sem compressão, como solicitado)
VCF_SAIDA="${VARIANT_DIR}/${SAMPLE_NAME}.vcf"

$GATK HaplotypeCaller \
    -R $GENOME_REF \
    -I $BAM_MARCADO \
    -O $VCF_SAIDA

echo ">>> [PASSO 4] Chamada de variantes concluída. Saída: $VCF_SAIDA"

