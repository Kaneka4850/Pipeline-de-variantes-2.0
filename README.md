# Pipeline de Chamada de Variantes Germinativas (FASTQ-para-VCF)

Este projeto é um pipeline de bioinformática robusto e reutilizável, escrito em Bash, que transforma dados brutos de sequenciamento (arquivos .fastq) em um arquivo de variantes (.vcf).

Pipeline de Chamada de Variantes Germinativas (FASTQ-para-VCF)
Este projeto é um pipeline de bioinformática robusto e reutilizável, escrito em Bash, que transforma dados brutos de sequenciamento (arquivos .fastq) em um arquivo de variantes (.vcf).

O fluxo de trabalho segue as Melhores Práticas (Best Practices) do GATK4 para a chamada de variantes germinativas, utilizando ferramentas padrão da indústria como BWA, Samtools e GATK.

Este pipeline foi desenhado para ser:

Generalista: Processa qualquer amostra (WGS, WES ou painel) através de argumentos de linha de comando.

Robusto: Inclui checagens de segurança (sanity checks) para verificar a existência de arquivos e pastas antes de iniciar o processamento.

Eficiente: Utiliza pipes (|) para evitar a criação de arquivos .sam intermediários gigantes, economizando espaço em disco e tempo de I/O.

Organizado: Mantém uma estrutura de diretórios limpa, separando dados brutos, referências, scripts e resultados.

Estrutura do Projeto

O pipeline espera e gera a seguinte estrutura de diretórios:

pipeline_variantes/
├── 01_raw_data/            # (Coloque seus arquivos .fastq/.fastq.gz aqui)
├── 02_genome_reference/    # (Coloque seu genoma de referência .fa e índices aqui)
│   └── hg19/
│       ├── hg19.fa
│       └── (outros índices .amb, .ann, .pac, .fai, .dict)
├── 03_scripts/
│   └── run_pipeline.sh     # (O script principal do pipeline)
├── 04_quality_control/     # (Saída dos relatórios do FastQC)
├── 05_alignment_results/   # (Saída dos arquivos .bam e .bai)
├── 06_variant_calling/     # (Saída do arquivo final .vcf)
├── 07_annotated_variants/  # (Reservado para anotação)
├── 08_final_reports/       # (Reservado para relatórios finais)
├── logs/                     # (Reservado para logs de execução)
└── README.md
Instalação e Configuração
Siga estes passos cuidadosamente para configurar o ambiente e as dependências.

## Instalação e Configuração

Siga estes passos cuidadosamente para configurar o ambiente e as dependências.

Pré-requisitos
Um sistema operacional Linux (ex: Ubuntu, CentOS) ou WSL2 no Windows.

O gerenciador de pacotes Miniconda instalado.

Bash

git clone [URL_DO_SEU_REPOSITORIO_GITHUB]
cd pipeline_variantes
2. Instale as Dependências (Ambiente Conda)
Todas as ferramentas de bioinformática são gerenciadas através de um ambiente Conda dedicado. Isso garante 100% de reprodutibilidade.

O comando abaixo irá criar um novo ambiente chamado `pipeline_variantes` com todas as ferramentas necessárias (BWA, Samtools, GATK4 e FastQC):

Bash

# (Certifique-se de ter os canais bioconda e conda-forge configurados)
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Crie o ambiente (isso pode levar alguns minutos)
conda create -n pipeline_variantes bwa samtools gatk4 fastqc

**Atenção**: Como o `bcftools` e `picard` foram mencionados no setup inicial mas não são usados ativamente no script final, eles foram omitidos por simplicidade. Se desejar, adicione-os ao comando `create`.

3. Prepare o Genoma de Referência (Passo ÚNICO!)
Este é o passo mais crucial e demorado. Ele só precisa ser feito uma vez.

Baixe seu genoma de referência em formato FASTA (ex: `hg19.fa`) e coloque-o na pasta `02_genome_reference/hg19/`.
Verifique o caminho do genoma dentro do script. Abra 03_scripts/run_pipeline.sh com um editor (ex: nano) e confirme se a variável GENOME_REF aponta para o local correto do seu arquivo .fa:

Bash

# Linha ~10 dentro de run_pipeline.sh
GENOME_REF="/root/pipeline_variantes/02_genome_reference/hg19/hg19.fa"
(Altere /root/ para seu caminho de usuário, se for diferente).

**(Altere `/root/` para seu caminho de usuário, se for diferente).**

Bash

# 1. Ative o ambiente
conda activate pipeline_variantes

# 2. Indexe para o BWA (MUITO DEMORADO: 1-2 horas)
bwa index /root/pipeline_variantes/02_genome_reference/hg19/hg19.fa

# 3. Indexe para o Samtools (Rápido)
samtools faidx /root/pipeline_variantes/02_genome_reference/hg19/hg19.fa

# 4. Crie o Dicionário para o GATK (Rápido)
gatk CreateSequenceDictionary -R /root/pipeline_variantes/02_genome_reference/hg19/hg19.fa

Seu ambiente está pronto!

Como Usar o Pipeline
Siga os passos abaixo toda vez que quiser processar uma nova amostra.

## Como Usar o Pipeline
1. Ative o Ambiente Conda
Sempre que você abrir um novo terminal, ative o ambiente para ter acesso às ferramentas:

Bash

conda activate pipeline_variantes
2. Adicione seus Dados Brutos (FASTQ)
Copie seus arquivos de sequenciamento (`.fastq` ou `.fastq.gz`) para a pasta `01_raw_data/`.

3. Execute o Pipeline
O script `run_pipeline.sh` espera três argumentos na linha de comando:

*   Um nome único para sua amostra (ex: `Amostra_01`).
*
O caminho completo para o arquivo FASTQ Read 1 (R1).

O caminho completo para o arquivo FASTQ Read 2 (R2).

Sintaxe do Comando:
Bash

```bash
Exemplo de Execução Real:
Bash

# (Certifique-se de estar na pasta raiz 'pipeline_variantes')

./03_scripts/run_pipeline.sh \
    cap-ngse-a-2019_S1 \
    01_raw_data/cap-ngse-a-2019_S1_L001_R1_001.fastq \
    01_raw_data/cap-ngse-a-2019_S1_L001_R2_001.fastq
```
O script irá executar todo o fluxo de trabalho automaticamente, exibindo o progresso de cada passo.

O que o Script Faz (Passo a Passo)

Verifica se todos os diretórios de saída (04_..., 05_..., etc.) existem. Se não, eles são criados.

Verifica se os arquivos R1 e R2 fornecidos realmente existem no disco. Se não, o script para com um erro claro.

[PASSO 1] Controle de Qualidade (FastQC):

Gera relatórios de qualidade para os arquivos FASTQ brutos.
## O que o Script Faz (Passo a Passo)
Saída: Relatórios .html na pasta 04_quality_control/.

[PASSO 2] Alinhamento e Pós-Processamento (BWA-MEM & Samtools):

Alinha as reads (R1 e R2) contra o genoma de referência.

Adiciona o Read Group (RG), essencial para o GATK.

Converte o alinhamento (SAM) para o formato binário (BAM).

Ordena o arquivo BAM por coordenadas genômicas.

Cria um índice (.bai) para o BAM ordenado.

Saída: Arquivos `.sorted.bam` e `.sorted.bam.bai` na pasta `05_alignment_results/`.

[PASSO 3] Marcar Duplicatas (GATK MarkDuplicates):

Lê o BAM ordenado e identifica duplicatas de PCR.

Gera um novo arquivo BAM "limpo" com as duplicatas marcadas.

Cria um índice para o novo BAM marcado.

Saída: Arquivos .sorted.marked.bam e .sorted.marked.bam.bai na pasta 05_alignment_results/.

[PASSO 4] Chamada de Variantes (GATK HaplotypeCaller):

Lê o BAM final (limpo e marcado) e o genoma de referência.

Analisa cada posição coberta e identifica onde a amostra difere da referência (SNPs e INDELs).

Saída (O Produto Final): Um arquivo `.vcf` contendo todas as variantes da amostra, salvo em `06_variant_calling/`.

Próximos Passos (Anotação)
Este pipeline entrega um arquivo VCF "cru" (sem anotação). O próximo passo lógico é anotar este VCF para entender o significado biológico de cada variante.

## Próximos Passos (Anotação)

Este pipeline entrega um arquivo VCF "cru" (sem anotação). O próximo passo lógico é anotar este VCF para entender o significado biológico de cada variante.

Um script de anotação separado (ex: `03_scripts/run_annotation_vep.sh`), que utiliza o Ensembl VEP (Variant Effect Predictor), está em desenvolvimento e será adicionado a este projeto.