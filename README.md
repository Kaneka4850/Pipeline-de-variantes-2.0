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
# Pipeline de Chamada de Variantes Germinativas (FASTQ → VCF)

Um pipeline em Bash para transformar arquivos FASTQ (dados brutos de sequenciamento) em um VCF com variantes germinativas. Segue as Best Practices do GATK4 e usa ferramentas padrão da indústria: BWA, Samtools, GATK e FastQC.

Principais características
- Generalista: funciona com WGS, WES ou painéis — parâmetros via linha de comando.
- Robusto: checks para existência de arquivos/pastas antes de rodar.
- Eficiente: usa pipes onde possível para reduzir I/O e arquivos intermediários.
- Organizado: estrutura de pastas clara para dados, referências, scripts e resultados.

Índice
- Estrutura do projeto
- Requisitos
- Instalação (Conda)
- Preparar o genoma de referência (indexação)
- Uso (exemplos)
- O que o script faz (resumo dos passos)
- Próximos passos / Anotação

## Estrutura do projeto

O repositório espera a seguinte organização (nomes de pastas sugeridos):

pipeline_variantes/
├── 01_raw_data/            # arquivos .fastq / .fastq.gz
├── 02_genome_reference/    # referência .fa e índices
│   └── hg19/
│       ├── hg19.fa
│       └── (hg19.fa.amb, .ann, .pac, .fai, .dict, ...)
├── 03_scripts/
│   └── run_pipeline.sh     # script principal
├── 04_quality_control/     # saída do FastQC
├── 05_alignment_results/   # arquivos .bam e índices
├── 06_variant_calling/     # arquivo final .vcf
├── 07_annotated_variants/  # (opcional) arquivos anotados
├── 08_final_reports/       # relatórios finais/plots
├── logs/                   # logs de execução
└── README.md

## Requisitos
- Sistema: Linux (Ubuntu/CentOS) ou Windows + WSL2
- Miniconda/Conda
- Bash
- Git

Ferramentas usadas (instaladas via Conda): BWA, Samtools, GATK4, FastQC. Opcional: bcftools, picard.

## Instalação (ambiente Conda)

1. Configure canais Conda (uma vez):

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

2. Crie o ambiente com as ferramentas básicas:

```bash
conda create -n pipeline_variantes bwa samtools gatk4 fastqc
```

Observação: adicione `bcftools` e `picard` se desejar.

## Preparar o genoma de referência

Coloque o FASTA do genoma em `02_genome_reference/hg19/hg19.fa` (ou altere o caminho no script). Indexe o genoma uma vez antes de rodar o pipeline:

```bash
conda activate pipeline_variantes

# BWA (pode levar 1-2 horas para genomas grandes)
bwa index 02_genome_reference/hg19/hg19.fa

# Samtools
samtools faidx 02_genome_reference/hg19/hg19.fa

# GATK (cria dicionário)
gatk CreateSequenceDictionary -R 02_genome_reference/hg19/hg19.fa
```

Observação: se o script aponta para um caminho absoluto (ex: /root/...), atualize a variável `GENOME_REF` dentro de `03_scripts/run_pipeline.sh` para o caminho correto no seu sistema.

## Uso

1. Ative o ambiente Conda:

```bash
conda activate pipeline_variantes
```

2. Coloque seus arquivos FASTQ em `01_raw_data/`.

3. Execute o pipeline (exemplo):

```bash
./03_scripts/run_pipeline.sh \
    AMOSTRA_01 \
    01_raw_data/AMOSTRA_01_R1.fastq.gz \
    01_raw_data/AMOSTRA_01_R2.fastq.gz
```

O script espera: <sample_name> <path_R1> <path_R2>

## O que o script faz (resumo dos passos)

1) Verificações iniciais
- Confere existência de diretórios de saída e arquivos de entrada; cria pastas faltantes.

2) Controle de qualidade (FastQC)
- Gera relatórios HTML em `04_quality_control/`.

3) Alinhamento e pós-processamento (BWA-MEM + Samtools)
- Alinha R1/R2 contra a referência, adiciona Read Groups, converte para BAM, ordena e indexa.
- Saída: `05_alignment_results/*.sorted.bam` e `.bai`.

4) Marcar duplicatas (GATK MarkDuplicates)
- Produz `*.sorted.marked.bam` e índice correspondente.

5) Chamada de variantes (GATK HaplotypeCaller)
- Gera o VCF final em `06_variant_calling/`.

## Próximos passos / Anotação

O resultado atual é um VCF "cru". Para interpretar variantes biologicamente, recomenda-se anotar com ferramentas como Ensembl VEP, ANNOVAR ou SnpEff.

Esse pipeline também possui um script de anotação via VEP, para utiliza-lo, basta entrar por esse link e baixar a versão de linha de comando:
https://www.ensembl.org/info/docs/tools/vep/index.html

## Como usar o anotador
- O anotador é bem simples de utilizar, utilizando apenas um parametro, conforme explicado abaixo

```bash
./03_scripts/run_annotation.sh \
    AMOSTRA_01 \
    06_variant_calling/AMOSTRA.vcf
```

Com esse comando ele ira fazer toda a organização do pipeline, gerando um arquivo no formato html e outro txt para analise.
OBSERVAÇÕES: Por motivos de otimização, os arquivos referencia, como os .fasta e os plugins do VEP, como dbsnp, revel, spliceAI, dentre outros, não foram colocados nesse script, eles deverão ser instalado conforme a necessidade do script. E provavelmente sera necessario alterar os caminhos conforme sua maquina/ambiente de trabalho

## Dicas e notas
- Ajuste a variável `GENOME_REF` em `03_scripts/run_pipeline.sh` para apontar ao caminho do seu FASTA.
- Para grandes conjuntos de dados, execute em uma máquina com bastante RAM e CPU; BWA e GATK são CPU-intensivos.
- Mantenha backups dos arquivos originais e dos índices do genoma — a indexação é feita apenas uma vez.

## Licença e contato

Adicione aqui a licença do projeto (por exemplo, MIT) e informações de contato ou link para issues no GitHub.

---

Se quiser, já posso:
- ajustar o README para incluir badges (build/conda/status),
- traduzir comandos para PowerShell/WSL quando necessário,
- ou criar um exemplo mínimo de `run_pipeline.sh` para acompanhar o README.
