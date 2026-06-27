# Changelog

Todas as mudanças notáveis deste projeto serão documentadas neste arquivo.

O formato segue [Keep a Changelog](https://keepachangelog.com/pt-BR/1.0.0/).

## [2.0.0] - 2026-06-26

### Adicionado
- **Snakemake** como orquestrador de pipeline (substitui scripts Bash monolíticos)
- **Módulos independentes** para cada etapa: reference, quality, alignment, caller, annotator
- **Detecção automática** de amostras (single-end e paired-end)
- **Configuração centralizada** via YAML (`config/config.yaml`)
- **BQSR condicional** — BaseRecalibrator + ApplyBQSR quando `known_sites` disponíveis
- **HaplotypeCaller em modo GVCF** — permite joint calling futuro
- **GenotypeGVCFs** — genotipagem individual e joint (via GenomicsDBImport)
- **Sistema de logging profissional** — log por amostra + log geral
- **Banco SQLite de rastreabilidade** — execuções, etapas, arquivos, checksums, versões
- **Validações robustas** — FASTQ, referência, disco, permissões, ferramentas
- **Menu interativo Bash** com 10 opções
- **Docker** — Dockerfile multi-stage + docker-compose.yml
- **Testes automatizados** — pytest para validação, amostras, banco, config
- **Suporte a VEP plugins** — ClinVar, gnomAD, REVEL, SpliceAI, CADD (desabilitados por padrão)
- **Checksums MD5** para arquivos críticos
- **MultiQC** — relatório agregado de qualidade (condicional)
- **Flagstat + samtools stats** — estatísticas de alinhamento
- **Ambientes Conda** por módulo
- **Schemas JSON** para validação do config
- **Documentação completa** em português

### Removido
- Scripts Bash monolíticos (`03_scripts/run_pipeline.sh`, `03_scripts/run_annotation.sh`)
- Estrutura de diretórios numerada (`01_raw_data/`, `03_scripts/`, etc.)
- Caminhos absolutos hardcoded

### Corrigido
- FASTA_FILE vazio no script de anotação (era não funcional)
- Ausência de BQSR no pipeline original
- HaplotypeCaller sem modo GVCF
- Threads hardcoded (agora configurável)
- Read Groups incompletos (adicionado Platform Unit)

## [1.0.0] - 2024-XX-XX

### Versão Original
- Pipeline em 2 scripts Bash
- Suporte apenas a paired-end
- Caminhos hardcoded
- Sem logging, testes, Docker ou banco de dados
