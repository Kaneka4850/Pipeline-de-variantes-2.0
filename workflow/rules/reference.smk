# ============================================================
# Pipeline de Variantes Germinativas — Módulo: Reference
# ============================================================
# Responsável por indexar o genoma de referência:
#   - BWA index (gera .amb, .ann, .bwt, .pac, .sa)
#   - Samtools faidx (gera .fai)
#   - GATK CreateSequenceDictionary (gera .dict)
#
# O Snakemake só executa as regras se os outputs não existirem
# ou estiverem desatualizados em relação ao input.
# ============================================================


rule bwa_index:
    """Indexa o genoma de referência com BWA."""
    input:
        genome=GENOME,
    output:
        amb=f"{GENOME}.amb",
        ann=f"{GENOME}.ann",
        bwt=f"{GENOME}.bwt",
        pac=f"{GENOME}.pac",
        sa=f"{GENOME}.sa",
    log:
        f"{LOGS_DIR}/reference/bwa_index.log",
    params:
        bwa=BWA,
    threads: 1
    resources:
        mem_mb=MEMORY_MB,
    shell:
        """
        echo "[$(date)] Iniciando indexação BWA: {input.genome}" > {log}
        {params.bwa} index {input.genome} 2>> {log}
        echo "[$(date)] Indexação BWA concluída" >> {log}
        """


rule samtools_faidx:
    """Cria o índice .fai do genoma com samtools."""
    input:
        genome=GENOME,
    output:
        fai=f"{GENOME}.fai",
    log:
        f"{LOGS_DIR}/reference/samtools_faidx.log",
    params:
        samtools=SAMTOOLS,
    threads: 1
    shell:
        """
        echo "[$(date)] Criando índice fai: {input.genome}" > {log}
        {params.samtools} faidx {input.genome} 2>> {log}
        echo "[$(date)] Índice fai criado" >> {log}
        """


rule gatk_dict:
    """Cria o dicionário de sequências com GATK."""
    input:
        genome=GENOME,
    output:
        # GATK gera o .dict trocando a extensão do FASTA
        dict_file=GENOME.rsplit(".", 1)[0] + ".dict",
    log:
        f"{LOGS_DIR}/reference/gatk_dict.log",
    params:
        gatk=GATK,
    threads: 1
    resources:
        mem_mb=MEMORY_MB,
    shell:
        """
        echo "[$(date)] Criando dicionário: {output.dict_file}" > {log}
        {params.gatk} CreateSequenceDictionary \
            -R {input.genome} \
            -O {output.dict_file} \
            2>> {log}
        echo "[$(date)] Dicionário criado" >> {log}
        """
