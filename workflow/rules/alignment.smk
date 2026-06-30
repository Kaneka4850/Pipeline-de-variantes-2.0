# ============================================================
# Pipeline de Variantes Germinativas — Módulo: Alignment
# ============================================================
# Responsável por:
#   - BWA-MEM (alinhamento com pipe para samtools)
#   - Samtools sort + index
#   - GATK MarkDuplicates
#   - Flagstat + Stats (estatísticas)
#
# Suporta amostras single-end e paired-end dinamicamente.
# ============================================================


def _get_bwa_input_cmd(wildcards):
    """Monta o comando de input do BWA-MEM baseado no tipo da amostra."""
    info = SAMPLES_DICT[wildcards.sample]
    if info["type"] == "PE":
        return f"{info['r1']} {info['r2']}"
    return info["r1"]


rule bwa_mem:
    """Alinha reads contra a referência com BWA-MEM.

    Usa pipes para evitar arquivos SAM intermediários:
    BWA-MEM → samtools view → samtools sort → BAM ordenado.
    """
    input:
        fastqs=get_fastqs,
        genome=GENOME,
        # Garante que os índices existem
        bwa_idx=f"{GENOME}.amb",
    output:
        bam=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.bam",
    log:
        f"{LOGS_DIR}/{{sample}}/bwa_mem.log",
    params:
        bwa=BWA,
        samtools=SAMTOOLS,
        rg_id=lambda wc: f"{wc.sample}_RG",
        rg_sm=lambda wc: wc.sample,
        rg_pl=RG_PLATFORM,
        rg_lb=RG_LIBRARY,
        rg_pu=RG_PU,
        fastq_cmd=_get_bwa_input_cmd,
        bwa_threads=lambda wildcards, input, output, threads, resources: max(1, threads - 2),
        sort_mem=lambda wildcards, input, output, threads, resources: f"{max(768, resources.mem_mb // 4)}M",
    threads: THREADS
    resources:
        mem_mb=MEMORY_MB,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam})
        echo "[$(date)] Iniciando alinhamento BWA-MEM: {wildcards.sample}" > {log}
        echo "[$(date)] Read Group: @RG\tID:{params.rg_id}\tSM:{params.rg_sm}\tPL:{params.rg_pl}\tLB:{params.rg_lb}\tPU:{params.rg_pu}" >> {log}

        {params.bwa} mem \
            -t {params.bwa_threads} \
            -R "@RG\tID:{params.rg_id}\tSM:{params.rg_sm}\tPL:{params.rg_pl}\tLB:{params.rg_lb}\tPU:{params.rg_pu}" \
            {input.genome} \
            {params.fastq_cmd} \
            2>> {log} \
        | {params.samtools} view -@ 1 -b - 2>> {log} \
        | {params.samtools} sort -@ 1 -m {params.sort_mem} -o {output.bam} - 2>> {log}

        echo "[$(date)] Alinhamento concluído: {output.bam}" >> {log}
        """


rule samtools_index_sorted:
    """Indexa o BAM ordenado."""
    input:
        bam=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.bam",
    output:
        bai=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.bam.bai",
    log:
        f"{LOGS_DIR}/{{sample}}/samtools_index_sorted.log",
    params:
        samtools=SAMTOOLS,
    threads: THREADS
    shell:
        """
        echo "[$(date)] Indexando BAM: {input.bam}" > {log}
        {params.samtools} index -@ {threads} {input.bam} 2>> {log}
        echo "[$(date)] Indexação concluída" >> {log}
        """


rule mark_duplicates:
    """Marca duplicatas com GATK MarkDuplicates."""
    input:
        bam=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.bam",
        bai=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.bam.bai",
    output:
        bam=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.markdup.bam",
        bai=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.markdup.bam.bai",
        metrics=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.markdup_metrics.txt",
    log:
        f"{LOGS_DIR}/{{sample}}/mark_duplicates.log",
    params:
        gatk=GATK,
    threads: THREADS
    resources:
        mem_mb=MEMORY_MB,
    shell:
        """
        echo "[$(date)] Marcando duplicatas: {wildcards.sample}" > {log}

        {params.gatk} MarkDuplicates \
            --java-options "-Xmx{resources.mem_mb}m" \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --CREATE_INDEX true \
            2>> {log}
            
        mv {output.bam:s}.bai {output.bai} 2>> {log} || true

        echo "[$(date)] Duplicatas marcadas: {output.bam}" >> {log}
        """


rule samtools_flagstat:
    """Gera estatísticas de alinhamento com flagstat."""
    input:
        bam=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.markdup.bam",
        bai=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.markdup.bai",
    output:
        flagstat=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.flagstat.txt",
    log:
        f"{LOGS_DIR}/{{sample}}/flagstat.log",
    params:
        samtools=SAMTOOLS,
    threads: THREADS
    shell:
        """
        echo "[$(date)] Gerando flagstat: {wildcards.sample}" > {log}
        {params.samtools} flagstat -@ {threads} {input.bam} > {output.flagstat} 2>> {log}
        echo "[$(date)] Flagstat concluído" >> {log}
        """


rule samtools_stats:
    """Gera estatísticas detalhadas de alinhamento."""
    input:
        bam=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.markdup.bam",
        bai=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.markdup.bai",
        genome=GENOME,
    output:
        stats=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.stats.txt",
    log:
        f"{LOGS_DIR}/{{sample}}/samtools_stats.log",
    params:
        samtools=SAMTOOLS,
    threads: THREADS
    shell:
        """
        echo "[$(date)] Gerando estatísticas: {wildcards.sample}" > {log}
        {params.samtools} stats \
            -@ {threads} \
            --reference {input.genome} \
            {input.bam} > {output.stats} 2>> {log}
        echo "[$(date)] Estatísticas concluídas" >> {log}
        """
