# ============================================================
# Pipeline de Variantes Germinativas — Módulo: Caller
# ============================================================
# Pipeline GATK completo para chamada de variantes germinativas:
#   - BaseRecalibrator + ApplyBQSR (condicional: requer known_sites)
#   - HaplotypeCaller em modo GVCF
#   - GenotypeGVCFs (amostra individual ou joint calling)
#   - Indexação VCF com tabix
#
# Se known_sites não está configurado, o BQSR é automaticamente
# ignorado e o HaplotypeCaller opera diretamente sobre o BAM
# com duplicatas marcadas.
# ============================================================


# ============================================================
# BQSR (condicional)
# ============================================================

if KNOWN_SITES:

    rule base_recalibrator:
        """Calcula tabela de recalibração de bases (BQSR)."""
        input:
            bam=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.markdup.bam",
            bai=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.markdup.bai",
            genome=GENOME,
            genome_fai=f"{GENOME}.fai",
            genome_dict=GENOME.rsplit(".", 1)[0] + ".dict",
            known_sites=KNOWN_SITES,
        output:
            table=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.recal_data.table",
        log:
            f"{LOGS_DIR}/{{sample}}/base_recalibrator.log",
        params:
            gatk=GATK,
            known_sites_args=lambda wc: " ".join(
                f"--known-sites {ks}" for ks in KNOWN_SITES
            ),
        threads: THREADS
        resources:
            mem_mb=MEMORY_MB,
        shell:
            """
            echo "[$(date)] Iniciando BaseRecalibrator: {wildcards.sample}" > {log}

            {params.gatk} BaseRecalibrator \
                --java-options "-Xmx{resources.mem_mb}m" \
                -R {input.genome} \
                -I {input.bam} \
                {params.known_sites_args} \
                -O {output.table} \
                2>> {log}

            echo "[$(date)] BaseRecalibrator concluído" >> {log}
            """


    rule apply_bqsr:
        """Aplica a recalibração de bases ao BAM."""
        input:
            bam=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.markdup.bam",
            bai=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.sorted.markdup.bai",
            table=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.recal_data.table",
            genome=GENOME,
            genome_fai=f"{GENOME}.fai",
            genome_dict=GENOME.rsplit(".", 1)[0] + ".dict",
        output:
            bam=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.recal.bam",
            bai=f"{RESULTS_DIR}/{{sample}}/alignment/{{sample}}.recal.bai",
        log:
            f"{LOGS_DIR}/{{sample}}/apply_bqsr.log",
        params:
            gatk=GATK,
        threads: THREADS
        resources:
            mem_mb=MEMORY_MB,
        shell:
            """
            echo "[$(date)] Aplicando BQSR: {wildcards.sample}" > {log}

            {params.gatk} ApplyBQSR \
                --java-options "-Xmx{resources.mem_mb}m" \
                -R {input.genome} \
                -I {input.bam} \
                --bqsr-recal-file {input.table} \
                -O {output.bam} \
                2>> {log}

            echo "[$(date)] BQSR aplicado: {output.bam}" >> {log}
            """


# ============================================================
# HAPLOTYPE CALLER
# ============================================================

rule haplotype_caller:
    """Chama variantes com GATK HaplotypeCaller em modo GVCF."""
    input:
        bam=get_bqsr_bam,
        genome=GENOME,
        genome_fai=f"{GENOME}.fai",
        genome_dict=GENOME.rsplit(".", 1)[0] + ".dict",
        bed=f"{RESULTS_DIR}/intervals/targets.bed",
    output:
        gvcf=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.g.vcf.gz",
        tbi=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.g.vcf.gz.tbi",
    log:
        f"{LOGS_DIR}/{{sample}}/haplotype_caller.log",
    params:
        gatk=GATK,
    threads: THREADS
    resources:
        mem_mb=MEMORY_MB,
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.gvcf})
        echo "[$(date)] Iniciando HaplotypeCaller (GVCF): {wildcards.sample}" > {log}
        
        MAPPED=$(samtools view -c -F 4 {input.bam} || echo 0)
        if [ "$MAPPED" -lt 1000 ]; then
            echo "ERRO: BAM com menos de 1000 reads alinhados ($MAPPED reads)." >> {log}
            exit 1
        fi

        {params.gatk} HaplotypeCaller \
            --java-options "-Xmx{resources.mem_mb}m" \
            -R {input.genome} \
            -I {input.bam} \
            -L {input.bed} \
            -O {output.gvcf} \
            -ERC GVCF \
            --native-pair-hmm-threads {threads} \
            2>> {log}

        echo "[$(date)] HaplotypeCaller concluído: {output.gvcf}" >> {log}
        """


# ============================================================
# GENOTYPE GVCFs
# ============================================================

rule genotype_gvcfs:
    """Genotipa GVCFs para produzir o VCF final.

    Para amostra individual: GenotypeGVCFs direto.
    Para joint calling: seria precedido por GenomicsDBImport (futuro).
    """
    input:
        gvcf=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.g.vcf.gz",
        tbi=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.g.vcf.gz.tbi",
        genome=GENOME,
        genome_fai=f"{GENOME}.fai",
        genome_dict=GENOME.rsplit(".", 1)[0] + ".dict",
        bed=f"{RESULTS_DIR}/intervals/targets.bed",
    output:
        vcf=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.raw.vcf.gz",
        tbi=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.raw.vcf.gz.tbi",
    log:
        f"{LOGS_DIR}/{{sample}}/genotype_gvcfs.log",
    params:
        gatk=GATK,
    threads: THREADS
    resources:
        mem_mb=MEMORY_MB,
    shell:
        """
        echo "[$(date)] Iniciando GenotypeGVCFs: {wildcards.sample}" > {log}

        {params.gatk} GenotypeGVCFs \
            --java-options "-Xmx{resources.mem_mb}m" \
            -R {input.genome} \
            -V {input.gvcf} \
            -L {input.bed} \
            -O {output.vcf} \
            2>> {log}

        echo "[$(date)] GenotypeGVCFs concluído: {output.vcf}" >> {log}
        """


rule variant_filtration:
    """Aplica Hard Filtering (VariantFiltration) no VCF individual."""
    input:
        vcf=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.raw.vcf.gz",
        tbi=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.raw.vcf.gz.tbi",
        genome=GENOME,
    output:
        vcf=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.vcf.gz",
        tbi=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.vcf.gz.tbi",
    log:
        f"{LOGS_DIR}/{{sample}}/variant_filtration.log",
    params:
        gatk=GATK,
    threads: 1
    resources:
        mem_mb=MEMORY_MB,
    shell:
        """
        echo "[$(date)] Iniciando VariantFiltration: {wildcards.sample}" > {log}

        {params.gatk} VariantFiltration \
            --java-options "-Xmx{resources.mem_mb}m" \
            -R {input.genome} \
            -V {input.vcf} \
            -O {output.vcf} \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            2>> {log}

        echo "[$(date)] VariantFiltration concluído: {output.vcf}" >> {log}
        """


# ============================================================
# JOINT CALLING (Opcional — ativado via config)
# ============================================================

if JOINT_CALLING and len(SAMPLES) > 1:

    rule genomics_db_import:
        """Cria banco GenomicsDB para joint calling."""
        input:
            gvcfs=expand(
                f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.g.vcf.gz",
                sample=SAMPLES,
            ),
            tbis=expand(
                f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.g.vcf.gz.tbi",
                sample=SAMPLES,
            ),
            genome=GENOME,
            genome_fai=f"{GENOME}.fai",
            genome_dict=GENOME.rsplit(".", 1)[0] + ".dict",
            bed=f"{RESULTS_DIR}/intervals/targets.bed",
        output:
            db=directory(f"{RESULTS_DIR}/joint/genomicsdb"),
        log:
            f"{LOGS_DIR}/joint/genomics_db_import.log",
        params:
            gatk=GATK,
            variant_args=lambda wc: " ".join(
                f"-V {RESULTS_DIR}/{s}/variants/{s}.g.vcf.gz"
                for s in SAMPLES
            ),
        threads: THREADS
        resources:
            mem_mb=MEMORY_MB,
        shell:
            """
            mkdir -p $(dirname {log})
            echo "[$(date)] Iniciando GenomicsDBImport" > {log}

            {params.gatk} GenomicsDBImport \
                --java-options "-Xmx{resources.mem_mb}m" \
                {params.variant_args} \
                -L {input.bed} \
                --genomicsdb-workspace-path {output.db} \
                2>> {log}

            echo "[$(date)] GenomicsDBImport concluído" >> {log}
            """


    rule joint_genotype_gvcfs:
        """Genotipa GVCFs conjuntamente a partir do GenomicsDB."""
        input:
            db=f"{RESULTS_DIR}/joint/genomicsdb",
            genome=GENOME,
            genome_fai=f"{GENOME}.fai",
            genome_dict=GENOME.rsplit(".", 1)[0] + ".dict",
            bed=f"{RESULTS_DIR}/intervals/targets.bed",
        output:
            vcf=f"{RESULTS_DIR}/joint/all_samples.vcf.gz",
            tbi=f"{RESULTS_DIR}/joint/all_samples.vcf.gz.tbi",
        log:
            f"{LOGS_DIR}/joint/joint_genotype_gvcfs.log",
        params:
            gatk=GATK,
        threads: THREADS
        resources:
            mem_mb=MEMORY_MB,
        shell:
            """
            echo "[$(date)] Iniciando Joint GenotypeGVCFs" > {log}

            {params.gatk} GenotypeGVCFs \
                --java-options "-Xmx{resources.mem_mb}m" \
                -R {input.genome} \
                -V gendb://{input.db} \
                -O {output.vcf} \
                2>> {log}

            echo "[$(date)] Joint GenotypeGVCFs concluído: {output.vcf}" >> {log}
            """
