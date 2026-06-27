# ============================================================
# Pipeline de Variantes Germinativas — Módulo: Annotator
# ============================================================
# Anotação de variantes com Ensembl VEP CLI.
#
# Gera:
#   - VCF anotado (.annotated.vcf.gz)
#   - TSV (.annotated.tsv)
#   - Relatório HTML de estatísticas
#
# Compatível com futuros plugins (ClinVar, gnomAD, OMIM)
# via seção vep.plugins no config.yaml.
# ============================================================


def _build_vep_cmd(wildcards):
    """Monta flags adicionais do VEP baseadas no config."""
    flags = []

    # Cache
    if VEP_CACHE:
        flags.append(f"--dir_cache {VEP_CACHE}")

    # Plugins directory
    if VEP_PLUGINS:
        flags.append(f"--dir_plugins {VEP_PLUGINS}")

    # Assembly
    flags.append(f"--assembly {VEP_ASSEMBLY}")

    # Extra flags do usuário
    if VEP_EXTRA:
        flags.append(VEP_EXTRA)

    # Plugins configurados
    plugins_config = config.get("vep", {}).get("plugins", {})
    if isinstance(plugins_config, dict):
        for plugin_name, plugin_conf in plugins_config.items():
            if isinstance(plugin_conf, dict) and plugin_conf.get("enabled", False):
                if plugin_name == "clinvar" and plugin_conf.get("file"):
                    flags.append(f"--plugin ClinVar,{plugin_conf['file']}")
                elif plugin_name == "revel" and plugin_conf.get("file"):
                    flags.append(f"--plugin REVEL,{plugin_conf['file']}")
                elif plugin_name == "cadd":
                    snv = plugin_conf.get("snv_file", "")
                    indel = plugin_conf.get("indel_file", "")
                    if snv:
                        flags.append(f"--plugin CADD,{snv},{indel}")
                elif plugin_name == "spliceai":
                    snv = plugin_conf.get("snv_file", "")
                    indel = plugin_conf.get("indel_file", "")
                    if snv:
                        flags.append(f"--plugin SpliceAI,snv={snv},indel={indel}")

    return " ".join(flags)


rule vep_annotate_vcf:
    """Anota variantes com Ensembl VEP — output VCF."""
    input:
        vcf=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.vcf.gz",
        tbi=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.vcf.gz.tbi",
        genome=GENOME,
    output:
        vcf=f"{RESULTS_DIR}/{{sample}}/annotation/{{sample}}.annotated.vcf.gz",
        stats=f"{RESULTS_DIR}/{{sample}}/annotation/{{sample}}_vep_summary.html",
    log:
        f"{LOGS_DIR}/{{sample}}/vep_annotate_vcf.log",
    params:
        vep=VEP,
        extra_flags=_build_vep_cmd,
    threads: THREADS
    resources:
        mem_mb=MEMORY_MB,
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        echo "[$(date)] Iniciando anotação VEP (VCF): {wildcards.sample}" > {log}

        {params.vep} \
            --input_file {input.vcf} \
            --output_file {output.vcf} \
            --format vcf \
            --vcf \
            --compress_output bgzip \
            --force_overwrite \
            --cache \
            --offline \
            --fasta {input.genome} \
            --fork {threads} \
            --stats_file {output.stats} \
            --no_progress \
            {params.extra_flags} \
            2>> {log}

        echo "[$(date)] Anotação VEP (VCF) concluída" >> {log}
        """


rule vep_annotate_tsv:
    """Anota variantes com Ensembl VEP — output TSV."""
    input:
        vcf=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.vcf.gz",
        tbi=f"{RESULTS_DIR}/{{sample}}/variants/{{sample}}.vcf.gz.tbi",
        genome=GENOME,
    output:
        tsv=f"{RESULTS_DIR}/{{sample}}/annotation/{{sample}}.annotated.tsv",
    log:
        f"{LOGS_DIR}/{{sample}}/vep_annotate_tsv.log",
    params:
        vep=VEP,
        extra_flags=_build_vep_cmd,
    threads: THREADS
    resources:
        mem_mb=MEMORY_MB,
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        echo "[$(date)] Iniciando anotação VEP (TSV): {wildcards.sample}" > {log}

        {params.vep} \
            --input_file {input.vcf} \
            --output_file {output.tsv} \
            --format vcf \
            --tab \
            --force_overwrite \
            --cache \
            --offline \
            --fasta {input.genome} \
            --fork {threads} \
            --no_stats \
            --no_progress \
            {params.extra_flags} \
            2>> {log}

        echo "[$(date)] Anotação VEP (TSV) concluída" >> {log}
        """
