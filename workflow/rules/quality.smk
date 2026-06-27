# ============================================================
# Pipeline de Variantes Germinativas — Módulo: Quality
# ============================================================
# Responsável pelo controle de qualidade:
#   - FastQC para cada arquivo FASTQ
#   - MultiQC para agregar relatórios (se disponível)
# ============================================================


rule fastqc:
    """Executa FastQC nos arquivos FASTQ de uma amostra."""
    input:
        fastqs=get_fastqs,
    output:
        # Flag file indicando conclusão (FastQC gera nomes dinâmicos)
        flag=f"{RESULTS_DIR}/{{sample}}/qc/{{sample}}_fastqc_complete.flag",
    log:
        f"{LOGS_DIR}/{{sample}}/fastqc.log",
    params:
        fastqc=FASTQC,
        outdir=f"{RESULTS_DIR}/{{sample}}/qc",
    threads: 2
    resources:
        mem_mb=2000,
    shell:
        """
        mkdir -p {params.outdir}
        echo "[$(date)] Iniciando FastQC para: {wildcards.sample}" > {log}

        for fq in {input.fastqs}; do
            echo "[$(date)] Processando: $fq" >> {log}
            {params.fastqc} \
                --outdir {params.outdir} \
                --threads {threads} \
                --quiet \
                "$fq" \
                2>> {log}
        done

        echo "[$(date)] FastQC concluído para: {wildcards.sample}" >> {log}
        touch {output.flag}
        """


rule multiqc:
    """Agrega todos os relatórios FastQC com MultiQC."""
    input:
        # Espera todos os FastQC terminarem
        flags=expand(
            f"{RESULTS_DIR}/{{sample}}/qc/{{sample}}_fastqc_complete.flag",
            sample=SAMPLES,
        ),
    output:
        report=f"{RESULTS_DIR}/multiqc/multiqc_report.html",
    log:
        f"{LOGS_DIR}/multiqc.log",
    params:
        multiqc=MULTIQC,
        qc_dirs=" ".join(
            f"{RESULTS_DIR}/{s}/qc" for s in SAMPLES
        ),
        outdir=f"{RESULTS_DIR}/multiqc",
    threads: 1
    resources:
        mem_mb=2000,
    run:
        import shutil
        import subprocess

        # Verifica se MultiQC está disponível
        if not shutil.which(params.multiqc):
            # MultiQC não instalado — cria output placeholder
            os.makedirs(params.outdir, exist_ok=True)
            with open(output.report, "w") as f:
                f.write("<html><body>")
                f.write("<h1>MultiQC não disponível</h1>")
                f.write("<p>Instale o MultiQC para agregar relatórios de qualidade.</p>")
                f.write("<p>pip install multiqc</p>")
                f.write("</body></html>")
            with open(log[0], "w") as logf:
                logf.write("[AVISO] MultiQC não encontrado. Placeholder criado.\n")
        else:
            shell(
                """
                mkdir -p {params.outdir}
                echo "[$(date)] Executando MultiQC" > {log}
                {params.multiqc} \
                    {params.qc_dirs} \
                    --outdir {params.outdir} \
                    --force \
                    --quiet \
                    2>> {log}
                echo "[$(date)] MultiQC concluído" >> {log}
                """
            )
