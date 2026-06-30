# ============================================================
# Pipeline de Variantes Germinativas — Módulo: Intervals
# ============================================================
# Responsável por preparar o arquivo BED com as regiões alvo.
# Se for WES, usa o BED configurado pelo usuário.
# Se for WGS, gera janelas otimizadas a partir do FAI.
# ============================================================

INTERVALS_BED = config.get("pipeline", {}).get("intervals", "")
GENERATE_INTERVALS_PY = os.path.join(workflow.basedir, "scripts/generate_intervals.py")

rule get_intervals:
    """Obtém ou gera o arquivo BED de intervalos."""
    input:
        fai=f"{GENOME}.fai"
    output:
        bed=f"{RESULTS_DIR}/intervals/targets.bed"
    log:
        f"{LOGS_DIR}/intervals/get_intervals.log"
    params:
        user_bed=INTERVALS_BED,
        script=GENERATE_INTERVALS_PY
    shell:
        """
        mkdir -p $(dirname {output.bed})
        echo "[$(date)] Preparando intervalos BED" > {log}
        
        if [ -n "{params.user_bed}" ] && [ -f "{params.user_bed}" ]; then
            echo "[$(date)] Usando BED fornecido pelo usuário: {params.user_bed}" >> {log}
            cp {params.user_bed} {output.bed}
        else
            echo "[$(date)] BED não fornecido. Gerando janelas de 50MB a partir do FAI para WGS" >> {log}
            python {params.script} {input.fai} {output.bed} --window 50000000 2>> {log}
        fi
        echo "[$(date)] Concluído: {output.bed}" >> {log}
        """
