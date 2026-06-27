# ============================================================
# Pipeline de Variantes Germinativas — Módulo: Reporting
# ============================================================
# Gera os relatórios JSON e HTML de execução ao final do pipeline.
# ============================================================

rule execution_report:
    """Gera relatório de reprodutibilidade em JSON e HTML."""
    input:
        # Pega a conclusão de todas as outras regras do pipeline
        # Usamos uma flag touch para garantir que rode apenas no fim
        flag=f"{RESULTS_DIR}/pipeline_completed.flag"
    output:
        json=f"{RESULTS_DIR}/execution_report.json",
        html=f"{RESULTS_DIR}/execution_report.html"
    log:
        f"{LOGS_DIR}/reporting/execution_report.log"
    params:
        db_path=f"{DB_DIR}/pipeline.db",
        script=f"{workflow.basedir}/scripts/generate_report.py"
    threads: 1
    shell:
        """
        echo "[$(date)] Iniciando geração de relatório" > {log}
        
        python3 {params.script} \
            --db {params.db_path} \
            --output $(dirname {output.json}) \
            2>> {log}
            
        echo "[$(date)] Relatório gerado com sucesso" >> {log}
        """
