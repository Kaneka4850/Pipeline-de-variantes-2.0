# ============================================================
# Pipeline de Variantes Germinativas — Snakefile Principal
# ============================================================
# Ponto de entrada do Snakemake. Carrega a configuração,
# detecta amostras automaticamente e inclui os módulos de
# regras para cada etapa do pipeline.
#
# Uso:
#   snakemake --cores 4 --configfile config/config.yaml
#   snakemake --cores 4 -n  (dry-run)
#   snakemake --cores 4 --until quality  (até uma etapa)
# ============================================================

import os
import sys
import pandas as pd
from pathlib import Path

# Adiciona o diretório de scripts ao PATH do Python
sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))

from common.utils import get_tool_path, detect_samples
from common.environment import capture_environment
from common.preflight import PreflightRunner
from common.reference_manager import ReferenceManager
from common.database import PipelineDatabase


# ============================================================
# 1. CONFIGURAÇÃO
# ============================================================

# Carrega o config padrão (pode ser sobrescrito via --configfile)
configfile: "config/config.yaml"

# Validação do config contra o schema
validate(config, schema="workflow/schemas/config.schema.yaml")


# ============================================================
# 2. VARIÁVEIS GLOBAIS
# ============================================================

# Diretórios
SAMPLES_DIR = config["dirs"]["samples"]
RESULTS_DIR = config["dirs"]["results"]
LOGS_DIR = config["dirs"]["logs"]
DB_DIR = config["dirs"]["database"]

# Referência
GENOME = config["reference"]["genome"]
KNOWN_SITES = config["reference"].get("known_sites", [])

# Recursos
THREADS = config["resources"]["threads"]
MEMORY_MB = config["resources"]["memory_mb"]

# Ferramentas
BWA = get_tool_path(config, "bwa")
GATK = get_tool_path(config, "gatk")
SAMTOOLS = get_tool_path(config, "samtools")
FASTQC = get_tool_path(config, "fastqc")
BCFTOOLS = get_tool_path(config, "bcftools")
MULTIQC = get_tool_path(config, "multiqc")
VEP = get_tool_path(config, "vep")

# Pipeline
PIPELINE_VERSION = config["pipeline"]["version"]
PIPELINE_NAME = config["pipeline"]["name"]
JOINT_CALLING = config["pipeline"].get("joint_calling", False)
ANNOTATION_ENABLED = config["pipeline"].get("annotation", True)

# Read Groups
RG_PLATFORM = config["read_groups"]["platform"]
RG_LIBRARY = config["read_groups"]["library"]
RG_PU = config["read_groups"]["platform_unit"]

# VEP
VEP_CACHE = config["vep"].get("cache_dir", "")
VEP_PLUGINS = config["vep"].get("plugins_dir", "")
VEP_ASSEMBLY = config["vep"].get("assembly", "GRCh38")
VEP_EXTRA = config["vep"].get("extra_flags", "")


# ============================================================
# 3. DETECÇÃO AUTOMÁTICA DE AMOSTRAS
# ============================================================

# Detecta amostras no diretório configurado
SAMPLES_DICT = detect_samples(
    samples_dir=SAMPLES_DIR,
    r1_patterns=config["sample_detection"]["r1_patterns"],
    r2_patterns=config["sample_detection"]["r2_patterns"],
)

SAMPLES = list(SAMPLES_DICT.keys())

if not SAMPLES:
    print(
        f"[ERRO] Nenhuma amostra encontrada em '{SAMPLES_DIR}/'.\n"
        f"       Coloque seus arquivos FASTQ nesse diretório e tente novamente.",
        file=sys.stderr,
    )
    sys.exit(1)

print(f"[INFO] {len(SAMPLES)} amostra(s) detectada(s): {', '.join(SAMPLES)}")
for name, info in SAMPLES_DICT.items():
    seq_type = info["type"]
    print(f"       → {name} ({seq_type}): R1={info['r1']}", end="")
    if info.get("r2"):
        print(f", R2={info['r2']}")
    else:
        print()


# ============================================================
# 4. FUNÇÕES AUXILIARES
# ============================================================

def get_fastq_r1(wildcards):
    """Retorna o caminho do FASTQ R1 para uma amostra."""
    return SAMPLES_DICT[wildcards.sample]["r1"]


def get_fastq_r2(wildcards):
    """Retorna o caminho do FASTQ R2 para uma amostra (PE only)."""
    return SAMPLES_DICT[wildcards.sample].get("r2", "")


def get_fastqs(wildcards):
    """Retorna lista de FASTQs para uma amostra (1 para SE, 2 para PE)."""
    info = SAMPLES_DICT[wildcards.sample]
    if info["type"] == "PE":
        return [info["r1"], info["r2"]]
    return [info["r1"]]


def is_paired(wildcards):
    """Verifica se a amostra é paired-end."""
    return SAMPLES_DICT[wildcards.sample]["type"] == "PE"


def get_bqsr_bam(wildcards):
    """Retorna o BAM final para chamada de variantes.
    Se known_sites estão configurados, retorna o BAM pós-BQSR.
    Caso contrário, retorna o BAM com duplicatas marcadas.
    """
    if KNOWN_SITES:
        return f"{RESULTS_DIR}/{wildcards.sample}/alignment/{wildcards.sample}.recal.bam"
    return f"{RESULTS_DIR}/{wildcards.sample}/alignment/{wildcards.sample}.sorted.markdup.bam"


# ============================================================
# 5. DEFINIÇÃO DOS OUTPUTS FINAIS (REGRA ALL)
# ============================================================

def get_final_targets():
    """Define todos os outputs finais desejados pelo pipeline."""
    targets = []

    for sample in SAMPLES:
        # Referência indexada (sempre necessário)
        targets.append(f"{GENOME}.fai")

        # QC
        targets.append(f"{RESULTS_DIR}/{sample}/qc/{sample}_fastqc_complete.flag")

        # Alinhamento — estatísticas
        targets.append(f"{RESULTS_DIR}/{sample}/alignment/{sample}.flagstat.txt")
        targets.append(f"{RESULTS_DIR}/{sample}/alignment/{sample}.stats.txt")

        # Variantes
        targets.append(
            f"{RESULTS_DIR}/{sample}/variants/{sample}.vcf.gz"
        )
        targets.append(
            f"{RESULTS_DIR}/{sample}/variants/{sample}.vcf.gz.tbi"
        )

        # Anotação (se habilitada)
        if ANNOTATION_ENABLED and VEP:
            targets.append(
                f"{RESULTS_DIR}/{sample}/annotation/{sample}.annotated.vcf.gz"
            )
            targets.append(
                f"{RESULTS_DIR}/{sample}/annotation/{sample}.annotated.tsv"
            )
            targets.append(
                f"{RESULTS_DIR}/{sample}/annotation/{sample}_vep_summary.html"
            )

    # MultiQC (se disponível — verificado na regra)
    targets.append(f"{RESULTS_DIR}/multiqc/multiqc_report.html")

    # Flag de conclusão
    targets.append(f"{RESULTS_DIR}/pipeline_completed.flag")

    # Report (se habilitado)
    if config.get("report", {}).get("enabled", True):
        targets.append(f"{RESULTS_DIR}/execution_report.json")
        targets.append(f"{RESULTS_DIR}/execution_report.html")

    return targets

# Regra intermediária para sincronizar o relatório
rule pipeline_completed_flag:
    input:
        [t for t in get_final_targets() if not t.endswith(".flag") and not t.endswith(".json") and not t.endswith(".html")]
    output:
        touch(f"{RESULTS_DIR}/pipeline_completed.flag")


rule all:
    input:
        get_final_targets(),


# ============================================================
# 6. INCLUSÃO DOS MÓDULOS DE REGRAS
# ============================================================

include: "workflow/rules/reference.smk"
include: "workflow/rules/quality.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/caller.smk"
include: "workflow/rules/annotator.smk"
include: "workflow/rules/reporting.smk"


# ============================================================
# 7. HOOKS — EVENTOS DO PIPELINE
# ============================================================

onstart:
    print(f"\n{'=' * 60}")
    print(f"  {PIPELINE_NAME} v{PIPELINE_VERSION}")
    print(f"  Amostras: {len(SAMPLES)}")
    print(f"  Threads por tarefa: {THREADS}")
    print(f"  Referência: {GENOME}")
    print(f"  BQSR: {'Habilitado' if KNOWN_SITES else 'Desabilitado (sem known_sites)'}")
    print(f"  Joint Calling: {'Sim' if JOINT_CALLING else 'Não'}")
    print(f"  Anotação VEP: {'Sim' if ANNOTATION_ENABLED and VEP else 'Não'}")
    print(f"{'=' * 60}\n")

    # ---- PREFLIGHT ----
    if config.get("preflight", {}).get("enabled", True):
        print("[PREFLIGHT] Iniciando verificações...")
        runner = PreflightRunner(config, SAMPLES_DICT)
        results = runner.run_all()
        print(runner.summary())
        
        if runner.has_failures:
            print("\n[ERRO CRÍTICO] O preflight encontrou falhas que impedem a execução segura.")
            # Dependendo da política, podemos dar sys.exit(1) aqui.
            # Como foi pedido para não interromper tudo, apenas avisamos e deixamos o Snakemake lidar com as falhas das regras.
        
        # Salva resultados no banco
        db = PipelineDatabase(f"{DB_DIR}/pipeline.db")
        # Para isso precisamos do run_id, que foi criado mas está perdido no estado do Snakemake.
        # Vamos assumir que a última execução é a nossa (já que a db foi instanciada acima).
        try:
            row = db.conn.execute("SELECT id FROM pipeline_runs ORDER BY id DESC LIMIT 1").fetchone()
            if row:
                run_id = row["id"]
                db.register_preflight(run_id, [r.to_dict() for r in results])
                
                # Registra ambiente
                snapshot = capture_environment(config)
                db.register_environment(run_id, snapshot.to_dict())
                
                # Atualiza pipeline_runs
                sys_info = snapshot.system
                import json
                db.update_run_environment(
                    run_id,
                    f"{sys_info.os_name} {sys_info.os_version}",
                    sys_info.cpu_count,
                    sys_info.memory_total_mb,
                    json.dumps(snapshot.config_dump)
                )
                
                # Registra VEP
                if ANNOTATION_ENABLED and VEP:
                    db.register_vep_info(run_id, snapshot.vep.to_dict())
                    
                # Registra Referência
                ref_manager = ReferenceManager(GENOME, config)
                ref_info = ref_manager.collect_info()
                db.register_reference(run_id, ref_info.to_dict())
                
        except Exception as e:
            print(f"[AVISO] Não foi possível registrar dados no banco durante o preflight: {e}")
        finally:
            db.close()


onsuccess:
    print(f"\n{'=' * 60}")
    print(f"  ✅ Pipeline concluído com sucesso!")
    print(f"  Resultados em: {RESULTS_DIR}/")
    print(f"{'=' * 60}\n")


onerror:
    print(f"\n{'=' * 60}")
    print(f"  ❌ Pipeline encontrou erros.")
    print(f"  Verifique os logs em: {LOGS_DIR}/")
    print(f"  Use 'snakemake --cores N' para retomar.")
    print(f"{'=' * 60}\n")
