#!/usr/bin/env bash
# ============================================================
# Pipeline de Variantes Germinativas — Menu Interativo
# ============================================================
# Interface principal para execução do pipeline.
#
# Uso:
#   ./run_pipeline.sh                          # Menu interativo
#   ./run_pipeline.sh --full                   # Pipeline completo
#   ./run_pipeline.sh --config custom.yaml     # Config customizado
#   ./run_pipeline.sh --threads 8              # Override de threads
#   ./run_pipeline.sh --samples-dir /data/fq   # Override de diretório
# ============================================================

set -euo pipefail

# ---- Cores ANSI ----
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# ---- Variáveis ----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config/config.yaml"
PIPELINE_VERSION="2.0.0"
SNAKEMAKE_CORES=4
SNAKEMAKE_EXTRA=""
DRY_RUN=false

# ============================================================
# FUNÇÕES UTILITÁRIAS
# ============================================================

print_header() {
    clear
    echo -e "${CYAN}"
    echo "╔══════════════════════════════════════════════════════════╗"
    echo "║                                                          ║"
    echo "║ Pipeline de Variantes Germinativas  v${PIPELINE_VERSION} ║"
    echo "║                                                          ║"
    echo "║   GATK4 Best Practices • BWA-MEM • Ensembl VEP           ║"
    echo "║                                                          ║"
    echo "╚══════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

print_separator() {
    echo -e "${BLUE}──────────────────────────────────────────────────────────${NC}"
}

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[AVISO]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERRO]${NC} $1"
}

check_snakemake() {
    if ! command -v snakemake &>/dev/null; then
        log_error "Snakemake não encontrado no PATH."
        echo ""
        echo "  Instale com:"
        echo "    conda install -c conda-forge -c bioconda snakemake"
        echo ""
        echo "  Ou crie o ambiente completo:"
        echo "    conda env create -f environment.yaml"
        echo "    conda activate pipeline-variantes"
        echo ""
        return 1
    fi
    return 0
}

run_snakemake() {
    local target="$1"
    local description="$2"

    check_snakemake || return 1

    print_separator
    log_info "Executando: ${description}"
    echo ""

    local cmd="snakemake --cores ${SNAKEMAKE_CORES} --configfile ${CONFIG_FILE}"

    if [ "${DRY_RUN}" = true ]; then
        cmd="${cmd} -n -p"
        log_warn "Modo dry-run (simulação)"
    fi

    cmd="${cmd} --keep-going --printshellcmds --rerun-incomplete"

    if [ -n "${target}" ]; then
        cmd="${cmd} ${target}"
    fi

    if [ -n "${SNAKEMAKE_EXTRA}" ]; then
        cmd="${cmd} ${SNAKEMAKE_EXTRA}"
    fi

    echo -e "${CYAN}Comando: ${cmd}${NC}"
    echo ""

    eval "${cmd}"
    local exit_code=$?

    echo ""
    if [ ${exit_code} -eq 0 ]; then
        log_info "${description} — concluído com sucesso ✅"
    else
        log_error "${description} — falhou com código ${exit_code} ❌"
        log_info "Verifique os logs em: logs/"
        log_info "Para retomar: execute novamente a mesma opção."
    fi

    echo ""
    read -rp "Pressione ENTER para voltar ao menu..."
}

# ============================================================
# OPÇÕES DO MENU
# ============================================================

menu_reference() {
    run_snakemake \
        "--until bwa_index samtools_faidx gatk_dict" \
        "Configuração da Referência (indexação)"
}

menu_quality() {
    run_snakemake \
        "--until fastqc multiqc" \
        "Controle de Qualidade (FastQC + MultiQC)"
}

menu_alignment() {
    run_snakemake \
        "--until samtools_flagstat samtools_stats" \
        "Alinhamento (BWA-MEM + MarkDuplicates + Stats)"
}

menu_calling() {
    run_snakemake \
        "--until genotype_gvcfs" \
        "Chamada de Variantes (BQSR + HaplotypeCaller + GenotypeGVCFs)"
}

menu_annotation() {
    run_snakemake \
        "--until vep_annotate_vcf vep_annotate_tsv" \
        "Anotação (Ensembl VEP)"
}

menu_full() {
    run_snakemake \
        "" \
        "Pipeline Completo"
}

menu_database() {
    print_separator
    echo -e "${BOLD}Consulta ao Banco de Dados${NC}"
    echo ""

    local db_file="${SCRIPT_DIR}/database/pipeline.db"

    if [ ! -f "${db_file}" ]; then
        log_warn "Banco de dados não encontrado: ${db_file}"
        log_info "O banco é criado automaticamente durante a execução do pipeline."
        echo ""
        read -rp "Pressione ENTER para voltar ao menu..."
        return
    fi

    echo "  1 - Listar todas as amostras"
    echo "  2 - Buscar amostra por nome"
    echo "  3 - Resumo geral"
    echo "  4 - Últimas execuções"
    echo "  0 - Voltar"
    echo ""
    read -rp "Opção: " db_option

    case "${db_option}" in
        1)
            echo ""
            sqlite3 -header -column "${db_file}" \
                "SELECT sample_name, seq_type, status, start_time, total_duration FROM samples ORDER BY start_time DESC LIMIT 50;"
            ;;
        2)
            read -rp "Nome da amostra: " sample_query
            echo ""
            sqlite3 -header -column "${db_file}" \
                "SELECT s.sample_name, s.status, s.start_time, s.total_duration, st.step_name, st.status as step_status, st.duration
                 FROM samples s LEFT JOIN steps st ON s.id = st.sample_id
                 WHERE s.sample_name LIKE '%${sample_query}%'
                 ORDER BY s.start_time DESC, st.id;"
            ;;
        3)
            echo ""
            echo -e "${BOLD}Resumo do Pipeline:${NC}"
            sqlite3 -header -column "${db_file}" \
                "SELECT
                    (SELECT COUNT(*) FROM pipeline_runs) as total_execucoes,
                    (SELECT COUNT(*) FROM samples) as total_amostras,
                    (SELECT COUNT(*) FROM samples WHERE status='SUCCESS') as sucesso,
                    (SELECT COUNT(*) FROM samples WHERE status='ERROR') as erros;"
            ;;
        4)
            echo ""
            sqlite3 -header -column "${db_file}" \
                "SELECT id, version, hostname, username, start_time, end_time, status, total_samples
                 FROM pipeline_runs ORDER BY start_time DESC LIMIT 10;"
            ;;
        *)
            return
            ;;
    esac

    echo ""
    read -rp "Pressione ENTER para voltar ao menu..."
}

menu_logs() {
    print_separator
    echo -e "${BOLD}Logs do Pipeline${NC}"
    echo ""

    echo "  1 - Log geral do pipeline"
    echo "  2 - Log de uma amostra específica"
    echo "  3 - Últimas 50 linhas do log geral"
    echo "  0 - Voltar"
    echo ""
    read -rp "Opção: " log_option

    case "${log_option}" in
        1)
            if [ -f "${SCRIPT_DIR}/logs/pipeline.log" ]; then
                less "${SCRIPT_DIR}/logs/pipeline.log"
            else
                log_warn "Log geral não encontrado."
            fi
            ;;
        2)
            read -rp "Nome da amostra: " sample_name
            local sample_log="${SCRIPT_DIR}/results/${sample_name}/logs/${sample_name}.log"
            if [ -f "${sample_log}" ]; then
                less "${sample_log}"
            else
                log_warn "Log não encontrado: ${sample_log}"
                echo "Logs disponíveis:"
                find "${SCRIPT_DIR}/results" -name "*.log" 2>/dev/null | head -20
            fi
            ;;
        3)
            if [ -f "${SCRIPT_DIR}/logs/pipeline.log" ]; then
                tail -50 "${SCRIPT_DIR}/logs/pipeline.log"
            else
                log_warn "Log geral não encontrado."
            fi
            ;;
        *)
            return
            ;;
    esac

    echo ""
    read -rp "Pressione ENTER para voltar ao menu..."
}

menu_config() {
    print_separator
    echo -e "${BOLD}Configurações${NC}"
    echo ""
    echo "  Arquivo de configuração: ${CONFIG_FILE}"
    echo ""
    echo "  1 - Visualizar configuração atual"
    echo "  2 - Editar configuração (vim)"
    echo "  3 - Editar configuração (nano)"
    echo "  4 - Alterar número de cores/threads"
    echo "  5 - Validar configuração (dry-run)"
    echo "  6 - Detectar amostras"
    echo "  0 - Voltar"
    echo ""
    read -rp "Opção: " config_option

    case "${config_option}" in
        1)
            echo ""
            cat "${CONFIG_FILE}"
            ;;
        2)
            vim "${CONFIG_FILE}"
            ;;
        3)
            nano "${CONFIG_FILE}"
            ;;
        4)
            read -rp "Número de cores para o Snakemake (atual: ${SNAKEMAKE_CORES}): " new_cores
            if [[ "${new_cores}" =~ ^[0-9]+$ ]] && [ "${new_cores}" -gt 0 ]; then
                SNAKEMAKE_CORES="${new_cores}"
                log_info "Cores alterado para: ${SNAKEMAKE_CORES}"
            else
                log_error "Valor inválido."
            fi
            ;;
        5)
            DRY_RUN=true
            run_snakemake "" "Validação (dry-run)"
            DRY_RUN=false
            ;;
        6)
            echo ""
            check_snakemake || return
            python3 "${SCRIPT_DIR}/scripts/parse_samples.py" \
                --samples-dir "$(grep 'samples:' "${CONFIG_FILE}" | head -1 | awk '{print $2}' | tr -d '\"')"
            ;;
        *)
            return
            ;;
    esac

    echo ""
    read -rp "Pressione ENTER para voltar ao menu..."
}

# ============================================================
# PARSING DE ARGUMENTOS CLI
# ============================================================

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --config)
                CONFIG_FILE="$2"
                shift 2
                ;;
            --threads|--cores)
                SNAKEMAKE_CORES="$2"
                shift 2
                ;;
            --dry-run|-n)
                DRY_RUN=true
                shift
                ;;
            --full)
                check_snakemake || exit 1
                menu_full
                exit 0
                ;;
            --reference)
                check_snakemake || exit 1
                menu_reference
                exit 0
                ;;
            --quality)
                check_snakemake || exit 1
                menu_quality
                exit 0
                ;;
            --alignment)
                check_snakemake || exit 1
                menu_alignment
                exit 0
                ;;
            --calling)
                check_snakemake || exit 1
                menu_calling
                exit 0
                ;;
            --annotation)
                check_snakemake || exit 1
                menu_annotation
                exit 0
                ;;
            --help|-h)
                echo "Uso: $0 [OPÇÕES]"
                echo ""
                echo "Opções:"
                echo "  --config FILE     Arquivo de configuração (default: config/config.yaml)"
                echo "  --threads N       Número de cores para Snakemake"
                echo "  --dry-run, -n     Simulação sem executar"
                echo "  --full            Executar pipeline completo"
                echo "  --reference       Apenas indexação da referência"
                echo "  --quality         Apenas controle de qualidade"
                echo "  --alignment       Apenas alinhamento"
                echo "  --calling         Apenas chamada de variantes"
                echo "  --annotation      Apenas anotação"
                echo "  --help, -h        Mostrar esta mensagem"
                exit 0
                ;;
            *)
                log_error "Argumento desconhecido: $1"
                echo "Use --help para ver as opções disponíveis."
                exit 1
                ;;
        esac
    done
}

# ============================================================
# LOOP PRINCIPAL DO MENU
# ============================================================

main() {
    parse_args "$@"

    while true; do
        print_header

        echo -e "${BOLD}  Selecione uma opção:${NC}"
        echo ""
        echo -e "  ${GREEN}1${NC} — Configurar Referência (indexação do genoma)"
        echo -e "  ${GREEN}2${NC} — Controle de Qualidade (FastQC + MultiQC)"
        echo -e "  ${GREEN}3${NC} — Alinhamento (BWA-MEM + MarkDuplicates)"
        echo -e "  ${GREEN}4${NC} — Chamada de Variantes (BQSR + HaplotypeCaller)"
        echo -e "  ${GREEN}5${NC} — Anotação (Ensembl VEP)"
        echo ""
        echo -e "  ${CYAN}6${NC} — Executar Pipeline Completo"
        echo ""
        echo -e "  ${YELLOW}7${NC} — Consultar Banco de Dados"
        echo -e "  ${YELLOW}8${NC} — Ver Logs"
        echo -e "  ${YELLOW}9${NC} — Configurações"
        echo ""
        echo -e "  ${RED}0${NC} — Sair"
        echo ""
        print_separator
        echo ""
        read -rp "  Opção: " choice

        case "${choice}" in
            1) menu_reference ;;
            2) menu_quality ;;
            3) menu_alignment ;;
            4) menu_calling ;;
            5) menu_annotation ;;
            6) menu_full ;;
            7) menu_database ;;
            8) menu_logs ;;
            9) menu_config ;;
            0)
                echo ""
                log_info "Até a próxima! 🧬"
                exit 0
                ;;
            *)
                log_error "Opção inválida: ${choice}"
                sleep 1
                ;;
        esac
    done
}

main "$@"
