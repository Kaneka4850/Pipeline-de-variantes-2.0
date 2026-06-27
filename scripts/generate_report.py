# ============================================================
# Pipeline de Variantes Germinativas — Relatório de Execução
# ============================================================
"""
Gera relatório de reprodutibilidade em JSON e HTML.

Conteúdo:
- Data/hora início e fim, duração total
- Usuário, hostname, SO, CPUs, memória
- Versões de todas as ferramentas
- Referência (nome, versão, assembly, checksum, índices)
- VEP (versão, cache, assembly, plugins)
- Configuração completa utilizada
- Tempo por etapa e por amostra
- Status final, erros encontrados

Uso:
    python scripts/generate_report.py --db database/pipeline.db --output results/
"""

import argparse
import json
import sys
from datetime import datetime
from html import escape
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common.database import PipelineDatabase
from common.utils import format_duration


# ============================================================
# GERAÇÃO DO RELATÓRIO JSON
# ============================================================

def generate_report_data(db: PipelineDatabase, run_id: int | None = None) -> dict[str, Any]:
    """Coleta todos os dados para o relatório.

    Args:
        db: Instância do banco de dados.
        run_id: ID da execução (None = última execução).

    Returns:
        Dicionário com todos os dados do relatório.
    """
    # Obtém a execução
    if run_id is None:
        row = db.conn.execute(
            "SELECT * FROM pipeline_runs ORDER BY id DESC LIMIT 1"
        ).fetchone()
    else:
        row = db.conn.execute(
            "SELECT * FROM pipeline_runs WHERE id = ?", (run_id,)
        ).fetchone()

    if not row:
        return {"error": "Nenhuma execução encontrada"}

    run = dict(row)
    actual_run_id = run["id"]

    # Amostras
    samples = [
        dict(r) for r in db.conn.execute(
            "SELECT * FROM samples WHERE run_id = ? ORDER BY id",
            (actual_run_id,),
        ).fetchall()
    ]

    # Etapas por amostra
    for sample in samples:
        sample["steps"] = [
            dict(r) for r in db.conn.execute(
                "SELECT * FROM steps WHERE sample_id = ? ORDER BY id",
                (sample["id"],),
            ).fetchall()
        ]

    # Ferramentas
    tools = [
        dict(r) for r in db.conn.execute(
            "SELECT * FROM tools WHERE run_id = ?", (actual_run_id,)
        ).fetchall()
    ]

    # Referência (nova tabela)
    reference = None
    try:
        ref_row = db.conn.execute(
            "SELECT * FROM reference_metadata WHERE run_id = ?", (actual_run_id,)
        ).fetchone()
        if ref_row:
            reference = dict(ref_row)
    except Exception:
        pass

    # VEP info (nova tabela)
    vep_info = None
    try:
        vep_row = db.conn.execute(
            "SELECT * FROM vep_info WHERE run_id = ?", (actual_run_id,)
        ).fetchone()
        if vep_row:
            vep_info = dict(vep_row)
    except Exception:
        pass

    # Environment snapshot (nova tabela)
    env_snapshot = None
    try:
        env_row = db.conn.execute(
            "SELECT * FROM environment_snapshots WHERE run_id = ?", (actual_run_id,)
        ).fetchone()
        if env_row:
            env_snapshot = dict(env_row)
            if env_snapshot.get("snapshot_json"):
                env_snapshot["snapshot"] = json.loads(env_snapshot["snapshot_json"])
    except Exception:
        pass

    # Preflight results (nova tabela)
    preflight = []
    try:
        preflight = [
            dict(r) for r in db.conn.execute(
                "SELECT * FROM preflight_results WHERE run_id = ?", (actual_run_id,)
            ).fetchall()
        ]
    except Exception:
        pass

    # Calcula duração total
    duration = ""
    if run.get("start_time") and run.get("end_time"):
        try:
            start = datetime.fromisoformat(run["start_time"])
            end = datetime.fromisoformat(run["end_time"])
            seconds = (end - start).total_seconds()
            duration = format_duration(seconds)
        except ValueError:
            pass

    # Config (do environment snapshot)
    config_dump = {}
    if env_snapshot and env_snapshot.get("snapshot"):
        config_dump = env_snapshot["snapshot"].get("config", {})

    return {
        "report_generated_at": datetime.now().isoformat(),
        "pipeline": {
            "name": "Pipeline de Variantes Germinativas",
            "version": run.get("version", ""),
            "run_id": actual_run_id,
            "status": run.get("status", ""),
            "start_time": run.get("start_time", ""),
            "end_time": run.get("end_time", ""),
            "duration": duration,
            "total_samples": run.get("total_samples", 0),
        },
        "environment": {
            "hostname": run.get("hostname", ""),
            "username": run.get("username", ""),
            "os_info": run.get("os_info", ""),
            "cpu_count": run.get("cpu_count", 0),
            "memory_total_mb": run.get("memory_total_mb", 0),
        },
        "tools": tools,
        "reference": reference,
        "vep": vep_info,
        "preflight": preflight,
        "samples": samples,
        "config": config_dump,
    }


# ============================================================
# GERAÇÃO JSON
# ============================================================

def write_json_report(data: dict[str, Any], output_path: str) -> None:
    """Escreve o relatório em formato JSON."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False, default=str)


# ============================================================
# GERAÇÃO HTML
# ============================================================

def write_html_report(data: dict[str, Any], output_path: str) -> None:
    """Gera relatório HTML formatado para auditoria."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    pipeline = data.get("pipeline", {})
    env = data.get("environment", {})
    tools = data.get("tools", [])
    ref = data.get("reference", {}) or {}
    vep = data.get("vep", {}) or {}
    samples = data.get("samples", [])
    preflight = data.get("preflight", [])

    html = f"""<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Relatório de Execução — Pipeline de Variantes v{escape(str(pipeline.get('version', '')))}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ font-family: 'Segoe UI', system-ui, sans-serif; background: #0f172a; color: #e2e8f0; line-height: 1.6; padding: 2rem; }}
        .container {{ max-width: 1100px; margin: 0 auto; }}
        h1 {{ font-size: 1.8rem; color: #38bdf8; margin-bottom: 0.5rem; }}
        h2 {{ font-size: 1.3rem; color: #7dd3fc; margin: 2rem 0 1rem; border-bottom: 1px solid #334155; padding-bottom: 0.5rem; }}
        h3 {{ font-size: 1.1rem; color: #bae6fd; margin: 1.5rem 0 0.5rem; }}
        .badge {{ display: inline-block; padding: 0.25rem 0.75rem; border-radius: 9999px; font-size: 0.8rem; font-weight: 600; }}
        .badge-success {{ background: #065f46; color: #6ee7b7; }}
        .badge-error {{ background: #7f1d1d; color: #fca5a5; }}
        .badge-warn {{ background: #78350f; color: #fde68a; }}
        .badge-running {{ background: #1e3a5f; color: #93c5fd; }}
        table {{ width: 100%; border-collapse: collapse; margin: 1rem 0; background: #1e293b; border-radius: 8px; overflow: hidden; }}
        th {{ background: #334155; color: #94a3b8; text-align: left; padding: 0.75rem 1rem; font-size: 0.85rem; text-transform: uppercase; letter-spacing: 0.05em; }}
        td {{ padding: 0.6rem 1rem; border-bottom: 1px solid #334155; font-size: 0.9rem; }}
        tr:last-child td {{ border-bottom: none; }}
        .meta-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr)); gap: 1rem; }}
        .meta-card {{ background: #1e293b; border-radius: 8px; padding: 1rem 1.2rem; }}
        .meta-label {{ font-size: 0.8rem; color: #94a3b8; text-transform: uppercase; letter-spacing: 0.05em; }}
        .meta-value {{ font-size: 1.1rem; color: #f1f5f9; margin-top: 0.25rem; word-break: break-all; }}
        .section {{ margin-bottom: 2rem; }}
        .preflight-pass {{ color: #6ee7b7; }}
        .preflight-warn {{ color: #fde68a; }}
        .preflight-fail {{ color: #fca5a5; }}
        .footer {{ margin-top: 3rem; text-align: center; color: #64748b; font-size: 0.8rem; }}
    </style>
</head>
<body>
<div class="container">
    <h1>🧬 Relatório de Execução</h1>
    <p style="color: #94a3b8;">Pipeline de Variantes Germinativas v{escape(str(pipeline.get('version', '')))}</p>
    <p style="color: #64748b; font-size: 0.85rem;">Gerado em: {escape(str(data.get('report_generated_at', '')))}</p>

    <h2>📋 Resumo da Execução</h2>
    <div class="meta-grid">
        <div class="meta-card">
            <div class="meta-label">Status</div>
            <div class="meta-value">
                <span class="badge {'badge-success' if pipeline.get('status') == 'SUCCESS' else 'badge-error' if pipeline.get('status') == 'ERROR' else 'badge-running'}">
                    {escape(str(pipeline.get('status', 'N/A')))}
                </span>
            </div>
        </div>
        <div class="meta-card">
            <div class="meta-label">Duração Total</div>
            <div class="meta-value">{escape(str(pipeline.get('duration', 'N/A')))}</div>
        </div>
        <div class="meta-card">
            <div class="meta-label">Amostras</div>
            <div class="meta-value">{pipeline.get('total_samples', 0)}</div>
        </div>
        <div class="meta-card">
            <div class="meta-label">Início</div>
            <div class="meta-value">{escape(str(pipeline.get('start_time', '')))}</div>
        </div>
        <div class="meta-card">
            <div class="meta-label">Fim</div>
            <div class="meta-value">{escape(str(pipeline.get('end_time', '')))}</div>
        </div>
        <div class="meta-card">
            <div class="meta-label">Run ID</div>
            <div class="meta-value">{pipeline.get('run_id', '')}</div>
        </div>
    </div>

    <h2>🖥️ Ambiente</h2>
    <div class="meta-grid">
        <div class="meta-card">
            <div class="meta-label">Hostname</div>
            <div class="meta-value">{escape(str(env.get('hostname', '')))}</div>
        </div>
        <div class="meta-card">
            <div class="meta-label">Usuário</div>
            <div class="meta-value">{escape(str(env.get('username', '')))}</div>
        </div>
        <div class="meta-card">
            <div class="meta-label">Sistema Operacional</div>
            <div class="meta-value">{escape(str(env.get('os_info', '')))}</div>
        </div>
        <div class="meta-card">
            <div class="meta-label">CPUs</div>
            <div class="meta-value">{env.get('cpu_count', 'N/A')}</div>
        </div>
        <div class="meta-card">
            <div class="meta-label">Memória Total</div>
            <div class="meta-value">{env.get('memory_total_mb', 0)} MB</div>
        </div>
    </div>

    <h2>🔧 Ferramentas</h2>
    <table>
        <tr><th>Ferramenta</th><th>Versão</th><th>Caminho</th></tr>
"""

    for tool in tools:
        html += f"        <tr><td>{escape(str(tool.get('tool_name', '')))}</td>"
        html += f"<td>{escape(str(tool.get('tool_version', '')))}</td>"
        html += f"<td style='color:#64748b;font-size:0.8rem;'>{escape(str(tool.get('tool_path', '')))}</td></tr>\n"

    html += "    </table>\n"

    # Referência
    if ref:
        html += """
    <h2>🧬 Referência Genômica</h2>
    <div class="meta-grid">
"""
        for label, key in [("Nome", "name"), ("Versão", "version"), ("Assembly", "assembly"),
                           ("Organismo", "organism"), ("Fonte", "source"), ("SHA256", "checksum_sha256")]:
            val = ref.get(key, "")
            if val:
                html += f'        <div class="meta-card"><div class="meta-label">{label}</div><div class="meta-value">{escape(str(val))}</div></div>\n'
        html += "    </div>\n"

    # VEP
    if vep:
        html += """
    <h2>🏷️ Ensembl VEP</h2>
    <div class="meta-grid">
"""
        for label, key in [("Versão", "vep_version"), ("Cache", "cache_version"),
                           ("Assembly", "assembly"), ("FASTA", "fasta_version")]:
            val = vep.get(key, "")
            if val:
                html += f'        <div class="meta-card"><div class="meta-label">{label}</div><div class="meta-value">{escape(str(val))}</div></div>\n'
        html += "    </div>\n"

    # Preflight
    if preflight:
        html += """
    <h2>✅ Verificações Pre-flight</h2>
    <table>
        <tr><th>Verificação</th><th>Status</th><th>Mensagem</th></tr>
"""
        for pf in preflight:
            status = pf.get("status", "")
            css_class = {"PASS": "preflight-pass", "WARN": "preflight-warn", "FAIL": "preflight-fail"}.get(status, "")
            icon = {"PASS": "✅", "WARN": "⚠️", "FAIL": "❌"}.get(status, "")
            html += f'        <tr><td>{escape(str(pf.get("check_name", "")))}</td>'
            html += f'<td class="{css_class}">{icon} {escape(status)}</td>'
            html += f'<td>{escape(str(pf.get("message", "")))}</td></tr>\n'
        html += "    </table>\n"

    # Amostras
    if samples:
        html += """
    <h2>📊 Amostras Processadas</h2>
    <table>
        <tr><th>Amostra</th><th>Tipo</th><th>Status</th><th>Duração</th></tr>
"""
        for sample in samples:
            status = sample.get("status", "")
            badge = "badge-success" if status == "SUCCESS" else "badge-error" if status == "ERROR" else "badge-running"
            dur = format_duration(sample["total_duration"]) if sample.get("total_duration") else "N/A"
            html += f'        <tr><td>{escape(str(sample.get("sample_name", "")))}</td>'
            html += f'<td>{escape(str(sample.get("seq_type", "")))}</td>'
            html += f'<td><span class="badge {badge}">{escape(status)}</span></td>'
            html += f'<td>{escape(dur)}</td></tr>\n'

            # Etapas da amostra
            if sample.get("steps"):
                for step in sample["steps"]:
                    s_status = step.get("status", "")
                    s_dur = format_duration(step["duration"]) if step.get("duration") else "-"
                    html += f'        <tr><td style="padding-left:2rem;color:#94a3b8;">↳ {escape(str(step.get("step_name", "")))}</td>'
                    html += '<td></td>'
                    html += f'<td style="font-size:0.85rem;">{escape(s_status)}</td>'
                    html += f'<td style="font-size:0.85rem;">{escape(s_dur)}</td></tr>\n'

        html += "    </table>\n"

    html += f"""
    <div class="footer">
        <p>Pipeline de Variantes Germinativas v{escape(str(pipeline.get('version', '')))} — Relatório gerado automaticamente</p>
    </div>
</div>
</body>
</html>
"""

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)


# ============================================================
# CLI
# ============================================================

def main() -> None:
    """CLI para geração de relatórios."""
    parser = argparse.ArgumentParser(
        description="Gera relatório de execução do pipeline"
    )
    parser.add_argument(
        "--db", default="database/pipeline.db",
        help="Caminho do banco SQLite",
    )
    parser.add_argument(
        "--output", default="results",
        help="Diretório de saída",
    )
    parser.add_argument(
        "--run-id", type=int, default=None,
        help="ID da execução (default: última)",
    )
    args = parser.parse_args()

    if not Path(args.db).is_file():
        print(f"[ERRO] Banco não encontrado: {args.db}")
        sys.exit(1)

    db = PipelineDatabase(args.db)
    data = generate_report_data(db, args.run_id)

    json_path = Path(args.output) / "execution_report.json"
    html_path = Path(args.output) / "execution_report.html"

    write_json_report(data, str(json_path))
    print(f"[INFO] Relatório JSON: {json_path}")

    write_html_report(data, str(html_path))
    print(f"[INFO] Relatório HTML: {html_path}")

    db.close()


if __name__ == "__main__":
    main()
