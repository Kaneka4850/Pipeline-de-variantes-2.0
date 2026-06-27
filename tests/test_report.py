from scripts.generate_report import generate_report_data, write_json_report, write_html_report

def test_report_generation(temp_db, tmp_path):
    # Simula dados no banco
    run_id = temp_db.start_run("2.1.0")
    temp_db.update_run_environment(
        run_id, "Linux", 8, 16000, '{"test": "config"}'
    )
    
    sample_id = temp_db.start_sample(run_id, "S1", "PE", "r1.fq")
    temp_db.finish_sample(sample_id, "SUCCESS")
    
    temp_db.register_reference(run_id, {
        "name": "GRCh38", "version": "110", "assembly": "GRCh38", "fasta_path": "/mock/ref.fa"
    })
    
    temp_db.register_vep_info(run_id, {
        "version": "112", "assembly": "GRCh38"
    })
    
    temp_db.register_preflight(run_id, [
        {"check_name": "TEST", "status": "PASS", "message": "OK"}
    ])
    
    temp_db.finish_run(run_id, "SUCCESS")
    
    # Gera dados
    data = generate_report_data(temp_db, run_id)
    
    assert "report_generated_at" in data
    assert data["pipeline"]["version"] == "2.1.0"
    assert data["pipeline"]["status"] == "SUCCESS"
    assert len(data["samples"]) == 1
    assert data["reference"]["name"] == "GRCh38"
    assert data["vep"]["vep_version"] == "112"
    assert len(data["preflight"]) == 1
    
    # Testa exportação JSON e HTML
    json_out = tmp_path / "report.json"
    html_out = tmp_path / "report.html"
    
    write_json_report(data, str(json_out))
    assert json_out.is_file()
    
    write_html_report(data, str(html_out))
    assert html_out.is_file()
    assert "GRCh38" in html_out.read_text(encoding="utf-8")
