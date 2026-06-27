from scripts.common.preflight import PreflightRunner

def test_preflight_no_samples(sample_config):
    runner = PreflightRunner(sample_config, {})
    runner.run_all()
    
    assert runner.has_failures
    fail_messages = [r.message for r in runner.failures]
    assert any("Nenhuma amostra detectada" in m for m in fail_messages)

def test_preflight_missing_reference(sample_config):
    # O FASTA não existe no tmp_path ainda
    runner = PreflightRunner(sample_config, {"S1": {"r1": "s1_R1.fq.gz"}})
    runner.run_all()
    
    assert runner.has_failures
    fail_checks = [r.check_name for r in runner.failures]
    assert "REFERENCE_FASTA" in fail_checks
    
def test_preflight_warnings_only(sample_config, tmp_path):
    # Cria os arquivos necessários para passar nas validações críticas
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nACGT")
    
    fastq = tmp_path / "s1.fq"
    fastq.write_text("@seq1\nACGT\n+\nIIII\n")
    
    runner = PreflightRunner(sample_config, {"S1": {"r1": str(fastq)}})
    runner.run_all()
    
    # Provavelmente teremos WARN para as ferramentas não instaladas e cache do VEP
    # Mas não falhas para o FASTA ou FASTQ
    fail_checks = [r.check_name for r in runner.failures]
    assert "REFERENCE_FASTA" not in fail_checks
    
    # Os índices não estão lá, então deve ter um WARN
    warn_checks = [r.check_name for r in runner.warnings]
    assert "REFERENCE_INDICES" in warn_checks
