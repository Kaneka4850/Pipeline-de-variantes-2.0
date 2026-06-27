from scripts.common.reference_manager import ReferenceManager

def test_reference_manager_detects_name(tmp_path):
    fasta = tmp_path / "GRCh38.p14.fa"
    fasta.write_text(">chr1\nACGT")
    
    manager = ReferenceManager(str(fasta), {})
    info = manager.collect_info()
    
    assert info.name == "GRCh38"
    assert info.version == "p14"

def test_reference_manager_missing_indices(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nACGT")
    
    manager = ReferenceManager(str(fasta), {})
    info = manager.collect_info()
    
    assert not info.integrity_ok
    assert "bwa_amb" in info.missing_indices
    assert "samtools_fai" in info.missing_indices
    
def test_reference_manager_caching(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nACGT")
    
    manager = ReferenceManager(str(fasta), {})
    info1 = manager.collect_info()
    assert info1.checksum_sha256 != ""
    
    # Verifica se o arquivo de cache foi criado
    cache = tmp_path / ".reference_meta.json"
    assert cache.is_file()
    
    # Executa novamente para garantir que usa o cache
    info2 = manager.collect_info()
    assert info1.checksum_sha256 == info2.checksum_sha256
