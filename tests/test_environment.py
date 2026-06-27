from scripts.common.environment import capture_environment, capture_system_info, capture_tool_versions

def test_capture_system_info():
    info = capture_system_info()
    assert info.os_name != ""
    assert info.cpu_count > 0
    assert info.hostname != ""

def test_capture_tool_versions():
    # Testa se a captura roda sem quebrar (ferramentas podem ou não estar instaladas no CI)
    tools = capture_tool_versions()
    assert "python3" in tools
    assert "bash" in tools
    
def test_capture_environment():
    config = {"pipeline": {"version": "2.1.0"}}
    snapshot = capture_environment(config)
    
    assert snapshot.pipeline_version == "2.1.0"
    assert snapshot.system.os_name != ""
    
    # Verifica serialização JSON
    json_data = snapshot.to_json()
    assert "2.1.0" in json_data
