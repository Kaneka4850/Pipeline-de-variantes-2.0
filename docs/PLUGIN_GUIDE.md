# ============================================================
# Guia de Extensão e Escalabilidade (Plugins & Perfis)
# ============================================================
#
# O Pipeline de Variantes 2.0 foi arquitetado para ser
# estendido sem a necessidade de reescrever seu core.

## 1. Adicionando Novos Annotators ou Callers (Plugins)

A estrutura modular em `workflow/rules/` permite que você crie 
seus próprios plugins de anotação ou de chamada de variantes.

### Exemplo: Adicionando SnpEff

1. Crie o arquivo `workflow/rules/annotator_snpeff.smk`.
2. Adicione as regras para o SnpEff, consumindo o `{sample}.vcf.gz` gerado pelo HaplotypeCaller.
3. No `Snakefile`, em vez de incluir `annotator.smk`, inclua seu arquivo.
4. (Recomendado) Adicione uma chave em `config.yaml` para alternar entre VEP e SnpEff.

## 2. Executando em HPC / Cluster SLURM

O pipeline é 100% compatível com gerenciadores de recursos (HPC) nativamente através do Snakemake.

Fornecemos um perfil base em `workflow/profiles/slurm/`.
Para utilizá-lo:

```bash
snakemake --profile workflow/profiles/slurm/
```

### Configuração do SLURM (config.yaml do profile)
O arquivo `workflow/profiles/slurm/config.yaml` contém mapeamentos de 
recursos. Exemplo de customização por regra:

```yaml
default-resources:
  - slurm_partition=normal
  - mem_mb=4000
  - runtime=60

set-resources:
  bwa_mem:
    - mem_mb=16000
    - runtime=240
  haplotype_caller:
    - mem_mb=8000
```

## 3. Ambientes em Nuvem (Kubernetes / AWS Batch)

O pipeline pode ser executado em nuvem usando a integração do Snakemake
com o **Tibanna** (AWS) ou nativamente via **Kubernetes**.

1. Crie um perfil Kubernetes (`workflow/profiles/k8s/config.yaml`).
2. Adicione as chaves `kubernetes:` ao perfil.
3. Garanta que os dados estão acessíveis em um volume montado ou via S3 (usando os remotes do Snakemake).
