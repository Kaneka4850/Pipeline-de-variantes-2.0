#!/bin/bash

# ==============================================================================
# Script de Download de Atualizações do GitHub (Pull)
# ==============================================================================
# Uso: ./git_pull.sh

echo "[1/2] Buscando atualizações no repositório remoto (Github)..."
git fetch origin

echo "[2/2] Aplicando atualizações (pull)..."
git pull origin main

if [ $? -eq 0 ]; then
    echo "✅ Sucesso! Seu repositório local está sincronizado e atualizado."
else
    echo "❌ Erro ao puxar código do GitHub. Pode haver conflitos no seu código local."
fi
