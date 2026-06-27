#!/bin/bash

# ==============================================================================
# Script de Atualização para o GitHub (Push)
# ==============================================================================
# Uso: ./git_push.sh "Sua mensagem de commit aqui"

if [ -z "$1" ]; then
    echo "Erro: Forneça uma mensagem de commit."
    echo "Uso correto: ./git_push.sh \"Mensagem explicativa das mudanças\""
    exit 1
fi

COMMIT_MESSAGE="$1"

echo "[1/3] Adicionando todos os arquivos modificados..."
git add .

echo "[2/3] Criando commit: '${COMMIT_MESSAGE}'"
git commit -m "${COMMIT_MESSAGE}"

echo "[3/3] Enviando para o repositório remoto (Github)..."
git push origin main

if [ $? -eq 0 ]; then
    echo "✅ Sucesso! Seu código foi atualizado no GitHub."
else
    echo "❌ Erro ao enviar para o GitHub. Verifique se há conflitos, sua conexão ou permissões."
fi
