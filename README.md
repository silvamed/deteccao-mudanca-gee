# Detecção de mudança de cobertura em propriedades militares da Mata Atlântica

App para o Google Earth Engine que mapeia perda de cobertura vegetal e exposição de solo
dentro de uma área delimitada pelo usuário, com imagens Sentinel-2 de 2019 a 2025.

Desenvolvido no Projeto de Fim de Curso de Engenharia Cartográfica do Instituto Militar de
Engenharia (IME), por Pedro Silva de Medeiros e Pablo dos Santos Gomes Artte.

## Arquivos

- `app_degradacao.js`: o aplicativo (versão 2.7).
- `extrair_endmembers.js`: script auxiliar que calcula os endmembers a partir de pixels puros
  da própria imagem.

## Como usar

1. Abra o [Code Editor do Earth Engine](https://code.earthengine.google.com/) com uma conta
   habilitada.
2. Crie um script novo, cole o conteúdo de `app_degradacao.js` e clique em **Run**.
3. Desenhe o polígono da área no mapa da esquerda ou informe o id de um asset seu e clique
   em **Usar asset**.
4. Clique em **Verificar bioma** e depois em **Processar**.
5. Use os botões de **Saídas** para exportar os produtos para o Google Drive (aba Tasks).

## O que o app faz

- Composição mediana de cada ano com imagens de junho a setembro (Sentinel-2 SR Harmonized,
  máscara de nuvem Cloud Score+).
- Modelo linear de mistura espectral (vegetação, solo e sombra) e NDVI.
- Classificação anual: fração de solo ≥ 0,30 e NDVI ≤ 0,45.
- Persistência (número de anos classificado) e redução contínua do NDVI entre 2019 e 2025.
- Amostra aleatória estratificada para verificação da acurácia, exportada em SHP.

## Limitações

Os parâmetros foram definidos e verificados apenas na Mata Atlântica. Fora desse bioma o app
pede confirmação antes de processar, e o resultado não tem verificação. A classificação indica
baixa cobertura vegetal e substrato aparente; não confirma degradação do solo, que depende de
diagnóstico de campo.

## Como citar

MEDEIROS, P. S. de; ARTTE, P. dos S. G. **Detecção de mudança de cobertura em propriedades
militares da Mata Atlântica**: aplicativo para o Google Earth Engine. Versão 2.7. 2026.
Disponível em: <ENDEREÇO DO REPOSITÓRIO>.
