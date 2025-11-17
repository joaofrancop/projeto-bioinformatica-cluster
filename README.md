# Projeto: Análise de Cluster de Proteínas (SCOPe)

Este projeto realiza uma análise de clusterização não-supervisionada em sequências de proteínas da base ASTRAL SCOPe.

O foco principal é a **metodologia não-subjetiva** para validar e escolher a melhor métrica de clusterização interna.

---

## 🎬 Vídeo de Apresentação

**O vídeo com a explicação completa da metodologia e dos resultados pode ser visto aqui:**

https://drive.google.com/file/d/13URtJYUO9IFK03tPHLqPLAh9sYTWlvti/view?usp=sharing

---

## 🚀 Como Executar

1.  **Baixar os Dados:**
    Os dados (arquivos `.fa`) podem ser baixados do site oficial: [ASTRAL SCOPe](https://scop.berkeley.edu/astral/ver=2.08)

2.  **Instalar Dependências:**
    ```bash
    pip install -r requirements.txt
    ```

3.  **Executar o Script:**
    (Lembre-se de atualizar os caminhos dos arquivos `.fa` dentro do script)
    ```bash
    python analise_cluster.py
    ```

---

## 📄 Resultado da Execução

<details>
<summary>Clique para ver o log completo da execução</summary>

```text
================================================================================
INICIANDO ANÁLISE PARA O ARQUIVO: astral-scopedom-seqres-gd-sel-gs-bib-95-2.08.fa
================================================================================

--- ETAPA A: CARREGAMENTO DE DADOS ---
-> Dados carregados: 35373 sequências.
-> Classes SCOPe (Gabarito): 7 classes (a, b, c, d, e, f, g).

--- ETAPA B: EXTRAÇÃO DE ATRIBUTOS (K-MER BINÁRIO) ---
  -> Total de 400 atributos (k-mers) gerados.
-> Matriz Binária gerada: (35373, 400)

--- ETAPA C: REDUÇÃO DE DIMENSIONALIDADE (PCA) ---
-> PCA finalizado. Matriz de atributos X: (35373, 300)
-> Variância explicada com 300 componentes: 85.71%
-> Matrizes intermediárias removidas da memória.

--- ETAPA D: TESTE DE ALGORITMOS (BASELINE) ---
  -> Rodando: KMeans (K=7)...
  -> Rodando: MiniBatchKMeans (K=7)...
  -> Rodando: Birch (K=7)...
  -> Rodando: DBSCAN (Base)...

--- RESULTADOS DA AVALIAÇÃO BASELINE ---
| Algoritmo             |    Nº Clusters |    Silhouette |    Davies-Bouldin |    V-Measure (Externo) |    Tempo (s) |
|:----------------------|---------------:|--------------:|------------------:|-----------------------:|-------------:|
| MiniBatchKMeans (K=7) |              7 |   -0.00317315 |           5.46644 |               0.22427  |     0.143215 |
| Birch (K=7)           |              7 |   -0.00589482 |           5.35167 |               0.22426  |    75.4981   |
| KMeans (K=7)          |              7 |    0.0117725  |           5.36169 |               0.188603 |     1.4318   |
| DBSCAN (Base)         |              0 |   -1          |         999       |               0        |     2.84153  |

================================================================================
ETAPA F: CORRELAÇÃO DE MÉTRICAS (A ESCOLHA NÃO-SUBJETIVA)
================================================================================

--- Matriz de Correlação (Internas vs Externa) ---
|                     |    Silhouette |    Davies-Bouldin |    V-Measure (Externo) |
|:--------------------|--------------:|------------------:|-----------------------:|
| Silhouette          |      1        |          -0.296571 |              -0.989679 |
| Davies-Bouldin      |     -0.296571 |           1        |               0.430365 |
| V-Measure (Externo) |     -0.989679 |           0.430365 |               1        |

--- Correlação com a Métrica Externa (V-Measure) ---
|                |    V-Measure (Externo) |
|:---------------|-----------------------:|
| Silhouette     |              -0.989679 |
| Davies-Bouldin |               0.430365 |

--- CONCLUSÃO DA METODOLOGIA ---
A métrica interna (não-supervisionada) mais confiável é: Silhouette
Justificativa: Ela possui a maior correlação absoluta (-0.9897) com a métrica externa (V-Measure).
A partir de agora, usaremos 'Silhouette' como nossa 'bússola' para otimizar os parâmetros.

--- ETAPA G: OTIMIZAÇÃO (ex: K-Means) USANDO A MÉTRICA VALIDADA ---
  -> Tunando K-Means com K variando de 2 a 12...
  -> Objetivo: Encontrar o 'k' com o melhor score de 'Silhouette'

--- MELHOR CONFIGURAÇÃO ENCONTRADA (Baseado em Silhouette) ---
|                     | 0                   |
|:--------------------|:--------------------|
| Algoritmo           | KMeans (k=2)        |
| Nº Clusters         | 2                   |
| Silhouette          | 0.09340622644301422 |
| Davies-Bouldin      | 3.496695085424174   |
| V-Measure (Externo) | 0.10162658494584167 |
| Tempo (s)           | 0.10062527656555176 |

Sugestão: O melhor K foi 2, que atingiu o score de 0.0934 (usando Silhouette).

================================================================================
INICIANDO ANÁLISE PARA O ARQUIVO: astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa
================================================================================

--- ETAPA A: CARREGAMENTO DE DADOS ---
-> Dados carregados: 15129 sequências.
-> Classes SCOPe (Gabarito): 7 classes (a, b, c, d, e, f, g).

--- ETAPA B: EXTRAÇÃO DE ATRIBUTOS (K-MER BINÁRIO) ---
  -> Total de 400 atributos (k-mers) gerados.
-> Matriz Binária gerada: (15129, 400)

--- ETAPA C: REDUÇÃO DE DIMENSIONALIDADE (PCA) ---
-> PCA finalizado. Matriz de atributos X: (15129, 300)
-> Variância explicada com 300 componentes: 85.43%
-> Matrizes intermediárias removidas da memória.

--- ETAPA D: TESTE DE ALGORITMOS (BASELINE) ---
  -> Rodando: KMeans (K=7)...
  -> Rodando: MiniBatchKMeans (K=7)...
  -> Rodando: Birch (K=7)...
  -> Rodando: DBSCAN (Base)...

--- RESULTADOS DA AVALIAÇÃO BASELINE ---
| Algoritmo             |    Nº Clusters |    Silhouette |    Davies-Bouldin |    V-Measure (Externo) |    Tempo (s) |
|:----------------------|---------------:|--------------:|------------------:|-----------------------:|-------------:|
| KMeans (K=7)          |              7 |   -0.0145484  |           6.39602 |               0.150725 |     0.169633 |
| Birch (K=7)           |              7 |   -0.0195497  |           7.78718 |               0.144299 |    14.3475   |
| MiniBatchKMeans (K=7) |              7 |   -0.00583201 |           6.87747 |               0.12874  |     0.0996928|
| DBSCAN (Base)         |              0 |   -1          |         999       |               0        |     0.555377 |

================================================================================
ETAPA F: CORRELAÇÃO DE MÉTRICAS (A ESCOLHA NÃO-SUBJETIVA)
================================================================================

--- Matriz de Correlação (Internas vs Externa) ---
|                     |    Silhouette |    Davies-Bouldin |    V-Measure (Externo) |
|:--------------------|--------------:|------------------:|-----------------------:|
| Silhouette          |      1        |          -0.51789  |              -0.792011 |
| Davies-Bouldin      |     -0.51789  |           1        |              -0.112083 |
| V-Measure (Externo) |     -0.792011 |          -0.112083 |               1        |

--- Correlação com a Métrica Externa (V-Measure) ---
|                |    V-Measure (Externo) |
|:---------------|-----------------------:|
| Silhouette     |              -0.792011 |
| Davies-Bouldin |              -0.112083 |

--- CONCLUSÃO DA METODOLOGIA ---
A métrica interna (não-supervisionada) mais confiável é: Silhouette
Justificativa: Ela possui a maior correlação absoluta (-0.7920) com a métrica externa (V-Measure).
A partir de agora, usaremos 'Silhouette' como nossa 'bússola' para otimizar os parâmetros.

--- ETAPA G: OTIMIZAÇÃO (ex: K-Means) USANDO A MÉTRICA VALIDADA ---
  -> Tunando K-Means com K variando de 2 a 12...
  -> Objetivo: Encontrar o 'k' com o melhor score de 'Silhouette'

--- MELHOR CONFIGURAÇÃO ENCONTRADA (Baseado em Silhouette) ---
|                     | 0                   |
|:--------------------|:--------------------|
| Algoritmo           | KMeans (k=2)        |
| Nº Clusters         | 2                   |
| Silhouette          | 0.09725861217534103 |
| Davies-Bouldin      | 3.446004392677131   |
| V-Measure (Externo) | 0.10727089398427898 |
| Tempo (s)           | 0.03355145454406738 |

Sugestão: O melhor K foi 2, que atingiu o score de 0.0973 (usando Silhouette).
