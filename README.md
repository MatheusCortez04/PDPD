# PDPD — Predição e Priorização de Drogas para Transtornos Psiquiátricos via Medicina de Rede

> **Nota:** Este documento é uma versão atualizada do README original, refletindo o estado atual do pipeline implementado em `src/PDPD.R` e nos módulos de `src/Utils/`. O `README.md` original foi mantido sem alterações.

Este projeto de Iniciação Científica utiliza princípios de **Medicina de Rede** (Network Medicine) sobre a rede de interação proteína-proteína (PPI) humana para investigar o reposicionamento de fármacos (*drug repurposing*) no contexto de dois transtornos psiquiátricos: o **Transtorno Depressivo Maior (TDM/MDD)** e o **Transtorno Bipolar (TB/BD)**.

A ideia central é medir a **proximidade topológica** entre os alvos proteicos de fármacos conhecidos e os módulos de doença (genes associados a MDD/BD) dentro da rede PPI, gerando um ranking de candidatos a reposicionamento, validando esse ranking contra uma base de evidências reais (RepoDB) e, por fim, interpretando biologicamente os resultados via análises de enriquecimento (GSEA/ORA).

## 1. Visão Geral do Pipeline

O ponto de entrada da aplicação é `src/PDPD.R`, que carrega os módulos abaixo e expõe um menu interativo em console (via `main()`, disparado automaticamente pelo `.Rprofile`):

```
--- Main Menu ---

 [1] Kernel Menu
 [2] Drug Menu
 [3] Evaluation Menu
 [4] Graph Menu
 [5] Enrichment Menu
 [Q] Exit
```

O fluxo de trabalho recomendado, em ordem, é:

1. **Kernel Menu** → gera as matrizes de similaridade (kernels) a partir da rede PPI.
2. **Drug Menu** → constrói a matriz droga‑proteína e calcula o *score* de cada droga por doença/kernel, consolidando em um ranking médio.
3. **Evaluation Menu** → valida o ranking contra o gold-standard (RepoDB), gera curvas ROC, Recall@K, e enriquece o top-N de candidatos com metadados clínicos via OpenTargets.
4. **Graph Menu** → visualiza o módulo de doença (subgrafo PPI) e sua maior componente conectada (LCC).
5. **Enrichment Menu** → realiza análises de enriquecimento (GSEA sobre o ranking de drogas e ORA sobre os genes de doença).

## 2. Metodologia

### 2.1 Kernels de rede

Os kernels são implementados através do pacote Bioconductor `diffuStats` e quantificam a proximidade/similaridade entre nós (proteínas) da rede PPI:

* **Diffusion Kernel** — difusão de calor sobre o grafo.
* **P-Step Kernel** (*Random Walk*) — passeio aleatório de `p` passos (padrão: 5).
* **Regularised Laplacian Kernel** — laplaciano regularizado do grafo.
* **Commute Time Kernel** — tempo de deslocamento (comute) entre nós.
* **Inverse Cosine Kernel** — similaridade de cosseno inverso.

### 2.2 Score droga-doença

Para cada kernel `K`, o score de uma droga é calculado como:

```
score = matriz_droga_proteina %*% K %*% vetor_binario_doenca
```

onde `vetor_binario_doenca` indica quais proteínas da rede pertencem ao módulo da doença (score de associação gene-doença ≥ 0.6, filtro configurável em `get_score_disease_gene_association()`).

Os scores de cada um dos 5 kernels são então convertidos em *rank* e combinados em um **rank médio** (`average_rank`) por droga, por doença (MDD e BD).

### 2.3 Validação

O ranking médio é validado contra o **RepoDB** (base de indicações droga-doença conhecidas/aprovadas), rotulando cada droga como validada ou não. A partir disso são geradas:

* **Curva ROC** (via `reportROC`), usando `-average_rank` como preditor.
* **Recall@K**, mostrando a fração de drogas validadas recuperadas nos top-K do ranking.

### 2.4 Enriquecimento de candidatos (OpenTargets)

Para o top-N de drogas do ranking, o pipeline consulta a API GraphQL pública do **OpenTargets** (`https://api.platform.opentargets.org/api/v4/graphql`) para:

* Resolver o ID DrugBank para o ID ChEMBL/nome da droga.
* Buscar candidatos clínicos e estágio clínico máximo (`PRECLINICAL` → `APPROVED`) para a doença.
* Recuperar relatórios clínicos (fonte, URL, resumo de evidência) associados.

### 2.5 Análises de enriquecimento biológico

* **DSEA (Drug Set Enrichment Analysis)** — um GSEA (`clusterProfiler::GSEA`) é executado sobre o ranking de drogas, usando como *TERM2GENE* o mapeamento pathway→droga construído a partir das bases Enrichr (`KEGG_2026`, `GO_Biological_Process_2026`, `Reactome_Pathways_2024`), interpretando quais vias biológicas estão associadas às drogas mais bem ranqueadas.
* **ORA (Over-Representation Analysis)** — análise de sobre-representação de vias KEGG para o conjunto de genes associados a cada doença (dotplot + tabela de sumário).

## 3. Instalação e Ambiente

Para garantir a reprodutibilidade, o projeto é totalmente encapsulado em um ambiente Docker (imagem base `bioconductor/bioconductor_docker:devel-r-4.6.0`).

### Pré-requisitos

* **Docker:** [Docker Docs](https://docs.docker.com/engine/install/)
* **Docker Compose:** [Docker Docs](https://docs.docker.com/compose/install/)
* **Git**

### Passos de Instalação

```bash
git clone git@github.com:MatheusCortez04/reposicionamento-farmacos-tdm-tb.gi
cd PDPD
docker-compose build --no-cache
```

O `Dockerfile` instala, via `pak`, as dependências R do projeto: `tidyverse`, `here`, `reportROC`, `httr2`, `diffuStats`, `igraphdata`, `igraph`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `ReactomePA`, `AnnotationDbi`, `enrichR`.

## 4. Como Executar a Análise

```bash
docker compose run --rm --service-ports rstudio
```

Ao iniciar a sessão R, o `.Rprofile` executa `source("src/PDPD.R")` e chama `main()`, abrindo o menu principal descrito na Seção 1.

## 5. Estrutura do Projeto

```
src/
├── PDPD.R                       # Entrypoint: carrega módulos e roda o menu principal
├── Utils/
│   ├── utils.R                  # Funções compartilhadas (I/O, genes de doença, filtros)
│   ├── kernels.R                # Geração dos 5 kernels de rede (Kernel Menu)
│   ├── drug.R                   # Matriz droga-proteína e score droga-doença (Drug Menu)
│   ├── evaluation.R             # ROC, Recall@K, validação RepoDB, top-N drugs (Evaluation Menu)
│   ├── graph_metrics.R          # Subgrafo de doença + LCC (Graph Menu)
│   ├── enrichment.R             # GSEA (drogas) e ORA (genes de doença) (Enrichment Menu)
│   └── open_targets_requests.R  # Cliente GraphQL para a API do OpenTargets
├── Data/
│   ├── PPI_gysi.csv             # Rede PPI de entrada
│   ├── Disease/                 # Vetores binários gene-doença (MDD/BD)
│   ├── Drug/                    # Matriz droga-alvo, scores por kernel/doença e rank médio
│   ├── Kernels/                 # Kernels calculados (.csv e .Rdata)
│   └── REPODB/                  # Gold-standard de indicações droga-doença (MDD/BIPOLAR)
├── Evaluation/                  # Predições, ROC e Recall@K por doença
├── Enrichment/                  # Resultados de GSEA/ORA por doença
└── Graph/                       # Plots do módulo de doença (LCC) por doença
```

## 6. Menus Detalhados

### 6.1 Kernel Menu (`src/Utils/kernels.R`)

Lê `src/Data/PPI_gysi.csv`, monta o grafo PPI não-direcionado (`generate_graph_from_dataframe`) e apresenta:

```
[1] Diffusion Kernel
[2] P-Step Kernel
[3] Regularised Laplacian Kernel
[4] Commute Time Kernel
[5] Inverse Cosine Kernel
[B] Back
```

Cada opção solicita interativamente parâmetros (ex.: número de passos do P-Step, se deve salvar `.RData`) e persiste o resultado em `src/Data/Kernels/` (`.csv` + `.Rdata`).

### 6.2 Drug Menu (`src/Utils/drug.R`)

```
[1] Generate Drug Protein Matrix
[2] Generate Drug Score
[B] Back
```

* **[1]** Lê `src/Data/Drug/drug_targets_DrugBank_Gysi.csv`, filtra apenas alvos presentes na rede PPI e monta a matriz binária Droga × Proteína, salva em `src/Data/Drug/`.
* **[2]** Para cada um dos 5 kernels e cada doença (MDD/BD), calcula o score de cada droga (`process_scores_for_disease`) e, em seguida, consolida os 5 rankings em um **rank médio** por droga (`generate_drug_rank`), salvo em `src/Data/Drug/Score/<DOENÇA>/average_kernel_rank.csv`.

### 6.3 Evaluation Menu (`src/Utils/evaluation.R`)

```
[1] Generate MDD ROC
[2] Generate MDD Recall@K
[3] Generate BD ROC
[4] Generate BD Recall@K
[5] Generate Top N drug rank to MDD
[6] Generate Top N drug rank to BD
[B] Back
```

* **[1]/[3]** Cruzam o rank médio com o RepoDB (gold-standard) e geram a curva ROC (`reportROC`) em `src/Evaluation/<DOENÇA>/ROC/`.
* **[2]/[4]** Geram o gráfico de Recall@K em `src/Evaluation/<DOENÇA>/Recall/`.
* **[5]/[6]** Enriquecem o ranking completo com metadados via OpenTargets (ChEMBL ID, estágio clínico, candidatos, relatórios clínicos) e salvam `top_drugs_<DOENÇA>_score_filter_<X>_.csv`.

### 6.4 Graph Menu (`src/Utils/graph_metrics.R`)

```
[1] Plot MDD Disease Module with LCC
[2] Plot BD Disease Module with LCC
[B] Back
```

Extrai o subgrafo induzido pelos genes de cada doença dentro da rede PPI, calcula a **Maior Componente Conectada (LCC)** e plota o módulo destacando genes comuns entre MDD/BD (azul) vs. genes específicos da doença (laranja), salvando o PDF em `src/Graph/<DOENÇA>/`.

### 6.5 Enrichment Menu (`src/Utils/enrichment.R`)

```
[1] Drug Set Enrichment Analysis MDD
[2] Drug Set Enrichment Analysis BD
[3] Over Representation Analysis MDD
[4] Over Representation Analysis BD
[B] Back
```

* **[1]/[2]** Rodam o pipeline de **DSEA/GSEA**: mapeiam alvos das drogas para símbolos gênicos, anotam vias via Enrichr (KEGG, GO-BP, Reactome), constroem o mapeamento via→droga e executam `clusterProfiler::GSEA` sobre o ranking de drogas, salvando tabelas (`*_gsea_<db>_.csv`) e gráficos (`gseaplot2`) em `src/Enrichment/Drugs/<DOENÇA>/<BASE>/`.
* **[3]/[4]** Rodam **ORA** sobre os genes associados à doença (base KEGG), gerando um dotplot (Gene Ratio × -log10 p-valor ajustado, tamanho = contagem de genes) e uma tabela de sumário em `src/Enrichment/<DOENÇA>/`.

## 7. Fontes de Dados

| Dado | Arquivo | Origem |
|---|---|---|
| Rede PPI | `src/Data/PPI_gysi.csv` | Rede de interação proteína-proteína humana |
| Alvos de fármacos | `src/Data/Drug/drug_targets_DrugBank_Gysi.csv` | DrugBank |
| Genes de doença | `src/Data/Disease/disease_genes.csv` | Associação gene-doença (score ≥ 0.6) |
| Gold-standard de validação | `src/Data/REPODB/MDD_REPODB.tsv`, `BIPOLAR_REPODB.tsv` | RepoDB |
| Metadados clínicos de drogas | — (requisição em tempo real) | API GraphQL do [OpenTargets](https://platform.opentargets.org/) |
| Vias biológicas | — (requisição em tempo real) | Enrichr (`KEGG_2026`, `GO_Biological_Process_2026`, `Reactome_Pathways_2024`) |

## 8. Referências

- [An Introduction to R (Manual) — seção The R Profile](https://cran.r-project.org/doc/manuals/r-release/R-intro.html#The-R-profile)
- [RStudio Documentation: Managing R Sessions and Profiles](https://docs.posit.co/ide/user/ide/guide/environments/r/managing-r.html#rprofile)
- [diffuStats (Bioconductor)](https://bioconductor.org/packages/release/bioc/html/diffuStats.html)
- [clusterProfiler (Bioconductor)](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
- [OpenTargets Platform API](https://platform.opentargets.org/)
- [RepoDB](http://apps.chiragjpgroup.org/repoDB/)
