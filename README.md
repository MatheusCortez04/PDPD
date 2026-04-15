# Análise de Reposicionamento de Fármacos para Transtornos Psiquiátricos


Este projeto utiliza redes de interação proteína-proteína (PPI) e a proximidade de alvos de medicamentos na rede para identificar potenciais novos usos para medicamentos existentes, focando especificamente no Transtorno Depressivo Maior (TDM) e no Transtorno Bipolar (TB).

## 1. Visão Geral

O objetivo é analisar a relação entre os alvos de medicamentos existentes e os genes associados a doenças específicas (TDM e TB) dentro da rede de interação proteína-proteína humana. Isso nos permite identificar medicamentos que potencialmente perturbam as redes moleculares associadas a essas doenças, fornecendo uma base para futuras investigações sobre seu potencial terapêutico ou de reposicionamento.

## 2. Metodologia

A análise central baseia-se no cálculo de proximidade entre conjuntos de proteínas na rede. Utilizamos diferentes métricas de "kernel" para quantificar a relação entre os nós (proteínas) na rede PPI.

Os kernels disponíveis são implementados através do pacote Bioconductor `diffuStats`:
* Difusão (Diffusion Kernel)
* Passeio Aleatório (P-Step Kernel)
* Laplaciano Regularizado (Regularised Laplacian Kernel)
* Tempo de Deslocamento (Commute Time Kernel)
* Cosseno Inverso (Inverse Cosine)

A pipeline de análise consiste em:
1.  **Construção do Kernel:** Gerar uma matriz de similaridade (kernel) a partir da rede PPI.
2.  **Mapeamento Droga-Alvo:** Gerar uma matriz binária de fármacos e seus alvos proteicos.
3.  **Análise de Proximidade:** Utilizar os dados gerados para analisar as relações entre os alvos dos fármacos e os genes associados às doenças.

## 3. Instalação e Ambiente

Para garantir a reprodutibilidade, o projeto é totalmente encapsulado em um ambiente Docker.

### Pré-requisitos

* **Docker:** Siga as instruções oficiais em [Docker Docs](https://docs.docker.com/engine/install/) para instalar o Docker Desktop.
* **Docker Compose:** Siga as instruções em [Docker Docs](https://docs.docker.com/compose/install/) para instalar o Docker Compose.
* **Git:** Necessário para clonar o repositório.

### Passos de Instalação

1.  **Clone o Repositório**
    Clone este repositório para a sua máquina local:
    ```bash
    git clone git@github.com:MatheusCortez04/PDPD.git
    cd PDPD
    ```

2.  **Construa a Imagem Docker**
       Este comando irá baixar a imagem base, instalar todas as dependências do sistema e pacotes R definidos no `Dockerfile`.
    ```bash
    docker-compose build --no-cache
    ```

## 4. Como Executar a Análise

1.  **Inicie o Contêiner**
    Use o Docker Compose para iniciar o serviço do RStudio. O script `.Rprofile` no diretório do projeto será executado automaticamente, iniciando o menu da aplicação no console.
    ```bash
    docker compose run --rm --service-ports rstudio
    ```

2.  **Siga o Menu Interativo**
    Após iniciar, o console do R apresentará o menu principal. O fluxo de trabalho recomendado é executar o Menu 1 primeiro para gerar os kernels antes de usar o Menu 2.

    ```
    [1] Kernel Menu
    [2] Drug-Disease Menu
    [Q] Exit
    ```

### 4.1. Menu 1: Kernel Menu

Esta opção permite calcular e salvar as matrizes de similaridade (kernels) baseadas na rede PPI.

* **Funcionalidade:** Calcula os kernels listados na seção de Metodologia.
* **Requisito:** Requer o arquivo `PPI_gysi.csv` localizado em `src/Data/`.
* **Saída:** Salva os kernels calculados em `src/Data/Kernels/` nos formatos `.csv` e `.Rdata`.

### 4.2. Menu 2: Drug-Disease Menu

Esta opção utiliza os kernels e os dados de fármacos para gerar a matriz de relação Droga-Proteína.

* **Funcionalidade:** Gera a matriz Droga (Linha) x Proteína (Coluna).
* **Requisito:** Requer o arquivo `drug_targets_DrugBank_Gysi.csv` localizado em `src/Data/Drug/`.
* **Saída:** Salva a matriz final como `drug_gene` (em `.csv` e `.Rdata`) dentro da pasta `src/Data/Drug/`.


### Detalhes dos Scripts

### `.Rprofile`



O arquivo `.Rprofile` é um script que é automaticamente executado toda vez que uma sessão R é iniciada no diretório do projeto. Ele tem como função preparar o ambiente para executar sua análise, criando pastas necessárias e carregando os scripts principais.


<!-- docker compose run --rm --service-ports rstudio  -->


#### Funcionalidades do `.Rprofile`



- **Importação do script principal:**  

  Executa o arquivo `src/PDPD.R` via `source()`, que contém as funções e o fluxo principal da aplicação.

- **Configuração de diretórios para dados:**  

  Define variáveis de caminho para as pastas onde os kernels e dados das drogas serão armazenados.

  Verifica se essas pastas existem e, caso não existam, cria os diretórios correspondentes (`src/Data/Kernels`, `src/Data/Drug`, e respectivas subpastas `RData`).


  - **Inicialização da aplicação:**  
  Após preparar o ambiente, chama a função `main()` (definida em `src/PDPD.R`), que dispara a execução da análise, iniciando o menu e fluxos interativos.


### Referências:

- [An Introduction to R (Manual) - seção The R Profile](https://cran.r-project.org/doc/manuals/r-release/R-intro.html#The-R-profile)  

- [RStudio Documentation: Managing R Sessions and Profiles](https://docs.posit.co/ide/user/ide/guide/environments/r/managing-r.html#rprofile)  

- [Personalizando a inicialização no R](https://kevinushey.github.io/blog/2015/02/02/rprofile-essentials/](https://www.datacamp.com/pt/doc/r/customizing))



---
### `src/PDPD.R`

Este é o script principal e o ponto de entrada (entrypoint) de toda a aplicação. Ele é o primeiro script a ser executado (via `source()` no `.Rprofile`) e é responsável por carregar as dependências e iniciar a interface interativa do usuário.

#### Funcionalidades do `src/PDPD.R`

* **Carregamento de Bibliotecas:** Importa as bibliotecas R essenciais para a operação, como `here` (para gerenciamento de caminhos), `tidyverse` (para manipulação de dados) e `utils` (funções básicas do R).
* **Importação de Módulos:** Carrega os scripts auxiliares (módulos) que contêm a lógica de negócio específica, usando `source()` e `here()` para garantir que os caminhos sejam encontrados corretamente:
    * `src/Utils/graph.R`
    * `src/Utils/kernels.R`
    * `src/Utils/drug.R`
* **Definição da Função `main()`:** Define a função `main()` como o "coração" da aplicação. Esta é a única função que o `.Rprofile` precisa chamar para iniciar todo o fluxo.
  
* **Loop de Aplicação (Menu Principal):** Implementa um loop `while(TRUE)` que apresenta continuamente o menu principal ao usuário, garantindo que o programa não termine até que o usuário decida sair.
  
* **Interface do Usuário (UI):** Utiliza `cat()` para exibir as opções de menu (`[1] Kernel Menu`, `[2] Drug-Disease Menu`) e `readline()` para capturar a entrada do usuário.
  
* **Delegação de Tarefas:** Atua como um "roteador". Ele não realiza os cálculos em si, mas sim delega as tarefas para as funções de submenu (`generate_kernel_menu()` e `scoring_drug_disease_menu()`), que foram carregadas a partir dos scripts em `src/Utils/`.
* **Saída Limpa:** Ao receber a opção 'Q' (e normalizá-la com `toupper()` e `trimws()`), o script sai do loop e chama `q(save="no")`, encerrando a sessão R de forma limpa, sem salvar o *workspace*.

####  Essa abordagem centraliza o fluxo da aplicação em um único script, que atua como um orquestrador, chamando os módulos especializados (em `src/Utils/`) para executar o trabalho pesado.
---
### `src/Utils/kernels.R` 
### generate_kernel_menu e Funções de Geração de Kernels

Este módulo é responsável por gerar diferentes matrizes Kernel baseadas na rede de interação proteína-proteína (PPI), que serão usadas para análises de similaridade e proximidade entre proteínas no contexto de doenças psiquiátricas.

#### `generate_kernel_menu()`

- Lê o arquivo CSV `src/Data/PPI_gysi.csv` contendo as interações proteína-proteína.
- Seleciona apenas as colunas relevantes que indicam os pares de proteínas (`proteinA_entrezid`, `proteinB_entrezid`).
- Gera um grafo não-direcionado com essas proteínas usando a função `generate_graph_from_dataframe`.
- Apresenta um menu interativo no console para o usuário escolher o tipo de kernel a ser gerado:
  - [1] Diffusion Kernel
  - [2] P-Step Kernel
  - [3] Regularised Laplacian Kernel
  - [4] Commute Time Kernel
  - [5] Inverse Cosine Kernel
  - [B] Para voltar ao menu anterior
- Valida a entrada do usuário e chama a função específica do kernel selecionado usando o mapa `kernel_function_mapper`.

#### `kernel_function_mapper`

- É uma lista que mapeia cada opção de kernel para uma função anônima que:
  - Solicita as opções específicas necessárias como o número de passos no P-Step Kernel ou se salvar o resultado como um arquivo `.RData`.
  - Chama a função correspondente de geração de kernel, passando o grafo e esses parâmetros.

#### Função `generate_graph_from_dataframe(dataframe)`

- Cria um grafo não-direcionado a partir do dataframe de pares proteína-proteína usando `igraph::graph_from_data_frame`.
- Retorna o grafo criado para uso posterior nos cálculos de kernels.

#### Funções de geração de kernel

Cada função calcula um tipo específico de kernel baseado no grafo de proteínas e salva os resultados em CSV e opcionalmente como objeto `.RData`.
Ao final, a função retorna  o kernel calculado para uso programático.

- `generate_difussion_kernel(graph, normalized=TRUE, save_rdata=TRUE)`
- `generate_pstep_kernel(graph, step=5, save_rdata=TRUE)`
- `generate_regularised_laplacian_kernel(graph, normalized=TRUE, save_rdata=TRUE)`
- `generate_commute_time_kernel(graph, normalized=TRUE, save_rdata=TRUE)`
- `generate_inverse_cosine_kernel(graph, save_rdata=TRUE)`

**Parâmetros comuns:**

- `graph`: objeto grafo criado a partir da rede PPI.
- `save_rdata`: lógico, indica se o resultado deve ser salvo como arquivo `.RData`.

**Armazenamento:**

- Os kernels são salvos em `src/Data/Kernels` como arquivos `.csv` para visualização e `.Rdata` para uso em R.
- Mensagens são exibidas no console para informar o progresso e o local dos arquivos gerados.

---

