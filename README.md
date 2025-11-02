# Análise de Reposicionamento de Fármacos para Transtornos Psiquiátricos

Este projeto utiliza redes de interação proteína-proteína (PPI) e a proximidade de alvos de medicamentos na
rede para identificar potenciais novos usos para medicamentos existentes, focando especificamente no Transtorno Depressivo Maior (TDM) e no Transtorno Bipolar (TB).

## Visão Geral

O objetivo é analisar a relação entre os alvos de medicamentos existentes e os genes associados a doenças específicas (TDM e TB) dentro da rede de interação proteína-proteína humana. Isso nos permite identificar medicamentos que potencialmente perturbam as redes moleculares associadas a essas doenças, fornecendo uma base para futuras investigações sobre seu potencial terapêutico ou de reposicionamento.
## Pré-requisitos

Para executar este projeto, você precisará ter o Docker e o Docker Compose instalados em seu sistema. Eles são usados para criar um ambiente consistente e isolado que contém todas as dependências necessárias.

* **Docker:** Siga as instruções oficiais em [Docker Docs](https://docs.docker.com/engine/install/) para instalar o Docker Desktop no seu sistema operacional (Windows, macOS ou Linux).
* **Docker Compose:** Siga as instruções em [Docker Docs](https://docs.docker.com/compose/install/) para instalar o Docker Compose.

## Como Executar o Projeto

1.  **Clone o Repositório**
    Clone este repositório para a sua máquina local usando o seguinte comando:
    ```
    git clone git@github.com:MatheusCortez04/PDPD.git
    ```

2.  **Execute o Docker Compose**
    Uma vez que você tenha o Docker instalado e em execução, navegue até o diretório raiz do projeto (onde o arquivo `docker-compose.yml` está localizado) e execute o seguinte comando:

    ```
    docker-compose build --no-cache
    ```.

3.  **Acesse o Ambiente R**
    Após a conclusão do processo de build, para inicializar a aplicacão deverá utilizar o comando:
    ```
     docker compose run --rm --service-ports rstudio       
    ```

4.  **Execute a Análise**
    Uma vez dentro do shell do contêiner, você podera utilizar os menus para gerar as analises desejadas.



5.  **Realizando  a Análise**
Após iniciar o programa, você será apresentado ao menu principal, que oferece duas opções principais e a saída do programa:

Kernel Menu: Esta opção permite calcular e salvar diferentes tipos de matrizes de similaridade (kernels) baseadas nas redes de proteínas. Os resultados são salvos em formatos .csv e .Rdata.

Drug-Disease Menu: Esta opção permite utilizar os kernels gerados para analisar as relações entre medicamentos e doenças.

Quit (Q): Sair do programa.

 [1] Kernel Menu
 [2] Drug-Disease Menu
 [Q] Exit


5.1  **Kernel Menu**

Este subMenu tem como objetivo gerar o calculo dos kernels através das funções da biblioteca 'diffuStats'[https://www.bioconductor.org/packages/release/bioc/html/diffuStats.html]

Os kernels disponiveis são:

* Difussão (Diffusion Kernel)
* Passeio Aleatorio (P-Step Kernel)
* Laplaciano Regularizado (Regularised Laplacian Kernel)
* Tempo de Deslocamento (Commute Time Kernel)
* Cosseno Inverso(Inverse Cosine)

Para que os kernels sejam calculados corretamente é necessário que o arquivo com as ligacoes entre as proteinas com o seguinte nome 'PPI_gysi.csv' e tal arquivo deve estar dentro da pasta 'Data'[src/Data]



5.2  **Kernel Menu**
    Este subMenu tem como objetivo gerar a matrix Droga X proteina  tanto como objeto RData quanto arquivo CSV

5.2.1 A primeira opcao deste menu 'Generate Drug Protein Matrix' tem como objetivo gerar a matriz de Droga(Linha)XProteina(Coluna).Para que este arquivo seja criado é necessario deixar o arquivo que contem a relacao droga proteina dentro da pasta 'Drug' com o nome 'drug_targets_DrugBank_Gysi.csv' e o arquivo de saida terá o nome 'drug_gene' dentro da pasta 'Drug' 

   


## Estrutura do Projeto


* `/src`: Contém os scripts R principais.
    * `main.R`: O script principal que orquestra a análise.
    * `utils.R`:  Funções auxiliares usadas no script principal.
    * `...`: (Outros scripts R, se houver).
    * * `/Data`: Contém os conjuntos de dados brutos, como `PPI_network.csv` e `drug_targets_DrugBank_Gysi.csv`.
* `Dockerfile`: Define a imagem Docker com o pacote 'tidyverse', instalando o R, as bibliotecas do sistema  R (como `igraph`, `dplyr`, `tidyr`).
* `docker-compose.yml`: Orquestra a inicialização do container.

## Observações

O uso do Docker garante que o ambiente (versões de pacotes, dependências do sistema) seja idêntico para todos os usuários, eliminando problemas de compatibilidade e facilitando a reprodutibilidade dos resultados.