# Exposure and sensitivity of terrestrial vertebrates to biological invasions worldwide.

[https://doi.org/10.5061/dryad.31zcrjdvh](https://doi.org/10.5061/dryad.31zcrjdvh)

Data for reproducing the results obtained by Marino et al. in the paper **Exposure and sensitivity of terrestrial vertebrates to biological invasions worldwide**.

## Description of the data and file structure

The data for analyses are divided into four main files related to the exposure and sensitivity of birds, mammals, and reptiles to biological invasions, as well as data completeness associated with these two metrics. They are all in [`data_for_analyses/`](https://github.com/claramarino/BIVA_framework/main/Data/data-for-analyses).  We also provide one additional spatial file for mapping the results and get the associated cells with all metric values, located in [`derived_data/`](https://github.com/claramarino/BIVA_framework/main/Data/derived_data).

### 1. Exposure to 304 IAS

Data file 10_Exposure_normalized_110_km.rds contains the exposure components and final values for each cell with at least one established alien species occurrence, for three normalization methods  (4 sheets x 6,629 lines x 16 columns).

| Column name            | Description                                                                                |
  | :--------------------- | :----------------------------------------------------------------------------------------- |
  | grid\_id               | ID of the cell in the spatial file of 110km land cell units                                |
  | SR\_tot\_ias\_max\_min | Alien species richness  normalized with max-min method                                     |
  | med\_range\_max\_min   | Median alien range normalized with max-min method                                          |
  | med\_ib\_max\_min      | Median alien impact breadth normalized with max-min method                                 |
  | expo\_max\_min         | Final exposure as sum of components with max-min normalization method                      |
  | expo\_prod\_max\_min   | Alternative exposure as product of components with max-min normalization method            |
  | SR\_tot\_ias\_log      | Alien species richness  normalized with log-transformed method                             |
  | med\_range\_log        | Median alien range normalized with log-transformed method                                  |
  | med\_ib\_log           | Median alien impact breadth normalized with log-transformed method                         |
  | expo\_log              | Alternative exposure as sum of components with log-transformed normalization method        |
  | expo\_prod\_log        | Alternative exposure as product of components with log-transformed normalization method    |
  | SR\_tot\_ias\_rank     | Alien species richness  normalized with cumulative ranking method                          |
  | med\_range\_rank       | Median alien range normalized with cumulative ranking method                               |
  | med\_ib\_rank          | Median alien impact breadth normalized with cumulative ranking method                      |
  | expo\_rank             | Alternative exposure as sum of components with cumulative ranking normalization method     |
  | expo\_prod\_rank       | Alternative exposure as product of components with cumulative ranking normalization method |
  
  ### 2. Sensitivity of birds, mammals, and reptiles
  
  The data file 13_Sensitivity_normalized_110_km.rds contains the sensitivity values for birds, mammals, and reptiles for each cell for three normalization methods (3 sheets x 16,052 lines x 4 columns).

| **Column name**      | **Description**                                                   |
  | :------------------- | :---------------------------------------------------------------- |
  | cell\_id             | ID of the cell in the spatial file of 110km land cell units       |
  | SR\_ias\_a\_max\_min | Final sensitivity  metric normalized with max-min method          |
  | SR\_ias\_a\_log      | Alternative sensitivity normalized with log-transformation method |
  | SR\_ias\_a\_rank     | Alternative sensitivity normalized with cumulative ranking method |
  
  ### 3. Completeness for exposure
  
  The data file 10_Completeness_Exposure_110_km.rds contains the completeness values associated with the exposure of all taxa for each cell (17,258 lines x 7 columns).

| **Column name**    | **Description**                                                               |
  | :----------------- | :---------------------------------------------------------------------------- |
  | grid\_id           | ID of the grid cell in the spatial file of 110km land cell units              |
  | log\_dens          | Log-transformed density of GBIF records for vertebrates from Meyer et al 2015 |
  | DB\_and\_GRIIS     | Proportion of alien species in both our database and the GRIIS checklist (ed) |
  | DB\_GRIIS\_or\_not | Proportion of alien species in  our database compared to total alien pool     |
  | norm\_dens         | Normalized density of GBIF records for vertebrates from Meyer et al 2015 (se) |
  | comp\_se\_sum      | Alternative completeness (sum of se and ed)                                   |
  | comp\_se\_product  | Final completeness for exposure (product of se and ed)                        |
  
  ### 4. Completeness for sensitivity
  
  The data file 14_Completeness_Sensitivity_110_km.rds contains the completeness values associated with the sensitivity of birds, mammals, and reptiles for each cell (3 sheets x 16,052 lines x 10 columns).

| **Column name** | **Description**                                                                     |
  | :-------------- | :---------------------------------------------------------------------------------- |
  | cell\_id        | ID of the cell in the spatial file of 110km land cell units                         |
  | SR\_tot         | Total species richness of the given taxon                                           |
  | SR\_ias\_a      | Number of species threatened by IAS in the IUCN RedList                             |
  | SR\_ias\_t      | Number of species threatened by IAS and at high extinction risk in the IUCN RedList |
  | SR\_ias\_nt     | Number of species threatened by IAS and at low extinction risk in the IUCN RedList  |
  | SR\_dd\_iasa    | Number of species threatened by IAS and Data Deficient in the IUCN RedList          |
  | SR\_dd\_all     | Number of species Data Deficient in the IUCN RedList                                |
  | se              | Sampling effort for sensitivity                                                     |
  | ed              | Effective detection for sensitivity                                                 |
  | comp\_prod      | Completeness for sensitivity (product of se and ed)                                 |
  
  ### 5. Land cell units (spatial file)
  
  The file 04_Grid_110km.rds contains the spatial information associated with all the cells for which we calculated exposure, sensitivity, and completeness. The associated dataset contains a column "grid_id" for linking with exposure data, and a column "cell_id" for linking with sensitivity data.

## Sharing/Access information

Raw data was derived from the following sources:

* Alien occurrences data were downloaded from the Global Biodiversity Information Facility (www.gbif.org). All DOIs associated with download requests on GBIF can be found in the file GBIF_download_doi.xlsx.
* The exotic and native ranges of alien species used to filter the occurrences were derived from multiple sources (see Table S1 of the associated paper) including the IUCN Red List of threatened species (Version 2022-1,  https://www.iucnredlist.org/resources/spatial-data-download), the Global Register of Introduced and Invasive Alien Species (www.griis.org), the Global Assessment of Reptile Distributions (Version 1.1, https://doi.org/10.5061/dryad.83s7k), Birdlife International (Version 2020.1, https://datazone.birdlife.org/species/requestdis), and the CABI compendium in Invasive Species (https://www.cabidigitallibrary.org/product/qi). Details on the sources associated with each alien species are provided in the file Alien_species_range_sources.xlsx.
* Native species ranges were accessed using the IUCN Red List of threatened species (Version 2022-1,  https://www.iucnredlist.org/resources/spatial-data-download) for mammals and reptiles, and Birdlife International (Version 2020.1, https://datazone.birdlife.org/species/requestdis) for birds.

## Code/Software

All analyses were conducted using R software version 4.2.2 (R Core Team, 2022). The associated scripts are provided in the [GitHub repository](https://github.com/claramarino/BIVA_framework).
