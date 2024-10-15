# WNS risk southern hemisphere
[![license](https://img.shields.io/badge/license-MIT%20+%20file%20LICENSE-lightgrey.svg)](https://choosealicense.com/)
![Open Source
Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)

This repository contains code and data needed to reproduce the article:

**Wu N. C., Welbergen J. A., Villada-Cadavid, T., Lumsden, L. F., & Turbill C.** (2024) Vulnerability of southern hemisphere bats to white-nose syndrome based on global analysis of fungal host specificity and cave temperatures. e14390. DOI: [![DOI](https://zenodo.org/badge/DOI/10.1111/cobi.14390.svg)](https://doi.org/10.1111/cobi.14390)


**Raw data**
- `bat_comp_data.csv` - Raw data used for the analysis.
- `hibernation_study.csv` - Data used for comparing hibernation studies between the northern and southern hemisphere.
- shp - Contains shape files from the [World Karst Aquifer Map](https://produktcenter.bgr.de/terraCatalog/OpenSearch.do?search=ab3b15cb-a6c3-42ea-ae0c-0b417d698949&type=/Query/OpenSearch.do). 

***R*** **code**
- `bat_comp_analysis_final.R` - Data cleaning, analysis and figure production.

## Abstract
White-nose syndrome (WNS), a disease affecting hibernating bats, is caused by the fungal pathogen *Pseudogymnoascus destructans* (*Pd*). Since the initial introduction of Pd from Eurasia to the United States in 2006, WNS has killed millions of bats throughout the temperate parts of North America. There is concern that if Pd is accidentally introduced to the Southern Hemisphere, WNS could pose similar threats to the bat fauna of the Southern Hemisphere's more temperate regions. Efforts are required to better understand the vulnerability of bats globally to WNS. We examined phylogenetic distances among cave roosting bat species globally to estimate the probability of infection by *Pd*. We predicted cave thermal suitability for *Pd* for 441 cave-roosting bat species across the globe via spatial analysis. We used host specificity models based on 65 species tested for Pd to determine phylogenetic specificity of *Pd*. Phylogenetic distance was not an important predictor of *Pd* infection, confirming that *Pd* has low host specificity. We found extensive areas (i.e., South America, Africa, and Australia) in the Southern Hemisphere with caves that were suitable for cave-roosting bat species and for *Pd* growth. Hence, if *Pd* spreads to the Southern Hemisphere, the risk of exposure is widespread for cave-roosting bats, and infection is possible regardless of relatedness to infected species in the Northern Hemisphere. Predicting the consequences of infection remains difficult due to lack of species-specific information about bat winter biology. Nevertheless, WNS is an important threat to naive Southern Hemisphere bat populations. Hence, biosecurity measures and planning of management responses that can help prevent or minimize a potential WNS outbreak in the Southern Hemisphere are urgently needed.


**Keywords:** Chiroptera, hibernation, *Pseudogymnoascus destructans*, disease, vulnerability ecological naivety.


## License
This repository is provided by the authors under the MIT License ([MIT](http://opensource.org/licenses/MIT)).
