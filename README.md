# Data and code for Bauer et al. (2023) J Veg Sci

_Markus Bauer <a href="https://orcid.org/0000-0001-5372-4174"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height = "16"></a>, Jakob K. Huber, and Johannes Kollmann <a href="https://orcid.org/0000-0002-4990-3636"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height = "16"></a>_  

Data and code for:

Bauer M, Huber JK, & Kollmann J (2024) __Beta diversity of restored dike grasslands is strongly influenced by uncontrolled spatio-temporal variability.__ &ndash; _Journal of Vegetation Science_ 35, e13293.

[![DOI:10.1111/jvs.13293](http://img.shields.io/badge/DOI-10.1111/jvs.13293-informational.svg)](https://doi.org/10.1111/jvs.13293)

**Study region**: [River Danube around Deggendorf](https://www.openstreetmap.org/#map=11/48.8127/12.9790)
<br>
<br>
## Content of the repository

1. The folder `data` contains  
    * `Raw` and `processed` data of the sites variables (.csv) 
    * `Raw` and `processed` data of the species' abundances (.csv) 
    * `Raw` and `processed` data of the species' traits (.csv)
    * Raw data of monthly `temperature` and `precipiation` data (.csv)
    * Raw and processed `spatial` data (.shp)
    * `photos` of the plots (.jpg)
    
2. The folder `outputs` contains  
    * The `figures` generated (.tiff)
    * The `statistics` tables from the principal component analyses (.csv)
    * The `tables` generated (.png)
    * The `models` calculated (.Rdata)
    
3. The folder `R` contains  
    * Scripts to calculate all models (.R)
    * Scripts to generate all figures and tables (.R)
    * Metadata script for creating EML file (.R)
    * Folder for calculating habitat types (ESY)
    
4. The folder `markdown` contains model checks to scroll through

#### Package versioning

The used versions of R and the packages are saved in `2023_danube_dike_survey/renv.lock`.

You can restore this state by executing `renv::restore()` in the console.

## Citation

[![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg


When using the __data available__ in this repository, please cite the original publication and the dataset.  

__Publication__

> Bauer M, Huber JK, & Kollmann J (2024) Beta diversity of restored dike grasslands is strongly influenced by uncontrolled spatio-temporal variability. &ndash; *Journal of Vegetation Science* 35, e13293. [https://doi.org/10.1111/jvs.13293](https://doi.org/10.1111/jvs.13293)

__Dataset__

> Bauer M, Huber JK, & Kollmann J (2024) Data and code for Bauer et al. (2024) J Veg Sci [Data set]. &ndash; *Zenodo*. [https://doi.org/10.5281/zenodo.6107806](https://doi.org/10.5281/zenodo.6107806)

This dataset is also linked to PANGAEA
> Bauer M, Huber JK & Kollmann J (2026) Vegetation surveys on dike grasslands between Straubing and Vilshofen at the river Danube &ndash; Sites [dataset]. &ndash; *PANGAEA*. [https://doi.org/10.1594/PANGAEA.974754](https://doi.org/10.1594/PANGAEA.974754)

> Bauer M, Huber JK & Kollmann J (2026) Vegetation surveys on dike grasslands between Straubing and Vilshofen at the river Danube &ndash; Vegetation [dataset]. &ndash; *PANGAEA*. [https://doi.org/10.1594/PANGAEA.973510](https://doi.org/10.1594/PANGAEA.973510)

Contact markus1.bauer@tum.de for any further information.  
