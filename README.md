# FrontFinder - an OHW22 project
Project proposed by Sophie Clayton. Repo started by Felipe and Maya

---
## One line description
---
### Collaborators
- Sophie Clayton (mentor)
- Maya Jakes
- Felipe Vilela da Silva
- Alessio Arena
- Mackenzie Blanusa
- Chaojiao Sun

---
### Background

### Goals

- Load and preprocess saildrone data
- Load and preprocess SSH data
- combine both datasets
- develop gradient based edge detection

### Datasets

- [Saildrone](https://data.saildrone.com/data/sets/antarctica-circumnavigation-2019), also available on [ERDDAP](https://erddap.ifremer.fr/erddap/info/index.html?page=1&itemsPerPage=1000)
  Current variables used:
   - TEMP_CTD_RBR_STDDEV
   - TEMP_CTD_RBR_MEAN
   - SAL_RBR_MEAN
   - SAL_RBR_STDDEV
   - O2_CONC_AANDERAA_MEAN
   - O2_CONC_AANDERAA_STDDEV
   - CHLOR_RBR_MEAN
   - CHLOR_RBR_STDDEV
- [Sea Surface Height](https://resources.marine.copernicus.eu/product-detail/SEALEVEL_GLO_PHY_L4_MY_008_047/INFORMATION)



### Workflow

### References

### Progress notes:
- Day 1
As an initial exercise, Maya and Felipe explored some Saildrone data (hydrography) and SSH from altimeters in the Southern Ocean between the beginning of 2019 and Sep 2019. We got stuck trying to figure out how to upload CMENS data with xr.open_dataset(). We are trying to add the user info (username and password) into the code.

- Day 2
    - Maya: reviewed and improved Sophie's code to calculate and classify gradients on time series (1D data)
    - Alessio: developed a preprocessing routine for saildrone data. This performs a spatial interpolation for missing lat/lon, then resample to 5 minutes and finally convert to spatial dataset (geodataframe)
        - see saildrone_processing.py, and view_data_AA.ipynb
    - Felipe: processed SSH data to calculate gradient and threshold. Compare and merge saildrone measurement with SSH in space and time
        - see gradient_peaks.ipynb
