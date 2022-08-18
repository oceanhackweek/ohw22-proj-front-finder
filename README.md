# FrontFinder - an OHW22 project
Project proposed by Sophie Clayton. Repo started by Felipe and Maya

---
## One line description
---
### Collaborators and roles
- Sophie Clayton (mentor)
- Maya Jakes
- Felipe Vilela da Silva (conceptualist and coder for the approach of exploring the in situ observations at the position of peaks in the SSH gradient available from altimeters)
- Alessio Arena
- Mackenzie Blanusa
- Chaojiao Sun

---
### Background
Ocean fronts play a key role in ocean dynamics and biogeochemical processes. There are different ways of defining fronts such as gradient thresholding, statistical approaches, and edge-detection. Ocean fronts may occur at various spatial scales and are often studied at sub-mesoscale (0-10km) and mesoscale (10-100km) resolutons. At large scales satellite altimeter data may be very useful for identifying fronts, but at small scales, finer resolution in-situ data is needed. This OHW project utliizes saildrone and satellite sea surface height (SSH) data to identify fronts along a 1D time series and compare them spatially. 

### Goals

- Load and preprocess saildrone data
- Load and preprocess SSH data
- combine both datasets
- develop gradient based edge detection
- develop functions to identify fronts along a 1D time series using gradient thresholding 
- compare fronts identified in the saildrone data with adcp data 

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
- Day 3
    - Felipe: combined the two versions of the frontal detection algorithms based on (i) peaks in the SSH gradient and 
    (ii) high gradients in surface data provided by in situ observations. Also, I started a notebook of a suggestion for the OHW presentation. 

- Day 2
    - Maya: reviewed and improved Sophie's code to calculate and classify gradients on time series (1D data)
    - Alessio: developed a preprocessing routine for saildrone data. This performs a spatial interpolation for missing lat/lon, then resample to 5 minutes and finally convert to spatial dataset (geodataframe)
        - see saildrone_processing.py, and view_data_AA.ipynb
    - Felipe: processed SSH data to calculate gradient and threshold. Also, compared and merged saildrone measurement with SSH in space and time
        - see gradient_peaks.ipynb
    - Mackenzie: Reviewed progress and started to convert Maya's code to xarray 

- Day 1
As an initial exercise, Maya and Felipe explored some Saildrone data (hydrography) and SSH from altimeters in the Southern Ocean between the beginning of 2019 and Sep 2019. We got stuck trying to figure out how to upload CMENS data with xr.open_dataset(). We are trying to add the user info (username and password) into the code.

Mackenzie accessed the saildrone CTD (1 min resolution) data and adcp (5 min) resolution data. She began preprocessing the data which includes removing Nans and  resampling the saildrone 1 min data to 5 min.
