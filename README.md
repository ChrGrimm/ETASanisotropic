# ETASanisotropic
This project contains Matlab code to estimate an Epidemic Type Aftershock Sequence (ETAS) model. The implementation advances the standard ETAS model by enabling anisotropic spatial kernels and accounting for short-term aftershock incompleteness.

Please refer to the user manuals for precise instructions on how to use the code.

You may start using the code with the exemplary earthquake dataset for Southern California, which has been downloaded from the Southern California Earthquake Data Center (https://scedc.caltech.edu/data/alt-2011-dd-hauksson-yang-shearer.html, last accessed on October 25, 2021). 
You can run the model either for the long-term dataset (1981-2019), or for the short extract of the 2019 M6.4 Ridgecrest foreshock sequence. You find the corresponding polygons and ini files prepared in the data folder as well.

Further citations:

•	Christian Grimm implemented the code in Matlab R2019a.

•	The standard functionality of the model goes back to the implementation of Jalilian (2019) in the R package ETAS:
Jalilian, A. (2019). ETAS: An R package for fitting the space-time ETAS model to earthquake data, J. Stat. Softw. 88, no. 1, 1–39, doi: 10.18637/jss.v088.c01.

•	Additional features such as anisotropic spatial kernels and short-term aftershock incompleteness refer to the model theory described in detail in the following papers:

C. Grimm, M. Käser, S. Hainzl, M. Pagani, and H. Küchenhoff (2021). Improving Earthquake Doublet Frequency Predictions by Modified Spatial Trigger Kernels in the Epidemic-Type Aftershock Sequence (ETAS) Model, Bull. Seismol. Soc. Am. XX, 1–20, doi: 10.1785/0120210097

Hainzl, S. (2021). ETAS-Approach Accounting for Short-Term Incompleteness of Earthquake 753 Catalogs. Bulletin of the Seismological Society of America. doi: 10.1785/0120210146.

C. Grimm, S. Hainzl, M. Käser, and H. Küchenhoff (2022). Solving three major biases of the ETAS model to improve forecasts of the 2019 Ridgecrest sequence., doi: https://doi.org/10.21203/rs.3.rs-1128731/v1 (accepted)

