1.) local_density (MATLAB) clusters into groups the GLCM parameters from diffraction images of cells

2.) DINet-R2 (Python) trains a deep learning network on the clusters and builds a model

3.) DINet-R2_testFeat (Python) tests the model, extracts features, and identifies support cells for each class

4.) reclass_corr (MATLAB) calculates a correlation density from the support cell features and reclassifies the non-support cells
