CSG toolbox for EEG data by Dorothée Coppieters
% ---------------------------------------------
To install the CSG toolbox:

1) install SPM 12
2) insert the 'fasst' and 'CSG_012017' folders in the Matlab setpath (CSG_012017 has to be on the top of the set path to have the priority).


To preprocess EEG data:
1) use the function CSG_preprocessing and select the raw file.
It computes the mat and dat files from SPM,
filters data (it creates file starting by 'P')
detects artefacts (file starting by 'IP_').

To display data:
use the CSG command in Matlab and click on 'display one file'. Select the 'IP' file created after the preprocessing.