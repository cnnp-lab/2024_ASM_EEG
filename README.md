# ASM tapering associated with delta band power reduction with dose, region and time-dependency

The code included on this repositoty generates figures and statistical analysis described on **ASM tapering associated with delta band power reduction with dose, region and time-dependency** by Guillermo M. Besne et. al. This work studies the effects of Anti-Seizure Medication tapering on delta band (1-4 HZ) power estimated from intracranial EEG, comparing 24h fragments at baseline ASM levels and reduced ASM levels.

## Provided Data
Under _data_ folder data used to obtain the published results is provided includinng: 
1) **iEEG_vs_Tap**: containintg a _.mat_ file for each subject with the two 24h windows of delta band power estimations and related information. These subjects have been subjected to ASM tapering
2) **MDL_Taper.mat**: File containtng the condensated informaition from *iEEG_vs_Tap*.
3) **MDL_nTaper.mat**: File containing similar information to *MDL_Taper.mat*, obtained from non-tapered patients, simulating the ASM tapering process

## Provided code
Three MatLab script are provided to generate the resutls for this publication. *gmb_ASM_EEG_P1* and *P2* are used to generate the results, while *gmb_ShadowPlot* assists on figure generation. Both scripts are configured to generate the expected resutls and not trequired configuration. 
### gmb_ASM_EEG_P1
Part 1 on the resutls generation workflow. It is not required to run this section to use *gmb_ASM_EEG_P2*. The scritp is configured to take as input the files on *iEEG_vs_Tap* and generates Figure 2 and Figure 5 from the Main text on the publication. Additionally, by uncomenting the last 2 lines on this script *MDL_Taper.mat* can be overwritten.
### gmb_ASM_EEG_P2
Part 2 on the resutls generation workflow. It is not required to run *gmb_ASM_EEG_P1* to use this section. This script takes as an input *MDL_Taper.mat* and *MDL_nTaper.mat* files to generate Figures 3 and 4 form the main text, Figure S3.1 from Supplementary material. Additionally, the followgin variables will contain results from statistical analysis and hierarchical modeling from both Main text and Supplemmentary materials: 
1) **RG_p_val & RG_ef_sz**: P-values and effect sizes studing if there is a statisticaly significant band power reduction (first row) and its correlation with ASM reduction (second row) considering different cortical segmentations.
2) **SJ_p_val & SJ_ef_sz**: P-values and effect sizes studing if there is a statisticaly significant band power reduction (first row) and its correlation with ASM reduction (second row) after averaging all RIOs fisrt and all subjects second.
3) **FB_p_val & FB_ef_sz**: P-values and effect sizes studing if there is a statisticaly significant band power reduction (first row) and relative band power reduction (second row).
4) **hier_mdl**: Variable contaitng the multple hierarchical models tested during the study to confirm our findings.
