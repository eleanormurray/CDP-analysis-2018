# CDP-analysis-2018
The repository CDP-analysis-2018 archives the SAS code for updated adherence-adjusted analyses in the placebo arm of the Coronary Drug Project trial using a survival analysis approach, published in Trials in 2018. [1] This code is also included in the SAS appendix of the Trials paper.

This repository contains the SAS programs used for the updated analysis - unadjusted, baseline standardized, and inverse probabiltiy weighted survival analyses comparing adherers and non-adherers in the placebo arm under a range of modeling assumptions about the appropriate functional form of adherence. SAS 9.4 was used for all analyses. The code appendix contains the following programs:

Program 1: Data management This program takes as input the cleaned CDP dataset used in our previous analyses (see CDP-analysis-2016), and available through application to the National Heart, Lung, and Blood Institute (NHLBI).

Program 2: Analysis This program contains a macro that can produce all output in the Table of the paper, by setting macro options for &adjust (0 = unadjusted, 1 = baseline standardized only, 2 = with inverse probability weighting), and &model (see code for the possible values corresponding to the Table analyses). 

The SciFigure is a visual summary comparing the original 1980 NEJM placebo arm analysis with our reproduction and replication efforts. The 2016 update uses inverse probability weighting and baseline standardization to obtain a null estimate of the effect of placebo adherence on the cumulative incidence of mortality, while the 2018 update uses the same techniques to obtain a null estimate of the effect of placebo adherence on survival time.

The datasets for the Coronary Drug Project will be available through application to the NHLBI data center. If you have any questions, comments or discover an error, please create an issue or contact me at emurray@mail.harvard.edu. For other similar programs, please visit www.hsph.harvard.edu/causal/.

Reference:

    Murray EJ, Hernan MA. Improved adherence adjustment in the Coronary Drug Project. Trials, 2018; 19:158. https://doi.org/10.1186/s13063-018-2519-5
    Murray EJ, Hernan MA. Adherence adjustment in the Coronary Drug Project: A call for better per-protocol effect estimates in randomized trials. Clinical Trials, 2016; 13(4):372-8.



