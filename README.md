The article (TUO sh et al. (2019) MP-HS-DHSI: Multi-group Harmony Search Algorithm for Detecting High-Order SNP Interactions) has been submitted to Bioinformatics, and its current status is a minor revision.

the flowchart as :
![Flowchart of MP-HS-DHSI](https://github.com/shouhengtuo/MP-HS-DHSI/blob/master/outline.jpg)

 ***Figure . Flow chart of the MP-HS-DHSI algorithm. ***
This figure shows a detection process of 4-order SNP interaction, in which CIRCLE POINT  represents an SNP marker. [o o o o]  represents a 4-order SNP combination. A. is the set of all SNP loci. B. is the set of all feasible 4-order SNP combinations. C. aims to find 4-order SNP combinations with the strongest association with disease status using the harmony search (HS) algorithm; in the HS, four harmony memories respectively employ four evaluation criteria (K2-score, JS-score, LR-score and ND-JE) for calculating the association. D. After the HS is terminated, four harmony memories (HM1, HM2, HM3 and HM4) are merged into candidate set CS1. E. Each candidate SNP combination in CS is tested using the G-test statistical method; it will be included in candidate set CS2 if the SNP combination meets the significance level. F. The SNP combinations are examined using MDR; it will pass the examination if the classification accuracy meets the threshold.


Please kindly cite the paper if you use the code, datasets or any results in this repo or in the paper:

Feel free to contact Shouheng TUO (tuo_sh@126.com) <tuo_sh AT 126 DOT com; tuo_sh AT xupt DOT edu DOT cn>  for any questions about the paper, datsaets, code and results.

code list:

 **main program**
  ============
     |-- MP_HS_DHSI_FOR_DME.m  : the main program for DME models 
   
     |-- MP_HS_DHSI_FOR_NDME.m : the main program for NDME mdoels
   
     |-- HS_2019_multiCRITERIA5.m : the search program of the 1st stage.
   
     |-- multi_criteriaEvaluationFuns3.m : the score functions for evaluating the association
   
     |-- Gtest_score.m : the G-test function 
   
     |-- MDR2.m  : The MDR program for the 3rd stage.


Since the paper is still under review, the final code will be uploaded after the paper is received.
