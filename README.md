
# Code review

## Fertility



1759 = cdg384
2360B5 = cdg225
2072 = cdg384_del
2301 = cdg338-384
D251 = cd225R
QA383P = cd384R

Transhet figures: The D251:2360 samples are actually D251HET and Trans-het and not WT and HET – see picture below:

Figure 2 survival after Smurf
qPCR


## Survival after blood feeding - 
colours and order of genotypes?



# Fig 2. Done! 
### This is weird and confusing!!! Not sure what to make of it? 

# Re-do with new data? 

# Fig 5 Transhet survival - is in blood_feeding.R Transehet

# Fig 1 and 3 are done?

# Fig 4 supplementation done

## 
I will resend the updated summary of what we agreed for each figure just to check we are all on the same page.

Figure 1 (FertilityData.xlsx): fecundity and fertility data of 1759, QA383P and 2360B5. We normally present the engorgement rate/hatching rate as a table, but if you have a different way of presenting it feel free to add it in the figure.
Figure 2 (SurvivalAfterBFData.xlsx): survival after blood-feeding of 1759, 2360B5, 2072, 2301, d251.
Figure 3 (SmurfAssayData.xlsx): There is an example of the figure in the draft. Unlike what we discussed, I would not include XA supplementation because in this experiment is where I saw rescue with NaOH and I don’t think the pH is good enough to argue it…
Figure 4 (XASupplementationData.xlsx): The 26mM was a typo in the headers, I have fixed it in the new document. For reference, the mM C- is the NaOH control in this dataset.
Figure 5 (TransHetSurvivalData.xlsx and qPCR.xlsx): There is an example of the figure in the draft. I was planning to add only one of the reference genes in the main manuscript, but I have added all the data.
Supplementary data

Hom survival (SurvivalData.xlsx): hom survival of 1759, cd KO or QA383P and 2072 (done and added to the file).


I will resend the updated summary of what we agreed for each figure just to check we are all on the same page.

Figure 1 (FertilityData.xlsx): fecundity and fertility data of 1759, QA383P and 2360B5. We normally present the engorgement rate/hatching rate as a table, but if you have a different way of presenting it feel free to add it in the figure.
Figure 2 (SurvivalAfterBFData.xlsx): survival after blood-feeding of 1759, 2360B5, 2072, 2301, d251.
Figure 3 (SmurfAssayData.xlsx): There is an example of the figure in the draft. Unlike what we discussed, I would not include XA supplementation because in this experiment is where I saw rescue with NaOH and I don’t think the pH is good enough to argue it…
Figure 4 (XASupplementationData.xlsx): The 26mM was a typo in the headers, I have fixed it in the new document. For reference, the mM C- is the NaOH control in this dataset.
Figure 5 (TransHetSurvivalData.xlsx and qPCR.xlsx): There is an example of the figure in the draft. I was planning to add only one of the reference genes in the main manuscript, but I have added all the data.
Supplementary data

Hom survival (SurvivalData.xlsx): hom survival of 1759, cd KO or QA383P and 2072 (done and added to the file).
 

Please do not hesitate to contact me if you have any questions or if I have missed anything.



#25.11.25
For the survival models - look at 6 and 5mM only???? Glu vs Xa? 
Remove Trp from all analyses. 

"To further explore survival patterns beyond the non-parametric Kaplan–Meier estimates, we fitted a series of parametric survival models to the data. Models were fitted separately for each supplementation group using the flexsurv package in  R. The following distributions were considered: exponential, Weibull, Gompertz, log-normal, and generalized gamma.
Each model included the covariates source, dosage, and genotype, as well as all two-way interaction terms (source × dosage, dosage × genotype, source × genotype). Models were fitted to the survival time variable (Hours) with censoring indicated by the event variable.
Model convergence was assessed, and any models that failed to converge were excluded from comparison. Model fit was evaluated using the Akaike Information Criterion (AIC), with lower AIC values indicating better relative fit. For each supplementation group, the best-fitting distribution was identified as the one with the minimum AIC.
To assess model adequacy, fitted survival curves from the best-fitting parametric model were overlaid with the corresponding Kaplan–Meier curves. This visual comparison provided an assessment of how well the parametric form captured the empirical survival patterns observed in the data"

Double check values - for dosage. 

# cardinal_fitness

1759: cdg384
2072: cdg384_del
2301: cdg338-384
2360: cdg225
QA383P or KO: cd384
D251: cd225


## Analysis notes

Done - Smurf after feeding - an effect of mechanism and genotype (no evidence of an interaction effect)

To-do Effect of smurf on survival? 

###

As a brief summary, all our KI cardinal lines present different degrees of death after blood feeding (which was also observed in kmo KI in An. stephensi). So, I did a bunch of experiments which they also did in a published manuscript (attaching below for reference) with a KI line that had approximately 50% of mortality after blood feeding. 

Experiment 1: XA and Trp supplementation. We hypothesized that supplementing sugar and blood with tryptophan (trp) (upstream of the pathway) would increase the phenotype and supplementing with xanthurenic acid (XA) (bottom of the pathway) would rescue it. Since the phenotype is survival after blood-feeding, all the results are shown in survival curves. The different tabs are for the four experiments that I did:  
•	Blood supplementation with XA 
•	Sugar supplementation with XA
•	Blood supplementation with Trp
•	Sugar supplementation with Trp
You will see that 6mM of XA was able to rescue and concentrations Trp above 10mM were able to rescue too. 

Experiment 2: smurf assay. Here I performed a smurf assay to assess midgut leakage with supplemented and non-supplemented mosquitoes and we saw that NaOH without XA rescued the phenotype better than XA and Trp. Therefore, we thought that maybe we were rescuing the phenotype by changing the pH in the midgut. 

Experiment 3: glycine supplementation. glycine is a neutral amino acid like tryptophan and we wanted to see if a neutral amino acid outside of the pathway was able to rescue the phenotype. In this case, some did and some did not. 


###
For the Smurf assay, I have added the total number of females that started in the cage. For those that were not included, I have added a column with the reason why. Furthermore, for the Smurf screening, I added the raw data using the survival sheet as a reference (each female number represents a female with the smurf group and when it died). 

Finally, in the supplementation document, I have added an extra sheet with a collection of assays that could be used to assess the variability of the 2360HEt non-supplemented control (highlighted in blue for all the experiments). 

Please let me know if there is anything that doesn't make sense or if there is anything else you may need. 


Yes, Smurf and XA data is unchanged and you should have it already.

## Phenotypes
In the phenotypes section we are going to include (@Michelle Anderson, correct me if I am wrong):
Homozygous viability of 1759, KO (QA383P), 2072, 2360 and D251. You should have the data for 1759, QA383P and 2072 in the Data.xlsx file. I think there has been some problems generating the cross for 2360  The hom viability will be the "survival" tabs.
We are no longer able to recover homozygotes for 2360, 100% of the females are dying after BF now so we will not have homozygous viability data for them beyond this observation - so no new data here.

## Fertility
Fertility assays of 1759, KO(QA383P), 2360 and potentially D251. Fertility for 1759 and KO should be in the Data.xlsx. This will be the "fertility" tabs
Fertility for 2360 hets vs WT has been set up this week - we'll have the data in ~2 weeks time. 
Are we doing d251? I don't think this experiment is currently planned.

## Survival
Survival after blood-feeding: you should have the survival after BFing of all the lines (2301 unchanged, 2072 unchanged, 1759 unchanged, 2360 unchanged and D251, which is new and has been added to the same Survival after BFing.xlsx file)
You should have access to the data, but I will send you everything later today.
