% save the normalized data for each receptor seperately and named them as hampa/mampa etc.
% check the data distribution
[h,p]=chi2gof(hampa)
hist(hampa)
% normalization
% zscore normalization
znor_mdpat=zscore(mdpat)
% minmax normalization
[m, n]  = size(mampa);
normalized_mampa = zeros(m, n);
    for i = 1:n
        ma = max( mampa(:, i) ); 
        mi = min( mampa(:, i) );
        normalized_mampa(:, i) = (mampa(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(mcgp5);
normalized_mcgp5 = zeros(m, n);
    for i = 1:n
        ma = max( mcgp5(:, i) ); 
        mi = min( mcgp5(:, i) );
        normalized_mcgp5(:, i) = (mcgp5(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(mdamp);
	normalized_mdamp = zeros(m, n);
    for i = 1:n
        ma = max( mdamp(:, i) ); 
        mi = min( mdamp(:, i) );
        normalized_mdamp(:, i) = (mdamp(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(mdpat);
	normalized_mdpat = zeros(m, n);
    for i = 1:n
        ma = max( mdpat(:, i) ); 
        mi = min( mdpat(:, i) );
        normalized_mdpat(:, i) = (mdpat(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(mflum);
	normalized_mflum = zeros(m, n);
    for i = 1:n
        ma = max( mflum(:, i) ); 
        mi = min( mflum(:, i) );
        normalized_mflum(:, i) = (mflum(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(mkain);
	normalized_mkain = zeros(m, n);
    for i = 1:n
        ma = max( mkain(:, i) ); 
        mi = min( mkain(:, i) );
        normalized_mkain(:, i) = (mkain(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(mketa);
	normalized_mketa = zeros(m, n);
    for i = 1:n
        ma = max( mketa(:, i) ); 
        mi = min( mketa(:, i) );
        normalized_mketa(:, i) = (mketa(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(mmk80);
	normalized_mmk80 = zeros(m, n);
    for i = 1:n
        ma = max( mmk80(:, i) ); 
        mi = min( mmk80(:, i) );
        normalized_mmk80(:, i) = (mmk80(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(mmusc);
	normalized_mmusc = zeros(m, n);
    for i = 1:n
        ma = max( mmusc(:, i) ); 
        mi = min( mmusc(:, i) );
        normalized_mmusc(:, i) = (mmusc(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(moxot);
	normalized_moxot = zeros(m, n);
    for i = 1:n
        ma = max( moxot(:, i) ); 
        mi = min( moxot(:, i) );
        normalized_moxot(:, i) = (moxot(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(mpire);
	normalized_mpire = zeros(m, n);
    for i = 1:n
        ma = max( mpire(:, i) ); 
        mi = min( mpire(:, i) );
        normalized_mpire(:, i) = (mpire(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(mpraz);
	normalized_mpraz = zeros(m, n);
    for i = 1:n
        ma = max( mpraz(:, i) ); 
        mi = min( mpraz(:, i) );
        normalized_mpraz(:, i) = (mpraz(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(msch2);
	normalized_msch2 = zeros(m, n);
    for i = 1:n
        ma = max( msch2(:, i) ); 
        mi = min( msch2(:, i) );
        normalized_msch2(:, i) = (msch2(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(muk14);
	normalized_muk14 = zeros(m, n);
    for i = 1:n
        ma = max( muk14(:, i) ); 
        mi = min( muk14(:, i) );
        normalized_muk14(:, i) = (muk14(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(hampa);
	normalized_hampa = zeros(m, n);
    for i = 1:n
        ma = max( hampa(:, i) ); 
        mi = min( hampa(:, i) );
        normalized_hampa(:, i) = (hampa(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(hcgp5);
	normalized_hcgp5 = zeros(m, n);
    for i = 1:n
        ma = max( hcgp5(:, i) ); 
        mi = min( hcgp5(:, i) );
        normalized_hcgp5(:, i) = (hcgp5(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(hdamp);
	normalized_hdamp = zeros(m, n);
    for i = 1:n
        ma = max( hdamp(:, i) ); 
        mi = min( hdamp(:, i) );
        normalized_hdamp(:, i) = (hdamp(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(hdpat);
	normalized_hdpat = zeros(m, n);
    for i = 1:n
        ma = max( hdpat(:, i) ); 
        mi = min( hdpat(:, i) );
        normalized_hdpat(:, i) = (hdpat(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(hflum);
	normalized_hflum = zeros(m, n);
    for i = 1:n
        ma = max( hflum(:, i) ); 
        mi = min( hflum(:, i) );
        normalized_hflum(:, i) = (hflum(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(hkain);
	normalized_hkain = zeros(m, n);
    for i = 1:n
        ma = max( hkain(:, i) ); 
        mi = min( hkain(:, i) );
        normalized_hkain(:, i) = (hkain(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(hketa);
	normalized_hketa = zeros(m, n);
    for i = 1:n
        ma = max( hketa(:, i) ); 
        mi = min( hketa(:, i) );
        normalized_hketa(:, i) = (hketa(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(hmk80);
	normalized_hmk80 = zeros(m, n);
    for i = 1:n
        ma = max( hmk80(:, i) ); 
        mi = min( hmk80(:, i) );
        normalized_hmk80(:, i) = (hmk80(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(hmusc);
	normalized_hmusc = zeros(m, n);
    for i = 1:n
        ma = max( hmusc(:, i) ); 
        mi = min( hmusc(:, i) );
        normalized_hmusc(:, i) = (hmusc(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(hoxot);
	normalized_hoxot = zeros(m, n);
    for i = 1:n
        ma = max( hoxot(:, i) ); 
        mi = min( hoxot(:, i) );
        normalized_hoxot(:, i) = (hoxot(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(hpire);
	normalized_hpire = zeros(m, n);
    for i = 1:n
        ma = max( hpire(:, i) ); 
        mi = min( hpire(:, i) );
        normalized_hpire(:, i) = (hpire(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(hpraz);
	normalized_hpraz = zeros(m, n);
    for i = 1:n
        ma = max( hpraz(:, i) ); 
        mi = min( hpraz(:, i) );
        normalized_hpraz(:, i) = (hpraz(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(hsch2);
	normalized_hsch2 = zeros(m, n);
    for i = 1:n
        ma = max( hsch2(:, i) ); 
        mi = min( hsch2(:, i) );
        normalized_hsch2(:, i) = (hsch2(:, i)-mi ) / ( ma-mi );
    end
[m, n]  = size(huk14);
	normalized_huk14 = zeros(m, n);
    for i = 1:n
        ma = max( huk14(:, i) ); 
        mi = min( huk14(:, i) );
        normalized_huk14(:, i) = (huk14(:, i)-mi ) / ( ma-mi );
    end
	
% save normalized human data into human_all.csv
% save normalized macaque data macaque_all.csv

% stepwise linear mixed-effects models for human
library(lmerTest)
library(tidyverse)
library(emmeans)
Data=read.csv('.../human_all.csv', header = T, sep = ',')
head(Data)
Data$brain=factor(Data$brain)
Data$area=factor(Data$area)
Data$receptor=factor(Data$receptor)

% first level omnibus test 
Modelmax=lmer(data=Data,
              formula=density_mm_nor~area*receptor+(1|brain),
              control=lmerControl(optimizer = 'bobyqa'))
isSingular(Modelmax)
summary(Modelmax)
anova(Modelmax)

% second level simple effect tests
sink("secend_level_results.txt")
joint_tests(Modelmax,by='receptor')
sink()
% save second_level_pvalue as pvalue_second
% FDR-corrected p-values for the second level simple effect tests 
p.adjust(pvalue_second,method = "fdr")

% third level post hoc tests with FDR-correction
emmeans(Modelmax, pairwise~area|receptor,adjust='fdr')

% stepwise linear mixed-effects models for macaque
library(lmerTest)
library(tidyverse)
library(emmeans)
Data=read.csv('.../macaque_all.csv', header = T, sep = ',')
head(Data)
Data$brain=factor(Data$brain)
Data$area=factor(Data$area)
Data$receptor=factor(Data$receptor)

% first level omnibus test 
Modelmax=lmer(data=Data,
              formula=density_mm_nor~area*receptor+(1|brain),
              control=lmerControl(optimizer = 'bobyqa'))
isSingular(Modelmax)
summary(Modelmax)
anova(Modelmax)

% second level simple effect tests
sink("secend_level_results.txt")
joint_tests(Modelmax,by='receptor')
sink()
% save second_level_pvalue as pvalue_second
% FDR-corrected p-values for the second level simple effect tests 
p.adjust(pvalue_second,method = "fdr")

% third level post hoc tests with FDR-correction
emmeans(Modelmax, pairwise~area|receptor,adjust='fdr')
