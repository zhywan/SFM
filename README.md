# SFM
Strategic Flipping Method for 2016 iDASH Challenge Track 1

I assume that the attacker can query SNPs in any order.

build_beacon_delta_rare_mitigate.py is the file contains the mitigation algorithm. Other files may have changed accordingly.

The mitigation algorithm includes four components:

1) Sorting strategy;

2) Mitigation strategy;

3) Noise-adding strategy; and

4) Optimization according to the evaluation metric.

These strategies should be applied in order.

1) Sorting strategy:

   All the SNPs are sorted according to a certain metric: 

   i) minor allele frequency, 

   ii) absolute difference between the pool and the underlying reference, 

   iii) difference between the pool and the underlying reference, 

   iv) log-likelihood ratio corresponding to each SNP, or

   v) default order.

2) Mitigation Strategy:

   The result from the server will be flip from "true" to "false" for a certain SNP if 

   i) it is within top K SNPs,

   ii) its sorting metric is within a threshold, 

   iii) it is within top k% SNPs, or 

   iv) same as iii but the value of "k" is determined automatically. For example, if the power without any mitigation strategy is 0, then k will be 0. If the power without any mitigation strategy is 1, then k will be kappa. If the power without any mitigation strategy is rho, then k should be kappa*rho. Note that to use this strategy, the whole process should be run beforehand once with a setting where K=0. The limitation for this choice is that it requires a provided case dataset.

3) Noise-adding strategy:

   This step brings some uncertainty to the result of the server, which would increase the accessing cost and false positive rate of the attacker.

4) Optimization according to the evaluation metric:

   As long as the defender has an evaluation metric, and he/she has enough computational budget, he/she can explorer the search space to obtain the global (local) optimum through various search strategy like random walk which is similar to the noise-adding strategy. However, this strategy requires a provided case dataset.

The commmands need to be executed including:

  * Precomputation stage

   "python vcf_bit.py -f data/diybu/challengeData16/c1beacon/wholeRecords/vcfSize500/chr10_selectedInd.vcf -n 500"

   "python vcf_bit.py -f data/diybu/challengeData16/c1beacon/wholeRecords/500notInBeacon/chr10_selectedInd.vcf -n 500"

  * Iteration stage

   "python build_beacon_delta_rare_mitigate.py -t 1 -d 1e-06"

   "python Query_fromBeacon_delta_wholeRecords.py -t 1 -d 1e-06"

   "python power_lrt.py"

All the input/output file names are almost the same with the provided sample codes. To run the code, please first add the corresponding input files to the right locations:

   "data/diybu/challengeData16/c1beacon/wholeRecords/vcfSize500/chr10_selectedInd.vcf"

   "data/diybu/challengeData16/c1beacon/wholeRecords/500notInBeacon/chr10_selectedInd.vcf"

The following files are intermediate and final output files:

   "vcfSize500wholeRecords.pkl"

   "500notInBeaconwholeRecords.pkl"

   "chr10returnValue_500wholeR.txt"

   "chr10_lrtScoreByStep.txt"

   "power_delta.pdf"

From the resulting power_delta.pdf, we can see the privacy risk has been mitigated substantially. Even for a very rare case where the log-likelihood scores are very similiar for the case and test, only 5% SNPs need to be filpped to mitigate the power dramatically.

Copyright 2016 Zhiyu Wan

HIPLAB, Department of Biomedical Informatics, School of Medicine

Vanderbilt University

9/14/2016
