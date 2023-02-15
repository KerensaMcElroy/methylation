# methylation

2023-02-15

Issues with counts of number of windows not matching up between actual calculated windows for each animal and treatment across all chromosomes and files defining window start and stop positions.

e.g.

cat window_mean_HE3_chrm*_ia11.txt | wc -l
100448333
cat window_positions_chrm*.csv | wc -l
99039749

Note these numbers are the same for each animal and treatment. Examination of the code, counts for each chromosome, and original files (e.g. filtered_meth_HE3_chrm10.dat) revealed that this is due to file window_positions_chrm22.csv being truncated. 

tail window_positions_chrm22.csv
chrm_22.1048566,26217427,26217477
chrm_22.1048567,26217452,26217502
chrm_22.1048568,26217477,26217527
chrm_22.1048569,26217502,26217552
chrm_22.1048570,26217527,26217577
chrm_22.1048571,26217552,26217602
chrm_22.1048572,26217577,26217627
chrm_22.1048573,26217602,26217652
chrm_22.1048574,26217627,26217677
chrm_22.1048575,26217652,26217702

tail filtered_meth_HE2_chrm22.dat 
61433015 12 0 5 0 NA NA 5 2 6 0 16 1 8 1 8 2 8 1 5 0 7 0 16 0 7 0 7 1
61433021 11 1 6 2 9 5 NA NA NA NA 14 3 8 0 8 3 9 0 NA NA 10 2 14 3 NA NA 7 3
61433022 12 1 5 0 NA NA NA NA 6 1 15 0 8 1 8 1 8 1 NA NA 7 2 16 0 7 0 7 1
61433029 9 0 7 1 9 2 NA NA NA NA 12 0 8 0 8 2 9 0 NA NA 10 1 14 0 NA NA 6 1
61433030 11 0 5 0 NA NA 5 0 6 0 16 0 8 1 8 0 8 0 5 0 7 0 16 0 7 0 7 1
61433033 10 0 6 1 9 1 NA NA NA NA 13 0 8 0 8 5 9 0 NA NA 10 0 14 1 NA NA 6 1
61433034 12 1 5 0 NA NA 5 1 6 1 16 0 8 1 8 0 8 0 5 0 6 0 16 0 7 0 7 1
61433039 10 0 5 0 9 0 NA NA NA NA 13 1 8 0 8 3 7 0 NA NA 10 1 14 1 NA NA 6 1
61433040 12 1 5 0 NA NA 5 0 6 1 16 0 8 1 8 2 8 0 5 0 7 0 16 0 7 0 7 1
61433043 10 9 6 6 9 9 NA NA NA NA 13 13 8 8 8 8 9 8 NA NA 10 10 14 13 NA NA 7 7

Will attempt to regenerate file window_positions_chrm22.csv
