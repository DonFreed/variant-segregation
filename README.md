# variant-segregation
Add a variant segregation metric to a VCF file.

### Usage

```
python3 vcf_annotate_segregation.py [options] family.ped

positional arguments:
  family.ped            A file containing pedigree information

optional arguments:
  -h, --help            show this help message and exit
  --infile INFILE
  --infile_gzipped
  --outfile OUTFILE
  --gzip_outfile
  --mendelian_inconsistency MENDELIAN_INCONSISTENCY
                        The probability of encountering a mendelian
                        inconsistency (0.0 - 0.1)
  --debug
  --min_geno_qual MIN_GENO_QUAL
                        The minimum quality genotype to include in the metric
  --min_num_informative MIN_NUM_INFORMATIVE
                        The minimum number of informative segregation events
                        to output a score
```

### Method

Adds the phred-scaled probability of the observed variant segregation to the INFO field of a VCF file. Probabilities are normalized to the maximum likelihood so that a phred-scaled value of '10' indicates that the observed segregation was 1/10th as likely as the most likely segregation of variants.  

Observed and maximum likelihood segregations are calculated using the multinomial distribution with mendelian inconsistencies are incorporated as rare events. The currently implementation only works at biallelic sites.
