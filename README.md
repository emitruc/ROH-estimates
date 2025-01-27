# ROH-estimates

There are different approaches to estimate Runs-Of-Homozygosity (ROHs) from genome sequencing reads mapped to a reference genome (BAM format) or from variant and invariant called sites (GVCF format). In this tutorial for the **Workshop on Population and Speciation Genomics 2025** we will use the output of a software called [ROHan](https://github.com/grenaud/ROHan).

Take your time to answer the _**questions**_ we left for you along the tutorial. We will then discuss them during the final wrap-up.

ROHan is a Bayesian framework to estimate local rates of heterozygosity, infer runs of homozygosity (ROHs) and compute global rates of heterozygosity. 
ROHan can work on modern and ancient samples with signs of ancient DNA damage [Renaud G et al (2019)](doi.org/10.1534/genetics.119.302057).

Alternative approaches: 
[BCFtools/RoH](https://samtools.github.io/bcftools/howtos/roh-calling.html) using also allele frequency to improve ROH inference; 
[Darwindow](https://github.com/mennodejong1986/Darwindow).

## How to run ROHan: understanding the options without actually running it 
(as it takes too long for this lab session)

`ROHan --rohmu 1e-5 -t 12 -o $IND --size 100000 --tstv 1.55 --auto ../autosomes.list $FASTA_FILE $BAMfile/$IND.sorted.dedup.bam`

$IND: the IDs of the sample you are running  
$FASTA_FILE: the reference genome sequence  
$BAMfile: the folder containing the BAM files to be analyzed  

Options to think about:  
--rohmu : heterozygosity threshold to call a ROH  
**What do you think is a good value? How to decide?**

--size : size of the non-overlapping genmic windows where local heterozygosity will be estimated  
**Any idea about what could be a relevant ROH length to analyze?** Below you will find more on this  

--tstv : transition/transversion ratio. This can be estimated on the sample of individuals from the population/species of interest  
**How would you do it? Which software? Which input file?**  

--auto : list of the chromosomes/scaffolds to use  
**Why do we have to exclude sexual chrs and too short chrs/scaffolds?**  

The coalescent times of ROH of a certain length is estimated as _g = 100/(2rL)_, where _g_ is the expected time (in generations) back to the parental common ancestor where the IBD haplotypes forming a ROH coalesce, _r_ is the recombination rate in cM/Mb, and _L_ is the length of the ROH in megabases (lots of references for this, eg. [Kardos et al 2018](doi.org/10.1038/s41559-017-0375-4); [Khan et al. 2021](doi.org/10.1073/pnas.202301811)).

Examples using a (recombination rate of ~2.8 cM/Mb as in butterflies)  
ROH(100kb) -> g = 100/(2 x 5 x 0.1) ~ 180 generations  
ROH(200kb) -> g ~ 90 generations ago  
ROH(500kb) -> g ~ 35 generations ago  
#ROH(1Mb) -> g ~ 18 generations ago  


## Summarizing ROHan output

The full output of ROhan from 30 lizard samples (20 _Podarcis raffonei_ from the small islets of La Canna (LC) and Strobolicchio (ST) and 10 _Podarcis wagleriana_ from the main island Sicily (DB)) is provided (one folder per individual) inside the activity folder (**30_genomic_load_roh/lizard**) on the AWS instance. More details about ROHan output files can be found in [ROHan](https://github.com/grenaud/ROHan) github page. 

Check the pdf files in one of the individual folder.  

**What metrics are they showing?**  

**Can you find the summary text file of the ROHan analysis for each sample?**   

**What is the minimum ROH length we are recoding?**   

**How can we get longer ROHs?**  

We will now focus on HMM point estimates (files with extension .mid.hmmp) and use some simple bash command to merge adjacent ROHs and summarize their length distribution.

Get the genomic windows with Probability of being a ROH > 0.9 for the point estimates of the HMM (.mid.hmmp).  
Every time you see $IND in the command line, it stands for one of the individuals in the lizard sample.

```
zcat $IND.mid.hmmp.gz | grep -v "NA" | awk '$5 < 0.1' > $IND.roh100kb
```

Using the command above, we selected ROHs using the probability cutoff of 0.1 of not being a ROH, meaning 0.9 of being a ROH. This saves us a bit of sanity check of the output.  

**Do you agree with this 9:1 threshold?**  

**What is the meaning of a Bayes Factor 9:1?**  

(After running the whole pipeline the follows using this cutoff, you can come back, change the cutoff and run the pipeline again to check the effect of a different cutoff on your estimates) 

We ran ROHan using a rather small window size (100kb) but we can merge adjacent small regions into larger ones starting from the output in the hmmp file. There is also another ROHAN output file (extension .hmmrohl) summarizing the results. Have a look a that too if you are interested.

To be sure we are not introducing any bias by merging short ROHs into larger ones, we can go back to ROHan and run it with larger window sizes and compare the results with the current run after merging.

Merge adjacent regions with `bedtools merge`. 

First, we need to sort the regions by chr and position. 
```
sort -k1,1 -k2,2n $IND.roh100kb > $IND.roh100kb.sorted
```
Then, we can merge the ROHs
```
bedtools merge -i $IND.roh100kb.sorted > $IND.roh100kb.sorted.merged
```
And check the distribution of ROHs by length
```
awk -F'\t' 'BEGIN {OFS=FS} {print $0, $3-$2}' $IND.roh100kb.sorted.merged | cut -f 4 | sort -n | uniq -c

Execute the command above for a few individuals from the three different populations of lizards. 

```
**How long are the longest ROH in the individual from LC?**  

**And in individuals from the Sicilian wall lizard (DB population)?**  

**What is the cumulative length of ROHs larger than 5Mb in the LC01 individual?**

Try the following to answer the question above:

```
awk -F'\t' 'BEGIN {OFS=FS} {print $0, $3-$2}' LC01/LC01.roh100kb.sorted.merged | cut -f 4 | sort -n | uniq -c | awk '$2 > 5000000  {sum+=$2} END {print sum}'
```

**Can you recollect the patterns we saw earlier today during the lecture?**

Here is a loop to iterate the commands above across individuals.  
To be executed wihtin the `30_genomic_load_roh/lizard` folder.

Example for La Canna population:

```
for IND in ST*
do
cat $IND/$IND.mid.hmmp | grep -v "NA" | awk '$5 < 0.1' | sort -k1,1 -k2,2n > $IND/$IND.roh100kb.sorted
bedtools merge -i $IND/$IND.roh100kb.sorted > $IND/$IND.roh100kb.sorted.merged
rm $IND/$IND.roh100kb.sorted
done
```

## Estimating the inbreeding coefficient F(ROH)

The inbreeding coefficient F(ROH) is calculated as the average fraction of the autosomal genome in ROH across individuals of a population. 95% confidence intervals can be estimated by jackknife resampling among individuals within each population.

To estimate F(ROH), we first need to get the total length of the genome that has been checked for ROHs. 

(ROHan is expected to be robust to low coverage data (that is the risk of calling a ROH because we do not have data to call a SNP) classifying regions as "Unclassified" if the probability of being ROH/nonROH is not enough to take a position.)

We can use the following command on one of the individuals and assume that the genome length is the same for all individuals.

```
zcat LC01/LC01.mid.hmmp.gz | sort -k1,1 -k2,2n > LC01/LC01.mid.hmmp.sorted
bedtools merge -i LC01/LC01.mid.hmmp.sorted | awk '{sum+=$3} END {print sum}'
```
The result should be 1429961685 bp

Calculate F(ROH) using all ROHs (in our estimate they are equal/longer than 100kb).  

```
awk -F'\t' 'BEGIN {OFS=FS} {print $0, $3-$2}' LC01/LC01.roh100kb.sorted.merged | cut -f 4 | sort -n | uniq -c | awk '$2 > 99999 {sum+=$2} END {print sum/1429961685}'
```

**What is the F(ROH) for this individual?**

Make a loop to estimate the F(ROH) per population:

```
for IND in LC*
do
echo $IND
awk -F'\t' 'BEGIN {OFS=FS} {print $0, $3-$2}' $IND/$IND.roh100kb.sorted.merged | cut -f 4 | sort -n | uniq -c | awk '$2 > 99999 {sum+=$2} END {print sum/1429961685}'
done
```

F(ROH) estimates can also be limited to the fraction of ROHs which are longer than a certain length to estimate the most recent inbreeding event

**Can you check the proportion of the genome in ROHs per population for 99999 < ROHs < 500000,  499999 < ROHs < 5 Mb, and ROHs > 5 Mb?**

For example, to get the proportion of the genome in ROHs longer the 100kb and shorter than 500kb in the population of lizard from La Canna:

```
for IND in LC*
do
awk -F'\t' 'BEGIN {OFS=FS} {print $0, $3-$2}' $IND/$IND.roh100kb.sorted.merged | cut -f 4 | sort -n | uniq -c | awk '$2 > 99999 && $2 < 500000 {sum+=$2} END {print sum/ 1429961685}'
done
```

We can use these results to calculate the mean F(ROH) per population.

**What is the mean F(ROH) for the population from ST?**


## EXTRA 1: 95% confidence intervals

It would be better to add 95% confidence intervals for these estimates by jackknife resampling.

**What are the 95% confidence intervals of total F(ROH) > 5Mb in ST?**

Make one list per population (use the prefix of each population)
```
for IND in ST*
do
echo $IND >> ST_list
done
```

Resample each population by leaving-one-out approach (100 times here, but 10000 when it matters)
```
echo 'ST' >> ST_FROH_size5Mb_95perc

for i in {1..100}
do
cat ST_list | shuf | tail -n +2 > ST_list_shuffled

while read IND
do
printf $IND >> F_ROH_size5Mb_shuffled
printf '\t'>> F_ROH_size5Mb_shuffled
awk -F'\t' 'BEGIN {OFS=FS} {print $0, $3-$2}' $IND/$IND.roh100kb.sorted.merged | cut -f 4 | sort -n | uniq -c | awk '$2> 5000000 {sum+=$2} END {print sum/1429961685}' >> F_ROH_size5Mb_shuffled
done < ST_list_shuffled

awk '{sum+=$2} END {print sum/NR}' F_ROH_size5Mb_shuffled >> ST_FROH_size5Mb_95perc
rm F_ROH_size5Mb_shuffled
done
```

Get the mean of the distribution

grep -v "ST" ST_FROH_size5Mb_95perc | awk '{sum += $0} END{print sum/NR}'

Get the 95% confidence intervals as the 0.025 and 0.975 percentiles of the distribution.

grep -v "ST" ST_FROH_size5Mb_95perc | sort -n | awk '{all[NR]= $0} END{print all[int(NR*0.025)]}'

grep -v "ST" ST_FROH_size5Mb_95perc | sort -n | awk '{all[NR]= $0} END{print all[int(NR*0.975)]}'

## EXTRA 2: Estimating Effective Population Size from ROH (same references as above).

If the estimated coalescent times for ROH are unbiased, then the average FROH based on ROH with estimated maximum coalescent times less than t generations back in time is an estimator of the inbreeding accumulated in the population from the time of sampling back to t generations ago. Ne is estimated from the following expression of mean expected individual inbreeding as a function of Ne over t generations:

F(ROH),t = 1-(1-1/2Ne)^t
 
This approach assumes that the probability of inferring an ROH segment due to several smaller homozygous segments is very low. 

**What was the population size when the most recent inbreeding event occurred for our populations?**  

Do not forget the confidence intervals!

## EXTRA 3 - The invasive fish quiz!

In the folder `30_genomic_load_roh/fish` you can find 20 samples which have been already analyzed with ROHan.  
This is an invasive species of which we collected data from the source range and three populations in the invasive range.

We know the invasion happened just a few generations in the past and we expect there could have been quite an intense bottlenck at introduction.

**Which one is the population from the native range and which are the invasive populations?**

