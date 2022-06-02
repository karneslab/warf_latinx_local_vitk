# LOCAL ANCESTRY-INFORMED CANDIDATE PATHWAY ANALYSIS OF WARFARIN STABLE DOSE IN LATINO POPULATIONS

## Preface

This document will walk you through all the code used from raw genetic
data to analysis but will not entirely cover how helper documents were
made or the exact code of the analysis.

**NOTE:** In the case of the replication data, this is not the exact
code that was run because some file locations are different.

## Quality Filtering

We first removed Non-SNPs from the dataset for ease of downstream
analyses.

    plink2 --bfile data --snps-only 'just-acgt' --chr 1-22 --recode vcf --out data_snpsonly

**2 038 233** variants were loaded from rawdata and **1 908 435**
variants remain.

6.37% of variants were removed with –snps-only.

–snps-only excludes all variants with one or more multi-character allele
codes. With ‘just-acgt’, variants with single-character allele codes
outside of {‘A’, ‘C’, ‘G’, ‘T’, ‘a’, ‘c’, ‘g’, ‘t’, <missing code>} are
also excluded.

We next included some quality filtering thresholds to remove *really*
ugly data prior to imputation.

`fail_IDs.txt` was made by playing around in plinkQC in R
(src/QC\_preimputation.R), rather than employing QC via PLINK alone.
Hapmap samples (n = 6) and duplicates (and one triplicate, n = 22) were
included in this list along with sex check failures (n=3), one
individual per pair with cryptic relatedness (n=2), and those with
&gt;10% missingness (n=7).

    plink1.9 --bfile data_snpsonly --geno 0.1 --hwe 1e-10 --remove fail_IDs.txt --maf 0.0000001 --make-bed --out data_snpsonly_clean

**1 028 367** variants and 142 people pass filters and QC at R2 &lt; 0.2

–geno filters out all variants with missing call rates exceeding the
provided value (default 0.1) to be removed –hwe filters out all variants
which have Hardy-Weinberg equilibrium exact test p-value below the
provided threshold. –maf filters out all variants with minor allele
frequency below the provided threshold (default 0.01)

46.11% of variant sites have been removed and 11.25% of the samples were
removed.

## Imputation Pre-Processing

The TOPMed imputation server has [some
documentation](https://topmedimpute.readthedocs.io/en/latest/prepare-your-data.html)
about data preprocessing. I didn’t do a good job following those steps.

### checkVCF

Then do [checkVCF.py](https://github.com/zhanxw/checkVCF) because TOPMed
wont run without it (or something similar).

I did the below so that it was easier to remove the SNPs that
checkVCF.py was calling. There are other ways to do this.

First, I changed the names of the files I was working with so that they
wouldn’t contain a1, a2 in the SNP ID column

    bcftools annotate -Ov -x 'ID' -I +'%CHROM:%POS' data_snpsonly_clean.vcf > data_snpsonly_clean_chrpos.vcf

Next, I ran the checkvcf script

    python2 /home/steiner/GWAS_MAY2021/data/preimpute/checkVCF.py -r /home/steiner/GWAS_MAY2021/data/preimpute/hs37d5.fa -o test data_snpsonly_clean_chrpos.vcf; done 

I wrote the dups to a new file

    awk '{print $2}' test.check.dup > test.check.dup2

and excluded the duplicates this way

    plink2 --bfile data_snpsonly_clean_chrpos --exclude test.check.dup2 --make-bed --out data_snpsonly_clean_chrpos_nodups

**987613** variants remain (0.0396298% were removed)

### Will Rayner Toolbox

Next, I ran [this tool](https://www.well.ox.ac.uk/~wrayner/tools/)

    # first get a frequency file from plink
    plink1.9 --bfile /home/shared/PR/qcdir/data_snpsonly_clean_chrpos_nodups --freq --out pr

    perl /home/shared/TOPMed/HRC-1000G-check-bim.pl -b /home/shared/PR/qcdir/data_snpsonly_clean_chrpos_nodups.bim -f pr.frq -r PASS.Variantsbravo-dbsnp-all.tab.gz -h

0.1219992% were removed here

### Re-organize

Next, I created new bash script to 1) call plink1.9, 2) ouput vcfs, and
3) include first file’s location

    bash Run-plink_hs.sh

Finally, zip up the chromosomes to import to TOPMed.

    for i in {1..22};do bcftools sort chr${i}.vcf -Oz -o chr${i}.vcf.gz ;done

## Imputation

### TOPMed

upload to
[TOPMed](%5Bhttps://imputation.biodatacatalyst.nhlbi.nih.gov/#!run/imputationserver%401.5.7)
Reference Panel: TOPMed r2 Array Build: GRCh37 rsq Filter: 0.2 Phasing:
Eagle v2.4 QC Frequency Check: Skip Mode: QC & Imputation

download the data with the TOPMed provided curl command you’ll ALWAYS
need the password that was emailed to you to unzip the files so **don’t
delete the email**.

**837 042** variants and 142 people pass QC filters.

    for i in {1..22}; do unzip chr_${i}.zip;done

### Filtering

Prior to concatenating the dosage files into an `autosomes.bcf` file,
use other R2 cutoffs to filter the data at varying quality:

    for i in {1..22};do bcftools view -i 'R2>.9' -Oz chr${i}.80.dose.vcf.gz > chr${i}.90.dose.vcf.gz; done

Then concatenate all the chromosomes to one file of autosomes

    bcftools concat chr{1..22}.dose.vcf.gz -Ou -o ../autosomes.bcf

move out of the `chr/` directory

    # convert bcf to bim/bed/fam
    plink2 --bcf autosomes.bcf --make-bed --out autosomes 

    # convert bfiles to ped/map/fam
    plink1.9 --bfile autosomes --recode --out autosomes

**18 614 232** variants and 142 people pass QC filters with R2 &lt; 0.2

**15 468 768** variants and 142 people pass QC filters with R2 &lt; 0.1

**9 819 921** variants and 142 people pass QC filters with R2 &lt; 0.01

### LiftOver

Then, `LiftOver` back to GRCh37 with
[liftOverPlink](https://github.com/sritchie73/liftOverPlink). I had to
change the pathname of the LiftOver program to `liftover` from
`LiftOver` with `nano`.

    # first liftover
    for f in 80 90 99; do python2 liftOverPlink.py --map autosomes.$f.map --out lifted --chain hg38ToHg19.over.chain.gz; done 

    # find bad lifts
    for f in 80 90 99; do python2 rmBadLifts.py --map lifted.$f.map --out  good_lifted.$f.map --log bad_lifted.$f.dat; done 

    # clean up the liftover
    for f in 80 90 99; do cut -f 2 bad_lifted.$f.dat > to_exclude.$f.dat; done 
    for f in 80 90 99; do cut -f 4 lifted.$f.bed.unlifted | sed "/^#/d" >> to_exclude.$f.dat ; done 

    # Note: this will clobber the lifted MAP file generated by `liftOverPlink`:
    for f in 80 90 99; do plink1.9 --file autosomes.$f --recode --out lifted.$f --exclude to_exclude.$f.dat  ; done 
    for f in 80 90 99; do plink1.9 --ped lifted.$f.ped --map good_lifted.$f.map --recode --out final.$f; done 

Use plink to create one VCF/chromosome

    for f in 80 90 99; do for i in {1..22}; do plink1.9 --file final.$f --recode vcf bgz --chr ${i} --out chr${i}_lifted.$f; done; done

Then, zip and use vcf tools to index the new VCFs

    for f in 80 90 99; do for i in {1..22}; do bcftools index chr${i}_lifted.$f.vcf.gz; done; done

### Phasing

Re-phase with
[Eagle\_v2.4.1](https://alkesgroup.broadinstitute.org/Eagle/) because
that’s what RFMix and ExtractTracts.py need.

    for f in 80 90 99; do for i in {1..22}; do eagle --vcfTarget chr${i}_lifted.$f.vcf.gz --vcfRef /home/karnes/general/1000g_phase3_nomulti_allelic/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz --geneticMapFile /home/karnes/general/genetic_map_hg19_withX.txt --chrom ${i} --numThreads 8 --outPrefix phased_chr${i}_lifted.$f; done; done

Create a phased, lifted file of concatenated chromosomes.

    for f in 80 90 99; do bcftools concat phased_chr{1..22}_lifted.$f.vcf.gz -Oz -o phased_autosomes_lifted.$f.vcf.gz; done

double check your file is phased by looking for [`|` between
alleles](https://www.biostars.org/p/82671/) for each variant in the
`.vcf` file

## Local Ancestry Inference

Local ancestry inference with
[RFmix\_v2](https://github.com/slowkoni/rfmix)

### RFMix

    for f in 80 90 99; do for i in {1..22}; do rfmix -f /home/steiner/GWAS_MAY2021/data/imputed/topmed/chr/phased_chr${i}_lifted.$f.vcf.gz -r /home/karnes/general/HGDP_1000G_Merge_Results/chr/chr${i}.vcf.gz --chromosome=${i} -m ../ref_pops.txt -g /home/shared/reference_data/b37_map_files/b37_map_files/b37.gmap -e 1 -o chr${i}.$f;done;done

Created `chr/chr*` files with
`bash for i in {1..22};do plink1.9 --bfile hgdp1000ghg19 --recode vcf bgz --chr ${i} --out chr${i}`
from the end of the [Merging HGDP and 1000 Genomes
Pipeline](https://nbviewer.org/github/tomszar/HGDP_1000G_Merge/blob/master/Code/2018-05-MergeGenotypes.ipynb)
(performed by John Feng in the Karnes Lab in 2021) and phased with
`bash eagle --vcfTarget chr${i}.vcf.gz --vcfRef /home/karnes/general/1000g_phase3_nomulti_allelic/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz --geneticMapFile /home/karnes/general/genetic_map_hg19_withX.txt --chrom ${i} --numThreads 8 --outPrefix phased_chr${i}`

### Painted Karyograms

as per [Dr.
Atkinson](https://github.com/eatkinson/ancestry_pipeline/blob/master/plot_painted_karyograms.sh):

    for f in 80 90 99; do for i in {1..22}; do awk '{ if ($1!~"#") print $0 }' chr${i}.$f.msp.tsv | awk '{ if ($1!~"#") print $0 }'| cut -f2,4 > chr${i}.$f.rfmix.snp_loc; done; done

    for f in 80 90 99; do for i in {1..22}; do awk '{ if ($1!~"#") print $0 }' chr${i}.$f.msp.tsv | awk '{ if ($1!~"#") print $0 }' > chr${i}.$f.msp1.tsv; done ; done

    for f in 80 90 99; do for i in {1..22}; do cut -f7- chr${i}.$f.msp1.tsv > chr${i}.$f.Viterbi.tsv; done; done

via
[armartin/ancestry\_pipeline](https://github.com/armartin/ancestry_pipeline)

`mkdir plots` before this step but don’t go inside of it

    for f in 80 90 99; do for i in {1..22}; do cat ../study.sample | while read line; do python2 /home/karnes/general/ancestry_pipeline/collapse_ancestry1.py --rfmix chr${i}.$f.Viterbi.tsv --snp_locations chr${i}.$f.rfmix.snp_loc --ind $line --ind_info ../study.sample --pop_labels "IBS,NAT,YRI" --out plots/$line ;done;done;done

`cd plots` before this step

    for f in 80 90 99; do printf '%s\n' *.$f.A.bed > bed_list_a.$f.txt; done

    for f in 80 90 99; do printf '%s\n' *.$f.B.bed > bed_list_b.$f.txt; done

    for f in 80 90 99; do paste -d ' ' bed_list_a.$f.txt bed_list_b.$f.txt > bed_list.$f.txt;done

plot karyograms

    cat ../study.sample | while read line; do \
    python2 /home/karnes/general/ancestry_pipeline/plot_karyogram1.py \
    --bed_a plots/$line.A.bed \
    --bed_b plots/$line.B.bed \
    --ind $line \
    --centromeres /home/karnes/general/ancestry_pipeline/centromeres_hg19.bed \
    --pop_order IBS,NAT \
    --out plots/$line.png ;done

### Global Ancestry Estimates

and obtain global ancestry estimates

    for POP in IBS,NAT; do python2 /home/karnes/general/ancestry_pipeline/lai_global_f.py \
    --bed_list bed_list.txt \
    --ind_list ../../study.sample \
    --pops IBS,NAT \
    --out lai_global.txt; done

The lai\_global\_f.py script was updated by John Feng to fit the bed
files I produced above.

### Optional: Admixture

double check with
[admixture](http://dalexander.github.io/admixture/admixture-manual.pdf)

    admixture ibs_nat_warf_merge.bed 2 --supervised  --j4

## Tractor

### Extract Tracts

combine MSP files prior to ExtractTracts.py (9/11/21 email from
Elizabeth Atkinson)

    for f in 80 90 99; do head -2 chr1.$f.msp.tsv > autosomes.$f.msp.tsv; done
    for f in 80 90 99; do tail -n +3 -q chr*.$f.msp.tsv >> autosomes.$f.msp.tsv; done

as per [Dr. Atkinson](https://github.com/eatkinson/Tractor), again

    for f in 80 90 99; do python3 /home/steiner/GWAS_MAY2021/Tractor/ExtractTracts.py --msp autosomes --vcf /home/steiner/GWAS_MAY2021/data/imputed/topmed/phased_autosomes_lifted --zipped --num-ancs 3; done

    mv /home/steiner/GWAS_MAY2021/data/imputed/topmed/*anc* /home/steiner/GWAS_MAY2021/Tractor/runs/ibs_nat_yri/

### Hail

then do GWAS in hail with python script `hail.py` as per the rest of the
materials in the github! There is an example python script, jupyter
notebook and a wiki that describes the process.

**STOP**: Need to “fix” the hapcount datasets that come out of
ExtractTracts.py prior to running hail

    for i in {0..2}; do cut -f4,5 --complement phased_autosomes_lifted.anc${i}.hapcount.txt > anc${i}.hapcount.hail.txt; done

**STOP**: Need to filter vcf for vitamin K gene variants prior to
running hail

    plink1.9 --vcf /home/steiner/GWAS_MAY2021/data/imputed/topmed/phased_autosomes_lifted.vcf.gz --extract range /home/steiner/GWAS_MAY2021/candidate_genes/vitkgenes_plink.txt --recode vcf-iid bgz --double-id --out /home/steiner/GWAS_MAY2021/data/imputed/topmed/vitk_phased_autosomes_lifted

The non-local ancestry adjusted regressions can be run with `lm.R` in R
