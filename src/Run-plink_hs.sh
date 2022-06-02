plink=plink1.9 
$plink --bfile /home/steiner/GWAS_MAY2021/qcdir/preimputation/data_snpsonly_clean_chrpos_nodups --exclude Exclude-data_snpsonly_clean_chrpos_nodups-HRC.txt --make-bed --out TEMP1
$plink --bfile TEMP1 --update-map Chromosome-data_snpsonly_clean_chrpos_nodups-HRC.txt --update-chr --make-bed --out TEMP2
$plink --bfile TEMP2 --update-map Position-data_snpsonly_clean_chrpos_nodups-HRC.txt --make-bed --out TEMP3
$plink --bfile TEMP3 --flip Strand-Flip-data_snpsonly_clean_chrpos_nodups-HRC.txt --make-bed --out TEMP4
$plink --bfile TEMP4 --reference-allele Force-Allele1-data_snpsonly_clean_chrpos_nodups-HRC.txt --make-bed --out data_snpsonly_clean_chrpos_nodups-updated
for i in {1..22}; do $plink --bfile data_snpsonly_clean_chrpos_nodups-updated --reference-allele Force-Allele1-data_snpsonly_clean_chrpos_nodups-HRC.txt --recode vcf --chr ${i} --out chr${i};done
rm TEMP*

