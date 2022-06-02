#### script that reads in imputed data and returns linear regression summary stats
#### need to change "hapcount" datasets after ExtractTracts.py before this 
#### heidi steiner
#### heidiesteiner@email.arizona.edu
#### 2022-02-16 
#### Last Updated: 2022-03-01


#### load python packages 
import argparse, hail as hl, numpy as np

#### initiate hail 
hl.init()
import hail as hl # is this needed here? 
hl.plot.output_notebook()

#### load hail libraries
import hail as hl
from hail.plot import show
from pprint import pprint
from bokeh.io import save


####### SET PATHS HERE ########


#### set path to data relative to where the plot will be made 
data_path = 'vitk_phased_autosomes_lifted.vcf.gz'
covar_filename = 'tractor3a_covariates_pr.txt'
out_path = '2022-03-01_vitk_lmLAadj_pr.tsv'

#### load GWAS data
ds = hl.import_vcf(data_path, reference_genome='GRCh37', force_bgz = True)


#### load covariate data 
table = (hl.import_table(covar_filename, impute=True).key_by('plate_id')) 


####### add or remove rfmix/K3/ here for PR data ########


#### put together genotype information
#### by ancestry 
row_fields={'CHROM': hl.tstr, 'POS': hl.tint, 'ID': hl.tstr, 'REF': hl.tstr, 'ALT': hl.tstr} 
anc0dos = hl.import_matrix_table('../rfmix/K3/phased_autosomes_lifted.anc0.dosage.txt', row_fields=row_fields,row_key=[], min_partitions=32)
anc0dos = anc0dos.key_rows_by().drop('row_id')
anc0dos = anc0dos.key_rows_by(locus=hl.locus(anc0dos.CHROM, anc0dos.POS)) 

anc1dos = hl.import_matrix_table('../rfmix/K3/phased_autosomes_lifted.anc1.dosage.txt', row_fields=row_fields,row_key=[], min_partitions=32)
anc1dos = anc1dos.key_rows_by().drop('row_id')
anc1dos = anc1dos.key_rows_by(locus=hl.locus(anc1dos.CHROM, anc1dos.POS)) 

anc2dos = hl.import_matrix_table('../rfmix/K3/phased_autosomes_lifted.anc2.dosage.txt', row_fields=row_fields,row_key=[], min_partitions=32)
anc2dos = anc2dos.key_rows_by().drop('row_id')
anc2dos = anc2dos.key_rows_by(locus=hl.locus(anc2dos.CHROM, anc2dos.POS)) 

row_fields={'CHROM': hl.tstr, 'POS': hl.tint, 'ID': hl.tstr} 
hapcounts0 = hl.import_matrix_table('../rfmix/K3/anc0.hapcount.hail.txt',  row_fields=row_fields, row_key=[], min_partitions=32) 
hapcounts0 = hapcounts0.key_rows_by().drop('row_id')
hapcounts0 = hapcounts0.key_rows_by(locus=hl.locus(hapcounts0.CHROM, hapcounts0.POS)) 

hapcounts1 = hl.import_matrix_table('../rfmix/K3/anc1.hapcount.hail.txt',  row_fields=row_fields, row_key=[], min_partitions=32) 
hapcounts1 = hapcounts1.key_rows_by().drop('row_id')
hapcounts1 = hapcounts1.key_rows_by(locus=hl.locus(hapcounts1.CHROM, hapcounts1.POS)) 


#### merge ancestries in a matrix table 
mt = ds.annotate_entries(anc0dos = anc0dos[ds.locus, ds.s], anc1dos = anc1dos[ds.locus, ds.s],anc2dos = anc2dos[ds.locus, ds.s], hapcounts0 = hapcounts0[ds.locus, ds.s], hapcounts1 = hapcounts1[ds.locus, ds.s])

#### merge genotypes and phenotypes
mt = mt.annotate_cols(pheno = table[mt.s])

#### write out mt to speed up the regressions later
mt.write('phenos_genos_dosages.mt', overwrite = True)


#### read that mt back in 
mt = hl.read_matrix_table('phenos_genos_dosages.mt')


#### stop here and comment out if doing PR data #### 
#### maybe write and if/else? 


#### remove outlier ?? 
mt = mt.filter_cols(mt.s != '0_0_WARFER034_1767-JK_Karnes_MEGA_Plate_01_E05_WARFER034') # dose outlier
mt = mt.filter_entries(mt.s != '0_0_WARFER034_1767-JK_Karnes_MEGA_Plate_01_E05_WARFER034') # dose outlier


#### QC steps 
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 1e-6)
mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.01)



#### GWAS 
mt = mt.annotate_rows(unadj = hl.agg.linreg(mt.pheno.dose,[1.0, mt.hapcounts0.x, mt.anc0dos.x, mt.hapcounts1.x, mt.anc1dos.x, mt.anc2dos.x]))

#### remove later if find bug? 
#### write out mt to speed up  later
mt.write('phenos_genos_dosages_qc_vitk_lmunadj.mt', overwrite = True)

#### read that mt back in 
mt = hl.read_matrix_table('phenos_genos_dosages_qc_vitk_lmunadj.mt')

#### export steps 
results = mt.rows() 

results = results.select(**results.unadj)

results.write('results.ht', overwrite = True)

results = hl.read_table('results.ht')
results = results.key_by() 

results = results.select( 
    variant=hl.variant_str(results.locus, results.alleles), 
    **{field: results[field] 
       for field in results.row if field not in ('locus', 'alleles')} 
) 


results = results.key_by('variant') 
results.export(out_path)

