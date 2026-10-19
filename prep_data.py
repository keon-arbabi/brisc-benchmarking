import os
import gc
import sys
from pathlib import Path
import polars as pl
_HOME_DIR = Path.home()
sys.path.append(f'{_HOME_DIR}/brisc-benchmarking')
from brisc import SingleCell
from utils_local import run

# SEAAD
# https://doi.org/10.21203/rs.3.rs-2921860/v1
# https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html

dir_data = f'{_HOME_DIR}/single-cell/SEAAD'
os.makedirs(dir_data, exist_ok=True)

file_data = f'{dir_data}/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad'
file_metadata = f'{dir_data}/sea-ad_cohort_donor_metadata.xlsx'
if not os.path.exists(file_data):
    run(f'wget https://sea-ad-single-cell-profiling.s3.amazonaws.com/'
        f'MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad '
        f'-O {file_data}')
if not os.path.exists(file_metadata):
    run(f'wget https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.'
        f'divio-media.net/filer_public/b4/c7/b4c727e1-ede1-4c61-b2ee-bf1ae4a3ef68/'
        f'sea-ad_cohort_donor_metadata_072524.xlsx '
        f'-O {file_metadata}')

donor_metadata = pl.read_excel(file_metadata)
cols = ['sample', 'donor', 'cell_type', 'cell_type_broad', 'is_ref', 'cond',
        'apoe4_dosage', 'pmi', 'age_at_death', 'sex']

sc = SingleCell(file_data)
sc = sc\
    .cast_obs({'Donor ID': pl.String})\
    .join_obs(
        donor_metadata.select(['Donor ID'] +
        list(set(donor_metadata.columns).difference(sc.obs.columns))),
        on='Donor ID', validate='m:1')\
    .rename_obs({
        'exp_component_name': 'cell_id', 'Donor ID': 'sample',
        'Subclass': 'cell_type','Class': 'cell_type_broad',
        'Neurotypical reference': 'is_ref',
        'APOE Genotype': 'apoe4_dosage',
        'Continuous Pseudo-progression Score': 'cp_score',
        'Age at Death': 'age_at_death', 'PMI': 'pmi', 'Sex': 'sex'})\
    .with_columns_obs(
        pl.col('cell_type_broad').cast(pl.String).replace_strict({
            'Neuronal: Glutamatergic': 'Excitatory',
            'Neuronal: GABAergic': 'Inhibitory',
            'Non-neuronal and Non-neural': 'Non-neuronal'}),
        pl.when(pl.col('Consensus Clinical Dx (choice=Alzheimers disease)')
                .eq('Checked')).then(pl.lit('AD'))
            .when(pl.col('Consensus Clinical Dx (choice=Control)')
                .eq('Checked')).then(pl.lit('Control'))
            .otherwise(None)
            .alias('cond'),
        pl.col('apoe4_dosage')
            .cast(pl.String).str.count_matches('4')
            .fill_null(strategy='mean').round(),
        pl.col('pmi')
            .cast(pl.String).cast(pl.Float32, strict=False),
        pl.col('sample').alias('donor'),
        pl.col('is_ref').cast(pl.String).eq('True').cast(pl.UInt8))\
    .select_obs(cols)\
    .rename_var({'gene_ids': 'gene_symbol'})\
    .set_var_names('gene_symbol')\
    .drop_var('_index')\
    .drop_obsm(list(sc.obsm))\
    .drop_obsp(list(sc.obsp))\
    .drop_uns(list(sc.uns))\
    .qc_metrics(
        num_counts_column='nCount_RNA',
        num_genes_column='nFeature_RNA',
        mito_fraction_column='percent.mt',
        allow_float=True)

print(sc)
sc.peek_obs()
sc.peek_var()

'''
SingleCell dataset in CSR format with 1,378,211 cells (obs), 36,601 genes (var),
and 7,621,750,342 non-zero entries (X)
    obs: cell_id, sample, donor, cell_type, cell_type_broad, is_ref, cond,
         apoe4_dosage, pmi, age_at_death, sex, nCount_RNA, nFeature_RNA, percent.mt
    var: gene_symbol
    uns: QCed, normalized

column           value
 cell_id          GGTGATTAGGTCACTT-L8TX_210722_0…
 sample           H20.33.034
 donor            H20.33.034
 cell_type        Oligodendrocyte
 cell_type_broad  Non-neuronal
 is_ref           0
 cond             null
 apoe4_dosage     0
 pmi              10.016666
 age_at_death     85.0
 sex              Female
 nCount_RNA       4401
 nFeature_RNA     2300
 percent.mt       0.00045444217

column       value
 gene_symbol  MIR1302-2HG
'''

sc.save(f'{dir_data}/SEAAD_raw.h5ad', overwrite=True)
sc.subsample_obs(n=5000, by_column='cell_type', QC_column=None)\
    .save(f'{dir_data}/SEAAD_raw_5K.h5ad', overwrite=True)

del sc; gc.collect()

# Parse 10M PBMC

# https://www.biorxiv.org/content/10.64898/2025.12.12.693897v1
# https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/

dir_data = f'{_HOME_DIR}/single-cell/Parse'
os.makedirs(dir_data, exist_ok=True)

file_data = f'{dir_data}/Parse_10M_PBMC_cytokines.h5ad'
if not os.path.exists(file_data):
    run(f'wget https://parse-wget.s3.us-west-2.amazonaws.com/10m/'
        f'Parse_10M_PBMC_cytokines.h5ad')

cell_type_broad_map = {
    # T cells
    'CD4 Naive': 'T cell',
    'CD4 Memory': 'T cell',
    'CD8 Naive': 'T cell',
    'CD14 Mono': 'Myeloid',
    'CD16 Mono': 'Myeloid',
    'cDC': 'Myeloid',
    'pDC': 'Myeloid',
    # NK/ILC
    'CD8 Memory': 'T cell',
    'Treg': 'T cell',
    'MAIT': 'T cell',
    'NKT': 'T cell',
    'NK': 'NK/ILC',
    'NK CD56bright':
    'NK/ILC',
    'ILC': 'NK/ILC',
    # B cells
    'B Naive': 'B cell',
    'B Intermediate/Memory': 'B cell',
    'Plasmablast': 'B cell',
    # Progenitors
    'HSPC': 'Progenitor',
}
cols = [
    'sample', 'donor', 'cell_type', 'cell_type_broad',
    'is_ref', 'cond', 'cytokine'
]

sc = SingleCell(file_data)\
    .rename_obs({'_index': 'cell_id'})\
    .cast_obs({'cell_type': pl.String})\
    .with_columns_obs(
        pl.col('cell_type').replace_strict(cell_type_broad_map)
            .alias('cell_type_broad'),
        pl.col('treatment').eq('PBS').cast(pl.UInt8).alias('is_ref'),
        pl.when(pl.col('cytokine').eq('IFN-gamma'))
            .then(pl.lit('IFN-gamma'))
            .when(pl.col('treatment').eq('PBS'))
            .then(pl.lit('PBS'))
            .otherwise(None)
            .alias('cond'))\
    .select_obs(cols)\
    .rename_var({'_index': 'gene_symbol'})\
    .drop_var('n_cells')\
    .qc_metrics(
        num_counts_column='nCount_RNA',
        num_genes_column='nFeature_RNA',
        mito_fraction_column='percent.mt',
        allow_float=True)

print(sc)
sc.peek_obs()
sc.peek_var()

'''
SingleCell dataset in CSR format with 9,697,974 cells (obs), 40,352 genes (var),
and 18,830,591,942 non-zero entries (X)
    obs: cell_id, sample, donor, cell_type, cell_type_broad, is_ref, cond,
         cytokine, nCount_RNA, nFeature_RNA, percent.mt
    var: gene_symbol
    uns: QCed, normalized

column           value
 cell_id          89_103_005__s1
 sample           Donor10_4-1BBL
 donor            Donor10
 cell_type        CD8 Naive
 cell_type_broad  T cell
 is_ref           0
 cond             null
 cytokine         4-1BBL
 nCount_RNA       4700
 nFeature_RNA     2236
 percent.mt       0.011914894

column       value
 gene_symbol  TSPAN6
'''

sc.save(f'{dir_data}/Parse_raw.h5ad', overwrite=True)
sc.subsample_obs(n=200000, by_column='cell_type', QC_column=None)\
    .save(f'{dir_data}/Parse_raw_200K.h5ad', overwrite=True)

del sc; gc.collect()

# PanSci
# https://pmc.ncbi.nlm.nih.gov/articles/PMC11910726/
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247719

dir_data = f'{_HOME_DIR}/single-cell/PanSci'
os.makedirs(dir_data, exist_ok=True)

geo_base = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE247nnn/GSE247719/suppl'
pansci_organs = [
    '01_Kidney', '02_Lung', '03_Heart', '04_Liver', '05_Muscle',
    '06_Stomach', '07_BAT', '08_iWAT', '09_gWAT', '10_Ileum',
    '11_Colon', '12_Jejunum', '13_Duodenum',
]
for organ in pansci_organs:
    fname = f'GSE247719_PanSci_{organ}_adata.h5ad'
    fpath = f'{dir_data}/{fname}'
    if not os.path.exists(fpath):
        url = f'{geo_base}/{fname.replace("_", "%5F").replace(".", "%2E")}'
        run(f'wget {url} -O {fpath}')

cols = ['sample', 'donor', 'cell_type', 'cell_type_broad', 'cell_type_fine',
        'is_ref', 'cond', 'organ', 'sex', 'age_group']

sc = SingleCell(
        f'{dir_data}/GSE247719_PanSci_{pansci_organs[0]}_adata.h5ad')\
    .concat_obs([
        SingleCell(f'{dir_data}/GSE247719_PanSci_{organ}_adata.h5ad')
        for organ in pansci_organs[1:]], flexible=True)\
    .set_var_names('gene_name')\
    .make_var_names_unique(separator='~')\
    .cast_obs({'ID': pl.String})\
    .rename_obs({
        'sample': 'cell_id', 'Main_cell_type': 'cell_type_fine',
        'Lineage': 'cell_type', 'Organ_name': 'organ',
        'Age_group': 'age_group', 'Sex': 'sex'})\
    .with_columns_obs(
        (pl.col('ID') + '_' + pl.col('organ')).alias('sample'),
        pl.col('ID').alias('donor'),
        pl.col('cell_type').alias('cell_type_broad'),
        pl.col('Genotype').eq('WT').cast(pl.UInt8).alias('is_ref'),
        pl.when(pl.col('Genotype').eq('WT') &
                pl.col('age_group').eq('06_months'))
            .then(pl.lit('Young'))
            .when(pl.col('Genotype').eq('WT') &
                  pl.col('age_group').eq('23_months'))
            .then(pl.lit('Aged'))
            .otherwise(None)
            .alias('cond'))\
    .select_obs(cols)\
    .drop_var(['gene_id', 'gene_type'])\
    .qc_metrics(
        num_counts_column='nCount_RNA',
        num_genes_column='nFeature_RNA',
        mito_fraction_column='percent.mt',
        allow_float=True)

print(sc)
sc.peek_obs()
sc.peek_var()

'''
SingleCell dataset in CSR format with 20,317,820 cells (obs), 55,416 genes (var),
and 11,863,129,008 non-zero entries (X)
    obs: cell_id, sample, donor, cell_type, cell_type_broad, is_ref, cond,
         organ, sex, age_group, nCount_RNA,
         nFeature_RNA, percent.mt
    var: gene_name
    uns: QCed, normalized

column           value
 cell_id          20221205_EXP096_01_AACCGATTGC_…
 sample           11_Kidney
 donor            11
 cell_type        Epithelial
 cell_type_broad  Epithelial
 cell_type_fine   Proximal tubule cells
 is_ref           1
 cond             Young
 organ            Kidney
 sex              Male
 age_group        06_months
 nCount_RNA       1560
 nFeature_RNA     813
 percent.mt       0.0025641026

column     value
 gene_name  4933401J01Rik
'''

sc.save(f'{dir_data}/PanSci_raw.h5ad', overwrite=True)
sc.subsample_obs(n=200000, by_column='cell_type', QC_column=None)\
    .save(f'{dir_data}/PanSci_raw_200K.h5ad', overwrite=True)

del sc; gc.collect()
