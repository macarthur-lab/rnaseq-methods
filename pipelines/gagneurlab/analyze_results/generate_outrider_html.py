#%%
import glob
import math
import os
import pandas as pd

from annotations.get_omim_table import get_omim_table
from annotations.get_hgnc_table import get_hgnc_table
from annotations.get_panel_app_table import get_panel_app_table
from annotations.get_clingen_table import get_clingen_dosage_sensitivity_table, get_clingen_gene_disease_validity_table

from sample_metadata.rnaseq_metadata_utils import get_seqr_info_and_other_metadata_df

#%%
BASE_DIR = os.path.expanduser("~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/OUTRIDER_results3")

os.chdir(BASE_DIR)

#%%

hgnc_df = get_hgnc_table()
hgnc_df = hgnc_df[[
    'Ensembl gene ID',
    'HGNC ID',
    'RefSeq IDs',
    'Pubmed IDs',
    'Locus type',
    # 'Mouse genome database ID',
    # 'Approved symbol',
    # 'Approved name',
    # 'Status',
    # 'Previous symbols',
    # 'Alias symbols',
    # 'Chromosome',
    # 'Accession numbers',
]]

omim_df = get_omim_table()
omim_df = omim_df[[
    'gene_id',
    'phenotype_description',
    'gene_description',
    'mim_number',
    'phenotype_mim_number',
    'phenotypic_series_number',
    'phenotype_inheritance',
    'pLI',
    'mis_z',
    'oe_lof_upper',
    'text',
    'comments',
    # 'chrom',
    # 'start',
    # 'end',
    # 'gene_symbols',
    # 'date_created',
    # 'date_updated',
    # 'mouse_gene_id',
]]

panel_app = get_panel_app_table()
panel_app = panel_app[[
    'hgnc',
    #'gene name',
    #'biotype',
    'gene id',
    'confidence',
    'penetrance',
    'mode_of_pathogenicity',
    'mode_of_inheritance',
    'publications',
    'evidence',
    'phenotypes',
    'panel name',
]]

clingen_dosage_sensitivity_df = get_clingen_dosage_sensitivity_table()
clingen_dosage_sensitivity_df = clingen_dosage_sensitivity_df[[
    #'GENE SYMBOL',
    'HGNC ID',
    'HAPLOINSUFFICIENCY',
    'TRIPLOSENSITIVITY',
    #'ONLINE REPORT',
    #'DATE',
]]

clingen_gene_disease_validity_df = get_clingen_gene_disease_validity_table()
clingen_gene_disease_validity_df = clingen_gene_disease_validity_df[[
    #'GENE SYMBOL',
    'GENE ID (HGNC)',
    'DISEASE LABEL',
    'DISEASE ID (MONDO)',
    'MOI',
    'SOP',
    'CLASSIFICATION',
    #'ONLINE REPORT',
    #'CLASSIFICATION DATE',
    'GCEP',
]]


hgnc_df.drop_duplicates(subset="Ensembl gene ID", inplace=True)  # TODO combine rows with duplicate gene ids
omim_df.drop_duplicates(subset="gene_id", inplace=True)          # TODO combine rows with duplicate gene ids
panel_app.drop_duplicates(subset="hgnc", inplace=True)
clingen_dosage_sensitivity_df.drop_duplicates(subset="HGNC ID", inplace=True)
clingen_gene_disease_validity_df.drop_duplicates(subset="GENE ID (HGNC)", inplace=True)


annotation_table = hgnc_df
annotation_table = pd.merge(annotation_table, omim_df, left_on="Ensembl gene ID", right_on="gene_id", how="outer")
annotation_table = pd.merge(annotation_table, panel_app, left_on="HGNC ID", right_on="hgnc", how="outer")
annotation_table = pd.merge(annotation_table, clingen_dosage_sensitivity_df, left_on="HGNC ID", right_on="HGNC ID", how="outer")
annotation_table = pd.merge(annotation_table, clingen_gene_disease_validity_df, left_on="HGNC ID", right_on="GENE ID (HGNC)", how="outer")

annotation_table = annotation_table[[
    'HGNC ID',
    'hgnc',
    'GENE ID (HGNC)',
    'gene id',
    'Ensembl gene ID',
    'gene_id',
    'RefSeq IDs',
    'Pubmed IDs',
    'Locus type',
    'phenotype_description',
    'gene_description',
    'mim_number',
    'phenotype_mim_number',
    'phenotypic_series_number',
    'phenotype_inheritance',
    'pLI',
    'mis_z',
    'oe_lof_upper',
    'text',
    'comments',
    'confidence',
    'penetrance',
    'mode_of_pathogenicity',
    'mode_of_inheritance',
    'publications',
    'evidence',
    'phenotypes',
    'panel name',
    'HAPLOINSUFFICIENCY',
    'TRIPLOSENSITIVITY',
    'DISEASE LABEL',
    'DISEASE ID (MONDO)',
    'MOI',
    'SOP',
    'CLASSIFICATION',
    'GCEP',
]]

#%%

meta_df = get_seqr_info_and_other_metadata_df()
meta_df['total reads x 10^6 (rnaseqc)'] = meta_df['total reads x 10^6 (rnaseqc)'].astype("float").astype("int").astype("str")

meta_df = meta_df[[
    'sample_id',
    'batch_date_from_hg19_bam_header',
    'star_pipeline_batch',
    #'sex',
    'imputed sex',
    'imputed tissue',
    #'analysis batch',
    'qc notes',
    'stranded? (rnaseqc)',
    'read length (rnaseqc)',
    'total reads x 10^6 (rnaseqc)',
    'mapping rate (rnaseqc)',
    'indiv (seqr)',
    'proj (seqr)',
    'fam (seqr)',
    'proj2 (seqr)',
    'fam2 (seqr)',
    'genome (seqr)',
    'population (seqr)',
    #'sample id (seqr)',
    'sample type (seqr)',
    'analysis status (seqr)',
    'variant notes (seqr)',
    'variant tags (seqr)',
    'coded phenotype (seqr)',
    'analysis summary + notes (seqr)',
    'internal case review notes (seqr)',
    #'cram path (seqr)',
    #'Sample ID_x',
    'Alias (Beryl:Supp.)',
    'Clinical Diagnosis (Beryl:Supp.)',
    'Sex (Beryl:Supp.)',
    'Age at muscle biopsy (Beryl:Supp.)',
    'Site of biopsy (Beryl:Supp.)',
    'Notes (Beryl:Supp.)',
    #'Sample ID_y',
    '%Contamin  RNAseq (Beryl:Probands)',
    'Age at  Biopsy (Beryl:Probands)',
    'Biopsy type (Beryl:Probands)',
    'Candidate  Variants (Beryl:Probands)',
    'CanditateGenes (culprit,if solved) (Beryl:Probands)',
    'Collab PI (Beryl:Probands)',
    'Data details (Beryl:Probands)',
    'Data_type (Beryl:Probands)',
    'Ethnicity (Beryl:Probands)',
    'Include in manuscript? (Beryl:Probands)',
    'Phenotype (Beryl:Probands)',
    'Short Phenotype (Beryl:Probands)',
    'Phenotype comments (Beryl:Probands)',
    'Other comments (Beryl:Probands)',
    'Variant  type (Beryl:Probands)',
    'Genetic diagnosis Status (Beryl:Probands)',
    'Variant consequence (Beryl:Probands)',
    'Phenotype (Beryl:Seqr-data)',
    'Look at (Beryl:Seqr-data)',
    'Notes (Beryl:Seqr-data)',
    'RNA tissue (definitive source) (Beryl:Seqr-data)',
]]

meta_df = meta_df.rename(columns={
    "star_pipeline_batch": "TGG pipeline batch",
    "batch_date_from_hg19_bam_header": "sequencing date",
})


#%%

def read_OUTRIDER_results_table(path):
    results_df = pd.read_table(path)

    results_df = results_df[[
        "sampleID",
        "geneID",
        "symbol",
        "biotype",
        "pValue",
        "padjust",
        "zScore",
        "rawcounts",
        "normcounts",
    ]]

    def round_to_sigdig(series, n=2):
        """Copied from https://stackoverflow.com/questions/57590046/round-a-series-to-n-number-of-significant-figures"""
        return series.apply(lambda x: round(x, n - int(math.floor(math.log10(abs(x))))))

    results_df["pValue"] = round_to_sigdig(results_df["pValue"])
    results_df["padjust"] = round_to_sigdig(results_df["padjust"])

    results_df = pd.merge(results_df, annotation_table, left_on="geneID", right_on="gene_id", how="left")
    results_df = results_df.fillna(" ")

    return results_df


#%%

"""
OUTRIDER results table example:

$1        sampleID : 204H_AM_M1
$2          geneID : ENSG00000260260
$3          symbol : SNHG19
$4         biotype : lincRNA
$5          pValue : 2.78716012098559e-07
$6         padjust : 0.00964888158420464
$7          zScore : 4.25
$8       rawcounts : 1669
$9      normcounts : 1130.87
$10  meanCorrected : 391.02
$11    description : small nucleolar RNA host gene 19 [Source:HGNC Symbol;Acc:HGNC:49574]
$12            chr : 16
$13          start : 2154797
$14            end : 2155358
$15         strand : -1
$16              q : 34
"""


def generate_html(dir_name, results_df, qc_plots, volcano_plots):

    qc_plots_html = ""
    for qc_plot in qc_plots:
        qc_plots_html += f"""
<div class="{'four' if len(qc_plots) <= 4 else 'three'} wide column">
    <a href="{qc_plot}"><img style="max-width: 100%" src="{qc_plot}"></a><br />
</div>"""

    sample_rows = {}
    for sample_id in sorted(set(results_df.sampleID)):
        sample_results_df = results_df[(results_df.sampleID == sample_id) & (results_df.padjust < 0.05)]
        sample_rows[sample_id] = sample_results_df

    sample_rows_sorted = sorted(sample_rows.items(), key=lambda x: len(x[1]), reverse=True)
    sample_links = "<hr />\n"
    for i, (sample_id, _) in enumerate(sample_rows_sorted):
        if i > 0: sample_links += ", &nbsp; "
        sample_links += f"""<a href="#anchor-{sample_id}">{sample_id}</a>\n"""

    samples_html = ""
    for sample_id, sample_results_df in sample_rows_sorted:
        # sample header
        sample_header_html = f"""
<hr />
<a name="anchor-{sample_id}"></a>
<table width="100%">
    <tr>
        <td>{len(sample_results_df)} outliers in <b>{sample_id}</b></td>
        <td style="text-align:right"><a href="#top-of-page">Jump to top</a></td>
    </tr>
</table>
"""
        # volcano plot
        volcano_plot_html = ""
        volcano_plot_paths = [v for v in volcano_plots if v.endswith(f"{sample_id}.png")]
        if len(volcano_plot_paths) > 1:
            raise ValueError(f"{len(volcano_plot_paths)} volcano plots found for {sample_id}: {volcano_plot_paths}")
        elif len(volcano_plot_paths) == 1:
            volcano_plot_path = volcano_plot_paths[0]
            volcano_plot_html = f"""<a href="{volcano_plot_path}"><img style="max-width: 100%" src="{volcano_plot_path}"></a>"""

        # sample info
        sample_meta = meta_df[meta_df.sample_id == sample_id]
        if len(sample_meta) > 1:
            raise ValueError(f"{len(sample_meta)} sample metadata rows found for {sample_id}")
        if len(sample_meta) == 0:
            if not sample_id.startswith("GTEX"):
                print(f"ERROR: '{sample_id}' not found in sample metadata table")
            continue

        sample_info_html = "<table>"
        for key, value in sample_meta.iloc[0].to_dict().items():
            if value is None or not value or not str(value).strip():
                continue
            #key = key.replace("_", "")
            sample_info_html += "<tr>"
            sample_info_html += f"""<td style="white-space: nowrap; padding-right: 10px"><b>{key}</b></td><td>{value}</td>"""
            sample_info_html += "</tr>"
        sample_info_html += "</table>"

        # sample table
        sample_table_html = "<tr>" + "".join([f"""<td>{c}</td>""" for c in sample_results_df.columns]) + "</tr>\n"
        for _, row in sample_results_df.sort_values("padjust").iterrows():
            sample_table_html += "<tr>" + "".join([f"""<td style="white-space: nowrap">{row[c]}</td>""" for c in sample_results_df.columns]) + "</tr>\n"

        # combined html
        samples_html += f"""
<div class="sixteen wide column" style="overflow-x: scroll; padding-top: 50px">
    {sample_header_html}
    
    {volcano_plot_html}
    
    {sample_info_html}
    
    <table class="ui stackable table">
        <tbody>
            {sample_table_html}
        </tbody>
    </table>    
</div>
"""

    html = f"""<html>
<head>
    <meta content="width=device-width, initial-scale=1" charset="utf-8" name="viewport" />
    <title>{dir_name}</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/semantic-ui@2.4.2/dist/semantic.min.css" />
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
    <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/semantic-ui@2.4.2/dist/semantic.min.js"></script>
</head>
<body>
    <div><a name="top-of-page"></a></div> 
    <div class="ui stackable grid">
        <div class="row">
            <div class="one wide column"></div>
            <div class="fourteen wide column">
                <div class="ui stackable grid">
                    <div class="row">
                        {qc_plots_html}
                    </div>
                </div>
                <div class="ui stackable grid">
                    <div class="row">
                        <div class="sixteen wide column">{sample_links}</div>
                    </div>
                    <div class="row">
                        {samples_html}
                    </div>
                </div>
            </div>
            <div class="one wide column"></div>
        </div>
    </div>
</body>
</html>"""

    return html


for dir_name in os.listdir(BASE_DIR):
    if not os.path.isdir(dir_name) or dir_name.startswith("."):
        continue
    print(dir_name)

    # get the OUTRIDER results table
    #results_table = glob.glob(f"{dir_name}/*padj_0.05_results.tsv*")
    results_table = glob.glob(f"{dir_name}/*_all_results.tsv*")
    if not results_table:
        raise ValueError(f"results table not found in {BASE_DIR}/{dir_name}")
    if len(results_table) > 1:
        raise ValueError(f"{len(results_table)} different results tables found in {BASE_DIR}/{dir_name}")
    results_table = results_table[0]
    results_df = read_OUTRIDER_results_table(results_table)

    volcano_plots = []
    qc_plots = []
    for plot_path in glob.glob(f"{dir_name}/*.png"):
        if "_volcano_" in plot_path:
            volcano_plots.append(plot_path)
        else:
            qc_plots.append(plot_path)

    html = generate_html(dir_name, results_df, qc_plots, volcano_plots)

    with open(f"{dir_name}.html", "wt") as f:
        f.write(html)

    print("\n".join(qc_plots))


#%%
#%%
