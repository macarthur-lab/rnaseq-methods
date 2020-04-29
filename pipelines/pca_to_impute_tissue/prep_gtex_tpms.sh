cat GTEx_v7_Annotations_SampleAttributesDS.txt | cut -f1,6,7 > tissue_types.txt
tail -n +4 GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct > gtex_gene_tpm.gct 
awk 'NR == 3' GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct > gtex_ids
for line in $(cat gtex_ids); do echo ${line} >> gtex_ids_transposed; done
for id in $(cat gtex_ids_transposed); do TISSUE_TYPE=$(grep ${id} tissue_types_sorted.txt | cut -f2) ; echo -e ${id} '\t' ${TISSUE_TYPE} >> tissue_types_matched.txt; done

for id in $(cat gtex_ids_transposed) \
do \
TISSUE_TYPE=$(grep ${id} tissue_types.txt | cut -f2) \
TISSUE_SUB_TYPE=$(grep ${id} tissue_types.txt | cut -f3) \
echo -e ${id} '\t' ${TISSUE_TYPE} '\t' ${TISSUE_SUB_TYPE} >> tpms_tissue_types_matched.tsv \
done

for id in $(cat gtex_ids_transposed); do TISSUE_TYPE=$(grep ${id} tissue_types.txt | cut -f2); TISSUE_SUB_TYPE=$(grep ${id} tissue_types.txt | cut -f3); echo -e ${id} '\t' ${TISSUE_TYPE} '\t' ${TISSUE_SUB_TYPE} >> tissue_types_matched.tsv; done

for id in $(cat tpms_ids_transposed); do TISSUE_TYPE=$(grep ${id} tissue_types.txt | cut -f2); TISSUE_SUB_TYPE=$(grep ${id} tissue_types.txt | cut -f3); echo -e ${id} '\t' ${TISSUE_TYPE} '\t' ${TISSUE_SUB_TYPE} >> tpms_tissue_types_matched.tsv; done

