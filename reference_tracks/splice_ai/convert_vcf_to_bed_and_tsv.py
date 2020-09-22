import collections
import gzip
from pprint import pprint
import sys
import tqdm

# delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL).
score_header = "ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL".split("|")
ds_header = "DS_AG|DS_AL|DS_DG|DS_DL".split("|")
dp_header = "DP_AG|DP_AL|DP_DG|DP_DL".split("|")

def process_high_scores(chrom, start, end, alleles, allele_score_dicts, strand, symbol):
    """Combine alleles with same start, end, and scores >= 0.2"""

symbol_to_strand_lookup = {}
with open("./annotations/grch38.txt", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip().split("\t")
        symbol = fields[0]
        strand = fields[2]
        symbol_to_strand_lookup[symbol] = strand

score_threshold = 0.5
input_path = sys.argv[1]
output_bed_path = input_path.replace(".vcf", "").replace(".gz", "").replace(".bgz", "") + f".score_{score_threshold}.bed"
output_allele_tsv_path = input_path.replace(".vcf", "").replace(".gz", "").replace(".bgz", "") + f".score_{score_threshold}.alleles.tsv"
output_junction_tsv_path = input_path.replace(".vcf", "").replace(".gz", "").replace(".bgz", "") + f".score_{score_threshold}.junctions.tsv"

junction_tsv_header = ["chrom", "start", "end", "length", "# of alleles", "alleles", "strand", "gene", "donor/acceptor", "loss/gain", "max score", "min score"]
allele_tsv_header = ["chrom", "start", "end", "length", "allele_type", "allele", "strand", "gene", "donor/acceptor", "loss/gain", "score"]

prev_chrom = prev_start = None
allele_score_dict = collections.defaultdict(lambda: collections.defaultdict(list))
with gzip.open(input_path, "rt") as f, \
    open(output_bed_path, "wt") as bed, \
    open(output_allele_tsv_path, "wt") as allele_tsv, \
    open(output_junction_tsv_path, "wt") as junction_tsv:

    allele_tsv.write("\t".join(allele_tsv_header) + "\n")
    junction_tsv.write("\t".join(junction_tsv_header) + "\n")

    for i, line in enumerate(tqdm.tqdm(f, unit=" lines")):
        if line.startswith("#"):
            continue

        fields = line.split()
        chrom = fields[0]
        if len(chrom) > 2:
            continue  # filter contigs that aren't chromosomes

        start = int(fields[1])
        ref = fields[3]
        alt = fields[4]
        score_fields = fields[-1].replace("SpliceAI=", "").split("|")
        alt_allele = score_fields[0]
        assert alt == alt_allele

        symbol = score_fields[1]
        strand = symbol_to_strand_lookup.get(symbol, "?")

        ds_dict = dict(zip(ds_header, map(float, score_fields[2:6])))
        ds_keys_with_high_score = [k for k, v in ds_dict.items() if v >= score_threshold]
        if not ds_keys_with_high_score:
            continue
        dp_dict = dict(zip(dp_header, map(int, score_fields[6:10])))

        #  process (chrom, start, end, strand, symbol, allele_score_dicts) tuples
        if chrom != prev_chrom or start != prev_start:
            #pprint(dict(allele_score_dict))
            for (_chrom, _start, _end, _strand, _symbol), _values_dict in allele_score_dict.items():
                for _ds_key, _values_list in _values_dict.items():
                    _max_score = max([_score for _, _, _score in _values_list])
                    _min_score = min([_score for _, _, _score in _values_list])
                    _alleles = []
                    for (_ref, _alt_allele, _score) in _values_list:
                        _allele_type = ""
                        if len(_ref) == len(_alt_allele):
                            _allele_type = "SNP"
                            _allele_string = _alt_allele
                        elif len(_ref) < len(_alt_allele):
                            _allele_type = "INS"
                            _allele_string = f"ins{_alt_allele[1:]}"
                        elif len(_ref) > len(_alt_allele):
                            _allele_type = "DEL"
                            _allele_string = f"del{_ref[1:]}"
                        _alleles.append(_allele_string)

                        allele_tsv_row = [_chrom, min(_start, _end), max(_start, _end), abs(_end - _start), _allele_type, _allele_string, _strand, _symbol, _ds_key[-2], _ds_key[-1], _score]
                        allele_tsv.write("\t".join(map(str, allele_tsv_row)) + "\n")

                    junction_tsv_row = [_chrom, min(_start, _end), max(_start, _end), abs(_end - _start), len(_alleles), ",".join(_alleles), _strand, _symbol, _ds_key[-2], _ds_key[-1], _max_score, _min_score]
                    junction_tsv.write("\t".join(map(str, junction_tsv_row)) + "\n")

                    _alleles = ",".join(sorted(_alleles, key=lambda a: (len(a), a))) if len(_alleles) <= 3 else f"{len(_alleles)} alleles"

                    _shape = "|" if _ds_key.endswith("G") else "x"
                    _color = "#11AA11" if _ds_key.endswith("G") else "#AA1111"
                    _loss_or_gain = "gain" if _ds_key.endswith("G") else "loss"
                    #_label = f"{_ds_key[-2]}-{_loss_or_gain}:" + _alleles.replace(" ", "_") + f":{_max_score}"
                    _label = _alleles.replace(" ", "_") + f"__{_max_score}__{_ds_key[-2]}-{_loss_or_gain}"
                    _shape_key = "left_shape" if _start < _end else "right_shape"
                    _line_width = float((_max_score - score_threshold)/(1 - score_threshold))*7
                    _gff_tags = f"label={_label};line_width={_line_width:0.1f};color={_color};{_shape_key}={_shape};alleles={_alleles}"
                    bed_row = [_chrom, min(_start, _end) - 1, max(_start, _end), _gff_tags, '0', _strand]
                    bed.write("\t".join(map(str, bed_row)) + "\n")

            #print("clear allele_score_dict")
            allele_score_dict.clear()
            prev_chrom = chrom
            prev_start = start

        for ds_key in ds_keys_with_high_score:
            score = ds_dict[ds_key]
            pos = dp_dict[ds_key.replace("DS", "DP")]
            end = start + pos

            allele_score_dict[(chrom, start, end, strand, symbol)][ds_key].append((ref, alt_allele, score))

        if i > 500000:
            break

    print("Done")
