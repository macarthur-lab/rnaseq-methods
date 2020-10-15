import argparse
import collections
import gzip
from pprint import pprint
import sys
import tqdm

# delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL).
SCORE_HEADER = "ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL".split("|")
DS_HEADER = "DS_AG|DS_AL|DS_DG|DS_DL".split("|")
DP_HEADER = "DP_AG|DP_AL|DP_DG|DP_DL".split("|")

MIN_ALLOWED_SCORE_THRESHOLD = 0.2

SYMBOL_TO_STRAND_LOOKUP = {}
with open("./annotations/grch38.txt", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip().split("\t")
        _symbol = fields[0]
        _strand = fields[2]

        assert _strand in "+-"
        SYMBOL_TO_STRAND_LOOKUP[_symbol] = _strand

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--test", action="store_true", help="output only the 1st part of the vcf")
    p.add_argument("-s", "--score-threshold", type=float, help="min score", default=0.5)
    gain_or_loss = p.add_mutually_exclusive_group(required=True)
    gain_or_loss.add_argument("--gain", action="store_true", help="output only splice donor/acceptor gain events")
    gain_or_loss.add_argument("--loss", action="store_true", help="output only splice donor/acceptor loss events")
    p.add_argument("input_vcf", help="Illumina precomputed SpliceAI score vcf")
    args = p.parse_args()

    if args.score_threshold < MIN_ALLOWED_SCORE_THRESHOLD:
        p.error(f"{args.score_threshold} < {MIN_ALLOWED_SCORE_THRESHOLD}")

    return args


def main():
    args = parse_args()
    input_path = args.input_vcf
    output_gain_or_loss = "gain" if args.gain else "loss"
    output_bed_path = input_path.replace(".vcf", "").replace(".gz", "").replace(".bgz", "") + f".score_{args.score_threshold}.splice_{output_gain_or_loss}.bed"

    #output_allele_tsv_path = input_path.replace(".vcf", "").replace(".gz", "").replace(".bgz", "") + f".score_{args.score_threshold}.alleles.tsv"
    #output_junction_tsv_path = input_path.replace(".vcf", "").replace(".gz", "").replace(".bgz", "") + f".score_{args.score_threshold}.junctions.tsv"
    #junction_tsv_header = ["chrom", "start", "end", "length", "# of alleles", "alleles", "strand", "gene", "donor/acceptor", "loss/gain", "max score", "min score"]
    #allele_tsv_header = ["chrom", "start", "end", "length", "allele_type", "allele", "strand", "gene", "donor/acceptor", "loss/gain", "score"]

    # maps (chrom, start, end, strand, symbol) => ("DS_AG" or "DS_AL" or "DS_DG" or "DS_DL") => list of (ref, alt_allele, score)
    allele_score_dict = collections.defaultdict(lambda: collections.defaultdict(list))

    with gzip.open(input_path, "rt") as f, open(output_bed_path, "wt") as bed:
            #open(output_allele_tsv_path, "wt") as allele_tsv, \
            #open(output_junction_tsv_path, "wt") as junction_tsv:

        #allele_tsv.write("\t".join(allele_tsv_header) + "\n")
        #junction_tsv.write("\t".join(junction_tsv_header) + "\n")

        prev_chrom = prev_start = None
        for i, line in enumerate(tqdm.tqdm(f, unit=" lines")):
            if args.test and i > 1000000:
                break

            if line.startswith("#"):
                continue

            # parse line
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

            ds_dict = dict(zip(DS_HEADER, map(float, score_fields[2:6])))

            ds_keys_with_high_score = []
            for k, v in ds_dict.items():
                if (output_gain_or_loss == "gain" and k.endswith("G")) or (
                    output_gain_or_loss == "loss" and k.endswith("L")):
                    if v >= args.score_threshold:
                        ds_keys_with_high_score.append(k)

            if not ds_keys_with_high_score:
                continue

            dp_dict = dict(zip(DP_HEADER, map(int, score_fields[6:10])))

            #  process accumulated (chrom, start, end, strand, symbol) for the previous chrom, pos in the vcf
            if chrom != prev_chrom or start != prev_start:
                #pprint(dict(allele_score_dict))
                for (_chrom, _start, _end, _strand, _symbol), _values_dict in allele_score_dict.items():
                    # _values_dict is ("DS_AG" or "DS_AL" or "DS_DG" or "DS_DL") => list of (ref, alt_allele, score)
                    for _ds_key, _values_list in _values_dict.items():
                        _alleles_by_score = sorted(_values_list, key=lambda t: -t[2])
                        _allele_with_max_score = _alleles_by_score[0]
                        _max_score = max([_score for _, _, _score in _values_list])
                        assert _allele_with_max_score[2] == _max_score
                        _min_score = min([_score for _, _, _score in _values_list])
                        _alleles_by_type = collections.defaultdict(list)
                        _alleles = []
                        for (_ref, _alt_allele, _score) in _values_list:  # sorted(_values_list, key=lambda t: -t[2]):
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
                            _alleles_by_type[_allele_type].append(f"{_allele_string}:{_score}")
                            _alleles.append(f"{_allele_string}:{_score}")

                            #allele_tsv_row = [_chrom, min(_start, _end), max(_start, _end), abs(_end - _start), _allele_type, _allele_string, _strand, _symbol, _ds_key[-2], _ds_key[-1], _score]
                            #allele_tsv.write("\t".join(map(str, allele_tsv_row)) + "\n")

                        #junction_tsv_row = [_chrom, min(_start, _end), max(_start, _end), abs(_end - _start), len(_alleles), ",".join(_alleles), _strand, _symbol, _ds_key[-2], _ds_key[-1], _max_score, _min_score]
                        #junction_tsv.write("\t".join(map(str, junction_tsv_row)) + "\n")

                        #_alleles = ",".join(sorted(_alleles, key=lambda a: (len(a), a))) if len(_alleles) <= 3 else f"{len(_alleles)} alleles"

                        _D_or_A = _ds_key[-2]  # donor or acceptor
                        _G_or_L = _ds_key[-1]  # gain or loss
                        _donor_or_acceptor = "donor" if _D_or_A == "D" else "acceptor"
                        _gain_or_loss = "gain" if _G_or_L == "G" else "loss"
                        _shape = "v" if _G_or_L == "G" else "x"  # shape symbol corresponding to gain or loss
                        _color = "#11AA11" if _G_or_L == "G" else "#AA1111"   # green for gain, red for loss
                        #_label = f"{_ds_key[-2]}-{_gain_or_loss}:" + _alleles.replace(" ", "_") + f":{_max_score}"
                        _shape_key = "left_shape" if (_start < _end and _strand == "+") or (_end < _start and _strand == "-") else "right_shape"
                        _line_width = float((_max_score - (MIN_ALLOWED_SCORE_THRESHOLD - 0.05)) / (1 - MIN_ALLOWED_SCORE_THRESHOLD)) * 4.4

                        _effect = f"{_donor_or_acceptor}_{_gain_or_loss}"
                        #_alleles = ",".join([f"{_a}" for _a in sorted(_alleles, key=lambda a: (len(a), a))])

                        _gff_tags = f"label={_max_score}_[{len(_alleles)}];line_width={_line_width:0.1f};color={_color};{_shape_key}={_shape};effect={_effect};"
                        _gff_tags += f"summary={len(_alleles)}_of_11_alleles_scored_above_{args.score_threshold};"

                        if len(_alleles) > 1:
                            _gff_tags += f"max_score={_max_score};"
                            _gff_tags += f"allele_with_max_score={_chrom}-{_start}-{_allele_with_max_score[0]}-{_allele_with_max_score[1]};"
                            if len(_alleles_by_type['SNP']) > 0: _gff_tags += f"SNPs={',_'.join(sorted(_alleles_by_type['SNP'], key=lambda a: (len(a), a)))};"
                            if len(_alleles_by_type['DEL']) > 0: _gff_tags += f"DELs={',_'.join(sorted(_alleles_by_type['DEL'], key=lambda a: (len(a), a)))};"
                            if len(_alleles_by_type['INS']) > 0: _gff_tags += f"INS={',_'.join(sorted(_alleles_by_type['INS'], key=lambda a: (len(a), a)))};"
                        else:
                            _gff_tags += f"score={_max_score};"
                            _gff_tags += f"allele={_chrom}-{_start}-{_allele_with_max_score[0]}-{_allele_with_max_score[1]};"
                        # _alleles.replace(" ", "_") + f"__{_max_score}__{_ds_key[-2]}-{_gain_or_loss}"
                        bed_row = [_chrom, min(_start, _end) - 1, max(_start, _end), _gff_tags, '0', _strand]
                        bed.write("\t".join(map(str, bed_row)) + "\n")

                #print("clear allele_score_dict")
                allele_score_dict.clear()
                prev_chrom = chrom
                prev_start = start

            strand = SYMBOL_TO_STRAND_LOOKUP.get(symbol, "?")
            for ds_key in ds_keys_with_high_score:  # ds_dict:  #
                score = ds_dict[ds_key]
                dp_key = ds_key.replace("DS", "DP")
                pos = dp_dict[dp_key]
                end = start + pos

                allele_score_dict[(chrom, start, end, strand, symbol)][ds_key].append((ref, alt_allele, score))

        print("Done")


if __name__ == "__main__":
    main()
