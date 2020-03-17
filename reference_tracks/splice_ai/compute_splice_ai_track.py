import argparse
import gzip
import os
import pyBigWig
import tqdm

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-r", "--fasta", help="Reference fasta file. This script reads the .fai index to get chromosome sizes.", default="hg38.fa")
    p.add_argument("-a", "--only", help="only output scores for this alt allele", choices="ACGT")
    p.add_argument("splice_ai_vcf", help="Path of spliceAI precomputed scores vcf")
    args = p.parse_args()

    args.fasta = os.path.realpath(os.path.expanduser(args.fasta))
    args.splice_ai_vcf = os.path.realpath(os.path.expanduser(args.splice_ai_vcf))

    fai_path = f"{args.fasta}.fai"
    if not os.path.isfile(fai_path):
        p.error(f"{fai_path} not found")
    if not os.path.isfile(args.splice_ai_vcf):
        p.error(f"{args.splice_ai_vcf} not found")

    return args


def get_chrom_sizes(fai_path):
    chrom_sizes = []
    with open(fai_path) as f:
        for line in f:
            fields = line.split("\t")
            chrom = fields[0]
            if "_" in chrom or "-" in chrom or "ebv" in chrom.lower():
                continue
            size = int(fields[1])
            chrom_sizes.append((chrom, size))
    
    return chrom_sizes


def main():
    args = parse_args()
    
    chrom_sizes = get_chrom_sizes(f"{args.fasta}.fai")

    output_path_prefix = args.splice_ai_vcf.replace(".vcf", "").replace(".gz", "") 
    if args.only:
        output_path_prefix += f".alt-allele-{args.only}"
    else:
        output_path_prefix += f".all_alleles"

    output_path = f"{output_path_prefix}.bigWig"
    print(f"Writing to {output_path}")

    info_splice_ai_fields = "ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL".split("|")

    bins = {
        0: 0,
        0.2: 0,
        0.3: 0,
        0.5: 0,
        0.8: 0,
    }

    bw = pyBigWig.open(output_path, "w")
    bw.addHeader(chrom_sizes)

    chroms_in_header = set([t[0] for t in chrom_sizes])

    use_chr_prefix = list(chroms_in_header)[0].startswith("chr")

    total_num_lines = 3_433_386_166 
    if "indel" in args.splice_ai_vcf:
        total_num_lines * 4/3.

    with gzip.open(args.splice_ai_vcf, "rt") as f:
        line = ''
        while not line.startswith("#CHROM"):
            line = next(f)

        skipped_counter = 0
        prev_chrom = None
        prev_pos = None
        current_scores = []
        for i, line in enumerate(tqdm.tqdm(f, unit=" lines", total=total_num_lines)):
            line = line.rstrip()
            fields = line.split("\t")
            chrom = fields[0].replace('chr', '')
            if use_chr_prefix:
                chrom = f"chr{chrom}" 
            if chrom not in chroms_in_header:
                skipped_counter += 1
                continue
            pos = int(fields[1])
            if not args.only and pos != prev_pos and current_scores:
                max_score = max(current_scores)
                if max_score >= 0.2:
                    bw.addEntries(chrom, [prev_pos], values=[max_score], span=1)
                    #bw.write(prev_chrom, prev_pos - 1, prev_pos, max_score)
        
                current_scores.clear()
                prev_chrom = chrom
                prev_pos = pos

                for k in bins:
                    if max_delta >= k:
                        bins[k] += 1
            
            #ref = fields[3]
            #alt = fields[4]
            idx = fields[7].index("=") + 1
            delta_scores = dict(zip(info_splice_ai_fields, fields[7][idx:].split("|")))
            max_delta = max(float(delta_scores[k]) for k in ('DS_AG', 'DS_AL', 'DS_DG', 'DS_DL'))
            if args.only:
                alt = fields[4]
                if alt == args.only and max_delta >= 0.2:
                    bw.addEntries(chrom, [pos], values=[max_delta], span=1)

                for k in bins:
                    if max_delta >= k:
                        bins[k] += 1
            else:
                current_scores.append(max_delta)

            #if i > 100_000_0000:
            #if i > total_num_lines:
            #    break

    bw.close()
    
    print("Done")
    for k, v in bins.items():
        print(f"{v:10d} rows >= {k:0.2f} {(100.0*v)/bins[0]:3.1f}%")    
    print(f"Skipped {skipped_counter} rows")

if __name__ == "__main__":
    main()
