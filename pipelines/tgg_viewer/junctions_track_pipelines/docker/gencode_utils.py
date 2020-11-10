import collections
import gzip


def _compute_introns_for_chrom(chrom, transcript_to_exons):
    introns = []
    for (transcript_id, gene_id, gene_name), exons_in_transcript in transcript_to_exons.items():
        prev_exon_3prime = None
        for exon_number, start_1based, end_1based, strand in sorted(exons_in_transcript):
            if strand == "+":
                current_exon_5prime = int(start_1based)
                if exon_number > 1:
                    assert prev_exon_3prime < current_exon_5prime
                    intron_start_1based = prev_exon_3prime + 1
                    intron_end_1based = current_exon_5prime - 1
                prev_exon_3prime = int(end_1based)
            elif strand == "-":
                current_exon_5prime = int(end_1based)
                if exon_number > 1:
                    assert prev_exon_3prime > current_exon_5prime
                    intron_start_1based = current_exon_5prime + 1
                    intron_end_1based = prev_exon_3prime - 1
                prev_exon_3prime = int(start_1based)
            else:
                print(f"ERROR: strand == {strand}")
                continue

            if exon_number > 1:
                introns.append((chrom, intron_start_1based, intron_end_1based))

                #introns.append({
                #"transcript_id": transcript_id,
                #"chrom": chrom,
                #"start_1based": str(intron_start_1based - 0),
                #"end_1based": str(intron_end_1based),
                #"strand": strand,
                #"feature_number": str(exon_number - 1),
                #"feature": "intron",
                #"num_exons_in_transcript": len(exons_in_transcript),
                #"gene_id": gene_id,
                #"gene_name": gene_name,
                #"key": f"{chrom}:{intron_start_1based}-{intron_end_1based} ({strand})"
                #})

    return introns


def parse_introns_from_gencode_gff(gencode_gff_path):
    print(f"Parsing: {gencode_gff_path}")
    # example line: chr1	HAVANA	exon	11869	12227	.	+	.	ID=exon:ENST00000456328.2:1;Parent=ENST00000456328.2;gene_id=ENSG00000223972.5;transcript_id=ENST00000456328.2;gene_type=transcribed_unprocessed_pseudogene;gene_name=DDX11L1;transcript_type=processed_transcript;transcript_name=DDX11L1-002;exon_number=1;exon_id=ENSE00002234944.1;level=2;transcript_support_level=1;tag=basic;havana_gene=OTTHUMG00000000961.2;havana_transcript=OTTHUMT00000362751.1
    introns = []
    with (gzip.open if gencode_gff_path.endswith("gz") else open)(gencode_gff_path, "rt") as f:
        prev_chrom = None
        transcript_to_exons = collections.defaultdict(list)
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "exon":
                continue
            chrom = fields[0]
            if chrom != prev_chrom:
                if prev_chrom is not None:
                    introns_for_prev_chrom = _compute_introns_for_chrom(prev_chrom, transcript_to_exons)
                    print(f"Found {len(introns_for_prev_chrom)} introns on {prev_chrom}")
                    introns.extend(introns_for_prev_chrom)
                transcript_to_exons = collections.defaultdict(list)
                prev_chrom = chrom

            start_1based = fields[3]
            end_1based = fields[4]
            strand = fields[6]
            annot = dict([key_value.split("=") for key_value in fields[8].split(";")])
            exon_number = annot['exon_number']
            transcript_id = annot['transcript_id'].split(".")[0]
            gene_id = annot['gene_id'].split(".")[0]
            gene_name = annot['gene_name']
            transcript_to_exons[(transcript_id, gene_id, gene_name)].append((int(exon_number), start_1based, end_1based, strand))

        introns_for_prev_chrom = _compute_introns_for_chrom(chrom, transcript_to_exons)
        print(f"Found {len(introns_for_prev_chrom)} introns on {chrom}")
        introns.extend(introns_for_prev_chrom)

    #for intron in introns:
    #    intron['feature_size (bp)'] = str(int(intron['end_1based']) - int(intron['start_1based']) + 1)
    #    intron["is_known"] = True

    return set(introns)

