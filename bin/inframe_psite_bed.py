#!/usr/bin/env python3
import sys
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option(
    "-i", "--input",
    dest="gtff",
    help="Input GTF file with ORF definitions"
)
parser.add_option(
    "-a", "--annot",
    dest="annot",
    default="yes",
    help="Is the input an annotation file (yes/no)? Default: yes"
)
parser.add_option(
    "-o", "--outdir",
    dest="outdir",
    help="Output directory for generated BED file"
)
parser.add_option(
    "-t", "--id-type",
    dest="id_type",
    default="transcript_id",
    help="Identifier field to use when parsing the GTF ('transcript_id' or 'ORF_id'). Default: transcript_id"
)
options, _ = parser.parse_args()

gtff    = options.gtff
annot   = options.annot
outdir  = options.outdir
id_type = options.id_type

# Validate required arguments
if not gtff or not outdir:
    parser.error("Both -i (input GTF) and -o (output directory) are required.")
if not gtff.endswith(".gtf"):
    sys.exit("Error: Only GTF input is supported (must end with .gtf).")

# Ensure output directory exists
os.makedirs(outdir, exist_ok=True)

# Open the output BED file
output_path = os.path.join(outdir, os.path.basename(gtff) + "_psites_plus_partial.bed")
out = open(output_path, "w+")

# ----------------------------
# Define a simple transcript container
# ----------------------------
class TransObject:
    def __init__(self, chrom, gene, strand):
        self.chrom = chrom
        self.gene = gene
        self.strand = strand
        self.start_coords = []
        self.end_coords = []


# ----------------------------
# parse_gtf: collect CDS coordinates per transcript
# ----------------------------
def parse_gtf(gtf_path, feature, id_field):
    transcripts = {}
    with open(gtf_path) as fh:
        for line in fh:
            if f"\t{feature}\t" not in line:
                continue

            fields = line.split("\t")
            chrom, source, feat, start, end, _, strand, _, attr = fields

            # Extract transcript or ORF ID
            try:
                if id_field == "transcript_id":
                    tid = attr.split('transcript_id "')[1].split('"')[0]
                else:  # ORF_id
                    tid = attr.split('ORF_id "')[1].split('"')[0]
                gid = attr.split('gene_id "')[1].split('"')[0]
            except IndexError:
                continue

            if tid not in transcripts:
                transcripts[tid] = TransObject(chrom, gid, strand)

            transcripts[tid].start_coords.append(int(start))
            transcripts[tid].end_coords.append(int(end))

    # Sort each transcript's start/end lists
    for tx in transcripts.values():
        paired = sorted(zip(tx.start_coords, tx.end_coords))
        tx.start_coords, tx.end_coords = zip(*paired)

    return transcripts


# ----------------------------
# collect_status: determine coding completeness per transcript
# ----------------------------
def collect_status(gtf_path):
    status_dict = {}
    with open(gtf_path) as fh:
        for line in fh:
            if not any(tag in line for tag in ["\tCDS\t", "\tstart_codon\t", "\tstop_codon\t"]):
                continue

            try:
                tid = line.split('transcript_id "')[1].split('"')[0]
            except IndexError:
                continue

            if "\tCDS\t" in line:
                status_dict.setdefault(tid, 0)
            elif "\tstart_codon\t" in line:
                if status_dict.get(tid, 0) == 2:
                    status_dict[tid] = 3
                else:
                    status_dict[tid] = 1
            elif "\tstop_codon\t" in line:
                if status_dict.get(tid, 0) == 1:
                    status_dict[tid] = 3
                else:
                    status_dict[tid] = 2

    return status_dict


# ----------------------------
# process_transcript: write P-site lines for one transcript
# ----------------------------
def process_transcript(tid, tx, coding_status):
    """
    For a given transcript 'tx', iterate its CDS coordinates,
    assign frames (p0, p1, p2) according to 'coding_status' and 'strand',
    and write to the global 'out' file handle.
    """
    # Flatten all CDS positions into a single list
    coords = []
    for s, e in zip(tx.start_coords, tx.end_coords):
        coords.extend(range(s, e + 1))

    # Reverse for minus strand so coords[0] is biologic start codon
    if tx.strand == "-":
        coords = coords[::-1]

    # Compute offset for partial codons
    offset = len(coords) % 3

    for i, pos in enumerate(coords):
        if coding_status in (1, 3):
            # Has start codon: frame anchors at coords[0]
            frame = i % 3
        elif coding_status == 2:
            # Stop-codon only: shift frame forward by offset (plus strand logic)
            if tx.strand == "+":
                frame = (i + offset) % 3
            else:
                # For minus strand, coords reversed: shift backward by offset
                frame = (i - offset) % 3
        else:
            # status == 0: no start or stop codon → exclude transcript entirely
            return

        # Write the P-site line
        out.write(f"{tx.chrom}\t{pos}\t{pos}\t{tid}\tp{frame}\t{tx.strand}\n")


# ----------------------------
# Main execution flow
# ----------------------------
gtf_data = parse_gtf(gtff, "CDS", id_type)
status   = collect_status(gtff) if annot.lower() == "yes" else {}

for tid, tx in gtf_data.items():
    st = status.get(tid, 3)
    if st == 0:
        print(f"{tid} excluded")
        continue

    process_transcript(tid, tx, st)

out.close()