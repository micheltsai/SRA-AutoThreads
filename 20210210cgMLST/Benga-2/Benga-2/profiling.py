import os
import hashlib
import argparse
from tempfile import TemporaryDirectory
import pandas as pd
from Bio import SeqIO
from utils import run_cmd, make_blast_database, run_blastp


def make_diamond_database(input_file, out):
    run_cmd(["diamond", "makedb", "--in", input_file, "-d", out])


def encode_sequence_id(record):
    """
    Convert nucleotide sequence to hash codes by hash function 'sha256'.
    """
    return hashlib.sha256(str(record.seq).encode("ascii")).hexdigest()


def predict_open_reading_frame(infile, outfile, prodigaltf=''):
    """
    Predict open reading frame and export nucleotide sequence
    """
    cmd = ['prodigal', '-i', infile, '-d', outfile, '-g', '11', '-c', '-m', '-q']
    if prodigaltf:
        cmd += ['-t', prodigaltf]
    run_cmd(cmd)


def run_diamond(query, db, out, threads=2):
    cmd = [
        'diamond', 'blastp',
        '-q', query,
        '-d', db,
        '-f', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
        '-o', out,
        '-p', str(threads), '-e', '1e-6',
        '--fast',
    ]
    run_cmd(cmd)


def filter_blast_result(blast_out, identity=95, min_cov=.75, max_cov=1.25):
    """
    Filter identity >= 95, ratio of query sequence length to subject sequence length between 0.75 and 1.25
    """
    handle = open(blast_out)
    for line in handle:
        qseqid, sseqid, pident, length, qlen, slen = line.strip().split()
        pident, length, qlen, slen = float(pident), float(length), float(qlen), float(slen)
        if qseqid != sseqid and pident >= identity and min_cov <= qlen/slen < max_cov and min_cov <= qlen/length < max_cov:
            yield sseqid, qseqid
    handle.close()


def profiling(infile, scheme, outfile, prodigaltf, fast, threads=2, match_output=False):
    locus_ids = set(record.id for record in SeqIO.parse(scheme, 'fasta'))

    with TemporaryDirectory(dir="/tmp") as tmpdir:
        nucl_orfs = os.path.join(tmpdir, 'orfs.fna')
        prot_orfs = os.path.join(tmpdir, 'orfs.faa')
        blast_db = os.path.join(tmpdir, 'scheme')
        blast_output = os.path.join(tmpdir, 'blast.out')
        
        predict_open_reading_frame(infile, nucl_orfs, prodigaltf=prodigaltf)
        nucl_records = []
        for record in SeqIO.parse(nucl_orfs, 'fasta'):
            record.id = encode_sequence_id(record)
            nucl_records.append(record)
        prot_records = [record.translate(table=11, id=True) for record in nucl_records]
        SeqIO.write(prot_records, prot_orfs, 'fasta')

        if fast:
            make_diamond_database(scheme, blast_db)
            run_diamond(prot_orfs, blast_db, blast_output, threads)
        else:
            make_blast_database(scheme, blast_db)
            run_blastp(prot_orfs, blast_db, blast_output, threads)

        df = pd.DataFrame(filter_blast_result(blast_output), columns=['locus_id', 'allele_id'])
        df = df.sort_values('allele_id', kind='mergesort').drop_duplicates('locus_id')
        df = df.set_index('locus_id').reindex(locus_ids).sort_values('locus_id')
        df.to_csv(outfile, sep='\t', index=True)

    if match_output:
        hits = df['allele_id'].dropna().to_dict()
        f = open(allele_nucleotide_output, 'w')
        for record in nucl_records:
            if record.id in hits:
                record.description = hits[record.id]
                f.write(record.format('fasta'))
        f.close()


def main():
    parser = argparse.ArgumentParser('Benga')
    parser.add_argument("-i", "--input",
                        required=True,
                        help="Path of query genome.")
    parser.add_argument("-o", "--output",
                        required=True,
                        help="Path of output file.")
    parser.add_argument("-s", "--scheme",
                        required=True,
                        help="Core-genome MLST scheme, it should is peptide sequence.")
    parser.add_argument("--prodigaltf",
                        default='',
                        help="Path of prodigal training file. default: None")
    parser.add_argument("-t", "--threads",
                        type=int,
                        default=2,
                        help="Number of threads. default: 2")
    parser.add_argument("--fast",
                        action="store_true",
                        help="Use diamond to alignment sequence. default: False")
    parser.add_argument("--match_output",
                        default="",
                        help="Output nucleotide FASTA file of allele nucleotide sequences. default: ''")
    args = parser.parse_args()

    profiling(
        infile=args.input,
        scheme=args.scheme,
        outfile=args.output,
        prodigaltf=args.prodigaltf,
        threads=args.threads,
        fast=args.fast,
        match_output=args.match_output
    )


if __name__ == '__main__':
    main()
