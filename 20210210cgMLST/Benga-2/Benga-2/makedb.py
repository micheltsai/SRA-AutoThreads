import os
import sys
import logging
import argparse
from pathlib import Path
from multiprocessing import Pool
from collections import defaultdict, Counter

import pandas as pd
from matplotlib import pyplot as plt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from utils import run_cmd, run_blastp, make_blast_database
from profiling import profiling, filter_blast_result
plt.style.use('ggplot')

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
console = logging.StreamHandler()
console.setLevel(logging.INFO)


def remove_non_cds_feature(in_gff, out_gff):
    with open(in_gff) as in_f, open(out_gff, 'w') as out_f:
        for line in in_f:
            if line != '##FASTA\n':
                if line.startswith('##'):
                    out_f.write(line)
                else:
                    if line.split()[2] == 'CDS':
                        out_f.write(line)
            else:
                out_f.write(line)
                out_f.write(in_f.read())


def annotate_genome(genome, out_path, prodigaltf):
    if os.path.exists(out_path):
        print(f"{out_path} is already exists, will use old result.", file=sys.stderr)
        return
    cmd = ['prokka', '--prefix', 'prokka', '--cpus', '2', '--outdir', out_path, genome]
    if prodigaltf:
        cmd += ["--prodigaltf", prodigaltf]
    run_cmd(cmd)


def plot_genome_coverage(data, figure_path, core=95):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(data, bins=range(102), histtype='step', lw=2)
    ax.set_title("Genome coverage distribution", fontsize=22)
    ax.set_xlabel("Percentage of genomes covered by loci (%)", fontsize=18)
    ax.set_ylabel("Number of locus", fontsize=18)
    ax.set_xticks(list(range(0, 101, 10)) + [core])
    fig.savefig(figure_path, dpi=300, facecolor='w', bbox_inches='tight')


def workflow(input_path, output_path, threads, prodigaltf):
    input_path = Path(input_path)
    output_path = Path(output_path)
    annot_dir = output_path/'Annotated'

    logging.info("Annotated")
    all_args = []
    for i in input_path.iterdir():
        if i.suffix in {'.fa', '.fna', '.fasta'}:
            args = (i, annot_dir/i.stem, prodigaltf)
            all_args.append(args)
    jobs_executor = Pool(threads//2 if threads//2 >= 2 else 1)
    for args in all_args:
        jobs_executor.apply_async(annotate_genome, args=args)
    jobs_executor.close()
    jobs_executor.join()

    gff_path = output_path/'GFF'
    gff_path.mkdir(exist_ok=True)
    for dirpath in annot_dir.iterdir():
        src_gff = dirpath/'prokka.gff'
        dst_gff = gff_path/(dirpath.name + '.gff')
        remove_non_cds_feature(src_gff, dst_gff)

    logging.info("Pan genome analysis")
    roary_path = output_path/'roary'
    run_cmd(["roary", '-p', str(threads), "-i", "95", "-s", "-f", roary_path] + [i for i in gff_path.iterdir()])

    clustered_proteins = roary_path / "clustered_proteins"
    seqid2gene = dict()
    f = open(clustered_proteins)
    for line in f:
        gene, cluster = line.strip().split(': ')
        for seqid in cluster.split():
            seqid2gene[seqid] = gene.replace(' ', '_')
    f.close()

    seq_occr = defaultdict(Counter)
    for i in annot_dir.iterdir():
        for record in SeqIO.parse(i/'prokka.ffn', 'fasta'):
            if record.id in seqid2gene:
                gene = seqid2gene[record.id]
                seq_occr[gene].update([record.seq])

    logging.info("Define locus reference sequence")
    pg_nucl_seqs = output_path/'pan_genome.fna'
    pg_prot_seqs = output_path/'pan_genome.faa'
    with open(pg_nucl_seqs, 'w') as nf, open(pg_prot_seqs, 'w') as pf:
        for gene, seq_count in seq_occr.items():
            rep_seq = seq_count.most_common(1)[0][0]  # select most occurrence sequence is representative
            nucl_rec = SeqRecord(rep_seq, id=gene, description='')
            prot_rec = SeqRecord(rep_seq.translate(table=11), id=gene, description='')
            SeqIO.write(nucl_rec, nf, 'fasta')
            SeqIO.write(prot_rec, pf, 'fasta')

    logging.info("Profiling")
    profile_path = output_path/'Profile'
    profile_path.mkdir(exist_ok=True)

    jobs_executor = Pool(threads)
    for i in annot_dir.iterdir():
        outfile = profile_path/(i.name + '.tsv')
        if outfile.exists():
            print(f"{outfile} is already exists.", file=sys.stderr)
        else:
            jobs_executor.apply_async(
                profiling,
                args=(i/'prokka.fna', pg_prot_seqs, outfile, prodigaltf),
            )
    jobs_executor.close()
    jobs_executor.join()

    logging.info("Recalculate loci occurrence")
    loci_occr = Counter()
    num_isolate = 0
    for filepath in profile_path.iterdir():
        profile = pd.read_csv(filepath, sep='\t')
        loci_occr.update(profile.dropna()['locus_id'])
        num_isolate += 1

    gene_presence_absence = pd.read_csv(
        roary_path/'gene_presence_absence.csv',
        usecols=['Gene', 'Annotation', 'No. isolates'],
        index_col=0
    )
    gene_presence_absence.index = gene_presence_absence.index.str.replace(' ', '_')
    new_gene_presence_absence = gene_presence_absence.filter(loci_occr.keys(), axis=0)  # keep loci occurrence not zero
    new_gene_presence_absence['No. isolates'] = new_gene_presence_absence.index.map(loci_occr)
    new_gene_presence_absence['occurrence'] = new_gene_presence_absence.index.map(
        {locus: count/num_isolate*100 for
         locus, count in loci_occr.items()}
    )

    logging.info('Drop duplicate loci')
    make_blast_database(pg_prot_seqs, pg_prot_seqs)

    blast_out = output_path/'self_align'
    run_blastp(pg_prot_seqs, pg_prot_seqs, blast_out, threads)

    drop_locus = set()
    for qseqid, sseqid in filter_blast_result(blast_out):
        if loci_occr[qseqid] > loci_occr[sseqid]:
            drop_locus.add(sseqid)
        else:
            drop_locus.add(qseqid)
    blast_out.unlink()

    pan_genome_results = new_gene_presence_absence.drop(index=drop_locus)
    pan_genome_results['Nucleotide'] = pan_genome_results.index.map(
        {record.id: str(record.seq) for record in SeqIO.parse(pg_nucl_seqs, 'fasta')}
    )
    pan_genome_results['Peptide'] = pan_genome_results.index.map(
        {record.id: str(record.seq) for record in SeqIO.parse(pg_prot_seqs, 'fasta')}
    )
    for file in output_path.glob('pan_genome*'):  # remove old pan genome sequence and blast database files
        file.unlink()
    pan_genome_results.to_csv(output_path/'pan_genome_info.txt', sep='\t')

    logging.info('Plot genome coverage')
    plot_genome_coverage(
        pan_genome_results['occurrence'],
        output_path/'genome_coverage.png',
    )
    plot_genome_coverage(
        pan_genome_results[pan_genome_results['occurrence'] > 5]['occurrence'],
        output_path/'genome_coverage_5_prec.png',
    )
    logging.info('Done')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        required=True,
                        help="Path of query genomes.")
    parser.add_argument("-o", "--output_path",
                        required=True,
                        help="Path of output directory.")
    parser.add_argument("--prodigaltf",
                        default='',
                        help="Path of prodigal training file. default:''")
    parser.add_argument("-t", "--threads",
                        type=int,
                        default=2,
                        help="Number of threads. default: 2")
    args = parser.parse_args()

    workflow(args.input, args.output_path, args.threads, args.prodigaltf)


if __name__ == '__main__':
    main()
