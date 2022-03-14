import subprocess
from Bio.Blast.Applications import NcbiblastpCommandline, NcbimakeblastdbCommandline


def run_cmd(cmd, stdout=False, stderr=False):
    if not stdout:
        stdout_arg = subprocess.DEVNULL
    else:
        stdout_arg = subprocess.PIPE
    if not stderr:
        stderr_arg = subprocess.DEVNULL
    else:
        stderr_arg = subprocess.PIPE
    p = subprocess.run(cmd, stdout=stdout_arg, stderr=stderr_arg, check=True, shell=False)
    return p


def make_blast_database(input_file, out, dbtype='prot'):
    cline = NcbimakeblastdbCommandline(input_file=input_file, dbtype=dbtype, out=out)
    cline()


def run_blastp(query, db, out, threads=2):
    cline = NcbiblastpCommandline(
        query=query,
        out=out,
        db=db,
        evalue=1e-6,
        num_threads=threads,
        outfmt='6 qseqid sseqid pident length qlen slen',
    )
    cline()
