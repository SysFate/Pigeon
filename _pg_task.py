#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Matthias Blum
# Contact: mat.blum@gmail.com

import os
import shutil
import tempfile
import zipfile
from subprocess import Popen, PIPE

from _pg_model import DT_CHIP, DT_HIC, DT_RNA

ERRORS = {
    'download': 1,
    'extract': 2,
    'align': 3,
    'post-align': 4,
    'merge': 5,
    'ratio': 6,
    'rsem': 7
}

SAM2BED = os.path.join(os.path.dirname(__file__), 'utils', 'sam2bed')
JOINBED = os.path.join(os.path.dirname(__file__), 'utils', 'joinbed')


def download_run(run, workdir):
    filename = os.path.join(workdir, os.path.basename(run.url))
    log_file = filename + '.log'

    with open(log_file, 'w') as fh:
        pop = Popen(['wget', '-T', '120', '-t', '5', '-O', filename, run.url], stderr=fh, stdout=fh)
        pop.wait()

    os.unlink(log_file)

    if pop.returncode == 0 and os.path.isfile(filename):
        run.sra_file = filename
        run.is_downloaded = True
        return True
    else:
        try:
            os.unlink(filename)
        except OSError:
            pass

        run.status = ERRORS['download']
        return False


def extract_sra(run, workdir, **kwargs):
    fastqdump = kwargs.get('fastqdump', 'fastq-dump')
    prefix = os.path.join(workdir, os.path.basename(run.sra_file)[:-4])

    if run.singleread:
        cmd = [fastqdump, '-F', '--gzip', '-O', workdir, run.sra_file]
        fq1 = prefix + '.fastq.gz'
        fq2 = fq1
    else:
        cmd = [fastqdump, '-F', '--gzip', '-O', workdir, '--split-files', run.sra_file]
        fq1 = prefix + '_1.fastq.gz'
        fq2 = prefix + '_2.fastq.gz'

    with open(os.devnull, 'w') as fh:
        pop = Popen(cmd, stdout=fh, stderr=fh)
        pop.wait()

    os.unlink(run.sra_file)

    if pop.returncode == 0 and os.path.isfile(fq1) and os.path.isfile(fq2):
        run.fq1_file = fq1
        run.fq2_file = None if run.singleread else fq2
        run.is_extracted = True
        return True
    else:
        try:
            os.unlink(fq1)
        except OSError:
            pass

        try:
            os.unlink(fq2)
        except OSError:
            pass

        run.status = ERRORS['extract']
        return False


def align_chip(run, workdir, ref_dir, **kwargs):
    bowtie2 = kwargs.get('bowtie2', 'bowtie2')
    flags = kwargs.get('bowtie2flags', '')
    mapq = kwargs.get('mapq', 0)
    threads = kwargs.get('threads', 1)

    ref = os.path.join(ref_dir, run.assembly)
    prefix = os.path.join(workdir, os.path.basename(run.url)[:-4])
    reads_file = prefix + '.bed.gz'

    cmd = [bowtie2, '-x', ref, '-p', str(threads)] + flags.split()

    if run.singleread:
        cmd += ['-U', run.fq1_file]
    else:
        cmd += ['-1', run.fq1_file, '-2', run.fq2_file]

    with open(reads_file, 'wb') as fh, open(os.devnull, 'w') as devnull:
        pop1 = Popen(cmd, stdout=PIPE, stderr=devnull)
        pop2 = Popen([SAM2BED, '-', str(mapq), '--no-tag'], stdin=pop1.stdout, stdout=PIPE, stderr=PIPE)
        pop3 = Popen(['gzip', '-c'], stdin=pop2.stdout, stdout=fh, stderr=devnull)
        pop1.stdout.close()
        pop1.wait()
        out, err = pop2.communicate()
        pop3.wait()

    os.unlink(run.fq1_file)

    try:
        # Paired-end data
        os.unlink(run.fq2_file)
    except (OSError, TypeError):
        pass

    if pop1.returncode == pop2.returncode == pop3.returncode == 0 and os.path.isfile(reads_file):
        reads, mapped_reads = err.decode('utf8').strip().split('\t')
        run.reads = int(reads) if run.singleread else int(reads) / 2
        run.mapped_reads = int(mapped_reads)
        run.is_aligned = True
        run.reads_file = reads_file
        return True
    else:
        try:
            os.unlink(reads_file)
        except OSError:
            pass

        run.status = ERRORS['align']
        return False


def align_hic(run, workdir, ref_dir, **kwargs):
    bowtie2 = kwargs.get('bowtie2', 'bowtie2')
    flags = kwargs.get('bowtie2_flags', '')
    mapq = kwargs.get('mapq', 0)
    threads = kwargs.get('threads', 1)

    ref = os.path.join(ref_dir, run.assembly)
    prefix = os.path.join(workdir, os.path.basename(run.url)[:-4])

    fw_reads_file = prefix + '_1.bed.gz'
    rv_reads_file = prefix + '_2.bed.gz'
    reads_file = prefix + '.tsv.gz'

    cmd = [bowtie2, '-x', ref, '-p', str(threads), '--reorder'] + flags.split()

    with open(fw_reads_file, 'wb') as fh, open(os.devnull, 'w') as devnull:
        pop1 = Popen(cmd + ['-U', run.fq1_file], stdout=PIPE, stderr=devnull)
        pop2 = Popen([SAM2BED, '-', str(mapq), '--report-all'], stdin=pop1.stdout, stdout=PIPE, stderr=PIPE)
        pop3 = Popen(['gzip', '-c'], stdin=pop2.stdout, stdout=fh, stderr=devnull)
        pop1.stdout.close()
        pop1.wait()
        pop2.wait()
        pop3.wait()

    os.unlink(run.fq1_file)

    if pop1.returncode != 0 or pop2.returncode != 0 or pop3.returncode != 0 or not os.path.isfile(fw_reads_file):
        os.unlink(run.fq2_file)

        try:
            os.unlink(fw_reads_file)
        except OSError:
            pass
        finally:
            run.status = ERRORS['align']
            return False

    with open(rv_reads_file, 'wb') as fh, open(os.devnull, 'w') as devnull:
        pop1 = Popen(cmd + ['-U', run.fq2_file], stdout=PIPE, stderr=devnull)
        pop2 = Popen([SAM2BED, '-', str(mapq), '--report-all'], stdin=pop1.stdout, stdout=PIPE, stderr=PIPE)
        pop3 = Popen(['gzip', '-c'], stdin=pop2.stdout, stdout=fh, stderr=devnull)
        pop1.stdout.close()
        pop1.wait()
        pop2.wait()
        pop3.wait()

    os.unlink(run.fq2_file)

    if pop1.returncode != 0 or pop2.returncode != 0 or pop3.returncode != 0 or not os.path.isfile(rv_reads_file):
        try:
            os.unlink(rv_reads_file)
        except OSError:
            pass
        finally:
            run.status = ERRORS['align']
            return False

    # Create named pipe (forward)
    fd, fw_tmp = tempfile.mkstemp()
    os.close(fd)
    os.unlink(fw_tmp)
    os.mkfifo(fw_tmp)

    # Create named pipe (forward)
    fd, fw_tmp = tempfile.mkstemp()
    os.close(fd)
    os.unlink(fw_tmp)
    os.mkfifo(fw_tmp)

    # Create named pipe (reverse)
    fd, rv_tmp = tempfile.mkstemp()
    os.close(fd)
    os.unlink(rv_tmp)
    os.mkfifo(rv_tmp)

    # Listen on named pipes
    fh = open(reads_file, 'wb')
    devnull = open(os.devnull, 'w')
    pop1 = Popen([JOINBED, fw_tmp, rv_tmp], stdout=PIPE, stderr=PIPE)
    pop2 = Popen(['gzip', '-c'], stdin=pop1.stdout, stdout=fh, stderr=devnull)

    # Write to named pipes
    fh_fw = open(fw_tmp, 'w')
    fh_rv = open(rv_tmp, 'w')
    pop3 = Popen(['gzip', '-cd', fw_reads_file], stdout=fh_fw, stderr=devnull)
    pop4 = Popen(['gzip', '-cd', rv_reads_file], stdout=fh_rv, stderr=devnull)

    # Close file handlers
    fh_fw.close()
    fh_rv.close()
    devnull.close()

    # Wait for writing processes to finish
    pop3.wait()
    pop4.wait()

    # Wait for joinbed/gzip to finish
    pop2.wait()
    out, err = pop1.communicate()

    # Delete named pipes
    os.unlink(fw_tmp)
    os.unlink(rv_tmp)

    # Delete independent aligned files
    os.unlink(fw_reads_file)
    os.unlink(rv_reads_file)

    if pop1.returncode == pop2.returncode == pop3.returncode == pop4.returncode and os.path.isfile(reads_file):
        reads, mapped_reads = err.decode('utf8').strip().split('\t')
        run.reads = int(reads) if run.singleread else int(reads) / 2
        run.mapped_reads = int(mapped_reads)
        run.is_aligned = True
        run.reads_file = reads_file
        return True
    else:
        try:
            os.unlink(reads_file)
        except OSError:
            pass
        finally:
            run.status = ERRORS['post-align']
            return False


def align_rna(run, workdir, ref_dir, gtf_dir=None, **kwargs):
    star = kwargs.get('star', 'STAR')
    flags = kwargs.get('star_flags', '')
    mapq = kwargs.get('mapq', 0)
    threads = kwargs.get('threads', 1)

    # Add a trailing '/' to force STAR to put files INSIDE the directory
    outdir = tempfile.mkdtemp(dir=workdir) + '/'

    ref = os.path.join(ref_dir, run.assembly)

    if gtf_dir:
        gtf = os.path.join(gtf_dir, run.assembly + '.gtf')
    else:
        gtf = None

    prefix = os.path.join(workdir, os.path.basename(run.url)[:-4])
    reads_file = prefix + '.bed.gz'
    reads_file_2 = os.path.join(outdir, 'Aligned.toTranscriptome.out.bam')

    cmd = [
        star, '--genomeDir', ref,
        '--outFileNamePrefix', outdir,
        '--outStd', 'SAM',
        # Report alignments to transcriptome (needed by RSEM)
        '--quantMode', 'TranscriptomeSAM',
        # Report unmapped reads
        '--outSAMunmapped', 'Within', 'KeepPairs',
        '--runThreadN', str(threads)
    ] + flags.split()

    if gtf and os.path.isfile(gtf):
        cmd += ['--sjdbGTFfile', gtf]

    if run.singleread:
        cmd += ['--readFilesIn', run.fq1_file]
    else:
        cmd += ['--readFilesIn', run.fq1_file, run.fq2_file]

    cmd += ['--readFilesCommand', 'gzip', '-cd']

    with open(reads_file, 'wb') as fh, open(os.devnull, 'w') as devnull:
        pop1 = Popen(cmd, stdout=PIPE, stderr=devnull)
        pop2 = Popen([SAM2BED, '-', str(mapq)], stdin=pop1.stdout, stdout=PIPE, stderr=PIPE)
        pop3 = Popen(['gzip', '-c'], stdin=pop2.stdout, stdout=fh, stderr=devnull)
        pop1.stdout.close()
        pop1.wait()
        out, err = pop2.communicate()
        pop3.wait()

    os.unlink(run.fq1_file)

    try:
        # Paired-end data
        os.unlink(run.fq2_file)
    except (OSError, TypeError):
        pass

    if pop1.returncode == pop2.returncode == pop3.returncode == 0 and os.path.isfile(reads_file) \
            and os.path.isfile(reads_file_2):
        reads, mapped_reads = err.decode('utf8').strip().split('\t')
        run.reads = int(reads) if run.singleread else int(reads) / 2
        run.mapped_reads = int(mapped_reads)
        run.is_aligned = True

        # Move BAM file
        shutil.move(reads_file_2, prefix + '.bam')
        run.reads_file_2 = prefix + '.bam'

        # Delete output directory (we don't use the other files)
        shutil.rmtree(outdir)

        # BED file
        run.reads_file = reads_file
        return True
    else:
        shutil.rmtree(outdir)

        try:
            os.unlink(reads_file)
        except OSError:
            pass
        finally:
            run.status = ERRORS['align']
            return False


def merge_runs(dset, workdir, ratio=0, **kwargs):
    maxmem = kwargs.get('maxmem', 0)
    samtools = kwargs.get('samtools', 'samtools')

    if dset.sample_id:
        basename = dset.sample_id + '_' + dset.exp_id
    else:
        basename = dset.exp_id

    if len(dset.runs) == 1:
        basename += '_' + os.path.basename(dset.runs[0].url[:-4])

    tot_reads = sum([r.reads for r in dset.runs])
    tot_mapped_reads = sum([r.mapped_reads for r in dset.runs])

    if tot_mapped_reads < ratio * tot_reads:
        for r in dset.runs:
            os.unlink(r.reads_file)

        if dset.data_type == DT_RNA:
            for r in dset.runs:
                os.unlink(r.reads_file_2)

        dset.status = ERRORS['ratio']
        return False
    else:
        tmp_dir = tempfile.mkdtemp(dir=workdir)

        if dset.data_type == DT_HIC:
            sort_cmd = ['sort', '-k1,1V', '-k2,2n', '-k4,4V', '-k5,5n', '-T', tmp_dir]
        else:
            sort_cmd = ['sort', '-k1,1V', '-k2,2n', '-k6,6', '-T', tmp_dir]

        if maxmem:
            sort_cmd.append('--buffer-size={}M'.format(maxmem))

        reads_files = [r.reads_file for r in dset.runs]
        reads_file = os.path.join(workdir, basename) + '.bed.gz'

        with open(reads_file, 'wb') as fh, open(os.devnull, 'w') as devnull:
            pop1 = Popen(['gzip', '-cd'] + reads_files, stdout=PIPE, stderr=devnull)
            pop2 = Popen(sort_cmd, stdin=pop1.stdout, stdout=PIPE, stderr=devnull)
            pop3 = Popen(['gzip', '-c'], stdin=pop2.stdout, stdout=fh, stderr=devnull)
            pop1.stdout.close()
            pop2.stdout.close()
            pop1.wait()
            pop2.wait()
            pop3.wait()

        shutil.rmtree(tmp_dir)

        for f in reads_files:
            os.unlink(f)

        if pop1.returncode == pop2.returncode == pop3.returncode == 0 and os.path.isfile(reads_file):
            if dset.data_type == DT_RNA:
                reads_files = [r.reads_file_2 for r in dset.runs]
                reads_file_2 = os.path.join(workdir, basename) + '.bam'

                if len(reads_files) > 1:
                    with open(os.devnull, 'w') as devnull:
                        pop = Popen([samtools, 'merge', reads_file_2] + reads_files, stdout=devnull, stderr=devnull)
                        pop.wait()

                    for f in reads_files:
                        os.unlink(f)

                    if pop.returncode == 0 and os.path.isfile(reads_file_2):
                        dset.reads_file_2 = reads_file_2
                    else:
                        os.unlink(reads_file)
                        try:
                            os.unlink(reads_file_2)
                        except OSError:
                            pass
                        finally:
                            return False
                else:
                    move(reads_files[0], reads_file_2)
                    dset.reads_file_2 = reads_file_2

            dset.is_merged = True
            dset.reads_file = reads_file
            return True
        else:
            try:
                os.unlink(reads_file)
            except OSError:
                pass
            finally:
                dset.status = ERRORS['merge']
                return False


def run_rsem(dset, workdir, ref_dir, **kwargs):
    rsem_calc_exp = kwargs.get('rsem', 'rsem-calculate-expression')
    threads = kwargs.get('threads', 1)

    cmd = [rsem_calc_exp, '--bam', '--no-bam-output', '-p', str(threads)]

    if not dset.singleread:
        cmd.append('--paired-end')

    # Prefix of files generated by rsem-prepare-reference
    reference_name = os.path.join(ref_dir, dset.assembly, 'STAR')

    outdir = tempfile.mkdtemp(dir=workdir)
    sample_name = os.path.join(outdir, 'rsem')

    cmd += [dset.reads_file_2, reference_name, sample_name]

    with open(os.devnull, 'w') as fh:
        pop = Popen(cmd, stdout=fh, stderr=fh)
        pop.wait()

    # Delete BAM file
    os.unlink(dset.reads_file_2)

    genes_file = sample_name + '.genes.results'
    isoforms_file = sample_name + '.isoforms.results'

    if pop.returncode == 0 and os.path.isfile(genes_file) and os.path.isfile(isoforms_file):
        zip_file = dset.reads_file_2[:-4] + '.zip'

        with zipfile.ZipFile(zip_file, 'w', compression=zipfile.ZIP_DEFLATED) as fh:
            fh.write(genes_file, arcname='genes.txt')
            fh.write(isoforms_file, arcname='isoforms.txt')

        dset.genes_file = zip_file
        shutil.rmtree(outdir)
        return True
    else:
        shutil.rmtree(outdir)
        dset.status = ERRORS['rsem']
        return False


def move(src, dest):
    if os.path.isdir(dest):
        dest = os.path.join(dest, os.path.basename(src))

    try:
        os.unlink(dest)
    except OSError:
        pass
    finally:
        shutil.move(src, dest)


def is_exec(cmd):
    return shutil.which(cmd) is not None


