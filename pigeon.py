#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Matthias Blum
# Contact: mat.blum@gmail.com

import argparse
import configparser
import logging
import os
import sys
import tempfile
import time

import _pg_pool
import _pg_model
import _pg_task

logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%y-%m-%d %H:%M:%S')

DATA_TYPES = ('ChIP', 'HiC', 'RNA')


def main():
    parser = argparse.ArgumentParser(description='PIGEON: a pipeline for GEO')
    parser.add_argument('-i', '--input', help='input JSON file')
    parser.add_argument('-c', '--config', help='configuration file', required=True)
    parser.add_argument('-d', '--datatypes', help='data types to handle', nargs='+', choices=DATA_TYPES, required=True)
    parser.add_argument('-a', '--assemblies', help='assemblies to handle', nargs='+', default=list())

    # Workers
    parser.add_argument('--download', help='number of download workers', type=int, default=1)
    parser.add_argument('--extract', help='number of extraction workers', type=int, default=1)
    parser.add_argument('--align', help='number of alignment workers', type=int, default=1)
    parser.add_argument('--merge', help='number of merging workers', type=int, default=1)
    parser.add_argument('--analyze', help='number of merging workers', type=int, default=1)

    # Systems
    parser.add_argument('-m', '--maxmem', help='memory (Mb) for sort', type=int)
    parser.add_argument('--threads-aln', dest='threads_aln',
                        help='number of threads for alignment (default: 1)', type=int, default=1)
    parser.add_argument('--threads-alz', dest='threads_alz',
                        help='number of threads for analyze (default: 1)', type=int, default=1)

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        sys.stderr.write('--config: {}: no such file or directory\n'.format(args.config))
        exit(1)
    elif args.input and not os.path.isfile(args.input):
        sys.stderr.write('--input: {}: no such file or directory\n'.format(args.input))
        exit(1)

    config = configparser.RawConfigParser()
    config.read(args.config)

    # Remove duplicates
    args.datatypes = list(set(args.datatypes))

    _bowtie2 = config.get('tools', 'bowtie2', fallback='bowtie2')
    _fastqdump = config.get('tools', 'fastq-dump', fallback='fastq-dump')
    _star = config.get('tools', 'STAR', fallback='STAR')
    _rsem = config.get('tools', 'RSEM_calc_exp', fallback='rsem-calculate-expression')
    _samtools = config.get('tools', 'samtools', fallback='samtools')

    for tool, path in zip(
            ['bowtie2', 'fastq-dump', 'STAR', 'RSEM', 'samtools', 'sam2bed', 'joinbed'],
            [_bowtie2, _fastqdump, _star, _rsem, _samtools, _pg_task.SAM2BED, _pg_task.JOINBED]):
        if not _pg_task.is_exec(path):
            sys.stderr.write('Missing dependency: {}\n'.format(tool))
            if tool in ('fastq-dump', 'samtools', 'sam2bed'):
                exit(1)
            elif tool == 'bowtie2':
                for dt in ('ChIP', 'HiC'):
                    if dt in args.datatypes:
                        args.datatypes.remove(dt)
            elif tool in ('STAR', 'RSEM'):
                try:
                    args.datatypes.remove('RNA')
                except ValueError:
                    pass
            elif tool == 'joinbed' and 'HiC' in args.datatypes:
                args.datatypes.remove('HiC')

    if not args.datatypes:
        sys.stderr.write('Cannot continue: please install missing dependencies\n')
        exit(1)

    workdir = os.path.abspath(config.get('paths', 'workdir'))
    genes_output_dir = os.path.abspath(config.get('paths', 'genes'))
    pairs_output_dir = os.path.abspath(config.get('paths', 'pairs'))
    reads_output_dir = os.path.abspath(config.get('paths', 'reads'))

    for _dir in (workdir, genes_output_dir, pairs_output_dir, reads_output_dir):
        if not os.path.isdir(_dir):
            os.makedirs(_dir)

    fd, lock_file = tempfile.mkstemp(dir=workdir, prefix='pigeon_', suffix='.lock')
    os.close(fd)
    sys.stderr.write('To stop Pigeon, delete the file: {}\n'.format(lock_file))

    dl_pool = _pg_pool.DownloadPool(workdir,
                                    processes=args.download)

    ex_pool = _pg_pool.ExtractPool(workdir,
                                   processes=args.extract,
                                   fastqdump=_fastqdump)

    aln_pool = _pg_pool.AlignPool(workdir,
                                  processes=args.align,
                                  bowtie2_ref_dir=config.get('alignment', 'bowtie2_ref'),
                                  star_ref_dir=config.get('alignment', 'STAR_ref'),
                                  star_gtf_dir=config.get('alignment', 'STAR_gtf'),
                                  bowtie2=_bowtie2,
                                  bowtie2_flags=config.get('alignment', 'bowtie2_flags'),
                                  star=_star,
                                  star_flags=config.get('alignment', 'STAR_flags'),
                                  mapq=config.getint('alignment', 'MAPQ'),
                                  threads=args.threads_aln)

    merge_pool = _pg_pool.MergePool(workdir,
                                    processes=args.merge,
                                    maxmem=args.maxmem,
                                    ratio=config.getfloat('alignment', 'ratio', fallback=0),
                                    samtools=_samtools)

    analyze_pool = _pg_pool.AnalyzePool(workdir,
                                        processes=args.analyze,
                                        rsem=_rsem,
                                        star_ref_dir=config.get('alignment', 'STAR_ref'),
                                        threads=args.threads_alz)

    if args.input:
        api = None
        dsets_in_queue = _pg_model.load_input(args.input)
    else:
        api = _pg_model.APIConnection(
            host=config.get('api', 'host'),
            port=config.getint('api', 'port'),
            secure_key=config['api']['secure_key']
        )
        dsets_in_queue = api.get(data_types=args.datatypes, assemblies=args.assemblies)

    """
    Counters: in order to not start to much downloads (which would fill the disk with SRA files),
    we need to keep track of how many workers are available
    """
    runs_wait_download = []
    runs_download = []
    runs_wait_extract = []
    runs_extract = []
    runs_wait_align = []
    runs_align = []
    dsets_wait_merge = []
    dsets_merge = []
    dsets_wait_analyze = []
    dsets_analyze = []

    while True:
        # Get completed jobs
        dled_runs = dl_pool.get()
        extracted_runs = ex_pool.get()
        aligned_runs = aln_pool.get()
        merged_dsets = merge_pool.get()
        analyzed_dsets = analyze_pool.get()

        # Update counters
        runs_wait_download = []
        for r in dled_runs.values():
            runs_download.remove(r.id)

            if r.status == 0:
                runs_wait_extract.append(r.id)

        for r in extracted_runs.values():
            runs_extract.remove(r.id)

            if r.status == 0:
                runs_wait_align.append(r.id)

        for r in aligned_runs.values():
            runs_align.remove(r.id)

        for dset in merged_dsets.values():
            dsets_merge.remove(dset.id)

            if dset.status == 0:
                dsets_wait_analyze.append(dset.id)

        for dset in analyzed_dsets.values():
            dsets_analyze.remove(dset.id)

        _dsets_in_queue = []
        dsets_out_queue = []

        for i, dset in enumerate(dsets_in_queue):
            for _queue in (merged_dsets, analyzed_dsets):
                if dset.id in _queue:
                    dset = _queue[dset.id]
                    dset.idle = True
                    dsets_in_queue[i] = dset
                    break

            if not dset.all_runs_ready:
                nruns_ready = 0

                for j, run in enumerate(dset.runs):
                    for _queue in (dled_runs, extracted_runs, aligned_runs):
                        if run.id in _queue:
                            run = _queue[run.id]
                            run.idle = True
                            dsets_in_queue[i].runs[j] = run
                            break

                    if run.status:
                        dset.status = run.status
                        nruns_ready += 1  # run is not really "ready" but won't make the dset wait any more
                    elif run.idle:
                        if dset.status:
                            """
                            Another run from the dataset failed:
                            there is not point to still process this run, so we have to
                                * delete its files
                                * remove it from waiting queues
                            """
                            try:
                                runs_wait_download.remove(run.id)
                            except ValueError:
                                pass

                            try:
                                runs_wait_extract.remove(run.id)
                            except ValueError:
                                pass

                            try:
                                runs_wait_align.remove(run.id)
                            except ValueError:
                                pass

                            run.destroy()
                            nruns_ready += 1
                        elif not run.is_downloaded:
                            if len(runs_download) < args.download and not runs_wait_extract:
                                dl_pool.put(run)
                                run.idle = False
                                runs_download.append(run.id)
                            else:
                                runs_wait_download.append(run.id)
                        elif not run.is_extracted:
                            if len(runs_extract) < args.extract and not runs_wait_align:
                                ex_pool.put(run)
                                run.idle = False
                                runs_wait_extract.remove(run.id)
                                runs_extract.append(run.id)
                        elif not run.is_aligned:
                            if len(runs_align) < args.align and len(dsets_wait_merge) <= 1:
                                aln_pool.put(run)
                                run.idle = False
                                runs_wait_align.remove(run.id)
                                runs_align.append(run.id)
                        else:
                            nruns_ready += 1

                if nruns_ready == len(dset.runs):
                    dset.all_runs_ready = True

            if dset.idle and dset.all_runs_ready:
                if dset.status:
                    dsets_out_queue.append(dset)
                    continue
                elif not dset.is_merged:
                    if len(dsets_merge) < args.merge:
                        merge_pool.put(dset)
                        dset.idle = False
                        dsets_merge.append(dset.id)
                        try:
                            dsets_wait_merge.remove(dset.id)
                        except ValueError:
                            # dset.id might not be in dsets_wait_merge
                            pass
                    elif dset.id not in dsets_wait_merge:
                        dsets_wait_merge.append(dset.id)
                elif not dset.is_analyzed:
                    if len(dsets_analyze) < args.analyze:
                        analyze_pool.put(dset)
                        dset.idle = False
                        dsets_wait_analyze.remove(dset.id)
                        dsets_analyze.append(dset.id)
                elif dset.is_analyzed:
                    # Move files
                    if dset.data_type == _pg_model.DT_CHIP:
                        _pg_task.move(dset.reads_file, reads_output_dir)
                    elif dset.data_type == _pg_model.DT_HIC:
                        _pg_task.move(dset.reads_file, pairs_output_dir)
                    elif dset.data_type == _pg_model.DT_RNA:
                        _pg_task.move(dset.genes_file, genes_output_dir)
                        _pg_task.move(dset.reads_file, reads_output_dir)

                    dsets_out_queue.append(dset)
                    continue

            _dsets_in_queue.append(dset)

        dsets_in_queue = _dsets_in_queue

        for dset in dsets_out_queue:
            print(dset.todict(as_json=True))

        if api:
            if os.path.isfile(lock_file) and not runs_wait_download:
                new_dsets = api.get(
                    data_types=args.datatypes,
                    assemblies=args.assemblies
                )

                for dset in new_dsets:
                    logging.info('Adding data set {}'.format(dset.exp_id))
                    dsets_in_queue.append(dset)

            if dsets_out_queue:
                api.set(dsets_out_queue)

        if dsets_in_queue:
            time.sleep(10)

            with open(os.path.join(workdir, 'stats.txt'), 'wt') as fh:
                lines = [
                    'datasets:           {}'.format(len(dsets_in_queue)),
                    'waiting downloads:  {}'.format(len(runs_wait_download)),
                    'downloading:        {}'.format(len(runs_download)),
                    'waiting extraction: {}'.format(len(runs_wait_extract)),
                    'extracting:         {}'.format(len(runs_extract)),
                    'waiting alignment:  {}'.format(len(runs_wait_align)),
                    'aligning:           {}'.format(len(runs_align)),
                    'waiting merge:      {}'.format(len(dsets_wait_merge)),
                    'merging:            {}'.format(len(dsets_merge)),
                    'waiting analyze:    {}'.format(len(dsets_wait_analyze)),
                    'analyzing:          {}'.format(len(dsets_analyze))
                ]
                fh.write('\n'.join(lines) + '\n')
        else:
            break


if __name__ == '__main__':
    main()
