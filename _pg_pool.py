#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Matthias Blum
# Contact: mat.blum@gmail.com

import logging
import os
import queue
import time
from multiprocessing import Pool, Queue

import _pg_task
from _pg_model import DT_CHIP, DT_HIC, DT_RNA

logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%y-%m-%d %H:%M:%S')


class _Pool:
    def __init__(self, workdir, processes=1, **kwargs):
        self.workdir = workdir
        self.kwargs = kwargs

        if not os.path.isdir(self.workdir):
            try:
                os.makedirs(self.workdir)
            except OSError:
                pass

        self._queue_in = Queue()
        self._queue_out = Queue()
        self._pool = Pool(processes, self.init_worker)

    def __del__(self):
        self.__exit__()

    def __exit__(self):
        self._pool.terminate()

    def put(self, obj):
        self._queue_in.put(obj, False)

    def get(self):
        objects = {}

        while True:
            try:
                obj = self._queue_out.get(False)
            except queue.Empty:
                break
            else:
                objects[obj.id] = obj

        return objects

    def init_worker(self):
        while True:
            try:
                obj = self._queue_in.get(False)
            except queue.Empty:
                pass
            else:
                self.process(obj)
                self._queue_out.put(obj)

            time.sleep(10)

    def process(self, obj):
        raise NotImplementedError


class DownloadPool(_Pool):
    def process(self, run):
        logging.info('downloading {}'.format(os.path.basename(run.url[:-4])))
        _pg_task.download_run(run, self.workdir)
        logging.info('{} downloaded (status: {})'.format(os.path.basename(run.url[:-4]), run.status))


class ExtractPool(_Pool):
    def process(self, run):
        logging.info('extracting {}'.format(os.path.basename(run.url[:-4])))
        _pg_task.extract_sra(run, self.workdir, **self.kwargs)
        logging.info('{} extracted (status: {})'.format(os.path.basename(run.url[:-4]), run.status))


class AlignPool(_Pool):
    def __init__(self, workdir, bowtie2_ref_dir, star_ref_dir, star_gtf_dir, processes=1, **kwargs):
        self.bowtie2_ref_dir = bowtie2_ref_dir
        self.star_ref_dir = star_ref_dir
        self.star_gtf_dir = star_gtf_dir
        super().__init__(workdir, processes, **kwargs)

    def process(self, run):
        logging.info('aligning {}'.format(os.path.basename(run.url[:-4])))

        if run.data_type == DT_CHIP:
            _pg_task.align_chip(run, self.workdir, self.bowtie2_ref_dir, **self.kwargs)
        elif run.data_type == DT_HIC:
            _pg_task.align_hic(run, self.workdir, self.bowtie2_ref_dir, **self.kwargs)
        elif run.data_type == DT_RNA:
            _pg_task.align_rna(run, self.workdir, self.star_ref_dir, self.star_gtf_dir, **self.kwargs)

        logging.info('{} aligned (status: {})'.format(os.path.basename(run.url[:-4]), run.status))


class MergePool(_Pool):
    def __init__(self, workdir, ratio=0, processes=1, **kwargs):
        self.ratio = ratio
        super().__init__(workdir, processes, **kwargs)

    def process(self, dset):
        logging.info('merging {}'.format(dset.exp_id))
        _pg_task.merge_runs(dset, self.workdir, self.ratio, **self.kwargs)
        logging.info('{} merged (status: {})'.format(dset.exp_id, dset.status))


class AnalyzePool(_Pool):
    def __init__(self, workdir, star_ref_dir, processes=1, **kwargs):
        self.star_ref_dir = star_ref_dir
        super().__init__(workdir, processes, **kwargs)

    def process(self, dset):
        logging.info('analyzing {}'.format(dset.exp_id))
        if dset.data_type == DT_HIC:
            # LOGIQA
            pass
        else:
            # NGS-QC Generator

            if dset.data_type == DT_RNA:
                # RSEM
                if _pg_task.run_rsem(dset, self.workdir, self.star_ref_dir, **self.kwargs):
                    pass
                else:
                    os.unlink(dset.reads_file)

        dset.is_analyzed = True
        logging.info('{} analyzed (status: {})'.format(dset.exp_id, dset.status))
