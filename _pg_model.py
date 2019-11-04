#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Matthias Blum
# Contact: mat.blum@gmail.com

import json
import os
import sys
import urllib.parse
import urllib.request

DT_CHIP = 'ChIP'
DT_HIC = 'HiC'
DT_RNA = 'RNA'


class _Run:
    def __init__(self, dset_id, assembly, singleread, data_type, **kwargs):
        self.dset_id = dset_id
        self.assembly = assembly
        self.singleread = singleread
        self.data_type = data_type
        self.id = kwargs['id']
        self.url = kwargs['url']

        # Files
        self.sra_file = None
        self.fq1_file = None
        self.fq2_file = None
        self.reads_file = None
        # Secondary alignments (to transcriptome for RNA-seq)
        self.reads_file_2 = None

        # Stats
        self.reads = None
        self.mapped_reads = None

        # Flags
        self.idle = True
        self.is_downloaded = False
        self.is_extracted = False
        self.is_aligned = False

        # Status
        self.status = 0

    def destroy(self):
        for f in (self.sra_file, self.fq1_file, self.fq2_file, self.reads_file, self.reads_file_2):
            try:
                os.unlink(f)
            except (OSError, TypeError):
                pass

    def todict(self, as_json=False):
        _dict = {
            'id': self.id,
            'reads': self.reads,
            'mapped_reads': self.mapped_reads
        }

        if as_json:
            return json.dumps(_dict)
        else:
            return _dict

    def update(self, run):
        self.sra_file = run.sra_file
        self.fq1_file = run.fq1_file
        self.fq2_file = run.fq2_file
        self.reads_file = run.reads_file
        self.reads_file_2 = run.reads_file_2
        self.reads = run.reads
        self.mapped_reads = run.mapped_reads
        self.is_downloaded = run.is_downloaded
        self.is_extracted = run.is_extracted
        self.is_aligned = run.is_aligned
        self.status = run.status


class _Dataset:
    def __init__(self, **kwargs):
        self.id = kwargs['id']
        self.sample_id = kwargs['sample']
        self.exp_id = kwargs['exp']
        self.assembly = kwargs['assembly']
        self.singleread = kwargs['singleread']
        self.data_type = kwargs['datatype']
        assert self.data_type in (DT_CHIP, DT_HIC, DT_RNA)

        # Compressed, sorted, BED file
        self.reads_file = None

        # BAM file (alignments to transcriptome, RNA-seq only)
        self.reads_file_2 = None

        # ZIP file containing results from RSEM (RNA-seq only)
        self.genes_file = None

        # Flags
        self.idle = True
        self.all_runs_ready = False
        self.is_merged = False
        self.is_analyzed = False
        self.runs = [_Run(self.id, self.assembly, self.singleread, self.data_type, **r) for r in kwargs['runs']]
        self.status = 0

    def todict(self, as_json=False):
        _dict = {
            'id': self.id,
            'sample': self.sample_id,
            'exp': self.exp_id,
            'assembly': self.assembly,
            'singleread': self.singleread,
            'reads_file': self.reads_file,
            'genes_file': self.genes_file,
            'status': self.status,
            'runs': [r.todict() for r in self.runs]
        }

        if as_json:
            return json.dumps(_dict)
        else:
            return _dict

    def update(self, dset):
        self.reads_file = dset.reads_file
        self.reads_file_2 = dset.reads_file_2
        self.genes_file = dset.genes_file
        self.is_merged = dset.is_merged
        self.is_analyzed = dset.is_analyzed
        self.status = dset.status


class APIConnection:
    def __init__(self, host, port, secure_key):
        self.host = host
        self.port = port
        self.secure_key = secure_key

    def get(self, data_types=list(), assemblies=list(), n=1):
        url = 'http://{}:{}/get'.format(self.host, self.port)

        data = urllib.parse.urlencode({
            'key': self.secure_key,
            'datatype': data_types,
            'assembly': assemblies,
            'n': n
        }, doseq=True).encode('utf8')

        try:
            response = urllib.request.urlopen(url, data)
        except (urllib.request.HTTPError, urllib.request.URLError):
            sys.stderr.write('API server unreachable\n')
            return []
        else:
            _obj = json.loads(response.read().decode('utf8'))
            if _obj.get('error'):
                sys.stderr.write('API: {}\n'.format(_obj['error']))
                return []
            else:
                return [_Dataset(**_dset) for _dset in _obj['datasets']]

    def set(self, dsets):
        url = 'http://{}:{}/set'.format(self.host, self.port)

        data = urllib.parse.urlencode({
            'key': self.secure_key,
            'datasets': json.dumps([dset.todict() for dset in dsets])
        }).encode('utf8')

        try:
            urllib.request.urlopen(url, data)
        except (urllib.request.HTTPError, urllib.request.URLError):
            sys.stderr.write('API server unreachable')


def load_input(filename):
    with open(filename, 'rt') as fh:
        _dict = json.load(fh)

    return [_Dataset(**_dset) for _dset in _dict['datasets']]
