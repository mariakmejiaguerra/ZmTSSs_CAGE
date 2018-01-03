#!/usr/bin/env python3

# ------------------------------------------------------------------------- 
# @author: Katherine Mejia-Guerra (mejia-guerra.1@osu.edu)
# Copyright (C) 2015 Katherine Mejia-Guerra
# -------------------------------------------------------------------------

import re
import sys
import time
import os
import subprocess
import logging
import dinopy
from collections import OrderedDict

def trimmed_filtered_tags(ids, adapter, reads_lst, sequence_qthreshold, sequence_lthreshold):
    import HTSeq
    reads_length = []
    good_reads = OrderedDict()
    pass_quality = 0; has_adapter = 0; pass_length = 0  
    for idread, dinopy_read in enumerate(reads_lst):
        dseq = dinopy_read[0].encode('utf-8')
        dname = dinopy_read[1]
        dquals = dinopy_read[2]
        allread = HTSeq.SequenceWithQualities(dseq, dname.decode("utf-8"), dquals)
        if mean_quals(allread.qual) > sequence_qthreshold and len(dinopy_read[0]) >= sequence_lthreshold:
            pass_quality += 1
            pos_tag = search_adapter(adapter, dseq)
            if len(pos_tag) == 2:
                has_adapter += 1
                cage_tag_seq = dseq[pos_tag[0]:pos_tag[1]] 
                cage_tag_quals = dquals[pos_tag[0]:pos_tag[1]]
                old_length = "length="+str(len(dseq))
                if len(cage_tag_seq) >= 27:
                    pass_length += 1
                    new_length = "length="+str(len(cage_tag_seq))
                    str_name  = dname.decode("utf-8").replace(old_length, new_length)
                    cutread = HTSeq.SequenceWithQualities(cage_tag_seq, str_name, cage_tag_quals)
                    if mean_quals(cutread.qual) > sequence_qthreshold and adapter.strip(".") not in cage_tag_seq.decode("utf-8"):
                        reads_length.append(len(cage_tag_seq))
                        good_reads[idread] = (cage_tag_seq.decode("utf-8"), str_name.encode('utf-8'), cage_tag_quals)
    trimmed_filtered_tags = list(good_reads.values())
    stats = [pass_quality, has_adapter, pass_length, mean_quals(reads_length)]
    return trimmed_filtered_tags, stats
    
def reads_any_adapter(lst_adapters, seq):
    sequence_ids = []
    for idx, adapter in enumerate(lst_adapters):
        matcher = re.compile('|'.join([adapter, adapter.translate(str.maketrans('TAGC', 'ATCG'))[::-1]]))
        for match in matcher.finditer(seq.decode("utf-8")):
            if len(match.group()) == len(adapter) :
                sequence_ids.append({"sample":idx})
    if len(sequence_ids) == 1:
        return sequence_ids[0]
    else:
        return dict()
    
def sample_adapter_index(filename):
    import csv
    reader = csv.reader(open(filename, 'r'))
    d = OrderedDict()
    for row in reader:
        sra, sample_id, adapter = row
        d[sra] = str(adapter).replace("N", ".")
    return d

def search_adapter(adapter, seq):
    tag_pos = []
    matcher = re.compile('|'.join([adapter, adapter.translate(str.maketrans('TAGC', 'ATCG'))[::-1]]))
    for match in matcher.finditer(seq.decode("utf-8")):
        if len(match.group()) >= 1 :
            strt_tag = int(match.end())
            end_tag = int(match.end())+28
            tag_pos = [strt_tag, end_tag]
    return tag_pos

def mean_quals(array):
    import warnings
    import numpy as np
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.mean(array)

if sys.argv[1] == "-help":
    print("Usage: python preprocessing.py [fastq file] [no header csv file with SRA number, samples, barcodes]")
    print("Example: python preprocessing.py multiplexed.fastq sample_barcode.csv")
    quit()
else:
    run_id = str(int(time.time()))
    reads_file = str(sys.argv[1])
    csv_file = str(sys.argv[2])
    pathname = os.path.dirname(sys.argv[0])
    write_demultiplex = False #should be arguments with flags
    WORKING_DIR = os.path.abspath(pathname)
    file_name = 'logging_preprocessing_reads_'+ run_id + '.txt'
    logging.basicConfig(level=logging.DEBUG, 
                        filename=file_name, 
                        filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info(run_id)
    logging.info("WORKING_DIR")
    logging.info(WORKING_DIR)
    logging.info("input: multiplex fastq file")
    logging.info(reads_file)
    logging.info("input: csv file sample:barcode")
    logging.info(csv_file)


lthreshold = 47 #should be arguments with flags
qthreshold = 39 #should be arguments with flags

logging.info("input: write_demultiplex")
logging.info(write_demultiplex)
logging.info("input: sequence length threshold")
logging.info(lthreshold)
logging.info("input: sequence quality threshold")
logging.info(qthreshold)

in_path_reads = os.path.join(WORKING_DIR,reads_file)
in_path_csv = os.path.join(WORKING_DIR,csv_file)

dict_sample_adapter = sample_adapter_index(in_path_csv) #sample,barcode (NO HEADERS)
samples_ids = list(dict_sample_adapter.keys())
in_adapters = list(dict_sample_adapter.values())

logging.info("input: total samples")
logging.info(len(samples_ids))
logging.info("input: total barcodes")
logging.info(len(in_adapters))

#demultiplexing fastq files
fqr = dinopy.FastqReader(in_path_reads)
reads = dict()
for seq, name, quals in fqr.reads(quality_values=True):
    adapterid = reads_any_adapter(in_adapters, seq)
    if adapterid:
        curread = (seq.decode("utf-8"), name, quals)
        sampleid = samples_ids[adapterid.get("sample")]
        if sampleid in reads:
            sample_reads = reads.get(sampleid)
            sample_reads.append(curread)
        else:
            sample_reads = []
            sample_reads.append(curread)
            reads[sampleid] = sample_reads

#output demultiplexed fastq files - intermediate files can be of use for QC analysis or in a different pipeline
if write_demultiplex:
    for ids in samples_ids:
        output_file_reads = ids + "_demultiplexed.fastq"
        out_path_reads = os.path.join(WORKING_DIR,output_file_reads)
        with dinopy.FastqWriter(out_path_reads) as fqw:
            fqw.write_reads(reads.get(ids), dtype=str)
        logging.info("output demultiplexed file")
        logging.info(out_path_reads)
            
#output trimmed and filtered fastq files - intermediate files can be of use for QC analysis or in a different pipeline
for idx, ids in enumerate(samples_ids):
    adapter = in_adapters[idx]
    demult_reads_list = reads.get(ids)
    print("Preprocessing %d reads for sample %s" % (len(demult_reads_list), ids) )
    print("Removing adapter %s" % (adapter.strip(".")) )
    cage_tags, stats = trimmed_filtered_tags(ids, adapter, demult_reads_list, qthreshold, lthreshold)
    print("Obtained %d tags for sample %s" % (len(cage_tags), ids) )
    print("Reads pass_quality: %d Tags after remove adapters: %d Tags after check length: %d Avg length: %f" % tuple(stats)  )
    logging.info("Preprocessing reads for sample %s with adapter %s " % (ids, adapter))
    logging.info("Total reads")
    logging.info(len(demult_reads_list))
    logging.info("Total tags")
    logging.info(len(cage_tags))
    logging.info("Avg length")
    logging.info(stats[3])
    output_file_tags = ids + "_demultiplexed_trimmed_filtered.fastq"
    out_path_tags =  os.path.join(WORKING_DIR,output_file_tags)
    with dinopy.FastqWriter(out_path_tags) as fqw:
        fqw.write_reads(cage_tags, dtype=str)
    logging.info("output demultiplexed/trimmed/filtered file")
    logging.info(out_path_tags)   
    print()

print("Output files are ready")
print('NEXT step: align using bowtie2')
completed = subprocess.run(['ls', '-1'])
print('returncode:', completed.returncode)
