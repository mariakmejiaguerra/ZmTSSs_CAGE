#!/usr/bin/env python3

# -------------------------------------------------------------------------
# @author: Katherine Mejia-Guerra (mejia-guerra.1@osu.edu)
# Copyright (C) 2015 Katherine Mejia-Guerra
# -------------------------------------------------------------------------

import re
import sys
import os
import time
import warnings
import logging
import dinopy
import numpy as np
import HTSeq
from collections import OrderedDict

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
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.mean(array)


if sys.argv[1] == "-help":
    print("Usage: python barcode_sequence_trimmer.py input_file_prefix")
    print("Example: python barcode_sequence_trimmer.py SRR2078288")
    quit()
else:
    run_id = str(int(time.time()))
    input_id = str(sys.argv[1])
    pathname = os.path.dirname(sys.argv[0])
    output_file = input_id + "_demultiplexed_trimmed_filtered.fastq"
    file_name = 'logging_barcode_sequence_trimmer_' + input_id + '_' + run_id + '.txt'
    logging.basicConfig(level=logging.DEBUG, filename=file_name, filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info(run_id)
    

WORKING_DIR = os.path.abspath(pathname)
logging.info("WORKING_DIR")
logging.info(WORKING_DIR)
logging.info("input: SRA fastq file")
logging.info(input_id)
logging.info("output: trimmed file")
logging.info(output_file)

input_path = input_id + ".fastq"

out_path_cbsu = os.path.join(WORKING_DIR,output_file)
fqr = dinopy.FastqReader(input_path)

good_reads = OrderedDict()
reads_length = []
pass_quality = 0
pass_length = 0
has_adapter = 0

sample2inadapter = {"SRR2078285":"GATCAGCAG",
                    "SRR2078286":"ACACAGCAG",
                    "SRR2078287":"ACTCAGCAG",
                    "SRR2078288":"ACGCAGCAG",
                    "SRR2078289":"AGACAGCAG",
                    "SRR2078290":"ATCCAGCAG",
                    "SRR2078291":"ATGCAGCAG",
                    "SRR2078292":"CTTCAGCAG"}

inadapter = sample2inadapter.get(input_id)
adapter = inadapter+"......"
sequence_lthreshold = int(len(adapter)) + 32
sequence_qthreshold = 39

idread = 0
for seq, name, quals in fqr.reads(quality_values=True):
    allread = HTSeq.SequenceWithQualities(seq, name.decode("utf-8"), quals)
    if mean_quals(allread.qual) > sequence_qthreshold and len(seq) >= sequence_lthreshold:
        pass_quality += 1
        pos_tag = search_adapter(adapter, seq)
        if len(pos_tag) == 2:
            has_adapter += 1
            cage_tag_seq = seq[pos_tag[0]:pos_tag[1]] 
            cage_tag_quals = quals[pos_tag[0]:pos_tag[1]]
            old_length = "length="+str(len(seq))
            if len(cage_tag_seq) >= 27:
                pass_length += 1
                new_length = "length="+str(len(cage_tag_seq))
                str_name  = name.decode("utf-8").replace(old_length, new_length)
                cutread = HTSeq.SequenceWithQualities(cage_tag_seq, str_name, cage_tag_quals)
                reads_length.append(len(cage_tag_seq))
                if mean_quals(cutread.qual) > sequence_qthreshold and inadapter not in cage_tag_seq.decode("utf-8"):
                    good_reads[idread] = (cage_tag_seq.decode("utf-8"), str_name.encode('utf-8'), cage_tag_quals)
    idread += 1

with dinopy.FastqWriter(out_path_cbsu) as fqw:
    fqw.write_reads(list(good_reads.values()), dtype=str)

print("Output file is ready")
print("NEXT step: align "+output_file+" using bowtie2")
logging.info("output: good read")
logging.info(pass_quality)  
logging.info("output: has adapter")
logging.info(has_adapter)    
logging.info("output: pass length")
logging.info(pass_length)
logging.info("output: pass quality")
logging.info(len(good_reads))    
logging.info("output: average read length")
logging.info(np.mean(reads_length))