"""Sequence encoding, FASTA I/O, peak detection, and data filters."""

import itertools
import numpy as np
from Bio import SeqIO
import os


def get_trids(length=5):
    nucle_com = []
    chars = ['A', 'C', 'G', 'U', 'N']
    for string in map(''.join, itertools.product(chars, repeat=length)):
        nucle_com.append(string)
    return nucle_com


KmerList = get_trids()


def get_4_nucleotide_composition(tris, seq):
    seq_len = len(seq)
    tri_feature = []
    k = len(tris[0])
    for x in range(len(seq) + 1 - k):
        kmer = seq[x:x + k]
        if kmer in tris:
            ind = tris.index(kmer)
            tri_feature.append(str(ind))
        else:
            print("kmer error")
    return tri_feature


def read_fasta_file(path, min_length=1, max_length=1000):
    all_records = SeqIO.parse(open(path), 'fasta')
    seqs = []
    names = []
    for seq_record in all_records:
        seq = str(seq_record.seq.upper()).replace('T', 'U')
        if 'N' not in seq and len(seq) >= min_length and len(seq) <= max_length:
            seqs.append(seq)
            names.append(seq_record.id)
    return seqs, names


def seq_length_distributon(seqs):
    seqLenList = [len(x) for x in seqs]
    seqLenList.sort()
    import matplotlib.pyplot as plt
    plt.figure()
    plt.scatter(range(len(seqLenList)), seqLenList, c='red')
    plt.xlabel('The number of sequence')
    plt.ylabel('Length(bp)')
    plt.legend()
    plt.show()


def onehot_encode1(inputs, vocab, dim=4):
    transformed = []
    mapper = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}
    for input in inputs:
        input = list(input)
        one_hot = np.zeros((len(input), len(vocab)))
        for i in range(len(input)):
            if input[i] in vocab:
                one_hot[i][mapper[input[i]]] = 1
            elif input[i] == 'N':
                pass
            else:
                print("Error: {message} not in vocab".format(message=input[i]))
                return False
        transformed.append(one_hot)
    transformed = np.asarray(transformed)
    return transformed


def onehot_encode_old(inputs, vocab, dim=4):
    transformed = []
    mapper = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}
    for sequence in inputs:
        one_hot = np.zeros((len(sequence), len(vocab)))
        for i, base in enumerate(sequence):
            if base in mapper:
                one_hot[i, mapper[base]] = 1
            elif base == 'N':
                pass
            else:
                print("Error: {message} not in vocab".format(message=base))
                return False
        transformed.append(one_hot)
    transformed = np.array(transformed)
    return transformed


def onehot_encode(inputs, vocab, dim=4):
    transformed = []
    mapper = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}

    for sequence in inputs:
        try:
            one_hot = np.zeros((len(sequence), len(vocab)))
            one_hot[np.arange(len(sequence)), [mapper[base] for base in sequence]] = 1
            transformed.append(one_hot)
        except KeyError as e:
            print(f"Warning: {e} not in vocab. Ignoring.", sequence)

    transformed = np.array(transformed)
    return transformed


base2int = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3, 'N': 4}


def sequence2int(sequence):
    return [base2int[base] for base in sequence]


def peakMark(signalList, signalMaxThr=3, signalSumThr=4, signalLen=1):
    subSignal = []
    peakStatistic = []
    for basePos in range(len(signalList)):
        if signalList[basePos] != 0:
            if basePos != (len(signalList) - 1):
                subSignal.append(signalList[basePos])
            else:
                subSignal.append(signalList[basePos])
                subSignalMax = max(subSignal)
                subSignalSum = sum(subSignal)
                peakStatistic.append([basePos - len(subSignal) + 1, basePos, subSignalMax, subSignalSum])
        elif (signalList[basePos] == 0) and len(subSignal) == 0:
            continue
        elif (signalList[basePos] == 0) and len(subSignal) != 0:
            subSignalMax = max(subSignal)
            subSignalSum = sum(subSignal)
            peakStatistic.append([basePos - len(subSignal), basePos - 1, subSignalMax, subSignalSum])
            subSignal = []

    textLabel = np.array([0] * len(signalList))
    peakCount = 0
    for peak in peakStatistic:
        if ((peak[1] - peak[0] + 1) >= signalLen and (peak[2] >= signalMaxThr or peak[3] >= signalSumThr)):
            textLabel[(peak[0]):peak[1] + 1] = 1
            peakCount += 1
    return (textLabel, peakCount)


def peakMarkBrief(signalList, signalMaxThr=3, signalSumThr=4, signalLen=1):
    subSignal = []
    peakStatistic = []
    for basePos in range(len(signalList)):
        if signalList[basePos] != 0:
            if basePos != (len(signalList) - 1):
                subSignal.append(signalList[basePos])
            else:
                subSignal.append(signalList[basePos])
                subSignalMax = max(subSignal)
                subSignalSum = sum(subSignal)
                peakStatistic.append([basePos - len(subSignal) + 1, basePos, subSignalMax, subSignalSum])
        elif (signalList[basePos] == 0) and len(subSignal) == 0:
            continue
        elif (signalList[basePos] == 0) and len(subSignal) != 0:
            subSignalMax = max(subSignal)
            subSignalSum = sum(subSignal)
            peakStatistic.append([basePos - len(subSignal), basePos - 1, subSignalMax, subSignalSum])
            subSignal = []

    textLabel = np.array([0] * len(signalList))
    peakCount = 0
    for peak in peakStatistic:
        if ((peak[1] - peak[0] + 1) >= signalLen and (peak[2] >= signalMaxThr or peak[3] >= signalSumThr)):
            textLabel[(peak[0]):peak[1] + 1] = 1
            peakCount += 1
    return textLabel


def find_peaks(signal_list, signal_max_threshold=2, signal_sum_threshold=4, signal_length=1):
    text_label = np.zeros(len(signal_list), dtype=int)
    peak_count = 0

    sub_signal_start = None
    sub_signal_max = 0
    sub_signal_sum = 0

    for base_pos, signal_value in enumerate(signal_list):
        if signal_value != 0:
            if sub_signal_start is None:
                sub_signal_start = base_pos
            sub_signal_max = max(sub_signal_max, signal_value)
            sub_signal_sum += signal_value
        elif sub_signal_start is not None:
            sub_signal_end = base_pos - 1
            if (sub_signal_end - sub_signal_start + 1) >= signal_length and (sub_signal_max >= signal_max_threshold or sub_signal_sum >= signal_sum_threshold):
                text_label[sub_signal_start:sub_signal_end + 1] = 1
                peak_count += 1
            sub_signal_start = None
            sub_signal_max = 0
            sub_signal_sum = 0

    if sub_signal_start is not None:
        sub_signal_end = len(signal_list) - 1
        if (sub_signal_end - sub_signal_start + 1) >= signal_length and (sub_signal_max >= signal_max_threshold or sub_signal_sum >= signal_sum_threshold):
            text_label[sub_signal_start:sub_signal_end + 1] = 1
            peak_count += 1

    return text_label, peak_count


def get_peak_signal(lst, start_index):
    result = []
    while start_index < len(lst) and lst[start_index] == 0:
        return result
    while start_index < len(lst) and lst[start_index] != 0:
        result.append(lst[start_index])
        start_index += 1
    return result


def FCSumFilter(encode_peak_DF, libraryRatio, treat_signal='treat_signal', ctl_signal='ctl_signal', FCThr=1, sumThr=1, FCThrOnly=True):
    encode_peak_DF[treat_signal] = encode_peak_DF[treat_signal].apply(eval)
    encode_peak_DF[ctl_signal] = encode_peak_DF[ctl_signal].apply(eval)

    def calculate_sum(lst):
        return sum(lst)

    def calculate_max(lst):
        return max(lst)

    encode_peak_DF['sum_treat_signal'] = (encode_peak_DF[treat_signal].apply(calculate_sum) + 1)
    encode_peak_DF['sum_ctl_signal'] = (encode_peak_DF[ctl_signal].apply(calculate_sum) + 1) * libraryRatio

    encode_peak_DF['max_treat_signal'] = (encode_peak_DF[treat_signal].apply(calculate_max) + 1)
    encode_peak_DF['max_ctl_signal'] = (encode_peak_DF[ctl_signal].apply(calculate_max) + 1) * libraryRatio

    if (FCThrOnly):
        encode_peak_DF_filter = encode_peak_DF[(encode_peak_DF['max_treat_signal'] / encode_peak_DF['max_ctl_signal'] >= FCThr)]
    else:
        encode_peak_DF_filter = encode_peak_DF[(encode_peak_DF['max_treat_signal'] / encode_peak_DF['max_ctl_signal'] >= FCThr) | (encode_peak_DF['sum_treat_signal'] / encode_peak_DF['sum_ctl_signal'] >= sumThr)]
    return encode_peak_DF_filter


def create_folder_if_not_exists(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Folder '{folder_path}' created.")
    else:
        print(f"Folder '{folder_path}' already exists.")
