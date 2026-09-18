"""Random reference sequence generation for baselines."""

import random


def randomRefSeq(sequence, seed):
    random.seed(seed)
    sequence_list = list(sequence)
    random.shuffle(sequence_list)
    shuffled_sequence = ''.join(sequence_list)
    return shuffled_sequence
