"""Parse numeric text with the same nearest-float contract as Python float().

This applies uniformly to original external measurements and derived CSVs.
It does not round values or merge near-equal measurements into biological ties.
"""
import pandas as pd


def read_csv_literal(*args, **kwargs):
    kwargs.setdefault('float_precision', 'round_trip')
    return pd.read_csv(*args, **kwargs)
