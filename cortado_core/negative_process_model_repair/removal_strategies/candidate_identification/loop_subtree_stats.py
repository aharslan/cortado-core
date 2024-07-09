from dataclasses import dataclass


@dataclass
class LoopSubtreeStats:
    loop_labels_annotated_keep = []
    do_frequency_keep = 0
    redo_frequency_keep = 0
    loop_labels_annotated_remove = []
    do_frequency_remove = 0
    redo_frequency_remove = 0
    loop_repitition_per_trace_keep = []
    count_positive_variants_highest_loop_repetitions = {}
    count_positive_variants_greater_repetitions_than_negative_repetitions = {}
    remove_redundant_redo_pm = 100
    optional_redo_mandatory_pm = 100
    repeat_exactly_n_pm = 100
    repeat_atleast_n_pm = 100
    min_pm = 100
    min_pm_approach = ""
    redo_frequency_trace_remove = 0
    redo_frequency_trace_keep = 0

    def __init__(self):
        return

