class Constants:
    # proportion of top loop repetitions to encode
    MAX_NUMBER_OF_LOOP_REPETITION_ENCODINGS = 4 # absolute number
    MAX_LENGTH_LOOP_REPETITION_ENCODING = 5  # absolute number
    MIN_THRESHOLD_POSITIVE_FITTING_VARIANTS = 75.0 # percentage.
    MAX_THRESHOLD_TREE_EDIT_DISTANCE = 5 # threshold of max allowed edit distance.
    STOP_WHEN_AN_UPDATE_MEETS_THRESHOLD = True
    USE_SUBTREE_BASED_CONFORMANCE_CHECKING = True

    def __init__(self):
        return
