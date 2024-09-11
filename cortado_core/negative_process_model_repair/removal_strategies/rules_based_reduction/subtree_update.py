import copy
from trace import Trace

from pm4py.objects.log.obj import Event
from pm4py.objects.process_tree.obj import Operator, ProcessTree
from cortado_core.freezing.reinsert_frozen_subtrees import post_process_tree
from cortado_core.models.infix_type import InfixType
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_activity import \
    CandidateActivity
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_subtree import \
    CandidateSubtree
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.removal_candidates_heuristics import \
    RemovalCandidatesHeuristics
from cortado_core.negative_process_model_repair.removal_strategies.rules_based_reduction.update_rule import \
    SequenceUpdateRule, ChoiceUpdateRule
from cortado_core.utils.alignment_utils import typed_trace_fits_process_tree
from cortado_core.negative_process_model_repair.temp_utils import calculate_process_tree_edit_distance, \
    find_tree_node_by_id
from cortado_core.negative_process_model_repair.temp_utils import (
    calculate_percentage_traces_conforming
)
import itertools
from collections import Counter
from cortado_core.negative_process_model_repair.constants import Constants
from cortado_core.utils.trace import TypedTrace


class SubtreeUpdate:
    removal_candidates_generator: RemovalCandidatesHeuristics
    removal_candidate_activities: list[CandidateActivity]
    removal_candidate_subtrees: list[CandidateSubtree]

    update_operations = 0

    def __init__(
        self,
        removal_candidates_generator: RemovalCandidatesHeuristics,
        removal_candidate_activities: list[CandidateActivity],
    ):
        self.removal_candidates_generator = removal_candidates_generator
        self.removal_candidate_activities = removal_candidate_activities
        self.removal_candidate_subtrees = (
            self.removal_candidates_generator.generate_removal_candidate_subtrees(
                removal_candidate_activities
            )
        )

    def calculate_candidate_tree_statistics(self, process_tree: ProcessTree,
                                            removal_candidate_subtree: CandidateSubtree,
                                            applied_rule: str):

        statistics = {
            'negative_trace_fits': True,
            'percentage_positive_traces_conforming': -1,
            'resulting_tree_edit_distance': 10000,
            'updated_tree': process_tree,
            'candidate_subtree': removal_candidate_subtree,
            'applied_rule': applied_rule,
            'thresholds_met': False
        }

        try:
            if (
                process_tree != self.removal_candidates_generator.process_tree
                and process_tree is not None
            ):
                process_tree = post_process_tree(process_tree, [])

                statistics['updated_tree'] = process_tree

                parent_sub_tree = process_tree
                if Constants.USE_SUBTREE_BASED_CONFORMANCE_CHECKING == True:
                    if removal_candidate_subtree.reference.parent is None:
                        parent_sub_tree = process_tree
                    else:
                        found = False
                        reference = removal_candidate_subtree.reference
                        while not found:
                            if reference.parent is not None:

                                if not hasattr(reference, 'id'):
                                    reference = reference.parent
                                    continue

                                tree_node = find_tree_node_by_id(reference.id, process_tree)
                                if tree_node is not None:
                                    parent_sub_tree = copy.deepcopy(tree_node)
                                    parent_sub_tree.parent = None
                                    found = True
                                else:
                                    reference = reference.parent
                            else:
                                parent_sub_tree = process_tree
                                found = True

                # apply_reduction_rules(candidate_updated_tree)

                tree_updated = True
                negative_sub_trace = self.removal_candidates_generator.negative_trace
                positive_sub_traces = self.removal_candidates_generator.fitting_traces

                if parent_sub_tree != process_tree and Constants.USE_SUBTREE_BASED_CONFORMANCE_CHECKING == True:
                    negative_sub_trace = get_sub_trace_using_sublog([self.removal_candidates_generator.negative_trace],
                                                                    self.removal_candidates_generator.remove_sublogs,
                                                                    parent_sub_tree.id)[0]
                    positive_sub_traces = get_sub_trace_using_sublog(self.removal_candidates_generator.fitting_traces,
                                                                     self.removal_candidates_generator.keep_sublogs,
                                                                     parent_sub_tree.id)

                # step 2: check if the negative variant conforms to the tree, stop removing candidates
                negative_trace_fits = typed_trace_fits_process_tree(
                    negative_sub_trace,
                    parent_sub_tree
                )

                percentage_positive_traces_conforming = (
                    calculate_percentage_traces_conforming(
                        positive_sub_traces,
                        parent_sub_tree,
                    )
                )

                resulting_tree_edit_distance = calculate_process_tree_edit_distance(
                    process_tree,
                    copy.deepcopy(self.removal_candidates_generator.process_tree)
                )

                statistics['negative_trace_fits'] = negative_trace_fits
                statistics['percentage_positive_traces_conforming'] = percentage_positive_traces_conforming
                statistics['resulting_tree_edit_distance'] = resulting_tree_edit_distance

                if (statistics['negative_trace_fits'] == False
                    and statistics[
                        'percentage_positive_traces_conforming'] >= Constants.MIN_THRESHOLD_POSITIVE_FITTING_VARIANTS
                    and statistics[
                        'resulting_tree_edit_distance'] <= Constants.MAX_THRESHOLD_TREE_EDIT_DISTANCE):
                    statistics['thresholds_met'] = True

            self.update_operations += 1

        except Exception as e:
            print(e)

        return statistics

    def handle_sequence_operator(
        self, removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
    ) -> ProcessTree:
        """Handle the sequence operator"""

        seq_update = SequenceUpdateRule(removal_candidate_subtree)
        try:
            result = self.calculate_candidate_tree_statistics(
                seq_update.apply_rule(tree_to_update),
                removal_candidate_subtree,
                'sequence'
            )
            return result
        except Exception as e:
            print("Rule application failed (sequence): ", e)

        return None

    def handle_choice_operator(
        self, removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
    ) -> ProcessTree:
        """Handle the choice operator"""

        choice_update = ChoiceUpdateRule(removal_candidate_subtree)
        try:
            return self.calculate_candidate_tree_statistics(
                choice_update.apply_rule(tree_to_update),
                removal_candidate_subtree,
                'choice'
            )
        except Exception as e:
            print("Rule application failed (choice): ", e)

        return None

    def get_all_possible_sequentializations_for_parallel_operator_using_neg_variant(
        self,
        removal_candidate_subtree: CandidateSubtree,
    ):
        (subsequent_sequences_excluding_repetitions,
         execution_sequence_of_child_subtrees,
         trace_frequencies_corresponding_to_sequences) = (
            get_subsequent_subtree_execution_sequences_excluding_repetitions(
                self.removal_candidates_generator.remove_sublogs,
                removal_candidate_subtree
            ))
        execution_sequence_of_child_subtrees = execution_sequence_of_child_subtrees[0]

        # evaluate pre, mid, post and dynamic sequences from the subsequent_sequences_excluding_repetitions
        pre_sequences = []
        mid_sequences = []
        post_sequences = []
        dynamic_sequences = []

        for i in range(len(subsequent_sequences_excluding_repetitions)):

            if i == 0 and len(set(execution_sequence_of_child_subtrees)) >= 3:  # first subsequence
                # only one item in first subsequence and there are atleast three unique elements in execution_sequence_of_child_subtrees
                if (len(subsequent_sequences_excluding_repetitions[0]) == 1):
                    pre_sequences.extend(subsequent_sequences_excluding_repetitions[i])

                else:  # more than one item in the last subsequence
                    pre_sequences.extend(generate_prefix_combinations_of_subsequent_items(
                        subsequent_sequences_excluding_repetitions[i]))

            if i == len(subsequent_sequences_excluding_repetitions) - 1 and len(
                set(execution_sequence_of_child_subtrees)) >= 3:  # last subsequence
                # only one item in last subsequence and there are atleast three unique elements in execution_sequence_of_child_subtrees
                if (len(subsequent_sequences_excluding_repetitions[
                            len(subsequent_sequences_excluding_repetitions) - 1]) == 1):
                    post_sequences.extend(subsequent_sequences_excluding_repetitions[i])
                else:  # more than one item in the first subsequence
                    post_sequences.extend(generate_postfix_combinations_of_subsequent_items(
                        subsequent_sequences_excluding_repetitions[i]))

            if len(subsequent_sequences_excluding_repetitions[i]) > 4:
                mid_sequences.extend(generate_mid_combinations_of_subsequent_items(
                    subsequent_sequences_excluding_repetitions[i]))

            # if length of current variant's executed subtrees is atleast 2. Or if this length is exactly 2 and there
            # are atleast 3 activities/subtrees (some optional activities/subtrees not executed by this variant)
            if (len(subsequent_sequences_excluding_repetitions[i]) > 2 or
                (len(removal_candidate_subtree.reference.children) and len(
                    subsequent_sequences_excluding_repetitions[i]) == 2)):
                dynamic_sequences.extend(generate_dynamic_permutations_of_subsequent_items(
                    subsequent_sequences_excluding_repetitions[i]))

        return pre_sequences, mid_sequences, post_sequences, dynamic_sequences, execution_sequence_of_child_subtrees


def get_sub_trace_using_sublog(traces: list[TypedTrace], sublog, subtree_id: int):
    temp_traces = []
    for item in sublog[subtree_id]:
        temp_trace = copy.deepcopy(traces[0])
        temp_trace.trace = item
        temp_traces.append(temp_trace)
    return temp_traces


def get_subsequent_subtree_execution_sequences_excluding_repetitions(
    sublogs,
    removal_candidate_subtree: ProcessTree
) -> (list[list[int]], list[list[int]]):
    execution_sequence_of_child_subtrees = []
    subsequent_sequences_excluding_repetitions = []
    ids_of_all_child_subtrees = get_ids_of_all_child_subtrees(removal_candidate_subtree.reference)
    trace_frequencies_corresponding_to_sequences = []

    # get the execution sequence of child subtrees
    for i in range(len(sublogs[removal_candidate_subtree.node_id])):
        sequence = []
        for activity in sublogs[removal_candidate_subtree.node_id][i]:
            child_subtree_id = -1
            for child_id in ids_of_all_child_subtrees:
                if is_activity_in_child_subtrees_sublog(sublogs,
                                                        activity['concept:id'],
                                                        child_id):
                    child_subtree_id = child_id
                    break
            if child_subtree_id != -1:
                sequence.append(child_subtree_id)
            else:
                sequence.append(activity['concept:id'])

        trace_frequencies_corresponding_to_sequences.append(
            sublogs[removal_candidate_subtree.node_id][i].attributes['frequency']
        )
        execution_sequence_of_child_subtrees.append(sequence)

        # remove duplicates from the execution sequence
        execution_sequence_of_child_subtrees[len(execution_sequence_of_child_subtrees) - 1] = \
            [k for k, g in
             itertools.groupby(execution_sequence_of_child_subtrees[len(execution_sequence_of_child_subtrees) - 1])]

    # get the subtrees which are repeating, to be used to break the sequence and get the subsequent sequences
    # excluding these repeating subtrees
    repeating_subtree_ids = []
    for i in range(len(execution_sequence_of_child_subtrees)):
        repeating_subtree_ids = {key: val for key, val in Counter(execution_sequence_of_child_subtrees[i]).items()
                                 if val > 1}

        # subsequent sequences excluding the repeating subtrees in repeating_subtree_ids
        temp_subsequent_sequences_excluding_repetitions = get_subsequent_sequences_excluding_repetitions(
            execution_sequence_of_child_subtrees[i],
            repeating_subtree_ids,
            0
        )
        if len(temp_subsequent_sequences_excluding_repetitions) == 1:
            subsequent_sequences_excluding_repetitions.extend(temp_subsequent_sequences_excluding_repetitions)
        else:
            for sequence in temp_subsequent_sequences_excluding_repetitions:
                subsequent_sequences_excluding_repetitions.append(sequence)

    return subsequent_sequences_excluding_repetitions, execution_sequence_of_child_subtrees, trace_frequencies_corresponding_to_sequences


def get_ids_of_all_child_subtrees(process_tree_node: ProcessTree) -> list[int]:
    ids_list = []
    for child in process_tree_node.children:
        if child.operator is not None:
            ids_list.append(child.id)
    return ids_list


def is_activity_in_child_subtrees_sublog(sublog, activity_id: int, child_subtree_id: int) -> bool:
    for i in range(len(sublog[child_subtree_id])):
        for activity in sublog[child_subtree_id][i]:
            if activity['concept:id'] == activity_id:
                return True
    return False


def get_subsequent_sequences_excluding_repetitions(
    execution_sequence_of_child_subtrees: list[int],
    repeating_subtree_ids: dict[int, int],
    min_length_of_sequence
):
    subsequent_sequences_excluding_repetitions = []
    temp_sequence = []
    for i in range(len(execution_sequence_of_child_subtrees)):
        if execution_sequence_of_child_subtrees[i] not in repeating_subtree_ids:
            temp_sequence.append(execution_sequence_of_child_subtrees[i])
        else:
            if len(temp_sequence) > min_length_of_sequence:
                subsequent_sequences_excluding_repetitions.append(temp_sequence)
            temp_sequence = []

    if len(temp_sequence) > min_length_of_sequence:
        subsequent_sequences_excluding_repetitions.append(temp_sequence)

    return subsequent_sequences_excluding_repetitions


def generate_combinations_of_subsequent_items(list_of_ids: list[int], min_length: int) -> list[list[int]]:
    combinations = []
    for i in range(len(list_of_ids)):
        if len([list_of_ids[i]]) >= min_length:
            combinations.extend([[list_of_ids[i]]])

        sequence = [list_of_ids[i]]
        for j in range(i + 1, len(list_of_ids)):
            sequence.append(list_of_ids[j])
            if len(sequence) >= min_length:
                combinations.extend(copy.deepcopy([sequence]))

    return sorted(
        combinations,
        key=lambda x: len(x),
        reverse=False,
    )


def generate_prefix_combinations_of_subsequent_items(list_of_ids: list[int]) -> list[list[int]]:
    permutations = []
    for i in range(1, len(list_of_ids)):
        remove_pairs = [list_of_ids[0:i]]
        perms = [list(t) for t in list(itertools.permutations(list_of_ids, i))]
        filtered_perms = [x for x in perms if x not in remove_pairs]
        permutations.extend(filtered_perms)

    return permutations


def generate_postfix_combinations_of_subsequent_items(list_of_ids: list[int]) -> list[list[int]]:
    permutations = []
    for i in range(1, len(list_of_ids)):
        remove_pairs = [list_of_ids[len(list_of_ids) - i:]]
        perms = [list(t) for t in list(itertools.permutations(list_of_ids, i))]
        filtered_perms = [x for x in perms if x not in remove_pairs]
        permutations.extend(filtered_perms)

    return permutations


def generate_mid_combinations_of_subsequent_items(list_of_ids: list[int]) -> list[list[int]]:
    mid_sequences = []
    mid_sequence_components = []

    for r in range(1, len(list_of_ids) - 3):
        perms = [list(t) for t in list(itertools.permutations(list_of_ids, r))]

        for p in range(len(perms)):
            middle_items = perms[p]

            for s in range(2, (len(list_of_ids) - r) - 1):
                remaining_items = [x for x in list_of_ids if x not in middle_items]
                left_combinations = [list(t) for t in list(itertools.combinations(remaining_items, s))]

                for q in range(len(left_combinations)):
                    left_combination = copy.deepcopy(left_combinations[q])
                    left_combination.extend(middle_items)

                    mid_sequence_components.append(
                        {
                            'left': copy.deepcopy(left_combinations[q]), 'middle': middle_items,
                            'right': [x for x in list_of_ids if x not in left_combination]
                        }
                    )

                    left_combination.extend([x for x in list_of_ids if x not in left_combination])
                    if left_combination == list_of_ids:
                        mid_sequence_components.pop()
                    else:
                        mid_sequences.append(left_combination)

    # permutations = sorted(permutations, key=lambda x: x)
    # [x for x in permutations if x[2]==16 and x[3]==19]

    # mid_sequences = np.unique(mid_sequences, axis=0).tolist()
    mid_sequences = [x for x in mid_sequences if x not in [list_of_ids]]

    return mid_sequence_components


def generate_dynamic_permutations_of_subsequent_items(list_of_ids: list[int]) -> list[list[int]]:
    permutations = []
    for i in range(2, len(list_of_ids) + 1):
        remove_pairs = [list(t) for t in list(itertools.combinations(list_of_ids, i))]
        perms = [list(t) for t in list(itertools.permutations(list_of_ids, i))]
        filtered_perms = [x for x in perms if x not in remove_pairs]
        permutations.extend(filtered_perms)

    return permutations


def generate_subsequent_pairs_of_given_length(list_of_ids: list[int], length: int) -> list[list[int]]:
    combinations = []
    for i in range(len(list_of_ids)):
        sequence = [list_of_ids[i]]
        for j in range(i + 1, min((i + length), len(list_of_ids))):
            sequence.append(list_of_ids[j])
        if (i + length) > len(list_of_ids):
            break
        combinations.extend([sequence])
    return combinations
    combinations
