import copy
import itertools
from collections import Counter
from pm4py.objects.process_tree.obj import Operator, ProcessTree
from cortado_core.negative_process_model_repair.constants import Constants
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_activity import \
    CandidateActivity
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_subtree import \
    CandidateSubtree
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.removal_candidates_heuristics import \
    RemovalCandidatesHeuristics
from cortado_core.negative_process_model_repair.removal_strategies.rules_based_reduction.subtree_update import \
    SubtreeUpdate
from cortado_core.negative_process_model_repair.removal_strategies.rules_based_reduction.update_rule import \
    SequenceUpdateRule, ChoiceUpdateRule, LoopUpdateRule, ParallelUpdateRule


class CompleteBruteForceSubtreeUpdate(SubtreeUpdate):

    def __init__(
        self,
        removal_candidates_generator: RemovalCandidatesHeuristics,
        removal_candidate_activities: list[CandidateActivity],
    ):
        super().__init__(removal_candidates_generator, removal_candidate_activities)

    def apply_complete_brute_force_subtree_update_based_reduction(self) -> (ProcessTree, bool, float):
        brute_force_results = []

        for removal_candidate_subtree in self.removal_candidate_subtrees:
            tree_to_update = copy.deepcopy(self.removal_candidates_generator.process_tree)
            result = None

            if removal_candidate_subtree.reference.operator == Operator.SEQUENCE:
                result = self.handle_sequence_operator(
                    removal_candidate_subtree, tree_to_update
                )

            elif removal_candidate_subtree.reference.operator == Operator.XOR:
                result = self.handle_choice_operator(
                    removal_candidate_subtree, tree_to_update
                )

            if (
                removal_candidate_subtree.reference.operator == Operator.PARALLEL
            ):
                result = self.handle_parallel_operator(
                    removal_candidate_subtree, tree_to_update
                )

            elif removal_candidate_subtree.reference.operator == Operator.LOOP:
                result = self.handle_loop_operator(
                    removal_candidate_subtree, tree_to_update
                )

            if result is not None:
                if isinstance(result, list):
                    brute_force_results.extend(result)
                else:
                    brute_force_results.append(result)

        successful_brute_Force_results = [token for token in brute_force_results if
                                          token['negative_trace_fits'] == False]
        successful_brute_Force_results = sorted(successful_brute_Force_results,
                                                key=lambda x: (
                                                    -x["percentage_positive_traces_conforming"],
                                                    x["resulting_tree_edit_distance"])
                                                )
        if len(successful_brute_Force_results) > 0:
            return (
                successful_brute_Force_results[0]['updated_tree'], True,
                successful_brute_Force_results[0]['percentage_positive_traces_conforming'],
                successful_brute_Force_results[0]['resulting_tree_edit_distance'],
                brute_force_results[0]['applied_rule'])
        else:
            return (None, False, None, None, None)

    def handle_sequence_operator(
        self, removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
    ) -> ProcessTree:
        """Handle the sequence operator"""

        seq_update = SequenceUpdateRule(removal_candidate_subtree)
        try:
            return self.calculate_candidate_tree_statistics(
                seq_update.apply_rule(tree_to_update),
                removal_candidate_subtree,
                'sequence'
            )
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

    def handle_loop_operator(
        self, removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
    ):
        loop_rules_results = []
        loop_update_rule = LoopUpdateRule(removal_candidate_subtree)

        removal_candidate_subtree.loop_subtree_stats = (
            self.removal_candidates_generator.calculate_trace_frequencies_based_heuristics_for_loop_operator(
                removal_candidate_subtree))

        try:
            loop_rules_results.append(
                self.calculate_candidate_tree_statistics(
                    loop_update_rule.apply_remove_redundant_redo(copy.deepcopy(tree_to_update)),
                    removal_candidate_subtree,
                    'loop: remove-redundant-redo'
                )
            )
        except Exception as e:
            print("Rule application failed (loop: remove-redundant-redo): ", e)

        try:
            loop_rules_results.append(
                self.calculate_candidate_tree_statistics(
                    loop_update_rule.apply_optional_redo_mandatory(copy.deepcopy(tree_to_update)),
                    removal_candidate_subtree,
                    'loop: optional_redo_mandatory'
                )
            )
        except Exception as e:
            print("Rule application failed (loop: optional_redo_mandatory): ", e)

        do_frequency_remove, redo_frequency_remove = (
            self.removal_candidates_generator.get_do_and_redo_frequencies_negative_variant(removal_candidate_subtree))

        combinations = self.get_all_loop_repetitions_encoding_combinations(
            do_frequency_remove,
            Constants.MAX_NUMBER_OF_LOOP_REPETITION_ENCODINGS,
            Constants.MAX_LENGTH_LOOP_REPETITION_ENCODING
        )

        for combination in combinations:
            try:
                loop_rules_results.append(
                    self.calculate_candidate_tree_statistics(
                        loop_update_rule.apply_repeat_exactly_n(copy.deepcopy(tree_to_update), combination),
                        removal_candidate_subtree,
                        'loop: repeat_exactly_n ' + str(combination)
                    )
                )
            except Exception as e:
                print("Rule application failed (loop: repeat_exactly_n): ", e)

        try:
            loop_rules_results.append(
                self.calculate_candidate_tree_statistics(
                    loop_update_rule.apply_repeat_at_least_n(
                        copy.deepcopy(tree_to_update), removal_candidate_subtree.loop_subtree_stats.do_frequency_remove
                    ),
                    removal_candidate_subtree,
                    'loop: repeat_at_least_n'
                )
            )
        except Exception as e:
            print("Rule application failed (loop: repeat_at_least_n): ", e)

        return loop_rules_results

    def handle_parallel_operator(
        self, removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
    ) -> ProcessTree:
        parallel_rules_results = []

        pre_sequences, mid_sequences, post_sequences, dynamic_sequences, execution_sequence_of_child_subtrees = self.get_all_possible_sequentializations_for_parallel_operator_using_neg_variant(
            removal_candidate_subtree)

        parallel_update_rule = ParallelUpdateRule(removal_candidate_subtree)

        for pre_sequence in pre_sequences:
            try:
                parallel_rules_results.append(
                    self.calculate_candidate_tree_statistics(
                        parallel_update_rule.apply_pre_sequeltialization_rule(copy.deepcopy(tree_to_update),
                                                                              pre_sequence),
                        removal_candidate_subtree,
                        'parallel: pre sequentialization ' + str(pre_sequence)
                    )
                )
            except Exception as e:
                print("Rule application failed (parallel: pre sequentialization): ", e)

        for post_sequence in post_sequences:
            try:
                parallel_rules_results.append(
                    self.calculate_candidate_tree_statistics(
                        parallel_update_rule.apply_post_sequeltialization_rule(copy.deepcopy(tree_to_update),
                                                                               post_sequence),
                        removal_candidate_subtree,
                        'parallel: post sequentialization ' + str(post_sequence)
                    )
                )
            except Exception as e:
                print("Rule application failed (parallel: post sequentialization): ", e)

        for mid_sequence in mid_sequences:
            try:
                parallel_rules_results.append(
                    self.calculate_candidate_tree_statistics(
                        parallel_update_rule.apply_mid_sequeltialization_rule(copy.deepcopy(tree_to_update),
                                                                              mid_sequence,
                                                                              execution_sequence_of_child_subtrees),
                        removal_candidate_subtree,
                        'parallel: mid sequentialization ' + str(mid_sequence)
                    )
                )
            except Exception as e:
                print("Rule application failed (parallel: mid sequentialization): ", e)

        for dynamic_sequence in dynamic_sequences:
            try:
                parallel_rules_results.append(
                    self.calculate_candidate_tree_statistics(
                        parallel_update_rule.apply_dynamic_sequeltialization_rule(copy.deepcopy(tree_to_update),
                                                                                  dynamic_sequence),
                        removal_candidate_subtree,
                        'parallel: dynamic sequentialization ' + str(dynamic_sequence)
                    )
                )
            except Exception as e:
                print("Rule application failed (parallel: dynamic sequentialization): ", e)

        return parallel_rules_results

    def get_all_loop_repetitions_encoding_combinations(
        self,
        do_frequency_remove: int,
        max_number_of_loop_repetition_encodings: int,
        max_length_loop_repetition_encoding: int
    ) -> list[list[int]]:
        repetitions = list(range(1, max_length_loop_repetition_encoding + 1))
        if do_frequency_remove in repetitions:
            repetitions.remove(do_frequency_remove)
        combinations = []

        if len(repetitions) > 2:

            max_num = max(max_number_of_loop_repetition_encodings, len(repetitions))
            for i in range(1, max_num + 1):
                if len(combinations) == 0:
                    combinations = (
                        [list(t) for t in list(itertools.combinations(repetitions, i))]
                    )
                else:
                    combinations.extend(
                        [list(t) for t in
                         list(itertools.combinations(repetitions, i))]
                    )

        return combinations

    def get_all_possible_sequentializations_for_parallel_operator_using_neg_variant(
        self,
        removal_candidate_subtree: CandidateSubtree,
    ):
        execution_sequence_of_child_subtrees: list[int] = []
        ids_of_all_child_subtrees = get_ids_of_all_child_subtrees(removal_candidate_subtree.reference)

        # get the execution sequence of child sub-trees
        for i in range(len(self.removal_candidates_generator.remove_sublogs[removal_candidate_subtree.node_id])):
            activity_sequences_per_subtree = {}
            sequence = []
            for activity in self.removal_candidates_generator.remove_sublogs[removal_candidate_subtree.node_id][i]:
                child_subtree_id = -1
                for child_id in ids_of_all_child_subtrees:
                    if is_activity_in_child_subtrees_sublog(self.removal_candidates_generator.remove_sublogs,
                                                            activity['concept:id'],
                                                            child_id):
                        child_subtree_id = child_id
                        break
                if child_subtree_id != -1:
                    execution_sequence_of_child_subtrees.append(child_subtree_id)
                else:
                    execution_sequence_of_child_subtrees.append(activity['concept:id'])

        # remove duplicates from the execution sequence
        execution_sequence_of_child_subtrees = [k for k, g in itertools.groupby(execution_sequence_of_child_subtrees)]

        # get the subtrees which are repeating, to be used to break the sequence and get the subsequent sequences excluding these repeating subtrees
        repeating_subtree_ids = {key: val for key, val in Counter(execution_sequence_of_child_subtrees).items() if
                                 val > 1}

        # subsequent sequences excluding the repeating subtrees in repeating_subtree_ids
        subsequent_sequences_excluding_repetitions = get_subsequent_sequences_excluding_repetitions(
            execution_sequence_of_child_subtrees, repeating_subtree_ids, 0)

        # evaluate pre, mid, post and dynamic sequences from the subsequent_sequences_excluding_repetitions
        pre_sequences = []
        mid_sequences = []
        post_sequences = []
        dynamic_sequences = []

        for i in range(len(subsequent_sequences_excluding_repetitions)):

            if i == 0:  # first subsequence
                # if len(subsequent_sequences_excluding_repetitions[i]) == 1:
                #     pre_sequences.extend([subsequent_sequences_excluding_repetitions[i]])
                # else:
                if len(subsequent_sequences_excluding_repetitions) == 1:  # only one single subsequence
                    pre_sequences.extend(generate_prefix_combinations_of_subsequent_items(
                        subsequent_sequences_excluding_repetitions[i]))

                    # mid_sequences.extend(generate_combinations_of_subsequent_items(
                    #     subsequent_sequences_excluding_repetitions[i][1:-1], 1))
                else:  # more than one subsequences
                    pre_sequences.extend(generate_prefix_combinations_of_subsequent_items(
                        subsequent_sequences_excluding_repetitions[i]))

                    # mid_sequences.extend(generate_combinations_of_subsequent_items(
                    #     subsequent_sequences_excluding_repetitions[i][1:], 1))

            if i == len(subsequent_sequences_excluding_repetitions) - 1:  # last subsequence
                # if len(subsequent_sequences_excluding_repetitions[i]) == 1:
                #     post_sequences.extend([subsequent_sequences_excluding_repetitions[i]])
                # else:
                if len(subsequent_sequences_excluding_repetitions) == 1:  # only one single subsequence
                    post_sequences.extend(generate_postfix_combinations_of_subsequent_items(
                        subsequent_sequences_excluding_repetitions[i]))
                else:  # more than one subsequences
                    post_sequences.extend(generate_postfix_combinations_of_subsequent_items(
                        subsequent_sequences_excluding_repetitions[i]))

                    # mid_sequences.extend(generate_combinations_of_subsequent_items(
                    #     subsequent_sequences_excluding_repetitions[i][:-1], 1))

            # if i != 0 and i != len(subsequent_sequences_excluding_repetitions) - 1: # middle subsequences
            if len(subsequent_sequences_excluding_repetitions[i]) > 4:
                mid_sequences.extend(generate_mid_combinations_of_subsequent_items(
                    subsequent_sequences_excluding_repetitions[i]))

            if len(subsequent_sequences_excluding_repetitions[i]) > 2:
                dynamic_sequences.extend(generate_dynamic_permutations_of_subsequent_items(
                    subsequent_sequences_excluding_repetitions[i]))

        return pre_sequences, mid_sequences, post_sequences, dynamic_sequences, execution_sequence_of_child_subtrees


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
    for i in range(1, len(list_of_ids) - 1):
        remove_pairs = [list_of_ids[0:i]]
        perms = [list(t) for t in list(itertools.permutations(list_of_ids, i))]
        filtered_perms = [x for x in perms if x not in remove_pairs]
        permutations.extend(filtered_perms)

    return permutations


def generate_postfix_combinations_of_subsequent_items(list_of_ids: list[int]) -> list[list[int]]:
    permutations = []
    for i in range(1, len(list_of_ids) - 1):
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
        remove_pairs = generate_subsequent_pairs_of_given_length(list_of_ids, i)
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

