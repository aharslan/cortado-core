import copy
from pm4py.objects.process_tree.obj import Operator, ProcessTree
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_activity import \
    CandidateActivity
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_subtree import \
    CandidateSubtree
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.removal_candidates_heuristics import \
    RemovalCandidatesHeuristics
from cortado_core.negative_process_model_repair.removal_strategies.rules_based_reduction.subtree_update import \
    SubtreeUpdate, get_subsequent_subtree_execution_sequences_excluding_repetitions
from cortado_core.negative_process_model_repair.removal_strategies.rules_based_reduction.update_rule import \
    SequenceUpdateRule, ChoiceUpdateRule, LoopUpdateRule, ParallelUpdateRule


class HeuristicBruteForceSubtreeUpdate(SubtreeUpdate):

    def __init__(
        self,
        removal_candidates_generator: RemovalCandidatesHeuristics,
        removal_candidate_activities: list[CandidateActivity],
    ):
        super().__init__(removal_candidates_generator, removal_candidate_activities)

    def apply_heuristic_brute_force_subtree_update_based_reduction(self) -> (ProcessTree, bool, float):
        brute_force_results = []

        for removal_candidate_subtree in self.removal_candidate_subtrees:
            tree_to_update = copy.deepcopy(self.removal_candidates_generator.process_tree)
            result = None

            # if removal_candidate_subtree.reference.operator == Operator.SEQUENCE:
            #     result = self.handle_sequence_operator(
            #         removal_candidate_subtree, tree_to_update
            #     )
            #
            # elif removal_candidate_subtree.reference.operator == Operator.XOR:
            #     result = self.handle_choice_operator(
            #         removal_candidate_subtree, tree_to_update
            #     )

            if (
                removal_candidate_subtree.reference.operator == Operator.PARALLEL
            ):
                result = self.handle_parallel_operator(
                    removal_candidate_subtree, tree_to_update
                )

            # elif removal_candidate_subtree.reference.operator == Operator.LOOP:
            #     result = self.handle_loop_operator(
            #         removal_candidate_subtree, tree_to_update
            #     )
            if result is not None:
                if isinstance(result, list):
                    brute_force_results.extend(result)
                else:
                    brute_force_results.append(result)

        brute_force_results = sorted(brute_force_results,
                                     key=lambda x: (
                                         -x["percentage_positive_traces_conforming"], x["resulting_tree_edit_distance"])
                                     )
        return (
            brute_force_results[0]['updated_tree'], True,
            brute_force_results[0]['percentage_positive_traces_conforming'],
            brute_force_results[0]['resulting_tree_edit_distance'], brute_force_results[0]['applied_rule'])

        # successful_brute_Force_results = [token for token in brute_force_results if
        #                                   token['negative_trace_fits'] == False]
        # successful_brute_Force_results = sorted(successful_brute_Force_results,
        #                                         key=lambda x: (
        #                                             -x["percentage_positive_traces_conforming"],
        #                                             x["resulting_tree_edit_distance"])
        #                                         )
        # if len(successful_brute_Force_results) > 0:
        #     return (
        #         successful_brute_Force_results[0]['updated_tree'], True,
        #         successful_brute_Force_results[0]['percentage_positive_traces_conforming'],
        #         successful_brute_Force_results[0]['resulting_tree_edit_distance'],
        #         brute_force_results[0]['applied_rule'])
        # else:
        #     return (None, False, None, None, None)

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

        try:
            repetitions_to_encode: list[int] = list(
                removal_candidate_subtree.loop_subtree_stats.count_positive_variants_highest_loop_repetitions.keys()
            )

            loop_rules_results.append(
                self.calculate_candidate_tree_statistics(
                    loop_update_rule.apply_repeat_exactly_n(copy.deepcopy(tree_to_update), repetitions_to_encode),
                    removal_candidate_subtree,
                    'loop: repeat_exactly_n'
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

        (subsequent_sequences_excluding_repetitions_positive,
         execution_sequence_of_child_subtrees_positive,
         trace_frequencies_corresponding_to_sequences_positive) = (
            get_subsequent_subtree_execution_sequences_excluding_repetitions(
                self.removal_candidates_generator.keep_sublogs,
                removal_candidate_subtree
            ))

        (pre_sequences_negative, mid_sequences_negative, post_sequences_negative,
         dynamic_sequences_negative, execution_sequence_of_child_subtrees_negative) = (
            self.get_all_possible_sequentializations_for_parallel_operator_using_neg_variant(
                removal_candidate_subtree))

        dynamic_sequences_negative_pruned = prune_dynamic_sequences_negative(
            subsequent_sequences_excluding_repetitions_positive,
            trace_frequencies_corresponding_to_sequences_positive,
            dynamic_sequences_negative
        )

        pre_sequences_negative_pruned = prune_pre_sequences_negative(
            subsequent_sequences_excluding_repetitions_positive,
            trace_frequencies_corresponding_to_sequences_positive,
            pre_sequences_negative
        )

        post_sequences_negative_pruned = prune_post_sequences_negative(
            subsequent_sequences_excluding_repetitions_positive,
            trace_frequencies_corresponding_to_sequences_positive,
            post_sequences_negative
        )

        mid_sequences_negative_pruned = prune_mid_sequences_negative(
            subsequent_sequences_excluding_repetitions_positive,
            trace_frequencies_corresponding_to_sequences_positive,
            mid_sequences_negative
        )

        parallel_update_rule = ParallelUpdateRule(removal_candidate_subtree)

        for key in pre_sequences_negative_pruned:
            try:
                parallel_rules_results.append(
                    self.calculate_candidate_tree_statistics(
                        parallel_update_rule.apply_pre_sequeltialization_rule(copy.deepcopy(tree_to_update),
                                                                              pre_sequences_negative_pruned[key][
                                                                                  'sequence']),
                        removal_candidate_subtree,
                        'parallel: pre sequentialization ' + str(pre_sequences_negative_pruned[key]['sequence'])
                    )
                )
            except Exception as e:
                print("Rule application failed (parallel: pre sequentialization): ", e)

        for key in post_sequences_negative_pruned:
            try:
                parallel_rules_results.append(
                    self.calculate_candidate_tree_statistics(
                        parallel_update_rule.apply_post_sequeltialization_rule(copy.deepcopy(tree_to_update),
                                                                               post_sequences_negative_pruned[key][
                                                                                   'sequence']),
                        removal_candidate_subtree,
                        'parallel: post sequentialization ' + str(post_sequences_negative_pruned[key]['sequence'])
                    )
                )
            except Exception as e:
                print("Rule application failed (parallel: post sequentialization): ", e)

        for key in mid_sequences_negative_pruned:
            try:
                parallel_rules_results.append(
                    self.calculate_candidate_tree_statistics(
                        parallel_update_rule.apply_mid_sequeltialization_rule(copy.deepcopy(tree_to_update),
                                                                              mid_sequences_negative_pruned[key][
                                                                                  'sequence'],
                                                                              execution_sequence_of_child_subtrees_negative),
                        removal_candidate_subtree,
                        'parallel: mid sequentialization ' + str(mid_sequences_negative_pruned[key]['sequence'])
                    )
                )
            except Exception as e:
                print("Rule application failed (parallel: mid sequentialization): ", e)

        for key in dynamic_sequences_negative_pruned:
            try:
                parallel_rules_results.append(
                    self.calculate_candidate_tree_statistics(
                        parallel_update_rule.apply_dynamic_sequeltialization_rule(copy.deepcopy(tree_to_update),
                                                                                  dynamic_sequences_negative_pruned[
                                                                                      key]['sequence']),
                        removal_candidate_subtree,
                        'parallel: dynamic sequentialization ' + str(dynamic_sequences_negative_pruned[key]['sequence'])
                    )
                )
            except Exception as e:
                print("Rule application failed (parallel: dynamic sequentialization): ", e)

        return parallel_rules_results


def prune_dynamic_sequences_negative(subsequent_sequences_excluding_repetitions_positive: list[list[int]],
                                     trace_frequencies_corresponding_to_sequences_positive: list[int],
                                     dynamic_sequences_negative: list[list[int]]) -> dict:
    dynamic_sequences_negative_pruned = {}
    for dynamic_sequence in dynamic_sequences_negative:
        for i in range(len(subsequent_sequences_excluding_repetitions_positive)):
            for j in range(len(subsequent_sequences_excluding_repetitions_positive[i])):
                if isinstance(subsequent_sequences_excluding_repetitions_positive[i][j], list):
                    if is_subarr_in_same_sequence_in_arr(dynamic_sequence,
                                                         subsequent_sequences_excluding_repetitions_positive[i][j]):
                        if str(dynamic_sequence) in dynamic_sequences_negative_pruned:
                            dynamic_sequences_negative_pruned[str(dynamic_sequence)]['trace_frequency'] += (
                                trace_frequencies_corresponding_to_sequences_positive[i])
                            break
                        else:
                            dynamic_sequences_negative_pruned[str(dynamic_sequence)] = \
                                {'sequence': dynamic_sequence,
                                 'trace_frequency': trace_frequencies_corresponding_to_sequences_positive[i]}
                            break
                else:
                    if is_subarr_in_same_sequence_in_arr(dynamic_sequence,
                                                         subsequent_sequences_excluding_repetitions_positive[i]):
                        if str(dynamic_sequence) in dynamic_sequences_negative_pruned:
                            dynamic_sequences_negative_pruned[str(dynamic_sequence)]['trace_frequency'] += (
                                trace_frequencies_corresponding_to_sequences_positive[i])
                            break
                        else:
                            dynamic_sequences_negative_pruned[str(dynamic_sequence)] = \
                                {'sequence': dynamic_sequence,
                                 'trace_frequency': trace_frequencies_corresponding_to_sequences_positive[i]}
                            break

    return dynamic_sequences_negative_pruned


def prune_pre_sequences_negative(subsequent_sequences_excluding_repetitions_positive: list[list[int]],
                                 trace_frequencies_corresponding_to_sequences_positive: list[int],
                                 pre_sequences_negative: list[list[int]]) -> dict:
    pre_sequences_negative_pruned = {}
    for dynamic_sequence in pre_sequences_negative:
        for i in range(len(subsequent_sequences_excluding_repetitions_positive)):

            if isinstance(subsequent_sequences_excluding_repetitions_positive[i][0], list):
                if is_left_subarray(dynamic_sequence, subsequent_sequences_excluding_repetitions_positive[i][0]):
                    if str(dynamic_sequence) in pre_sequences_negative_pruned:
                        pre_sequences_negative_pruned[str(dynamic_sequence)]['trace_frequency'] += (
                            trace_frequencies_corresponding_to_sequences_positive[i])
                    else:
                        pre_sequences_negative_pruned[str(dynamic_sequence)] = \
                            {'sequence': dynamic_sequence,
                             'trace_frequency': trace_frequencies_corresponding_to_sequences_positive[i]}
            else:
                if is_left_subarray(dynamic_sequence, subsequent_sequences_excluding_repetitions_positive[i]):
                    if str(dynamic_sequence) in pre_sequences_negative_pruned:
                        pre_sequences_negative_pruned[str(dynamic_sequence)]['trace_frequency'] += (
                            trace_frequencies_corresponding_to_sequences_positive[i])
                    else:
                        pre_sequences_negative_pruned[str(dynamic_sequence)] = \
                            {'sequence': dynamic_sequence,
                             'trace_frequency': trace_frequencies_corresponding_to_sequences_positive[i]}

    return pre_sequences_negative_pruned


def prune_post_sequences_negative(subsequent_sequences_excluding_repetitions_positive: list[list[int]],
                                  trace_frequencies_corresponding_to_sequences_positive: list[int],
                                  post_sequences_negative: list[list[int]]) -> dict:
    post_sequences_negative_pruned = {}
    for dynamic_sequence in post_sequences_negative:
        for i in range(len(subsequent_sequences_excluding_repetitions_positive)):

            if isinstance(subsequent_sequences_excluding_repetitions_positive[i][
                              len(subsequent_sequences_excluding_repetitions_positive[i]) - 1], list):
                if is_right_subarray(dynamic_sequence, subsequent_sequences_excluding_repetitions_positive[i][
                    len(subsequent_sequences_excluding_repetitions_positive[i]) - 1]):
                    if str(dynamic_sequence) in post_sequences_negative_pruned:
                        post_sequences_negative_pruned[str(dynamic_sequence)]['trace_frequency'] += (
                            trace_frequencies_corresponding_to_sequences_positive[i])
                    else:
                        post_sequences_negative_pruned[str(dynamic_sequence)] = \
                            {'sequence': dynamic_sequence,
                             'trace_frequency': trace_frequencies_corresponding_to_sequences_positive[i]}
            else:
                if is_right_subarray(dynamic_sequence, subsequent_sequences_excluding_repetitions_positive[i]):
                    if str(dynamic_sequence) in post_sequences_negative_pruned:
                        post_sequences_negative_pruned[str(dynamic_sequence)]['trace_frequency'] += (
                            trace_frequencies_corresponding_to_sequences_positive[i])
                    else:
                        post_sequences_negative_pruned[str(dynamic_sequence)] = \
                            {'sequence': dynamic_sequence,
                             'trace_frequency': trace_frequencies_corresponding_to_sequences_positive[i]}

    return post_sequences_negative_pruned


def prune_mid_sequences_negative(subsequent_sequences_excluding_repetitions_positive: list[list[int]],
                                 trace_frequencies_corresponding_to_sequences_positive: list[int],
                                 mid_sequences_negative: list[list[int]]) -> dict:
    mid_sequences_negative_pruned = {}
    for mid_sequence in mid_sequences_negative:
        for i in range(len(subsequent_sequences_excluding_repetitions_positive)):
            for j in range(len(subsequent_sequences_excluding_repetitions_positive[i])):
                if isinstance(subsequent_sequences_excluding_repetitions_positive[i][j], list):
                    if is_mid_sequence_in_execution_sequence(mid_sequence,
                                                             subsequent_sequences_excluding_repetitions_positive[i][j]):
                        if str(mid_sequence) in mid_sequences_negative_pruned:
                            mid_sequences_negative_pruned[str(mid_sequence)]['trace_frequency'] += (
                                trace_frequencies_corresponding_to_sequences_positive[i])
                            break
                        else:
                            mid_sequences_negative_pruned[str(mid_sequence)] = \
                                {'sequence': mid_sequence,
                                 'trace_frequency': trace_frequencies_corresponding_to_sequences_positive[i]}
                            break
                else:
                    if is_mid_sequence_in_execution_sequence(mid_sequence,
                                                             subsequent_sequences_excluding_repetitions_positive[i]):
                        if str(mid_sequence) in mid_sequences_negative_pruned:
                            mid_sequences_negative_pruned[str(mid_sequence)]['trace_frequency'] += (
                                trace_frequencies_corresponding_to_sequences_positive[i])
                            break
                        else:
                            mid_sequences_negative_pruned[str(mid_sequence)] = \
                                {'sequence': mid_sequence,
                                 'trace_frequency': trace_frequencies_corresponding_to_sequences_positive[i]}
                            break

    return mid_sequences_negative_pruned


def is_mid_sequence_in_execution_sequence(mid_sequence, execution_sequence):
    if len(execution_sequence) >= (len(mid_sequence['left']) + len(mid_sequence['middle'])):
        if (mid_sequence['middle'] ==
            execution_sequence[len(mid_sequence['left']):len(mid_sequence['left']) + len(mid_sequence['middle'])]):
            is_middle_equal = True
            if (set(mid_sequence['left']) ==
                set(execution_sequence[0:len(mid_sequence['left'])])):
                return True
    return False


def is_subarr_in_same_sequence_in_arr(subarr, arr):
    assert len(subarr) > 0
    assert len(arr) > 0

    if len(subarr) == len(arr) and subarr == arr:
        return True

    if len(subarr) < len(arr):
        sequential_matches = 0
        current_match = subarr[0]
        for i in range(len(arr)):
            if current_match == arr[i]:
                sequential_matches += 1
            if sequential_matches == len(subarr):
                return True
            current_match = subarr[sequential_matches]
    return False


def is_subarray(subarr, arr):
    for i in range(len(arr)):
        if arr[i:i + len(subarr)] == subarr:
            return True
    return False


def is_left_subarray(subarr, arr):
    if arr[0:0 + len(subarr)] == subarr:
        return True
    return False


def is_right_subarray(subarr, arr):
    if arr[-len(subarr):len(arr)] == subarr:
        return True
    return False
