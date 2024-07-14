import copy
import itertools
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
        if len(brute_force_results) > 0:
            return (
                brute_force_results[0]['updated_tree'], True,
                brute_force_results[0]['percentage_positive_traces_conforming'],
                brute_force_results[0]['resulting_tree_edit_distance'],
                brute_force_results[0]['applied_rule'])
        else:
            return (None, False, None, None, None)

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
        #         successful_brute_Force_results[0]['applied_rule'])
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

    def handle_parallel_operator(
        self, removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
    ) -> ProcessTree:
        parallel_rules_results = []

        pre_sequences, mid_sequences, post_sequences, dynamic_sequences, execution_sequence_of_child_subtrees = (
            self.get_all_possible_sequentializations_for_parallel_operator_using_neg_variant(
                removal_candidate_subtree))

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
