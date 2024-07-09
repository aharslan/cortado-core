import copy
from pm4py.objects.process_tree.obj import Operator, ProcessTree
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

            if removal_candidate_subtree.reference.operator == Operator.SEQUENCE:
                result = self.handle_sequence_operator(
                    removal_candidate_subtree, tree_to_update
                )

            elif removal_candidate_subtree.reference.operator == Operator.XOR:
                result = self.handle_choice_operator(
                    removal_candidate_subtree, tree_to_update
                )

            elif (
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

        brute_force_results = sorted(brute_force_results,
                                     key=lambda x: (
                                         -x["percentage_positive_traces_conforming"], x["resulting_tree_edit_distance"])
                                     )
        return (
            brute_force_results[0]['updated_tree'], True,
            brute_force_results[0]['percentage_positive_traces_conforming'],
            brute_force_results[0]['resulting_tree_edit_distance'], brute_force_results[0]['applied_rule'])

    def handle_sequence_operator(
            self, removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
    ) -> ProcessTree:
        """Handle the sequence operator"""

        seq_update = SequenceUpdateRule(removal_candidate_subtree)

        return self.calculate_candidate_tree_statistics(
            seq_update.apply_rule(tree_to_update),
            removal_candidate_subtree,
            'sequence'
        )

    def handle_choice_operator(
            self, removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
    ) -> ProcessTree:
        """Handle the choice operator"""

        choice_update = ChoiceUpdateRule(removal_candidate_subtree)

        return self.calculate_candidate_tree_statistics(
            choice_update.apply_rule(tree_to_update),
            removal_candidate_subtree,
            'choice'
        )

    def handle_loop_operator(
            self, removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
    ):
        loop_rules_results = []
        loop_update_rule = LoopUpdateRule(removal_candidate_subtree)

        removal_candidate_subtree.loop_subtree_stats = (
            self.removal_candidates_generator.calculate_trace_frequencies_based_heuristics_for_loop_operator(
                removal_candidate_subtree))

        loop_rules_results.append(
            self.calculate_candidate_tree_statistics(
                loop_update_rule.apply_remove_redundant_redo(copy.deepcopy(tree_to_update)),
                removal_candidate_subtree,
                'loop: remove-redundant-redo'
            )
        )

        loop_rules_results.append(
            self.calculate_candidate_tree_statistics(
                loop_update_rule.apply_optional_redo_mandatory(copy.deepcopy(tree_to_update)),
                removal_candidate_subtree,
                'loop: optional_redo_mandatory'
            )
        )

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

        loop_rules_results.append(
            self.calculate_candidate_tree_statistics(
                loop_update_rule.apply_repeat_at_least_n(
                    copy.deepcopy(tree_to_update), removal_candidate_subtree.loop_subtree_stats.do_frequency_remove
                ),
                removal_candidate_subtree,
                'loop: repeat_at_least_n'
            )
        )

        return loop_rules_results

    def handle_parallel_operator(
            self, removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
    ) -> ProcessTree:
        parallel_rules_results = []

        removal_candidate_subtree.candidate_patterns = (
            self.removal_candidates_generator.calculate_trace_frequencies_based_heuristics_for_parallel_operator(
                removal_candidate_subtree, False)
        )

        for i in range(len(removal_candidate_subtree.candidate_patterns)):

            if len(removal_candidate_subtree.candidate_patterns[i]["alternate-combinations"]) > 0:
                alternate_combination = removal_candidate_subtree.candidate_patterns[i][
                    "alternate-combinations"
                ][0]["combination"]

                parallel_update_rule = ParallelUpdateRule(removal_candidate_subtree)

                parallel_rules_results.append(
                    self.calculate_candidate_tree_statistics(
                        parallel_update_rule.apply_rule(copy.deepcopy(tree_to_update), alternate_combination),
                        removal_candidate_subtree,
                        'parallel'
                    )
                )

        return parallel_rules_results

