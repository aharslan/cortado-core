import copy
from pm4py.objects.process_tree.obj import Operator, ProcessTree
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_activity import \
    CandidateActivity
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_subtree import \
    CandidateSubtree, Candidate
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.removal_candidates_heuristics import \
    RemovalCandidatesHeuristics
from cortado_core.negative_process_model_repair.removal_strategies.rules_based_reduction.heuristic_brute_force_subtree_update import \
    HeuristicBruteForceSubtreeUpdate
from cortado_core.negative_process_model_repair.removal_strategies.rules_based_reduction.subtree_update import \
    SubtreeUpdate, get_subsequent_subtree_execution_sequences_excluding_repetitions
from cortado_core.negative_process_model_repair.removal_strategies.rules_based_reduction.update_rule import \
    LoopUpdateRule, ParallelUpdateRule
from cortado_core.negative_process_model_repair.constants import Constants
from cortado_core.negative_process_model_repair.temp_utils import find_tree_node_by_id


class RatingBasedSubtreeUpdate(HeuristicBruteForceSubtreeUpdate):

    def __init__(
        self,
        removal_candidates_generator: RemovalCandidatesHeuristics,
        removal_candidate_activities: list[CandidateActivity],
    ):
        super().__init__(removal_candidates_generator, removal_candidate_activities)

    def apply_rating_based_subtree_update_reduction(self) -> (
        ProcessTree, bool, bool, float, float, str):
        brute_force_results = []

        candidates: [Candidate] = []

        for removal_candidate_subtree in self.removal_candidate_subtrees:
            tree_to_update = copy.deepcopy(self.removal_candidates_generator.process_tree)
            result = None

            if removal_candidate_subtree.reference.operator == Operator.SEQUENCE:
                removal_candidate_subtree.rating = self.rate_sequence_sub_tree(removal_candidate_subtree)
                candidates.append(
                    Candidate(removal_candidate_subtree, 'sequence', removal_candidate_subtree.rating, {}))

            elif removal_candidate_subtree.reference.operator == Operator.XOR:
                removal_candidate_subtree.rating = self.rate_choice_sub_tree(removal_candidate_subtree)
                candidates.append(
                    Candidate(removal_candidate_subtree, 'choice', removal_candidate_subtree.rating, None))

            elif removal_candidate_subtree.reference.operator == Operator.LOOP:
                removal_candidate_subtree.rating, removal_candidate_subtree.rating_approach = (
                    self.rate_loop_sub_tree(removal_candidate_subtree))

                candidates.append(
                    Candidate(removal_candidate_subtree, 'remove_redundant_redo',
                              removal_candidate_subtree.loop_subtree_stats.remove_redundant_redo_pm,
                              None))
                candidates.append(
                    Candidate(removal_candidate_subtree, 'optional_redo_mandatory',
                              removal_candidate_subtree.loop_subtree_stats.optional_redo_mandatory_pm,
                              None))
                candidates.append(
                    Candidate(removal_candidate_subtree, 'repeat_exactly_n',
                              removal_candidate_subtree.loop_subtree_stats.repeat_exactly_n_pm,
                              None))
                candidates.append(
                    Candidate(removal_candidate_subtree, 'repeat_at_least_n',
                              removal_candidate_subtree.loop_subtree_stats.repeat_atleast_n_pm,
                              None))
            elif (
                removal_candidate_subtree.reference.operator == Operator.PARALLEL
            ):
                (dynamic_sequences_negative_pruned,
                 pre_sequences_negative_pruned,
                 post_sequences_negative_pruned,
                 mid_sequences_negative_pruned,
                 execution_sequence_of_child_subtrees_negative) = (
                    self.rate_parallel_sub_tree(removal_candidate_subtree))

                for seq in dynamic_sequences_negative_pruned:
                    candidates.append(
                        Candidate(removal_candidate_subtree, 'dynamic_sequence',
                                  dynamic_sequences_negative_pruned[seq]['pm'],
                                  {'sequence': dynamic_sequences_negative_pruned[seq]['sequence'],
                                   'execution_sequence_of_child_subtrees_negative': execution_sequence_of_child_subtrees_negative}))
                for seq in pre_sequences_negative_pruned:
                    candidates.append(
                        Candidate(removal_candidate_subtree, 'pre_sequence', pre_sequences_negative_pruned[seq]['pm'],
                                  {'sequence': pre_sequences_negative_pruned[seq]['sequence'],
                                   'execution_sequence_of_child_subtrees_negative': execution_sequence_of_child_subtrees_negative}))
                for seq in post_sequences_negative_pruned:
                    candidates.append(
                        Candidate(removal_candidate_subtree, 'post_sequence', post_sequences_negative_pruned[seq]['pm'],
                                  {'sequence': post_sequences_negative_pruned[seq]['sequence'],
                                   'execution_sequence_of_child_subtrees_negative': execution_sequence_of_child_subtrees_negative}))
                for seq in mid_sequences_negative_pruned:
                    candidates.append(
                        Candidate(removal_candidate_subtree, 'mid_sequence', mid_sequences_negative_pruned[seq]['pm'],
                                  {'sequence': mid_sequences_negative_pruned[seq]['sequence'],
                                   'execution_sequence_of_child_subtrees_negative': execution_sequence_of_child_subtrees_negative}))

        candidates = sorted(candidates, key=lambda x: x.rating)
        # removal_candidate_subtrees = sorted(self.removal_candidate_subtrees, key=lambda x: x.rating)

        successful_results = []
        failed_updates = []
        for candidate in candidates:

            if candidate.removal_candidate_subtree.reference.operator == Operator.SEQUENCE:
                result = self.handle_sequence_operator(
                    candidate.removal_candidate_subtree, copy.deepcopy(tree_to_update)
                )
            elif candidate.removal_candidate_subtree.reference.operator == Operator.XOR:
                result = self.handle_choice_operator(
                    candidate.removal_candidate_subtree, copy.deepcopy(tree_to_update)
                )
            elif candidate.removal_candidate_subtree.reference.operator == Operator.LOOP:
                result = self.handle_loop_operator(
                    candidate, copy.deepcopy(tree_to_update)
                )
            elif candidate.removal_candidate_subtree.reference.operator == Operator.PARALLEL:
                result = self.handle_parallel_operator(
                    candidate, copy.deepcopy(tree_to_update)
                )

            if result is not None and result['negative_trace_fits'] == False:
                result['rating'] = candidate.rating
                result['rating'] = candidate.rating
                if isinstance(result, list):
                    successful_results.extend(result)
                else:
                    successful_results.append(result)
            else:
                if isinstance(result, list):
                    failed_updates.extend({'candidate_sub_tree': removal_candidate_subtree, 'result': result})
                else:
                    failed_updates.append({'candidate_sub_tree': removal_candidate_subtree, 'result': result})

            if (Constants.STOP_WHEN_AN_UPDATE_MEETS_THRESHOLD == True and
                len(successful_results) > 0 and
                successful_results[len(successful_results) - 1]['thresholds_met'] == True and
                successful_results[len(successful_results) - 1]['negative_trace_fits'] == False):
                break

        successful_results = sorted(successful_results,
                                    key=lambda x: (
                                        -x["percentage_positive_traces_conforming"],
                                        x["resulting_tree_edit_distance"])
                                    )
        if len(successful_results) > 0:
            return (
                successful_results[0]['updated_tree'],
                True,
                successful_results[0]['thresholds_met'],
                successful_results[0]['percentage_positive_traces_conforming'],
                successful_results[0]['resulting_tree_edit_distance'],
                successful_results[0]['applied_rule'],
                failed_updates)
        else:
            return (None, False, False, None, None, None, [])

    def rate_sequence_sub_tree(self, removal_candidate_subtree: CandidateSubtree):
        # count total number of positive traces executing this subtree
        positive_trace_count = 0

        if removal_candidate_subtree.reference.id in self.removal_candidates_generator.keep_sublogs:
            for item in self.removal_candidates_generator.keep_sublogs[removal_candidate_subtree.reference.id]:
                positive_trace_count += item.attributes['frequency']

        rating = positive_trace_count / self.removal_candidates_generator.positive_traces_frequency

        return rating

    def rate_choice_sub_tree(self, removal_candidate_subtree: CandidateSubtree):
        # count total number of positive traces executing this subtree
        positive_trace_count = 0

        if len(removal_candidate_subtree.candidate_subtrees) > 0:
            reference = removal_candidate_subtree.candidate_subtrees[0]
            while True:
                if len(reference.candidate_activities) > 0:
                    positive_trace_count = reference.candidate_activities[0].frequency_in_positive_traces
                    break
                else:
                    reference = reference.candidate_subtrees[0]

        elif len(removal_candidate_subtree.candidate_activities) > 0:
            positive_trace_count = removal_candidate_subtree.candidate_activities[0].frequency_in_positive_traces

        rating = positive_trace_count / self.removal_candidates_generator.positive_traces_frequency

        return rating

    def rate_loop_sub_tree(self, removal_candidate_subtree: CandidateSubtree):
        try:
            removal_candidate_subtree.loop_subtree_stats = (
                self.removal_candidates_generator.calculate_trace_frequencies_based_heuristics_for_loop_operator(
                    removal_candidate_subtree))
        except Exception as e:
            print(e)

        return (removal_candidate_subtree.loop_subtree_stats.min_pm,
                removal_candidate_subtree.loop_subtree_stats.min_pm_approach)

    def rate_parallel_sub_tree(self, removal_candidate_subtree: CandidateSubtree):

        maximum_positive_traces_in_a_sequence = -1
        sequence_approach_with_maximum_traces = None
        sequence = None

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

        dynamic_sequences_negative_pruned = self.prune_dynamic_sequences_negative(
            subsequent_sequences_excluding_repetitions_positive,
            trace_frequencies_corresponding_to_sequences_positive,
            dynamic_sequences_negative
        )

        dynamic_sequences_negative_pruned = {k: v for k, v in
                                             sorted(dynamic_sequences_negative_pruned.items(),
                                                    key=lambda item: item[1]['trace_frequency'])}

        if len(dynamic_sequences_negative_pruned) > 0:

            trace_count = 0
            for value in dynamic_sequences_negative_pruned:
                seq = dynamic_sequences_negative_pruned[value]['sequence']
                f = dynamic_sequences_negative_pruned[value]['trace_frequency']
                for pos_child_sequence in execution_sequence_of_child_subtrees_positive:
                    if not str(seq)[1:-1] in str(pos_child_sequence)[1:-1]:
                        trace_count += f
                dynamic_sequences_negative_pruned[value][
                    'pm'] = trace_count / self.removal_candidates_generator.positive_traces_frequency

            maximum_positive_traces_in_a_sequence = \
                dynamic_sequences_negative_pruned[list(dynamic_sequences_negative_pruned.keys())[0]]['trace_frequency']
            sequence_approach_with_maximum_traces = 'dynamic_sequence'
            sequence = dynamic_sequences_negative_pruned[list(dynamic_sequences_negative_pruned.keys())[0]]['sequence']

        pre_sequences_negative_pruned = self.prune_pre_sequences_negative(
            subsequent_sequences_excluding_repetitions_positive,
            trace_frequencies_corresponding_to_sequences_positive,
            pre_sequences_negative
        )

        pre_sequences_negative_pruned = {k: v for k, v in
                                         sorted(pre_sequences_negative_pruned.items(),
                                                key=lambda item: item[1]['trace_frequency'])}

        if len(pre_sequences_negative_pruned) > 0:

            trace_count = 0
            for value in pre_sequences_negative_pruned:
                seq = pre_sequences_negative_pruned[value]['sequence']
                f = pre_sequences_negative_pruned[value]['trace_frequency']
                for pos_child_sequence in execution_sequence_of_child_subtrees_positive:
                    if not pos_child_sequence[:len(seq)] == seq:
                        trace_count += f
                pre_sequences_negative_pruned[value][
                    'pm'] = trace_count / self.removal_candidates_generator.positive_traces_frequency

            maximum_positive_traces_in_a_sequence = \
                pre_sequences_negative_pruned[list(pre_sequences_negative_pruned.keys())[0]]['trace_frequency']
            sequence_approach_with_maximum_traces = 'pre_sequence'
            sequence = pre_sequences_negative_pruned[list(pre_sequences_negative_pruned.keys())[0]]['sequence']

        post_sequences_negative_pruned = self.prune_post_sequences_negative(
            subsequent_sequences_excluding_repetitions_positive,
            trace_frequencies_corresponding_to_sequences_positive,
            post_sequences_negative
        )

        post_sequences_negative_pruned = {k: v for k, v in
                                          sorted(post_sequences_negative_pruned.items(),
                                                 key=lambda item: item[1]['trace_frequency'])}

        if len(post_sequences_negative_pruned) > 0:

            trace_count = 0
            for value in post_sequences_negative_pruned:
                seq = post_sequences_negative_pruned[value]['sequence']
                f = post_sequences_negative_pruned[value]['trace_frequency']
                for pos_child_sequence in execution_sequence_of_child_subtrees_positive:
                    if not pos_child_sequence[-len(seq):] == seq:
                        trace_count += f
                post_sequences_negative_pruned[value][
                    'pm'] = trace_count / self.removal_candidates_generator.positive_traces_frequency

            maximum_positive_traces_in_a_sequence = \
                post_sequences_negative_pruned[list(post_sequences_negative_pruned.keys())[0]]['trace_frequency']
            sequence_approach_with_maximum_traces = 'post_sequence'
            sequence = post_sequences_negative_pruned[list(post_sequences_negative_pruned.keys())[0]]['sequence']

        mid_sequences_negative_pruned = self.prune_mid_sequences_negative(
            subsequent_sequences_excluding_repetitions_positive,
            trace_frequencies_corresponding_to_sequences_positive,
            mid_sequences_negative
        )

        mid_sequences_negative_pruned = {k: v for k, v in
                                         sorted(mid_sequences_negative_pruned.items(),
                                                key=lambda item: item[1]['trace_frequency'])}

        if len(mid_sequences_negative_pruned) > 0:

            trace_count = 0
            for value in mid_sequences_negative_pruned:
                seq = mid_sequences_negative_pruned[value]['sequence']
                f = mid_sequences_negative_pruned[value]['trace_frequency']
                for pos_child_sequence in execution_sequence_of_child_subtrees_positive:
                    if (not set(seq['left']) == set(pos_child_sequence[0:len(seq['left'])]) and
                        not seq['middle'] == pos_child_sequence[
                                             len(seq['left']):len(seq['left']) + len(seq['middle'])]):
                        trace_count += f
                mid_sequences_negative_pruned[value][
                    'pm'] = trace_count / self.removal_candidates_generator.positive_traces_frequency

            maximum_positive_traces_in_a_sequence = \
                mid_sequences_negative_pruned[list(mid_sequences_negative_pruned.keys())[0]]['trace_frequency']
            sequence_approach_with_maximum_traces = 'mid_sequence'
            sequence = mid_sequences_negative_pruned[list(mid_sequences_negative_pruned.keys())[0]]['sequence']

        # rating = 1
        # if maximum_positive_traces_in_a_sequence != -1:
        #     rating = 1 - (
        #         maximum_positive_traces_in_a_sequence / self.removal_candidates_generator.positive_traces_frequency)
        return (dynamic_sequences_negative_pruned, pre_sequences_negative_pruned, post_sequences_negative_pruned,
                mid_sequences_negative_pruned, execution_sequence_of_child_subtrees_negative)

    def handle_loop_operator(
        self, candidate: Candidate, tree_to_update: ProcessTree
    ):
        loop_rule_result = None
        loop_update_rule = LoopUpdateRule(candidate.removal_candidate_subtree)

        do_frequency_remove, redo_frequency_remove = (
            self.removal_candidates_generator.get_do_and_redo_frequencies_negative_variant(
                candidate.removal_candidate_subtree))

        if candidate.update_rule == 'remove_redundant_redo':
            try:
                loop_rule_result = self.calculate_candidate_tree_statistics(
                    loop_update_rule.apply_remove_redundant_redo(copy.deepcopy(tree_to_update)),
                    candidate.removal_candidate_subtree,
                    'loop: remove-redundant-redo'
                )

            except Exception as e:
                print("Rule application failed (loop: remove-redundant-redo): ", e)

        elif candidate.update_rule == 'optional_redo_mandatory':
            try:
                loop_rule_result = self.calculate_candidate_tree_statistics(
                    loop_update_rule.apply_optional_redo_mandatory(copy.deepcopy(tree_to_update)),
                    candidate.removal_candidate_subtree,
                    'loop: optional_redo_mandatory'
                )

            except Exception as e:
                print("Rule application failed (loop: optional_redo_mandatory): ", e)

        elif candidate.update_rule == 'repeat_exactly_n':
            try:
                repetitions_to_encode: list[int] = list(
                    candidate.removal_candidate_subtree.loop_subtree_stats.count_positive_variants_highest_loop_repetitions.keys()
                )

                loop_rule_result = self.calculate_candidate_tree_statistics(
                    loop_update_rule.apply_repeat_exactly_n(copy.deepcopy(tree_to_update), repetitions_to_encode),
                    candidate.removal_candidate_subtree,
                    'loop: repeat_exactly_n ' + str(repetitions_to_encode)
                )

            except Exception as e:
                print("Rule application failed (loop: repeat_exactly_n): ", e)

        elif candidate.update_rule == 'repeat_at_least_n':
            try:
                loop_rule_result = self.calculate_candidate_tree_statistics(
                    loop_update_rule.apply_repeat_at_least_n(copy.deepcopy(tree_to_update),
                                                             candidate.removal_candidate_subtree.loop_subtree_stats.do_frequency_remove),
                    candidate.removal_candidate_subtree,
                    'loop: repeat_at_least_n ' + str(
                        candidate.removal_candidate_subtree.loop_subtree_stats.do_frequency_remove)
                )

            except Exception as e:
                print("Rule application failed (loop: repeat_at_least_n): ", e)

        return loop_rule_result

    def handle_parallel_operator(
        self, candidate: Candidate, tree_to_update: ProcessTree
    ) -> ProcessTree:
        parallel_rules_result = None
        parallel_update_rule = ParallelUpdateRule(candidate.removal_candidate_subtree)

        if candidate.update_rule == 'pre_sequence':
            try:
                parallel_rules_result = self.calculate_candidate_tree_statistics(
                    parallel_update_rule.apply_pre_sequeltialization_rule(copy.deepcopy(tree_to_update),
                                                                          candidate.parameters['sequence']),
                    candidate.removal_candidate_subtree,
                    'parallel: pre sequentialization ' + str(candidate.parameters['sequence'])
                )

            except Exception as e:
                print("Rule application failed (parallel: pre sequentialization): ", e)

        elif candidate.update_rule == 'post_sequence':
            try:
                parallel_rules_result = self.calculate_candidate_tree_statistics(
                    parallel_update_rule.apply_post_sequeltialization_rule(copy.deepcopy(tree_to_update),
                                                                           candidate.parameters['sequence']),
                    candidate.removal_candidate_subtree,
                    'parallel: post sequentialization ' + str(candidate.parameters['sequence'])
                )

            except Exception as e:
                print("Rule application failed (parallel: post sequentialization): ", e)

        elif candidate.update_rule == 'mid_sequence':
            try:
                parallel_rules_result = self.calculate_candidate_tree_statistics(
                    parallel_update_rule.apply_mid_sequeltialization_rule(copy.deepcopy(tree_to_update),
                                                                          candidate.parameters['sequence'],
                                                                          candidate.parameters[
                                                                              'execution_sequence_of_child_subtrees_negative']),
                    candidate.removal_candidate_subtree,
                    'parallel: mid sequentialization ' + str(candidate.parameters['sequence'])
                )

            except Exception as e:
                print("Rule application failed (parallel: mid sequentialization): ", e)

        elif candidate.update_rule == 'dynamic_sequence':
            try:
                parallel_rules_result = self.calculate_candidate_tree_statistics(
                    parallel_update_rule.apply_dynamic_sequeltialization_rule(copy.deepcopy(tree_to_update),
                                                                              candidate.parameters['sequence']),
                    candidate.removal_candidate_subtree,
                    'parallel: dynamic sequentialization ' + str(candidate.parameters['sequence'])
                )

            except Exception as e:
                print("Rule application failed (parallel: dynamic sequentialization): ", e)

        return parallel_rules_result

