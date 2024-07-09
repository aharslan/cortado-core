import copy
from pm4py.objects.process_tree.obj import Operator, ProcessTree
from cortado_core.freezing.reinsert_frozen_subtrees import post_process_tree
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_activity import \
    CandidateActivity
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_subtree import \
    CandidateSubtree
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.removal_candidates_heuristics import \
    RemovalCandidatesHeuristics
from cortado_core.utils.alignment_utils import typed_trace_fits_process_tree
from cortado_core.negative_process_model_repair.temp_utils import calculate_process_tree_edit_distance
from cortado_core.negative_process_model_repair.temp_utils import (
    calculate_percentage_traces_conforming
)


class SubtreeUpdate:

    removal_candidates_generator: RemovalCandidatesHeuristics
    removal_candidate_activities: list[CandidateActivity]
    removal_candidate_subtrees: list[CandidateSubtree]

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
        }

        try:
            if (
                process_tree != self.removal_candidates_generator.process_tree
                and process_tree is not None
            ):
                post_process_tree(process_tree, [])
                # apply_reduction_rules(candidate_updated_tree)

                tree_updated = True

                # step 2: check if the negative variant conforms to the tree, stop removing candidates
                negative_trace_fits = typed_trace_fits_process_tree(
                    self.removal_candidates_generator.negative_trace, process_tree
                )

                percentage_positive_traces_conforming = (
                    calculate_percentage_traces_conforming(
                        self.removal_candidates_generator.fitting_traces,
                        process_tree,
                    )
                )
                resulting_tree_edit_distance = calculate_process_tree_edit_distance(
                    process_tree,
                    copy.deepcopy(self.removal_candidates_generator.process_tree)
                )

                statistics['negative_trace_fits'] = negative_trace_fits
                statistics['percentage_positive_traces_conforming'] = percentage_positive_traces_conforming
                statistics['resulting_tree_edit_distance'] = resulting_tree_edit_distance

        except Exception as e:
            print(e)

        return statistics
