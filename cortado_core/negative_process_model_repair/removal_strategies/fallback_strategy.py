import copy
import math
from pm4py import ProcessTree
from cortado_core.freezing.reinsert_frozen_subtrees import post_process_tree
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_activity import \
    CandidateActivity
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.removal_candidates_generator import \
    RemovalCandidatesGenerator
from cortado_core.negative_process_model_repair.temp_utils import (
    calculate_percentage_traces_conforming, calculate_process_tree_edit_distance
)
from cortado_core.process_tree_utils.miscellaneous import is_tau_leaf, get_root
from cortado_core.utils.alignment_utils import typed_trace_fits_process_tree


class FallbackStrategy:
    """Class responsible for the fallback strategy for removal of
    a trace from the given sublogs including process tree and
    candidate activities"""

    removal_candidates_generator: RemovalCandidatesGenerator
    removal_candidate_activities: list[CandidateActivity]

    def __init__(
        self,
        removal_candidates_generator: RemovalCandidatesGenerator,
        removal_candidate_activities: list[CandidateActivity],
    ):
        self.removal_candidates_generator: RemovalCandidatesGenerator = (
            removal_candidates_generator
        )
        self.removal_candidate_activities: list[
            CandidateActivity
        ] = removal_candidate_activities

    # def apply_fallback_strategy(
    #     self, max_frequency_threshold=0, allow_breach_threshold=False
    # ):
    #     """Removes the least frequent candidate in the list from the tree. Starts by removing all
    #     having frequency 0, then removing all having frequency equal to max_frequency_threshold and
    #     if no candidates with frequency less than or equal to max_frequency_threshold, removes the
    #     first (least frequent) candidate from the list of candidates if allow_breach_threshold is True
    #     """
    #
    #     updated = False
    #     tree_to_update = copy.deepcopy(self.sublogs_projection.process_tree)
    #
    #     for candidate in self.removal_candidates:
    #         if candidate.frequency == 0:
    #             # step 1: remove all activities having frequency 0
    #             tree_to_update = remove_activity_from_tree_by_id(
    #                 candidate.parent_id, candidate.activity_name, tree_to_update
    #             )
    #             updated = True
    #             print("Removed candidates having frequency: 0")
    #
    #         elif candidate.frequency == max_frequency_threshold:
    #             # step 2: remove all activities having frequency == max_frequency_threshold
    #             tree_to_update = remove_activity_from_tree_by_id(
    #                 candidate.parent_id, candidate.activity_name, tree_to_update
    #             )
    #             updated = True
    #             print("Removed candidates having frequency: max_frequency_threshold")
    #
    #         elif candidate.frequency > max_frequency_threshold:
    #             break
    #
    #     if not updated and allow_breach_threshold:
    #         # step 2: removes first (least frequent) activity in the list of candidates if no updates were applied
    #         tree_to_update = self.remove_candidate_activity_from_tree(
    #             self.removal_candidates[0],
    #             tree_to_update,
    #         )
    #         updated = True
    #         print("Breaching threshold: removing the least frequent candidate with frequency: " + str(
    #             self.removal_candidates[0].frequency))
    #
    #     if updated:
    #         post_process_tree(tree_to_update, [])
    #     return tree_to_update

    def apply_fallback_strategy(self) -> (ProcessTree, bool, float):
        """
        Keep removing the least frequent activities in list of candidates from the tree
        until the variant to remove is no longer fitting the given process tree.
        """

        tree_to_update = copy.deepcopy(self.removal_candidates_generator.process_tree)

        top_candidate = {
            "candidate_updated_tree": copy.deepcopy(self.removal_candidates_generator.process_tree),
            "tree_updated": False,
            "percentage_positive_traces_conforming": 0,
            "resulting_tree_edit_distance": -1
        }

        # removal_candidates_grouped = group_list_to_dict(
        #     self.removal_candidate_activities
        # )
        # minimum_frequency_in_positive_traces = next(iter(removal_candidates_grouped))
        candidate_trees: [any] = []

        position_temp = 0
        top_removal_candidate_activities = self.removal_candidate_activities[
                                           0: math.ceil(len(self.removal_candidate_activities) / 3)
                                           ]
        for candidate in top_removal_candidate_activities:
            # step 1: remove first (least frequent) activity in the list of candidates
            candidate_updated_tree = copy.deepcopy(tree_to_update)
            candidate_updated_tree = self.remove_candidate_activity_from_tree(
                candidate,
                candidate_updated_tree,
            )
            tree_updated: bool = False
            post_process_tree(candidate_updated_tree, [])
            if (
                candidate_updated_tree != self.removal_candidates_generator.process_tree
                and candidate_updated_tree is not None
            ):
                # apply_reduction_rules(candidate_updated_tree)
                tree_updated = True

            # step 2: check if the negative variant conforms to the tree, stop removing candidates
            negative_trace_fits = True

            try:
                negative_trace_fits = typed_trace_fits_process_tree(
                    self.removal_candidates_generator.negative_trace, candidate_updated_tree
                )
            except Exception as e:
                print("Invalid tree after removing fallback candidate and post-processing-")
                print(e)

            if not negative_trace_fits:
                percentage_positive_traces_conforming = (
                    calculate_percentage_traces_conforming(
                        self.removal_candidates_generator.fitting_traces,
                        candidate_updated_tree,
                    )
                )
                resulting_tree_edit_distance = calculate_process_tree_edit_distance(
                    candidate_updated_tree,
                    copy.deepcopy(self.removal_candidates_generator.process_tree)
                )
                candidate_trees.append(
                    {
                        "position_temp": position_temp,
                        "candidate": candidate,
                        "negative_trace_fits": negative_trace_fits,
                        "percentage_positive_traces_conforming": percentage_positive_traces_conforming,
                        "resulting_tree_edit_distance": resulting_tree_edit_distance,
                        "candidate_updated_tree": copy.deepcopy(candidate_updated_tree),
                        "tree_updated": tree_updated,
                    }
                )

            position_temp += 1
        candidate_trees = sorted(candidate_trees,
                                 key=lambda x: (-x["percentage_positive_traces_conforming"], x["resulting_tree_edit_distance"])
                                 )

        if len(candidate_trees) > 0:
            top_candidate = candidate_trees[0]

        return (
            top_candidate["candidate_updated_tree"],
            top_candidate["tree_updated"],
            top_candidate["percentage_positive_traces_conforming"],
            top_candidate["resulting_tree_edit_distance"],
        )

    def remove_candidate_activity_from_tree(
        self, candidate: CandidateActivity, tree_to_update: ProcessTree
    ):
        tree_to_update = remove_activity_from_tree_by_id(
            candidate.parent_id,
            candidate.activity_name,
            tree_to_update,
        )
        return tree_to_update


def remove_activity_from_tree_by_id(
    parent_id: int, activity_name: str, pt: ProcessTree, replace_with_tau: bool = False
) -> ProcessTree:
    if pt.id == parent_id:
        return delete_activity(activity_name, pt, replace_with_tau)
    elif len(pt.children) > 0:
        i = 0
        for child in pt.children:
            if parent_id > child.id:
                i += 1
                if len(pt.children) == i:
                    return remove_activity_from_tree_by_id(
                        parent_id, activity_name, child
                    )
                else:
                    continue
            elif parent_id < child.id:
                return remove_activity_from_tree_by_id(
                    parent_id, activity_name, pt.children[i - 1]
                )
            elif child.id == parent_id:
                return remove_activity_from_tree_by_id(parent_id, activity_name, child)
            else:
                # operator node not found
                return get_root(pt)
    else:
        # operator node not found
        return get_root(pt)


def delete_activity(
    activity_name: str, pt: ProcessTree, replace_with_tau: bool = False
):
    for child in pt.children:
        if child.label == activity_name and not replace_with_tau:
            pt.children.remove(child)

        elif child.label == activity_name and replace_with_tau:
            child.children = []
            child.operator = None
            child.label = None

        elif is_tau_leaf(child) and activity_name == "tau":
            pt.children.remove(child)

    return get_root(pt)


def group_list_to_dict(lst):
    groups = {}
    for item in lst:
        if item.frequency_in_positive_traces in groups:
            groups[item.frequency_in_positive_traces].append(item)
        else:
            groups[item.frequency_in_positive_traces] = [item]
    return groups
