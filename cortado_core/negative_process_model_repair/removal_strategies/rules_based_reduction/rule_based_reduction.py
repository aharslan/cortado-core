import copy
from pm4py.objects.process_tree.obj import Operator, ProcessTree
from cortado_core.freezing.reinsert_frozen_subtrees import post_process_tree
from cortado_core.negative_process_model_repair.constants import Constants
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_activity import \
    CandidateActivity
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_subtree import \
    CandidateSubtree
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.removal_candidates_generator import \
    RemovalCandidatesGenerator
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.removal_candidates_heuristics import \
    RemovalCandidatesHeuristics
from cortado_core.negative_process_model_repair.removal_strategies.fallback_strategy import (
    remove_activity_from_tree_by_id,
)
from cortado_core.utils.alignment_utils import typed_trace_fits_process_tree
from cortado_core.negative_process_model_repair.temp_utils import is_tau, calculate_process_tree_edit_distance, \
    find_tree_node_by_id, remove_and_return_child_activity, remove_and_return_child_subtree
from cortado_core.negative_process_model_repair.temp_utils import (
    calculate_percentage_traces_conforming,
)


class RuleBasedReduction:
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
        removal_candidate_subtrees = (
            self.removal_candidates_generator.generate_removal_candidate_subtrees(
                removal_candidate_activities
            )
        )
        self.removal_candidate_subtrees = (
            self.removal_candidates_generator.apply_trace_frequencies_based_heuristics_to_rate_subtrees
            (removal_candidate_subtrees)
        )

    def apply_rule_based_reduction(self) -> (ProcessTree, bool, float):
        tree_to_update = copy.deepcopy(self.removal_candidates_generator.process_tree)
        tree_updated: bool = False
        applied_rules = []
        percentage_positive_traces_conforming = 100.0

        for removal_candidate_subtree in self.removal_candidate_subtrees:
            if removal_candidate_subtree.min_pm_approach == 'subtree':
                continue

            if removal_candidate_subtree.reference.operator == Operator.SEQUENCE:
                tree_to_update = handle_sequence_operator(
                    removal_candidate_subtree, tree_to_update
                )
                applied_rules.append('sequence')

            elif removal_candidate_subtree.reference.operator == Operator.XOR:
                tree_to_update = handle_choice_operator(
                    removal_candidate_subtree, tree_to_update
                )
                applied_rules.append('choice')

            elif (
                removal_candidate_subtree.reference.operator == Operator.PARALLEL
                and len(removal_candidate_subtree.candidate_patterns) > 0
            ):
                tree_to_update = handle_parallel_operator(
                    removal_candidate_subtree, tree_to_update
                )
                applied_rules.append('parallel')

            elif removal_candidate_subtree.reference.operator == Operator.LOOP:
                tree_to_update = handle_loop_operator(
                    removal_candidate_subtree, tree_to_update
                )
                applied_rules.append('loop')

            if (
                tree_to_update != self.removal_candidates_generator.process_tree
                and tree_to_update is not None
            ):
                post_process_tree(tree_to_update, [])
                # apply_reduction_rules(candidate_updated_tree)

                tree_updated = True

            # step 2: check if the negative variant conforms to the tree, stop removing candidates
            negative_trace_fits = typed_trace_fits_process_tree(
                self.removal_candidates_generator.negative_trace, tree_to_update
            )

            if not negative_trace_fits:
                break

        percentage_positive_traces_conforming = (
            calculate_percentage_traces_conforming(
                self.removal_candidates_generator.fitting_traces,
                tree_to_update,
            )
        )
        resulting_tree_edit_distance = calculate_process_tree_edit_distance(
            tree_to_update,
            copy.deepcopy(self.removal_candidates_generator.process_tree)
        )
        return (tree_to_update, tree_updated, percentage_positive_traces_conforming,
                resulting_tree_edit_distance, applied_rules)


def handle_sequence_operator(
    removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
) -> ProcessTree:
    """Handle the sequence operator. The function removes all
    candidate activities from the identified sequence operator
    Test: receipt log first 20 variants"""

    for candidate_activity in removal_candidate_subtree.candidate_activities:
        tree_to_update = remove_activity_from_tree_by_id(
            candidate_activity.parent_id,
            candidate_activity.activity_name,
            tree_to_update,
        )

    return tree_to_update


def handle_choice_operator(
    removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
) -> ProcessTree:
    """Handle the choice operator"""

    # desirable behavior or activity no longer optional case:
    if (
        len(removal_candidate_subtree.candidate_activities) == 1
        and removal_candidate_subtree.candidate_activities[0].activity_name == "tau"
    ):
        tree_to_update = remove_activity_from_tree_by_id(
            removal_candidate_subtree.candidate_activities[0].parent_id,
            removal_candidate_subtree.candidate_activities[0].activity_name,
            tree_to_update,
        )

    # undesirable behavior or activity replaced with tau case:
    elif (
        len(removal_candidate_subtree.candidate_activities) == 1
        and removal_candidate_subtree.candidate_activities[0].activity_name != "tau"
    ):
        tree_to_update = remove_activity_from_tree_by_id(
            removal_candidate_subtree.candidate_activities[0].parent_id,
            removal_candidate_subtree.candidate_activities[0].activity_name,
            tree_to_update,
            True,
        )

    return tree_to_update


def handle_loop_operator(
    removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
):
    if (
        removal_candidate_subtree.loop_subtree_stats.min_pm_approach
        == "remove_redundant_redo"
    ):
        remove_redundant_redo(removal_candidate_subtree, tree_to_update)

    elif (
        removal_candidate_subtree.loop_subtree_stats.min_pm_approach
        == "optional_redo_mandatory"
    ):
        optional_redo_mandatory(removal_candidate_subtree, tree_to_update)

    elif (
        removal_candidate_subtree.loop_subtree_stats.min_pm_approach
        == "repeat_exactly_n"
    ):
        repeat_exactly_n(removal_candidate_subtree, tree_to_update)

    elif (
        removal_candidate_subtree.loop_subtree_stats.min_pm_approach
        == "repeat_at_least_n"
    ):
        repeat_atleast_n(removal_candidate_subtree, tree_to_update)

    return tree_to_update


def remove_redundant_redo(removal_candidate_subtree, tree_to_update: ProcessTree):
    tree_node_to_update = find_tree_node_by_id(
        removal_candidate_subtree.node_id, tree_to_update
    )
    temp_node: ProcessTree = removal_candidate_subtree.reference.children[0]

    for i in range(len(tree_node_to_update.parent.children)):
        if tree_node_to_update.parent.children[i].id == tree_node_to_update.id:
            tree_node_to_update.parent.children[i] = temp_node
            temp_node.parent = tree_node_to_update.parent
            break

    return tree_to_update


def optional_redo_mandatory(
    removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
):
    tree_node_to_update = find_tree_node_by_id(
        removal_candidate_subtree.node_id, tree_to_update
    )

    if is_tau(tree_node_to_update.children[0]):
        left_child = copy.deepcopy(tree_node_to_update.children[0])

        tree_node_to_update.children[0] = copy.deepcopy(tree_node_to_update.children[1])
        tree_node_to_update.children[0].parent = tree_node_to_update

        tree_node_to_update.children[1] = left_child
        tree_node_to_update.children[1].parent = tree_node_to_update

    return tree_to_update


def repeat_exactly_n(
    removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
):
    tree_node_to_update = find_tree_node_by_id(
        removal_candidate_subtree.node_id, tree_to_update
    )

    choice_node = ProcessTree()
    choice_node.operator = Operator.XOR
    choice_node.children = []

    index_in_parents_children = 0
    for i in range(len(tree_node_to_update.parent.children)):
        if tree_node_to_update.parent.children[i].id == tree_node_to_update.id:
            index_in_parents_children = i
            break

    for key, value in removal_candidate_subtree.loop_subtree_stats.count_positive_variants_highest_loop_repetitions.items():
        seq_node = ProcessTree()
        seq_node.operator = Operator.SEQUENCE
        seq_node.parent = choice_node
        seq_node.children = []

        left_child_loop = copy.deepcopy(tree_node_to_update.children[0])
        left_child_loop.parent = seq_node
        right_child_loop = copy.deepcopy(tree_node_to_update.children[1])
        right_child_loop.parent = seq_node

        if key >= 1:
            if key == 1:
                seq_node.children.append(left_child_loop)
            elif key > 1:
                for i in range(key - 1):
                    seq_node.children.append(left_child_loop)
                    seq_node.children.append(right_child_loop)
                seq_node.children.append(left_child_loop)

            choice_node.children.append(seq_node)

    choice_node.parent = tree_node_to_update.parent
    tree_node_to_update.parent.children[index_in_parents_children] = choice_node

    return tree_to_update


def repeat_atleast_n(
    removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
):
    tree_node_to_update = find_tree_node_by_id(
        removal_candidate_subtree.node_id, tree_to_update
    )

    tree_node_to_update = find_tree_node_by_id(
        removal_candidate_subtree.node_id, tree_to_update
    )
    temp_node = ProcessTree()
    temp_node.operator = Operator.SEQUENCE
    temp_node.children = []

    index_in_parents_children = 0
    for i in range(len(tree_node_to_update.parent.children)):
        if tree_node_to_update.parent.children[i].id == tree_node_to_update.id:
            index_in_parents_children = i
            break

    left_child_loop = copy.deepcopy(tree_node_to_update.children[0])
    left_child_loop.parent = temp_node
    right_child_loop = copy.deepcopy(tree_node_to_update.children[1])
    right_child_loop.parent = temp_node

    removal_candidate_subtree.loop_subtree_stats.do_frequency_remove
    for i in range(removal_candidate_subtree.loop_subtree_stats.do_frequency_remove):
        temp_node.children.append(copy.deepcopy(left_child_loop))
        temp_node.children.append(copy.deepcopy(right_child_loop))

    loop_node = copy.deepcopy(tree_node_to_update)
    loop_node.parent = temp_node
    temp_node.children.append(loop_node)

    temp_node.parent = tree_node_to_update.parent
    tree_node_to_update.parent.children[index_in_parents_children] = temp_node

    return tree_to_update


def handle_parallel_operator(
    removal_candidate_subtree: CandidateSubtree, tree_to_update: ProcessTree
) -> ProcessTree:
    tree_updated = copy.deepcopy(tree_to_update)
    tree_update_reference_node = find_tree_node_by_id(
        removal_candidate_subtree.node_id, tree_updated
    )

    for i in range(len(removal_candidate_subtree.candidate_patterns)):
        alternate_combination = removal_candidate_subtree.candidate_patterns[i][
            "alternate-combinations"
        ][0]["combination"]
        sequence: ProcessTree = ProcessTree()
        sequence.operator = Operator.SEQUENCE
        sequence.parent = tree_update_reference_node
        sequence.label = None

        for j in range(len(alternate_combination)):
            child = None
            if (
                alternate_combination[j]["parent-id"]
                == removal_candidate_subtree.node_id
            ):
                child = remove_and_return_child_activity(
                    tree_update_reference_node, alternate_combination[j]["names"][0]
                )

            elif (
                alternate_combination[j]["parent-id"]
                != removal_candidate_subtree.node_id
            ):
                child = remove_and_return_child_subtree(
                    tree_update_reference_node, alternate_combination[j]["parent-id"]
                )
            if child is not None:
                child.parent = sequence
                sequence.children.append(child)

        tree_update_reference_node.children.append(sequence)
        break

    return tree_updated
