import operator
import copy
from pm4py import ProcessTree
from cortado_core.negative_process_model_repair.neg_repair_sublog_utils import \
    calculate_sublog_for_negative_process_model_repair_with_frequency
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_activity import \
    CandidateActivity
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_subtree import \
    CandidateSubtree
from cortado_core.negative_process_model_repair.temp_utils import is_tau, find_tree_node_by_id, \
    is_activity_in_children, get_activity_reference_in_children
from cortado_core.utils.trace import TypedTrace


class RemovalCandidatesGenerator:
    process_tree: ProcessTree
    fitting_traces: [TypedTrace]
    negative_trace: {TypedTrace, None}
    negative_trace_frequency: int
    positive_traces_frequency: int

    def __init__(self, p_t, f_t, n_t, n_t_f, p_t_f):
        self.keep_sublogs = None
        self.remove_sublogs = None
        self.process_tree = p_t
        self.fitting_traces = f_t
        self.negative_trace = n_t
        self.negative_trace_frequency = n_t_f
        self.positive_traces_frequency = p_t_f

    def generate_sublogs_for_traces(
        self,
        pool=None
    ):
        try:
            self.keep_sublogs = (
                calculate_sublog_for_negative_process_model_repair_with_frequency(
                    self.process_tree, self.fitting_traces, pool
                )
            )
            self.remove_sublogs = (
                calculate_sublog_for_negative_process_model_repair_with_frequency(
                    self.process_tree,
                    [self.negative_trace],
                    pool,
                )
            )
        except Exception as e:
            self.keep_sublogs = None
            self.remove_sublogs = None
            print("Error calculating sublogs: " + str(e))

    def generate_removal_candidate_activities(self):
        removal_candidates: list[CandidateActivity] = []
        removal_candidates_set = (
            set()
        )  # maintain unique names of activities in negative trace

        if self.keep_sublogs is not None and self.remove_sublogs is not None:
            for tree_node_id in self.remove_sublogs:
                for remove_sublog_traces in self.remove_sublogs[tree_node_id]:
                    tree_node = find_tree_node_by_id(tree_node_id, self.process_tree)

                    candidate_activity: str
                    if len(remove_sublog_traces) > 0:
                        for activity in remove_sublog_traces:
                            candidate_activity = activity["concept:name"]

                            if is_activity_in_children(candidate_activity, tree_node):
                                """Only add if activity leaf node is a direct child"""

                                if (
                                    candidate_activity + "_" + str(tree_node_id)
                                    not in removal_candidates_set
                                ):
                                    (
                                        frequency_in_positive_traces,
                                        number_of_positive_traces_with_activity,
                                    ) = get_activity_frequency_in_sublogs(
                                        candidate_activity,
                                        tree_node_id,
                                        self.keep_sublogs
                                    )

                                    (
                                        frequency_in_negative_trace,
                                        number_of_negative_traces_with_activity,
                                    ) = get_activity_frequency_in_sublogs(
                                        candidate_activity,
                                        tree_node_id,
                                        self.remove_sublogs
                                    )

                                    removal_candidates.append(
                                        CandidateActivity(
                                            tree_node_id,
                                            tree_node.operator,
                                            candidate_activity,
                                            frequency_in_positive_traces,
                                            frequency_in_negative_trace,
                                            number_of_positive_traces_with_activity,
                                            get_activity_reference_in_children(
                                                candidate_activity, tree_node
                                            ),
                                            self.positive_traces_frequency,
                                        )
                                    )
                                    removal_candidates_set.add(
                                        candidate_activity + "_" + str(tree_node_id)
                                    )

                                    # if operator is loop and sibling is tau, add it also to the set of candidates
                                    if tree_node.operator.name == "LOOP" and (
                                        is_tau(tree_node.children[0])
                                        or is_tau(tree_node.children[1])
                                    ):
                                        tau_candidate = self.get_tau_on_loop_to_be_considered_as_candidate(
                                            tree_node_id,
                                            tree_node,
                                            frequency_in_positive_traces,
                                            frequency_in_negative_trace,
                                            number_of_positive_traces_with_activity,
                                        )
                                        if (
                                            tau_candidate.frequency_in_negative_trace
                                            != 0
                                        ):
                                            removal_candidates.append(tau_candidate)
                                            removal_candidates_set.add(
                                                "tau" + "_" + str(tree_node_id)
                                            )

                    elif len(remove_sublog_traces) == 0:
                        candidate_activity = "tau"

                        if is_activity_in_children(candidate_activity, tree_node):
                            """Only add if activity leaf node is a direct child"""

                            if (
                                candidate_activity + "_" + str(tree_node_id)
                                not in removal_candidates_set
                            ):
                                (
                                    frequency_in_positive_traces,
                                    number_of_positive_traces_with_activity,
                                ) = get_activity_frequency_in_sublogs(
                                    candidate_activity, tree_node_id, self.keep_sublogs
                                )

                                (
                                    frequency_in_negative_trace,
                                    number_of_negative_traces_with_activity,
                                ) = get_activity_frequency_in_sublogs(
                                    candidate_activity,
                                    tree_node_id,
                                    self.remove_sublogs
                                )

                                # check if tau on left of loop, count redo executions in positive variants and add it to number_of_positive_traces_with_activity
                                if tree_node.operator.name == "LOOP" and not is_tau(tree_node.children[1]):
                                    right_label = tree_node.children[1].label

                                    redo_frequency_in_positive_traces, redo_number_of_positive_traces_with_activity = (
                                        get_activity_frequency_in_sublogs(right_label, tree_node_id, self.keep_sublogs)
                                    )
                                    frequency_in_positive_traces += redo_frequency_in_positive_traces
                                    number_of_positive_traces_with_activity += redo_number_of_positive_traces_with_activity

                                removal_candidates.append(
                                    CandidateActivity(
                                        tree_node_id,
                                        tree_node.operator,
                                        candidate_activity,
                                        frequency_in_positive_traces,
                                        frequency_in_negative_trace,
                                        number_of_positive_traces_with_activity,
                                        get_activity_reference_in_children(
                                            candidate_activity, tree_node
                                        ),
                                        self.positive_traces_frequency,
                                    )
                                )
                                removal_candidates_set.add(
                                    candidate_activity + "_" + str(tree_node_id)
                                )

        return sorted(
            removal_candidates, key=lambda x: x.presence_measure, reverse=False
        )

    def get_tau_on_loop_to_be_considered_as_candidate(
        self,
        tree_node_id,
        tree_node: ProcessTree,
        sibling_frequency_in_positive_traces,
        sibling_frequency_in_negative_trace,
        sibling_number_of_positive_traces_with_activity,
    ):
        tau_negative_frequency = 0
        tau_positive_frequency = 0
        tau_number_of_positive_traces_with_activity = (
            sibling_number_of_positive_traces_with_activity
        )
        tau_candidate: CandidateActivity

        if is_tau(tree_node.children[0]):
            tau_negative_frequency = sibling_frequency_in_negative_trace + 1
            tau_positive_frequency = (sibling_frequency_in_positive_traces + 1) + (
                sibling_number_of_positive_traces_with_activity - 1
            )
        elif is_tau(tree_node.children[1]):
            if sibling_frequency_in_negative_trace != 0:
                tau_negative_frequency = sibling_frequency_in_negative_trace - 1
                tau_positive_frequency = (sibling_frequency_in_positive_traces - 1) - (
                    sibling_number_of_positive_traces_with_activity - 1
                )

        tau_candidate = CandidateActivity(
            tree_node_id,
            tree_node.operator,
            "tau",
            tau_positive_frequency,
            tau_negative_frequency,
            tau_number_of_positive_traces_with_activity,
            get_activity_reference_in_children("tau", tree_node),
            self.positive_traces_frequency,
        )

        return tau_candidate

    def generate_removal_candidate_subtrees(
        self, removal_candidate_activities: list[CandidateActivity]
    ):
        removal_subtree_candidates: list[CandidateSubtree] = []
        removal_subtree_candidate_ids = set()

        for candidate_activity in removal_candidate_activities:
            if candidate_activity.parent_id not in removal_subtree_candidate_ids:
                candidate_subtree = CandidateSubtree()
                candidate_subtree.node_id = candidate_activity.parent_id
                candidate_subtree.reference = (
                    candidate_activity.activity_reference.parent
                )
                candidate_subtree.candidate_activities.append(candidate_activity)
                removal_subtree_candidates.append(candidate_subtree)
                removal_subtree_candidate_ids.add(candidate_activity.parent_id)

            elif candidate_activity.parent_id in removal_subtree_candidate_ids:
                for removal_subtree_candidate in removal_subtree_candidates:
                    if (
                        candidate_activity.parent_id
                        == removal_subtree_candidate.node_id
                    ):
                        removal_subtree_candidate.candidate_activities.append(
                            candidate_activity
                        )
                        break

        removal_subtree_candidates = sorted(
            removal_subtree_candidates,
            key=lambda x: x.node_id,
            reverse=True,
        )

        # 1. Add subtree candidates to parent candidates. If parent not a candidate, make it a candidate giving it
        # the average presence measure from the children subtrees and/or activities.
        for subtree_candidate in removal_subtree_candidates:
            if subtree_candidate.node_id != 0:
                parent_id = subtree_candidate.get_parent_id()
                if parent_id in removal_subtree_candidate_ids:
                    parent_subtree_candidate = get_subtree_candidate_by_id(
                        removal_subtree_candidates, parent_id
                    )
                    parent_subtree_candidate.candidate_subtrees.append(
                        subtree_candidate
                    )
                else:
                    candidate_subtree = CandidateSubtree()
                    candidate_subtree.node_id = parent_id
                    candidate_subtree.reference = subtree_candidate.reference.parent
                    candidate_subtree.candidate_subtrees.append(subtree_candidate)
                    removal_subtree_candidates.append(candidate_subtree)
                    removal_subtree_candidate_ids.add(parent_id)

        removal_subtree_candidates = sorted(
            removal_subtree_candidates,
            key=lambda x: x.node_id,
            reverse=True,
        )
        return copy.deepcopy(removal_subtree_candidates)


def get_activity_frequency_in_sublogs(activity_name: str, parent_id: int, sublogs):
    frequency: int = 0
    number_of_traces_with_activity = 0

    if parent_id in sublogs:
        for sublog_traces in sublogs[parent_id]:
            exists_in_trace = False

            if len(sublog_traces) > 0:
                for activity in sublog_traces:
                    if activity["concept:name"] == activity_name:
                        frequency += sublog_traces.attributes['frequency']
                        exists_in_trace = True
            elif len(sublog_traces) == 0 and activity_name == "tau":
                frequency += sublog_traces.attributes['frequency']
                exists_in_trace = True

            if exists_in_trace:
                number_of_traces_with_activity += sublog_traces.attributes['frequency']

    return frequency, number_of_traces_with_activity


def get_subtree_candidate_by_id(removal_subtree_candidates, node_id):
    for candidate in removal_subtree_candidates:
        if candidate.node_id == node_id:
            return candidate
    return None
