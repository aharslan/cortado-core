import copy
from dataclasses import dataclass
from typing import List
from pm4py import ProcessTree
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_activity import \
    CandidateActivity
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.loop_subtree_stats import \
    LoopSubtreeStats


@dataclass
class CandidateSubtree:
    node_id: int
    average_presence_measure: float
    candidate_activities: List[CandidateActivity]
    candidate_subtrees: List[any]
    candidate_patterns: List[any]
    loop_subtree_stats: LoopSubtreeStats
    reference: ProcessTree

    min_pm: float
    min_pm_approach: str

    rating: float = None
    rating_approach: str = None
    sequence_in_parallel: list[int] = None
    execution_sequence_of_child_subtrees_negative: list[int] = None

    def __init__(self):
        self.node_id = None
        self.candidate_activities = []
        self.candidate_subtrees = []
        self.candidate_patterns = []
        self.loop_subtree_stats = None

    def get_child_candidate_id(self):
        if len(self.candidate_activities) > 0:
            return self.candidate_activities[0].activity_reference.id
        elif len(self.candidate_subtrees) > 0:
            return self.candidate_subtrees[0].node_id

    def set_average_presence_measure(self):
        avg_activity_presence_measure = 0
        avg_subtree_presence_measure = 0
        avg_pattern_presence_measure = 0
        average_loop_presence_measure = 0

        self.min_pm = 100000000

        if self.loop_subtree_stats is not None:
            average_loop_presence_measure = self.loop_subtree_stats.min_pm
            if self.min_pm > average_loop_presence_measure:
                self.min_pm = average_loop_presence_measure
                self.min_pm_approach = 'loop'

        if len(self.candidate_patterns) > 0:
            avg_pattern_presence_measure = self.candidate_patterns[0]['presence-measure']
            if self.min_pm > avg_pattern_presence_measure:
                self.min_pm = avg_pattern_presence_measure
                self.min_pm_approach = 'pattern'

        removal_activities = copy.deepcopy(self.candidate_activities)
        # count activities if not present in the pattern
        remove_indexes = set()
        for i in range(len(self.candidate_patterns)):
            for j in range(len(self.candidate_patterns[i]["combination"])):
                for k in range(
                    len(self.candidate_patterns[i]["combination"][j]["names"])
                ):
                    for l in range(len(removal_activities)):
                        if (
                            self.candidate_patterns[i]["combination"][j]["names"][k]
                            == removal_activities[l].activity_name
                        ):
                            remove_indexes.add(l)
                            break

        for i in sorted(remove_indexes, reverse=True):
            del removal_activities[i]

        if len(removal_activities) > 0:
            avg_activity_presence_measure = min(
                c.presence_measure for c in removal_activities
            )
            if self.min_pm > avg_activity_presence_measure:
                self.min_pm = avg_activity_presence_measure
                self.min_pm_approach = 'activity'

        if len(self.candidate_subtrees) > 0:
            avg_subtree_presence_measure = min(
                c.min_pm for c in self.candidate_subtrees
            )
            if self.min_pm > avg_subtree_presence_measure:
                self.min_pm = avg_subtree_presence_measure
                self.min_pm_approach = 'subtree'

        self.min_pm

    def get_parent_id(self, counter=0):
        if len(self.candidate_activities) > 0:
            parent_reference = self.candidate_activities[
                0
            ].activity_reference.parent.parent
            for i in range(counter):
                parent_reference = parent_reference.parent
            return parent_reference.id
        else:
            return self.candidate_subtrees[0].get_parent_id(counter + 1)

    def get_tree_node_reference(self, counter=0):
        if len(self.candidate_activities) > 0:
            reference = self.candidate_activities[0].activity_reference.parent
            for i in range(counter):
                reference = reference.parent
            return reference
        else:
            return self.candidate_subtrees[0].get_tree_node_reference(counter + 1)
