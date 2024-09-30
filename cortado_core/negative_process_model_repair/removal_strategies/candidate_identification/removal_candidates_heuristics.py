from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.removal_candidates_generator import \
    RemovalCandidatesGenerator
import copy
import math
from pm4py.objects.process_tree.obj import Operator
from cortado_core.negative_process_model_repair.constants import Constants
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_subtree import \
    CandidateSubtree
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.loop_subtree_stats import \
    LoopSubtreeStats
from cortado_core.negative_process_model_repair.temp_utils import is_tau, find_tree_node_by_id
from cortado_core.process_tree_utils.miscellaneous import is_subtree
import collections


class RemovalCandidatesHeuristics(RemovalCandidatesGenerator):

    def __init__(self, p_t, f_t, n_t, n_t_f, p_t_f):
        super().__init__(p_t, f_t, n_t, n_t_f, p_t_f)

    def apply_trace_frequencies_based_heuristics_to_rate_subtrees(self,
                                                                  removal_subtree_candidates: list[CandidateSubtree]):

        for subtree_candidate in removal_subtree_candidates:
            pt_reference = subtree_candidate.get_tree_node_reference()
            if pt_reference.operator == Operator.PARALLEL:
                remove_combinations_parallel = (
                    self.calculate_trace_frequencies_based_heuristics_for_parallel_operator(
                        subtree_candidate
                    )
                )
                subtree_candidate.candidate_patterns = remove_combinations_parallel

            elif pt_reference.operator == Operator.LOOP:
                loop_subtree_stats = (
                    self.calculate_trace_frequencies_based_heuristics_for_loop_operator(
                        subtree_candidate
                    )
                )
                subtree_candidate.loop_subtree_stats = loop_subtree_stats

            subtree_candidate.set_average_presence_measure()

        return sorted(
            removal_subtree_candidates,
            key=lambda x: x.min_pm,
            reverse=False,
        )

    def calculate_trace_frequencies_based_heuristics_for_loop_operator(
        self, subtree_candidate: CandidateSubtree
    ):
        loop_subtree_stats: LoopSubtreeStats = LoopSubtreeStats()

        # given loop operator, evaluate if the activities seen under the loop are do or redo for remove_sublogs and keep sublogs..
        loop_subtree_stats.loop_labels_annotated_remove = (
            self.annotate_loop_activities_with_do_or_redo(
                subtree_candidate, self.remove_sublogs
            )
        )
        if self.keep_sublogs:
            loop_subtree_stats.loop_labels_annotated_keep = (
                self.annotate_loop_activities_with_do_or_redo(
                    subtree_candidate, self.keep_sublogs
                )
            )

        (
            loop_subtree_stats.do_frequency_remove,
            loop_subtree_stats.redo_frequency_remove,
            x,
            loop_subtree_stats.do_frequency_trace_remove,
            loop_subtree_stats.redo_frequency_trace_remove
        ) = self.calculate_loop_statistics_from_annotated_labels(
            loop_subtree_stats.loop_labels_annotated_remove
        )

        (
            loop_subtree_stats.do_frequency_keep,
            loop_subtree_stats.redo_frequency_keep,
            loop_subtree_stats.loop_repitition_per_trace_keep,
            loop_subtree_stats.do_frequency_trace_keep,
            loop_subtree_stats.redo_frequency_trace_keep
        ) = self.calculate_loop_statistics_from_annotated_labels(
            loop_subtree_stats.loop_labels_annotated_keep
        )

        self.calculate_presence_measures_for_loop_scenarios(
            loop_subtree_stats, subtree_candidate
        )

        loop_subtree_stats.min_pm = 100
        loop_subtree_stats.min_pm_approach = ""
        if loop_subtree_stats.min_pm > loop_subtree_stats.remove_redundant_redo_pm:
            loop_subtree_stats.min_pm = loop_subtree_stats.remove_redundant_redo_pm
            loop_subtree_stats.min_pm_approach = "remove_redundant_redo"

        if loop_subtree_stats.min_pm > loop_subtree_stats.optional_redo_mandatory_pm:
            loop_subtree_stats.min_pm = loop_subtree_stats.optional_redo_mandatory_pm
            loop_subtree_stats.min_pm_approach = "optional_redo_mandatory"

        if loop_subtree_stats.min_pm >= loop_subtree_stats.repeat_exactly_n_pm:
            loop_subtree_stats.min_pm = loop_subtree_stats.repeat_exactly_n_pm
            loop_subtree_stats.min_pm_approach = "repeat_exactly_n"

        if loop_subtree_stats.min_pm >= loop_subtree_stats.repeat_atleast_n_pm:
            loop_subtree_stats.min_pm = loop_subtree_stats.repeat_atleast_n_pm
            loop_subtree_stats.min_pm_approach = "repeat_at_least_n"

        return loop_subtree_stats

    def calculate_presence_measures_for_loop_scenarios(
        self, loop_subtree_stats: LoopSubtreeStats, subtree_candidate: CandidateSubtree
    ):

        if loop_subtree_stats.redo_frequency_trace_remove > 0:
            loop_subtree_stats.remove_redundant_redo_pm = (
                loop_subtree_stats.redo_frequency_trace_keep
                / self.positive_traces_frequency
            )

        count_positive_traces_never_executing_redo = 0
        for i in range(len(loop_subtree_stats.loop_labels_annotated_keep)):
            if len(loop_subtree_stats.loop_labels_annotated_keep[i]) == 1:
                count_positive_traces_never_executing_redo += loop_subtree_stats.loop_labels_annotated_keep[i][0][
                    'frequency']

        if is_tau(
            subtree_candidate.reference.children[0]
        ) and loop_subtree_stats.redo_frequency_trace_remove == 0:
            loop_subtree_stats.optional_redo_mandatory_pm = (
                count_positive_traces_never_executing_redo
                / self.positive_traces_frequency
            )

        tally_loop_repititions_keep_traces = {}
        # find frequency of most common repetition of the loop in keep_traces != repition of loop in remove_trace:
        for i in range(len(loop_subtree_stats.loop_repitition_per_trace_keep)):
            if (
                loop_subtree_stats.loop_repitition_per_trace_keep[i]['repititions']
                in tally_loop_repititions_keep_traces
            ):
                tally_loop_repititions_keep_traces[
                    loop_subtree_stats.loop_repitition_per_trace_keep[i]['repititions']
                ] += loop_subtree_stats.loop_repitition_per_trace_keep[i]['trace-count']
            else:
                tally_loop_repititions_keep_traces[
                    loop_subtree_stats.loop_repitition_per_trace_keep[i]['repititions']
                ] = loop_subtree_stats.loop_repitition_per_trace_keep[i]['trace-count']

        tally_loop_repititions_keep_traces = collections.OrderedDict(
            sorted(
                tally_loop_repititions_keep_traces.items(),
                key=lambda kv: kv[1],
                reverse=True,
            )
        )
        most_frequent_n = 0
        tally_most_frequent_n = 0
        tally_loop_repititions_keep_traces_original = copy.deepcopy(tally_loop_repititions_keep_traces)

        # remove the n that corresponds to the loop repetitions in neg variant ..
        if loop_subtree_stats.do_frequency_remove in tally_loop_repititions_keep_traces:
            tally_loop_repititions_keep_traces.pop(
                loop_subtree_stats.do_frequency_remove)

        loop_subtree_stats.count_positive_variants_highest_loop_repetitions = copy.deepcopy(
            tally_loop_repititions_keep_traces)

        different_than_encoding_tuples = {}
        count_different_than_encoding_tuples = 0
        # remove individual repetitions exceeding Constants.MAX_LENGTH_LOOP_REPETITION_ENCODING
        for key in list(loop_subtree_stats.count_positive_variants_highest_loop_repetitions.keys()):
            if key > Constants.MAX_LENGTH_LOOP_REPETITION_ENCODING:
                different_than_encoding_tuples[key] = \
                loop_subtree_stats.count_positive_variants_highest_loop_repetitions[key]
                del loop_subtree_stats.count_positive_variants_highest_loop_repetitions[key]

        if (len(
            loop_subtree_stats.count_positive_variants_highest_loop_repetitions) >
            Constants.MAX_NUMBER_OF_LOOP_REPETITION_ENCODINGS):
            # trim to only encode the top frequent n's
            for i in range(
                len(loop_subtree_stats.count_positive_variants_highest_loop_repetitions) -
                Constants.MAX_NUMBER_OF_LOOP_REPETITION_ENCODINGS
            ):
                loop_subtree_stats.count_positive_variants_highest_loop_repetitions.popitem()

        for key, value in loop_subtree_stats.count_positive_variants_highest_loop_repetitions.items():
            tally_most_frequent_n += value

        positive_variants_having_lesser_repetitions_than_negative_variant = 0
        for key, value in tally_loop_repititions_keep_traces_original.items():
            if key not in loop_subtree_stats.count_positive_variants_highest_loop_repetitions:
                count_different_than_encoding_tuples += value
            if key <= loop_subtree_stats.do_frequency_remove:
                positive_variants_having_lesser_repetitions_than_negative_variant += value

        loop_subtree_stats.repeat_exactly_n_pm = (
            count_different_than_encoding_tuples
            / self.positive_traces_frequency
        )

        # count number of positive variants having greater repetitions than the negative variant repetitions

        loop_subtree_stats.repeat_atleast_n_pm = (
            positive_variants_having_lesser_repetitions_than_negative_variant
            / self.positive_traces_frequency
        )

    def calculate_loop_statistics_from_annotated_labels(self, loop_labels_annotated):
        do_frequency = 0
        redo_frequency = 0
        loop_frequency = 0
        redo_frequency_trace = 0

        loop_repitition_per_trace = []

        last_label = "redo"

        for i in range(len(loop_labels_annotated)):
            do_frequency_trace = 0
            do_frequency = 0
            is_redo__trace = False
            for j in range(len(loop_labels_annotated[i])):
                if (
                    loop_labels_annotated[i][j]["loop-location"] == "do"
                    and last_label != "do"
                ):
                    do_frequency += 1
                    do_frequency_trace += loop_labels_annotated[i][j]['frequency']
                    last_label = "do"
                elif (
                    loop_labels_annotated[i][j]["loop-location"] == "redo"
                    and last_label != "redo"
                ):
                    redo_frequency += 1
                    is_redo__trace = True
                    last_label = "redo"

            if is_redo__trace:
                redo_frequency_trace += loop_labels_annotated[i][j]['frequency']
            loop_repitition_per_trace.append(
                {'repititions': do_frequency, 'trace-count': loop_labels_annotated[i][j]['frequency']})
            last_label = "redo"

        return do_frequency, redo_frequency, loop_repitition_per_trace, do_frequency_trace, redo_frequency_trace

    def annotate_loop_activities_with_do_or_redo(
        self, subtree_candidate: CandidateSubtree, sublogs
    ):
        loop_activities = []

        index_k_redo = 0
        index_l_redo = 0
        index_k_do = 0
        index_l_do = 0

        try:
            if not subtree_candidate.node_id in sublogs:
                return [[{"name": "tau", "loop-location": "do", "frequency": 0}]]

            for i in range(len(sublogs[subtree_candidate.node_id])):
                loop_activities.append([])
                loop_location = "do"

                if len(sublogs[subtree_candidate.node_id][i]) > 0:
                    # for activity in sublogs[subtree_candidate.node_id][i]:
                    for j in range(len(sublogs[subtree_candidate.node_id][i])):
                        success = False
                        while not success:
                            success = True
                            activity = sublogs[subtree_candidate.node_id][i][j]
                            loop_activity = {
                                "name": activity["concept:name"],
                                "loop-location": "none",
                                "frequency": sublogs[subtree_candidate.node_id][i].attributes['frequency']
                            }

                            if (
                                subtree_candidate.reference.children[0].operator is None
                                and (
                                activity["concept:name"]
                                == subtree_candidate.reference.children[0].label
                            )
                                and loop_location == "do"
                            ):
                                loop_activity["loop-location"] = loop_location
                                loop_location = "redo"

                            elif (
                                subtree_candidate.reference.children[1].operator is None
                                and (
                                    activity["concept:name"]
                                    == subtree_candidate.reference.children[1].label
                                )
                                and loop_location == "redo"
                            ):
                                loop_activity["loop-location"] = loop_location
                                loop_location = "do"

                            # check if activity is in the left subtree
                            elif loop_location == "do":
                                found_in_do = False
                                is_tau = False

                                if subtree_candidate.reference.children[0].id in sublogs:
                                    for k in range(
                                        index_k_do,
                                        len(sublogs[subtree_candidate.reference.children[0].id]),
                                    ):
                                        match = False

                                        if (len(sublogs[subtree_candidate.reference.children[0].id][k]) > 0):
                                            for l in range(index_l_do,
                                                           len(sublogs[subtree_candidate.reference.children[0].id][k])):
                                                if (
                                                    sublogs[
                                                        subtree_candidate.reference.children[
                                                            0
                                                        ].id
                                                    ][k][l]["concept:name"]
                                                    == activity["concept:name"]
                                                ):
                                                    match = True
                                                    index_l_do = l + 1
                                                    break
                                        else:
                                            is_tau = True
                                            break

                                        if match:
                                            found_in_do = True
                                            loop_activity["loop-location"] = loop_location
                                            index_k_do = k

                                            if index_l_do == len(
                                                sublogs[
                                                    subtree_candidate.reference.children[
                                                        0
                                                    ].id
                                                ][k]
                                            ):
                                                index_k_do = k + 1
                                                index_l_do = 0
                                                loop_location = "redo"
                                                if index_k_do == len(
                                                    sublogs[
                                                        subtree_candidate.reference.children[
                                                            0
                                                        ].id
                                                    ]
                                                ):
                                                    index_k_do = 0
                                            break
                                else:
                                    is_tau = True

                                if is_tau:
                                    loop_activity["name"] = "tau"
                                    loop_activity["loop-location"] = loop_location
                                    loop_activity["frequency"] = sublogs[subtree_candidate.node_id][i].attributes[
                                        'frequency']
                                    loop_activities[i].append(loop_activity)
                                    loop_location = "redo"
                                    found_in_do = True
                                    success = False
                                    index_k_do = index_k_do + 1

                                if not found_in_do:
                                    loop_location = "redo"
                                    success = False
                                    continue

                            # check if activity is in the right subtree
                            elif loop_location == "redo":
                                found_in_redo = False
                                is_tau = False

                                if subtree_candidate.reference.children[1].id in sublogs:
                                    for k in range(
                                        index_k_redo,
                                        len(
                                            sublogs[
                                                subtree_candidate.reference.children[1].id
                                            ]
                                        ),
                                    ):
                                        match = False

                                        if (
                                            len(
                                                sublogs[
                                                    subtree_candidate.reference.children[
                                                        1
                                                    ].id
                                                ][k]
                                            )
                                            > 0
                                        ):
                                            for l in range(
                                                index_l_redo,
                                                len(
                                                    sublogs[
                                                        subtree_candidate.reference.children[
                                                            1
                                                        ].id
                                                    ][k]
                                                ),
                                            ):
                                                if (
                                                    sublogs[
                                                        subtree_candidate.reference.children[
                                                            1
                                                        ].id
                                                    ][k][l]["concept:name"]
                                                    == activity["concept:name"]
                                                ):
                                                    match = True
                                                    index_l_redo = l + 1
                                                    break
                                        else:
                                            is_tau = True
                                            break

                                        if match:
                                            found_in_redo = True
                                            loop_activity["loop-location"] = loop_location
                                            index_k_redo = k

                                            if index_l_redo == len(
                                                sublogs[
                                                    subtree_candidate.reference.children[
                                                        1
                                                    ].id
                                                ][k]
                                            ):
                                                index_k_redo = k + 1
                                                index_l_redo = 0
                                                loop_location = "do"
                                                if index_k_redo == len(
                                                    sublogs[
                                                        subtree_candidate.reference.children[
                                                            1
                                                        ].id
                                                    ]
                                                ):
                                                    index_k_redo = 0
                                            break
                                else:
                                    is_tau = True

                                if is_tau:
                                    loop_activity["name"] = "tau"
                                    loop_activity["loop-location"] = loop_location
                                    loop_activity["frequency"] = sublogs[subtree_candidate.node_id][i].attributes[
                                        'frequency']
                                    loop_activities[i].append(loop_activity)
                                    loop_location = "do"
                                    found_in_redo = True
                                    success = False
                                    index_k_redo = index_k_redo + 1

                                if not found_in_redo:
                                    loop_location = "do"
                                    success = False
                                    continue

                        loop_activities[i].append(loop_activity)
                        if j == len(sublogs[subtree_candidate.node_id][i]) - 1:
                            # cater for last missed do with tau, sublog of the loop(parent) doesn't contain evidence of the
                            # last tau and the for loop(j) terminates before the case can be handled
                            if loop_activities[i][-1]["loop-location"] == "redo":
                                loop_activities[i].append(
                                    {"name": "tau", "loop-location": "do",
                                     "frequency": sublogs[subtree_candidate.node_id][i].attributes['frequency']}
                                )
                                index_k_do += 1

                else:
                    loop_activities[i].append({"name": "tau", "loop-location": "do",
                                               "frequency": sublogs[subtree_candidate.node_id][i].attributes[
                                                   'frequency']})
                    index_k_do += 1

        except Exception as e:
            print("Error in annotate_loop_activities_with_do_or_redo: " + str(e))

        return loop_activities

    def get_do_and_redo_frequencies_negative_variant(self, subtree_candidate) -> (int, int):

        loop_labels_annotated_remove = (
            self.annotate_loop_activities_with_do_or_redo(
                subtree_candidate, self.remove_sublogs
            )
        )

        (
            do_frequency_remove,
            redo_frequency_remove,
            x,
            do_frequency_trace_remove,
            redo_frequency_trace_remove
        ) = self.calculate_loop_statistics_from_annotated_labels(
            loop_labels_annotated_remove
        )

        return do_frequency_remove, redo_frequency_remove

    def calculate_trace_frequencies_based_heuristics_for_parallel_operator(
        self, subtree_candidate: CandidateSubtree, drop_least_frequent_patterns: bool = True
    ):
        # generation of candidate activity sequences for the candidate subtree

        # label combinations in tree used by the variant to remove
        remove_combinations = self.generate_candidate_activity_sequences_for_parallel(
            subtree_candidate, self.remove_sublogs
        )
        remove_combinations_stats = []
        if len(remove_combinations) > 0:
            remove_combinations_stats = copy.deepcopy(remove_combinations[0])

        for i in range(len(remove_combinations_stats)):
            remove_combinations_stats[i] = {
                "combination": remove_combinations_stats[i],
                # "pos-freq": 0,
                # "neg-freq": 1,
                "pos-traces-with-only-this-pattern": 0,
                "presence-measure": 0,
                "alternate-combinations": [],
            }

        # label combinations in tree used by the variants to keep (each index refers to a patterns in the corresponding variant)
        keep_combinations = self.generate_candidate_activity_sequences_for_parallel(
            subtree_candidate,
            self.keep_sublogs,
            remove_combinations,
            remove_combinations_stats,
        )

        top_remove_combinations_stats = []
        if len(keep_combinations) > 0:
            # count frequency of label patterns seen in the variants to keep
            frequent_combinations = get_frequent_combinations(
                keep_combinations, subtree_candidate
            )

            # calculate statisitcs for the patterns seen in the negative variant to potentially remove
            for i in range(len(remove_combinations_stats)):
                # calculate if there are any alternate patterns in the positive variants with the same activity labels as the current pattern in the negative variant
                remove_combinations_stats[i][
                    "alternate-combinations"
                ] = get_alternate_combinations(
                    remove_combinations_stats[i]["combination"],
                    frequent_combinations,
                    subtree_candidate,
                )

                frequency_best_alternate_combination = 0.1
                if len(remove_combinations_stats[i]["alternate-combinations"]) > 0:
                    frequency_best_alternate_combination = remove_combinations_stats[i][
                        "alternate-combinations"
                    ][0]["frequency"]

                remove_combinations_stats[i]["presence-measure"] = 1 - (
                    frequency_best_alternate_combination / self.positive_traces_frequency
                )

            remove_combinations_stats = sorted(
                remove_combinations_stats,
                key=lambda x: x["presence-measure"],
                reverse=False,
            )

            if drop_least_frequent_patterns:
                # drop 50% of the patterns having least presence-measure
                top_remove_combinations_stats = remove_combinations_stats[
                                                0: math.ceil(len(remove_combinations_stats) / 5)
                                                ]
            else:
                top_remove_combinations_stats = remove_combinations_stats

            # drop patterns that are sub-patterns of the others
            indexes_to_remove = set()
            for i in range(len(top_remove_combinations_stats)):
                for j in range(i, len(top_remove_combinations_stats)):
                    if i != j and is_combination_subsetOf(
                        top_remove_combinations_stats[i]["combination"],
                        top_remove_combinations_stats[j]["combination"],
                        subtree_candidate,
                    ):
                        # remove j
                        indexes_to_remove.add(j)
            for i in sorted(list(indexes_to_remove), reverse=True):
                del top_remove_combinations_stats[i]

        return top_remove_combinations_stats

    def generate_candidate_activity_sequences_for_parallel(
        self,
        subtree_candidate: CandidateSubtree,
        sublogs,
        remove_combinations=None,
        remove_combinations_stats=None,
    ):
        combinations = []
        combination_map = {}
        remove_trace_indexes = {}
        for i in range(len(sublogs[subtree_candidate.node_id])):
            activity_sequences_per_subtree = {}
            sequence = []
            for activity in sublogs[subtree_candidate.node_id][i]:
                matched = False
                for key in sublogs:
                    if key != subtree_candidate.node_id and is_subtree(
                        subtree_candidate.reference,
                        find_tree_node_by_id(key, self.process_tree),
                    ):
                        for j in range(len(sublogs[key])):
                            for matching_activity in sublogs[key][j]:
                                if (
                                    matching_activity["concept:name"]
                                    == activity["concept:name"]
                                ):
                                    matched = True
                                    sequence.append(
                                        {
                                            "name": activity["concept:name"],
                                            "parent-id": key,
                                            "trace-index": sublogs[subtree_candidate.node_id][i].attributes['index'],
                                            "frequency": sublogs[subtree_candidate.node_id][i].attributes['frequency']
                                        }
                                    )
                                    if activity_sequences_per_subtree.get(key) is None:
                                        activity_sequences_per_subtree[key] = [
                                            activity["concept:name"]
                                        ]
                                    else:
                                        activity_sequences_per_subtree[key].append(
                                            activity["concept:name"]
                                        )
                                    break
                            if matched:
                                break
                    else:
                        continue
                    if matched:
                        break

                if not matched:
                    sequence.append(
                        {
                            "name": activity["concept:name"],
                            "parent-id": subtree_candidate.node_id,
                            "trace-index": sublogs[subtree_candidate.node_id][i].attributes['index'],
                            "frequency": sublogs[subtree_candidate.node_id][i].attributes['frequency']
                        }
                    )
                    if (
                        activity_sequences_per_subtree.get(subtree_candidate.node_id)
                        is None
                    ):
                        activity_sequences_per_subtree[subtree_candidate.node_id] = [
                            activity["concept:name"]
                        ]
                    else:
                        activity_sequences_per_subtree[
                            subtree_candidate.node_id
                        ].append(activity["concept:name"])

            index = 0
            parent_id = sequence[0]["parent-id"]
            parallel_activities_in_sequence = []
            for activity in sequence:
                for activity_in_subtree_sequence in activity_sequences_per_subtree[
                    activity["parent-id"]
                ]:
                    if activity_in_subtree_sequence == activity["name"]:
                        if (
                            parent_id != activity["parent-id"]
                            or activity["parent-id"] == subtree_candidate.node_id
                        ):
                            parent_id = activity["parent-id"]
                            index += 1
                        if activity["parent-id"] == subtree_candidate.node_id:
                            parent_id = activity["parent-id"]
                            parallel_activities_in_sequence.append(
                                {"names": [activity["name"]], "parent-id": parent_id,
                                 "trace-index": activity["trace-index"], "frequency": activity["frequency"]}
                            )
                        else:
                            if len(parallel_activities_in_sequence) - 1 < index:
                                parallel_activities_in_sequence.append(
                                    {
                                        "names": [activity["name"]],
                                        "parent-id": parent_id,
                                        "trace-index": activity["trace-index"],
                                        "frequency": activity["frequency"]
                                    }
                                )
                            else:
                                parallel_activities_in_sequence[index]["names"].append(
                                    activity["name"]
                                )
                        break

            combination = generate_combinations_from_sequence(
                parallel_activities_in_sequence,
                subtree_candidate,
                remove_combinations,
                remove_combinations_stats,
            )

            # remove_trace_indexes to hold indexes of traces which have both a viable alternate combination as well as the negative combination

            if len(combination) > 0:
                if remove_combinations is None:
                    if combination_map.get(str(combination)) is None:
                        combinations.append(combination)
                        combination_map[str(combination)] = 1
                else:
                    # dont use the combination from a variant that also includes the remove combination in case parallel is under a loop..
                    is_comb_in_remove_combs = False
                    for rem_comb in remove_combinations:
                        if is_combination_in_remove_combinations(
                            rem_comb,
                            combination,
                            subtree_candidate,
                        ):
                            is_comb_in_remove_combs = True
                            break

                    if is_comb_in_remove_combs and remove_trace_indexes.get(combination[0][0]['trace-index']) is None:
                        remove_trace_indexes[combination[0][0]['trace-index']] = 1
                    else:
                        if combination_map.get(str(combination)) is None:
                            combinations.append(combination)
                            combination_map[str(combination)] = 1

        # prune the combinations having trace_index in remove_trace_indexes
        if remove_combinations is not None:
            for index in remove_trace_indexes:
                combinations = [k for k in combinations if k[0][0]['trace-index'] != index]

        return combinations


def generate_combinations_from_sequence(
    parallel_activities_in_sequence,
    subtree_candidate: CandidateSubtree,
    remove_combinations=None,
    remove_combinations_stats=None,
):
    combinations = []
    for i in range(len(parallel_activities_in_sequence)):
        prefix = [parallel_activities_in_sequence[i]]
        if (
            len(prefix[0]["names"]) > 1
            and subtree_candidate.node_id == prefix[0]["parent-id"]
        ):
            if remove_combinations is not None:
                if (len(remove_combinations) > 0):
                    combinations.append(copy.deepcopy(prefix))
            else:
                combinations.append(copy.deepcopy(prefix))

        for j in range(i + 1, len(parallel_activities_in_sequence)):
            prefix.append(parallel_activities_in_sequence[j])
            if remove_combinations is not None:
                if (len(remove_combinations) > 0):
                    combinations.append(copy.deepcopy(prefix))
            else:
                combinations.append(copy.deepcopy(prefix))

    return combinations


def is_combination_in_remove_combinations(
    remove_combinations,
    combination,
    subtree_candidate: CandidateSubtree
):
    if len(combination) == len(remove_combinations):
        for i in range(len(remove_combinations)):
            all_combination_nodes_match = False

            all_combination_nodes_match = are_combinations_equal(
                remove_combinations[i], combination[i], subtree_candidate
            )

            if all_combination_nodes_match:
                # remove_combinations_stats[i]["pos-freq"] += 1
                # remove_combinations_stats[i]["pos-traces-with-pattern"] += 1
                return True

    return False


def are_combinations_equal(
    combination1, combination2, subtree_candidate: CandidateSubtree
):
    all_combination_nodes_match = False
    if len(combination1) == len(combination2):
        for j in range(len(combination1)):
            all_names_match = False
            if (
                len(combination1[j]["names"]) == len(combination2[j]["names"])
                and combination1[j]["parent-id"] == combination2[j]["parent-id"]
                and subtree_candidate.node_id == combination2[j]["parent-id"]
            ):
                for k in range(len(combination1[j]["names"])):
                    if combination2[j]["names"][k] == combination1[j]["names"][k]:
                        all_names_match = True
                    else:
                        all_names_match = False
                        break
            elif (
                combination1[j]["parent-id"] == combination2[j]["parent-id"]
                and subtree_candidate.node_id != combination2[j]["parent-id"]
            ):
                all_names_match = True

            if all_names_match:
                all_combination_nodes_match = True
            else:
                all_combination_nodes_match = False
                break

    if all_combination_nodes_match:
        return True

    return False


def is_combination_subsetOf(combination1, combination2, subtree_candidate):
    is_subset = False
    for i in range((len(combination2) - len(combination1)) + 1):
        j = 0
        for j in range(len(combination1)):
            if (
                len(combination1[j]["names"]) == len(combination2[i + j]["names"])
                and combination1[j]["parent-id"] == combination2[i + j]["parent-id"]
                and subtree_candidate.node_id == combination2[i + j]["parent-id"]
            ):
                all_names_match = True
                for k in range(len(combination1[j]["names"])):
                    if combination2[i + j]["names"][k] != combination1[j]["names"][k]:
                        all_names_match = False
                        break
                if not all_names_match:
                    break

            elif (
                combination1[j]["parent-id"] != combination2[i + j]["parent-id"]
                and subtree_candidate.node_id != combination2[i + j]["parent-id"]
            ):
                break

        if j == len(combination1) - 1:
            return True

    return False


def combination_without_trace_index(combination):
    temp = copy.deepcopy(combination)
    combination_without_trace_index = []
    for i in range(len(temp)):
        combination_without_trace_index.append({'names': temp[i]['names'], 'parent-id': temp[i]['parent-id']})
    return combination_without_trace_index


def get_frequent_combinations(combinations, subtree_candidate):
    frequent_combinations = []
    all_combinations = {}
    all_combinations_with_trace_index = {}
    for i in range(len(combinations)):
        for j in range(len(combinations[i])):
            key = str(combination_without_trace_index(combinations[i][j]))
            if all_combinations.get(key) is None:
                all_combinations[key] = {
                    "combination": combinations[i][j],
                    "frequency": combinations[i][j][0]['frequency'],
                    "pos-traces-with-only-this-pattern": 1
                }
            else:
                all_combinations[key]["frequency"] += combinations[i][j][0]['frequency']

            if all_combinations_with_trace_index.get(str(combinations[i][j])) is None:
                all_combinations[key]["pos-traces-with-only-this-pattern"] += 1
                all_combinations_with_trace_index[str(combinations[i][j])] = 1

    for key in all_combinations:
        frequent_combinations.append(all_combinations[key])

    return sorted(frequent_combinations, key=lambda x: x["frequency"], reverse=True)


def get_alternate_combinations(
    combination_in_neg_trace: any,
    frequent_combinations_in_pos_traces: [any],
    subtree_candidate: CandidateSubtree,
):
    alternate_combinations = []

    for i in range(len(frequent_combinations_in_pos_traces)):
        frequent_combination_pos = frequent_combinations_in_pos_traces[i]["combination"]
        alternate_found = True

        if len(frequent_combination_pos) == len(combination_in_neg_trace) and not are_combinations_equal(
            frequent_combination_pos, combination_in_neg_trace, subtree_candidate):
            negCounter = collections.Counter()
            posCounter = collections.Counter()

            for j in range(len(combination_in_neg_trace)):
                single_name_match = False

                for k in range(len(frequent_combination_pos)):
                    for l in range(len(frequent_combination_pos[k]["names"])):
                        posCounter[frequent_combination_pos[k]["names"][l] + '-' + str(
                            frequent_combination_pos[k]["parent-id"])] += 1
                    for m in range(len(combination_in_neg_trace[j]["names"])):
                        negCounter[combination_in_neg_trace[j]["names"][m] + '-' + str(
                            combination_in_neg_trace[j]["parent-id"])] += 1

            if len(negCounter - posCounter) == 0:
                alternate_combinations.append(frequent_combinations_in_pos_traces[i])

    return sorted(alternate_combinations, key=lambda x: x["frequency"], reverse=True)
