import copy
from typing import List
import pm4py.objects.process_tree.importer.importer as ptml_importer
import pm4py.objects.log.importer.xes.importer as xes_importer
from cortado_core.negative_process_model_repair.constants import Constants
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.removal_candidates_generator import \
    RemovalCandidatesGenerator
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.removal_candidates_heuristics import \
    RemovalCandidatesHeuristics
from cortado_core.negative_process_model_repair.removal_strategies.rules_based_reduction.complete_brute_force_subtree_update import \
    CompleteBruteForceSubtreeUpdate
from cortado_core.negative_process_model_repair.removal_strategies.rules_based_reduction.heuristic_brute_force_subtree_update import \
    HeuristicBruteForceSubtreeUpdate
from cortado_core.negative_process_model_repair.removal_strategies.rules_based_reduction.rating_based_subtree_update import \
    RatingBasedSubtreeUpdate
from cortado_core.utils.trace import TypedTrace
from pm4py import ProcessTree
from cortado_core.negative_process_model_repair.removal_strategies.fallback_strategy import (
    FallbackStrategy,
)
from cortado_core.negative_process_model_repair.removal_strategies.rules_based_reduction.rule_based_reduction import (
    RuleBasedReduction,
)
from cortado_core.negative_process_model_repair.temp_utils import (
    get_traces_from_variant,
    get_c_variants,
    get_traces_from_variants,
)
from cortado_core.lca_approach import set_preorder_ids_in_tree
from cortado_core.utils.timestamp_utils import TimeUnit
import pm4py.objects.process_tree.exporter.exporter as ptml_exporter
from pm4py.objects.process_tree.utils.generic import parse as pt_parse
from cortado_core.negative_process_model_repair.PoolFactory import PoolFactory


def apply_complete_brute_force_based_negative_process_model_repair(
    pt: ProcessTree,
    negative_variant: List[any],
    positive_variants: List[any],
):
    approach_used = None
    applied_rules = None
    percentage_positive_variants_conforming = -1
    resulting_tree_edit_distance = -1
    update_operations = 0

    if type(negative_variant[0]) is not TypedTrace:
        negative_variant, negative_trace_frequency = get_traces_from_variants(negative_variant)
    if type(positive_variants[0]) is not TypedTrace:
        positive_variants, positive_traces_frequency = get_traces_from_variants(positive_variants)

    set_preorder_ids_in_tree(pt)

    removal_candidates_generator = RemovalCandidatesHeuristics(
        pt, positive_variants, negative_variant[0], negative_trace_frequency, positive_traces_frequency
    )

    removal_candidates_generator.generate_sublogs_for_traces()
    # removal_candidates_generator.generate_sublogs_for_traces(PoolFactory.instance().get_pool())

    removal_candidate_activities = (
        removal_candidates_generator.generate_removal_candidate_activities()
    )

    resulting_tree = copy.deepcopy(pt)

    repair_strategy = CompleteBruteForceSubtreeUpdate(
        removal_candidates_generator, removal_candidate_activities
    )
    (
        resulting_tree_com_brute_force_subtree_update_based,
        tree_updated_com_brute_force_subtree_update_based,
        thresholds_met_heu_brute_force_subtree_update_based,
        percentage_positive_traces_conforming_com_brute_force,
        resulting_tree_edit_distance_com_brute_force,
        applied_rules_com_brute_force
    ) = repair_strategy.apply_complete_brute_force_subtree_update_based_reduction()

    update_operations += repair_strategy.update_operations
    if tree_updated_com_brute_force_subtree_update_based:
        resulting_tree = resulting_tree_com_brute_force_subtree_update_based
        approach_used = 'Complete Brute Force'
        percentage_positive_variants_conforming = percentage_positive_traces_conforming_com_brute_force
        resulting_tree_edit_distance = resulting_tree_edit_distance_com_brute_force
        applied_rules = applied_rules_com_brute_force

    else:
        resulting_tree = pt

    return (tree_updated_com_brute_force_subtree_update_based, resulting_tree, approach_used,
            percentage_positive_variants_conforming, resulting_tree_edit_distance, applied_rules, [], update_operations)


def apply_heuristic_brute_force_based_negative_process_model_repair(
    pt: ProcessTree,
    negative_variant: List[any],
    positive_variants: List[any],
):
    approach_used = None
    applied_rules = None
    percentage_positive_variants_conforming = -1
    resulting_tree_edit_distance = -1
    update_operations = 0

    if type(negative_variant[0]) is not TypedTrace:
        negative_variant, negative_trace_frequency = get_traces_from_variants(negative_variant)
    if type(positive_variants[0]) is not TypedTrace:
        positive_variants, positive_traces_frequency = get_traces_from_variants(positive_variants)

    set_preorder_ids_in_tree(pt)

    removal_candidates_generator = RemovalCandidatesHeuristics(
        pt, positive_variants, negative_variant[0], negative_trace_frequency, positive_traces_frequency
    )

    removal_candidates_generator.generate_sublogs_for_traces()
    # removal_candidates_generator.generate_sublogs_for_traces(PoolFactory.instance().get_pool())

    removal_candidate_activities = (
        removal_candidates_generator.generate_removal_candidate_activities()
    )

    resulting_tree = copy.deepcopy(pt)

    repair_strategy = HeuristicBruteForceSubtreeUpdate(
        removal_candidates_generator, removal_candidate_activities
    )
    (
        resulting_tree_heu_brute_force_subtree_update_based,
        tree_updated_heu_brute_force_subtree_update_based,
        thresholds_met_heu_brute_force_subtree_update_based,
        percentage_positive_traces_conforming_heu_brute_force,
        resulting_tree_edit_distance_heu_brute_force,
        applied_rules_heu_brute_force
    ) = repair_strategy.apply_heuristic_brute_force_subtree_update_based_reduction()

    update_operations += repair_strategy.update_operations
    if tree_updated_heu_brute_force_subtree_update_based:
        resulting_tree = resulting_tree_heu_brute_force_subtree_update_based
        approach_used = 'Heuristic Brute Force'
        percentage_positive_variants_conforming = percentage_positive_traces_conforming_heu_brute_force
        resulting_tree_edit_distance = resulting_tree_edit_distance_heu_brute_force
        applied_rules = applied_rules_heu_brute_force

    else:
        resulting_tree = pt

    return (tree_updated_heu_brute_force_subtree_update_based, resulting_tree, approach_used,
            percentage_positive_variants_conforming, resulting_tree_edit_distance, applied_rules, [], update_operations)


def apply_sub_tree_rating_based_negative_process_model_repair(
    pt: ProcessTree,
    negative_variant: List[any],
    positive_variants: List[any],
):
    approach_used = None
    applied_rules = None
    percentage_positive_variants_conforming = -1
    resulting_tree_edit_distance = -1
    update_operations = 0

    if type(negative_variant[0]) is not TypedTrace:
        negative_variant, negative_trace_frequency = get_traces_from_variants(negative_variant)
    if type(positive_variants[0]) is not TypedTrace:
        positive_variants, positive_traces_frequency = get_traces_from_variants(positive_variants)

    set_preorder_ids_in_tree(pt)

    removal_candidates_generator = RemovalCandidatesHeuristics(
        pt, positive_variants, negative_variant[0], negative_trace_frequency, positive_traces_frequency
    )

    removal_candidates_generator.generate_sublogs_for_traces()
    # removal_candidates_generator.generate_sublogs_for_traces(PoolFactory.instance().get_pool())

    removal_candidate_activities = (
        removal_candidates_generator.generate_removal_candidate_activities()
    )

    resulting_tree = copy.deepcopy(pt)

    repair_strategy = RatingBasedSubtreeUpdate(
        removal_candidates_generator, removal_candidate_activities
    )
    (
        resulting_tree_rating_based_subtree_update,
        tree_updated_rating_based_subtree_update,
        thresholds_met_rating_based_subtree_update,
        percentage_positive_traces_conforming_rating_based,
        resulting_tree_edit_distance_rating_based,
        applied_rules_rating_based,
        failed_updates
    ) = repair_strategy.apply_rating_based_subtree_update_reduction()

    update_operations += repair_strategy.update_operations
    if tree_updated_rating_based_subtree_update:
        resulting_tree = resulting_tree_rating_based_subtree_update
        approach_used = 'Sub-tree Rating Based Update'
        percentage_positive_variants_conforming = percentage_positive_traces_conforming_rating_based
        resulting_tree_edit_distance = resulting_tree_edit_distance_rating_based
        applied_rules = applied_rules_rating_based

    else:
        resulting_tree = pt

    return (tree_updated_rating_based_subtree_update, resulting_tree, approach_used,
            percentage_positive_variants_conforming, resulting_tree_edit_distance, applied_rules, failed_updates,
            update_operations)


def apply_negative_process_model_repair(
    pt: ProcessTree,
    negative_variant: List[any],
    positive_variants: List[any],
):
    """Apply neg proc model repair using heuristic based subtree pruning and case identification
        with the fallback strategy that brute forces all possible update cases on all subtrees
    """
    assert len(positive_variants) > 0
    assert len(negative_variant) > 0

    approach_used = None
    applied_rules = None
    percentage_positive_variants_conforming = -1
    resulting_tree_edit_distance = -1
    failed_updates = []
    tree_updated = False
    update_operations = 0

    if type(negative_variant[0]) is not TypedTrace:
        negative_variant, negative_trace_frequency = get_traces_from_variants(negative_variant)
    if type(positive_variants[0]) is not TypedTrace:
        positive_variants, positive_traces_frequency = get_traces_from_variants(positive_variants)

    set_preorder_ids_in_tree(pt)

    removal_candidates_generator = RemovalCandidatesHeuristics(
        pt, positive_variants, negative_variant[0], negative_trace_frequency, positive_traces_frequency
    )

    removal_candidates_generator.generate_sublogs_for_traces()
    # removal_candidates_generator.generate_sublogs_for_traces(PoolFactory.instance().get_pool())

    removal_candidate_activities = (
        removal_candidates_generator.generate_removal_candidate_activities()
    )

    resulting_tree = copy.deepcopy(pt)

    repair_strategy = RatingBasedSubtreeUpdate(
        removal_candidates_generator, removal_candidate_activities
    )
    (
        resulting_tree_rating_based_subtree_update,
        tree_updated_rating_based_subtree_update,
        thresholds_met_rating_based_subtree_update,
        percentage_positive_traces_conforming_rating_based,
        resulting_tree_edit_distance_rating_based,
        applied_rules_rating_based,
        failed_updates_rating_based
    ) = repair_strategy.apply_rating_based_subtree_update_reduction()

    update_operations += repair_strategy.update_operations
    if tree_updated_rating_based_subtree_update:
        resulting_tree = resulting_tree_rating_based_subtree_update
        tree_updated = tree_updated_rating_based_subtree_update
        approach_used = 'Rating Based Update'
        percentage_positive_variants_conforming = percentage_positive_traces_conforming_rating_based
        resulting_tree_edit_distance = resulting_tree_edit_distance_rating_based
        applied_rules = applied_rules_rating_based
        failed_updates = failed_updates_rating_based

    if thresholds_met_rating_based_subtree_update == False:
        repair_strategy = HeuristicBruteForceSubtreeUpdate(
            removal_candidates_generator, removal_candidate_activities
        )
        (
            resulting_tree_heu_brute_force_subtree_update_based,
            tree_updated_heu_brute_force_subtree_update_based,
            thresholds_met_heu_brute_force_subtree_update_based,
            percentage_positive_traces_conforming_heu_brute_force,
            resulting_tree_edit_distance_heu_brute_force,
            applied_rules_heu_brute_force
        ) = repair_strategy.apply_heuristic_brute_force_subtree_update_based_reduction()

        update_operations += repair_strategy.update_operations
        if tree_updated_heu_brute_force_subtree_update_based:
            resulting_tree = resulting_tree_heu_brute_force_subtree_update_based
            tree_updated = tree_updated_heu_brute_force_subtree_update_based
            approach_used = 'Heuristic Brute Force'
            percentage_positive_variants_conforming = percentage_positive_traces_conforming_heu_brute_force
            resulting_tree_edit_distance = resulting_tree_edit_distance_heu_brute_force
            applied_rules = applied_rules_heu_brute_force
            failed_updates = []

        else:
            resulting_tree = pt

        # if thresholds_met_heu_brute_force_subtree_update_based == False:
        #     repair_strategy = CompleteBruteForceSubtreeUpdate(
        #         removal_candidates_generator, removal_candidate_activities
        #     )
        #     (
        #         resulting_tree_com_brute_force_subtree_update_based,
        #         tree_updated_com_brute_force_subtree_update_based,
        #         thresholds_met_com_brute_force_subtree_update_based,
        #         percentage_positive_traces_conforming_com_brute_force,
        #         resulting_tree_edit_distance_com_brute_force,
        #         applied_rules_com_brute_force,
        #     ) = repair_strategy.apply_complete_brute_force_subtree_update_based_reduction()
        #
        #     if tree_updated_com_brute_force_subtree_update_based:
        #         resulting_tree = resulting_tree_com_brute_force_subtree_update_based
        #         tree_updated = tree_updated_com_brute_force_subtree_update_based
        #         approach_used = 'Complete Brute Force'
        #         percentage_positive_variants_conforming = percentage_positive_traces_conforming_com_brute_force
        #         resulting_tree_edit_distance = resulting_tree_edit_distance_com_brute_force
        #         applied_rules = applied_rules_com_brute_force
        #         failed_updates = []
        #         update_operations += repair_strategy.update_operations
        #
        #     else:
        #         resulting_tree = pt

    return (tree_updated, resulting_tree, approach_used, percentage_positive_variants_conforming,
            resulting_tree_edit_distance, applied_rules, failed_updates, update_operations)


def apply_frequency_rating_based_negative_process_model_repair(
    pt: ProcessTree,
    negative_variant: List[any],
    positive_variants: List[any],
):
    """Apply neg proc model repair using the frequency based rating method for subtree identification
    with the fallback strategy that removes individual activities independently and returns the 'best'
    resulting tree
    """
    approach_used = 'Rule Based Strategy'
    percentage_positive_variants_conforming = 0
    resulting_tree_edit_distance = 0

    if type(negative_variant[0]) is not TypedTrace:
        negative_variant, negative_trace_frequency = get_traces_from_variants(negative_variant)
    if type(positive_variants[0]) is not TypedTrace:
        positive_variants, positive_traces_frequency = get_traces_from_variants(positive_variants)

    set_preorder_ids_in_tree(pt)

    removal_candidates_generator_heuristic = RemovalCandidatesHeuristics(
        pt, positive_variants, negative_variant[0], negative_trace_frequency, positive_traces_frequency
    )

    removal_candidates_generator_heuristic.generate_sublogs_for_traces()
    # removal_candidates_generator_heuristic.generate_sublogs_for_traces(PoolFactory.instance().get_pool())

    removal_candidate_activities = (
        removal_candidates_generator_heuristic.generate_removal_candidate_activities()
    )

    resulting_tree = copy.deepcopy(pt)

    rule_based_strategy = RuleBasedReduction(
        removal_candidates_generator_heuristic, removal_candidate_activities
    )
    (
        resulting_tree_rule_based,
        tree_updated_rule_based,
        percentage_positive_traces_conforming_rule_based,
        resulting_tree_edit_distance_rule_based,
        applied_rules
    ) = rule_based_strategy.apply_rule_based_reduction()

    if ((percentage_positive_traces_conforming_rule_based < Constants.MIN_THRESHOLD_POSITIVE_FITTING_VARIANTS or
         resulting_tree_edit_distance_rule_based > Constants.MAX_THRESHOLD_TREE_EDIT_DISTANCE) or
        not tree_updated_rule_based):
        fallback_strategy = FallbackStrategy(
            removal_candidates_generator_heuristic, removal_candidate_activities
        )
        (
            resulting_tree_fallback,
            tree_updated_fallback,
            percentage_positive_traces_conforming_fallback,
            resulting_tree_edit_distance_fallback
        ) = fallback_strategy.apply_fallback_strategy()

        if ((percentage_positive_traces_conforming_fallback - resulting_tree_edit_distance_fallback) >
            (percentage_positive_traces_conforming_rule_based - resulting_tree_edit_distance_rule_based)
        ):
            resulting_tree = resulting_tree_fallback
            approach_used = 'Fallback Strategy'
            percentage_positive_variants_conforming = percentage_positive_traces_conforming_fallback
            resulting_tree_edit_distance = resulting_tree_edit_distance_fallback
        else:
            resulting_tree = resulting_tree_rule_based
            percentage_positive_variants_conforming = percentage_positive_traces_conforming_rule_based
            resulting_tree_edit_distance = resulting_tree_edit_distance_rule_based
    else:
        resulting_tree = resulting_tree_rule_based
        percentage_positive_variants_conforming = percentage_positive_traces_conforming_rule_based
        resulting_tree_edit_distance = resulting_tree_edit_distance_rule_based

    return resulting_tree, approach_used, percentage_positive_variants_conforming, resulting_tree_edit_distance, applied_rules


if __name__ == "__main__":
    pt = ptml_importer.apply(
        r"C:\Users\ahmad\OneDrive\Desktop\Sem6\Thesis\Implementation\Sublogs Implementation Artifacts\reciept_pt_first_20.ptml"
    )

    pt_string = str(pt)
    pt_test = pt_parse("-> (*(X(->('A','B'), 'C')")
    set_preorder_ids_in_tree(pt_test)

    set_preorder_ids_in_tree(pt)

    variant = xes_importer.Variants.ITERPARSE
    parameters = {variant.value.Parameters.TIMESTAMP_SORT: True}
    event_log = xes_importer.apply(
        r"C:\Users\ahmad\OneDrive\Desktop\MS DS\Cortado HiWi\Tasks Artifacts\Processes Event Logs\receipt\log.xes",
        variant=variant,
        parameters=parameters,
    )

    # variants = get_experiment_variants(event_log)
    res_variants, cache_variants = get_c_variants(event_log, False, min(TimeUnit))

    keep_fitting_traces = []
    negative_trace = get_traces_from_variant(res_variants[19]["variant"], True)[0]

    i = 0
    for v in res_variants:
        traces = get_traces_from_variant(v["variant"])

        keep_fitting_traces.extend(traces)
        i += 1
        if i == 19:
            break

    removal_candidates_generator = RemovalCandidatesGenerator(
        pt, keep_fitting_traces, negative_trace
    )
    removal_candidates_generator.generate_sublogs_for_traces()
    removal_candidate_activities = (
        removal_candidates_generator.generate_removal_candidate_activities()
    )

    # sublogs.process_tree = pt_test

    fallback_strategy = FallbackStrategy(
        removal_candidates_generator, removal_candidate_activities
    )
    # resulting_tree, tree_updated = fallback_strategy.apply_fallback_strategy()

    rule_based_strategy = RuleBasedReduction(
        removal_candidates_generator, removal_candidate_activities
    )
    resulting_tree, tree_updated = rule_based_strategy.apply_rule_based_reduction()

    ptml_exporter.apply(
        resulting_tree,
        r"C:\Users\ahmad\OneDrive\Desktop\Sem6\Thesis\Implementation\Sublogs Implementation Artifacts\reciept_pt_first_20_reduced.ptml",
    )
