import copy
import dataclasses
from typing import Callable
from collections import Counter
from pm4py import ProcessTree
from pm4py.objects.log.obj import Trace, EventLog, Event
from pm4py.objects.petri_net.utils.align_utils import STD_MODEL_LOG_MOVE_COST
from pm4py.util.variants_util import variant_to_trace

from cortado_core.negative_process_model_repair.frequency_trace import FrequencyTypedTrace
from cortado_core.utils.alignment_utils import (
    calculate_alignment_typed_trace,
    alignment_contains_deviation,
)
from cortado_core.utils.cvariants import get_concurrency_variants, get_detailed_variants
from cortado_core.models.infix_type import InfixType
from cortado_core.utils.sequentializations import generate_sequentializations
from cortado_core.utils.split_graph import Group
from cortado_core.utils.timestamp_utils import TimeUnit
from cortado_core.utils.trace import TypedTrace
from cortado_core.utils.process_tree import LabelWithIndex
from pm4py.algo.conformance.alignments.process_tree.variants import (
    search_graph_pt as tree_alignment,
)
import zss


@dataclasses.dataclass
class VariantInformation:
    infix_type: InfixType
    is_user_defined: bool


def get_c_variants(
    event_log: EventLog,
    use_mp: bool = False,
    time_granularity: TimeUnit = min(TimeUnit),
):
    variants: dict[Group, list[Trace]] = get_concurrency_variants(
        event_log, use_mp, time_granularity
    )

    total_traces: int = len(event_log)
    info_generator: Callable[
        [list[Trace]], VariantInformation
    ] = lambda _: VariantInformation(
        infix_type=InfixType.NOT_AN_INFIX, is_user_defined=False
    )

    return variants_to_variant_objects(
        variants, time_granularity, total_traces, info_generator
    )


def variants_to_variant_objects(
    variants: dict[Group, list[Trace]],
    time_granularity: TimeUnit,
    total_traces: int,
    info_generator: Callable[[list[Trace]], VariantInformation],
):
    res_variants = []

    cache_variants = dict()

    for bid, (v, ts) in enumerate(
        sorted(list(variants.items()), key=lambda e: len(e[1]), reverse=True)
    ):
        info: VariantInformation = info_generator(ts)
        v.infix_type = info.infix_type
        variant, sub_vars = create_variant_object(
            time_granularity, total_traces, bid, v, ts, info
        )

        res_variants.append(variant)
        cache_variants[bid] = (v, ts, sub_vars, info)

    return (
        sorted(res_variants, key=lambda variant: variant["count"], reverse=True),
        cache_variants,
    )


def create_variant_object(
    time_granularity: TimeUnit,
    total_traces: int,
    bid: int,
    v: Group,
    ts: list[Trace],
    info: VariantInformation,
):
    sub_variants = create_subvariants(ts, time_granularity)

    # Default value of clusterId in a variant = -1
    variant = {
        "count": len(ts),
        "variant": v.serialize(),
        "bid": bid,
        "length": len(v),
        "number_of_activities": v.number_of_activities(),
        "percentage": round(len(ts) / total_traces * 100, 2),
        "nSubVariants": len(sub_variants.keys()),
        "userDefined": info.is_user_defined,
        "infixType": info.infix_type.value,
        "clusterId": -1,
    }

    # If the variant is only a single activity leaf, wrap it up as a sequence
    if "leaf" in variant["variant"].keys() or "parallel" in variant["variant"].keys():
        variant["variant"] = {"follows": [variant["variant"]]}

    return variant, sub_variants


def create_subvariants(ts: list[Trace], time_granularity: TimeUnit):
    sub_vars = get_detailed_variants(ts, time_granularity=time_granularity)

    return sub_vars


def get_traces_from_variants(variants):
    is_n_sequentialization_reduction_enabled = True
    number_of_sequentializations_per_variant = 1
    total_frequency = 0

    n_sequentializations = (
        -1
        if not is_n_sequentialization_reduction_enabled
        else number_of_sequentializations_per_variant
    )
    traces = []

    for cvariant, infix_type, frequency in variants:
        if frequency == 0:
            frequency = 1
        sequentializations = generate_sequentializations(
            Group.deserialize(cvariant), n_sequentializations=n_sequentializations
        )
        traces += [
            FrequencyTypedTrace(variant_to_trace(seq), InfixType(infix_type), frequency)
            for seq in sequentializations
        ]
        total_frequency += frequency

    return traces, total_frequency


def get_traces_from_variant(variant, is_negative_trace: bool = False):
    is_n_sequentialization_reduction_enabled = True
    number_of_sequentializations_per_variant = 1

    n_sequentializations = (
        -1
        if not is_n_sequentialization_reduction_enabled
        else number_of_sequentializations_per_variant
    )

    traces = []
    cvariant = {"follows": variant["follows"]}

    sequentializations = generate_sequentializations(
        Group.deserialize(cvariant), n_sequentializations=n_sequentializations
    )
    traces += [
        TypedTrace(variant_to_trace(seq), InfixType.NOT_AN_INFIX)
        for seq in sequentializations
    ]

    # if is_negative_trace:
    #     traces = {
    #         "first_sequentialization_trace": TypedTrace(
    #             variant_to_trace(sequentializations[0]), InfixType.NOT_AN_INFIX
    #         ),
    #         "c_variant": cvariant,
    #     }
    #
    # else:
    #     traces += [
    #         TypedTrace(variant_to_trace(sequentializations[0]), InfixType.NOT_AN_INFIX)
    #     ]

    return traces


def calculate_percentage_traces_conforming(
    traces: [TypedTrace], process_tree: ProcessTree
):
    deviations = 0

    try:
        for trace in traces:
            alignment = calculate_alignment_typed_trace(process_tree, trace)
            deviation = alignment_contains_deviation(alignment)
            if deviation:
                deviations += 1

        percentage_traces_conforming = ((len(traces) - deviations) / len(traces)) * 100

    except Exception as e:
        percentage_traces_conforming = 100
        print("Error calculating percentage traces conforming: " + str(e))

    return percentage_traces_conforming


# def calculate_variant_conformance(c_variant, process_tree: ProcessTree):
#     total_cost = 0
#     deviations = 0
#     all_variants_length = 1
#     c_variant_indexed = index_leafs(c_variant)
#
#     try:
#         n_sequentializations = (
#             -1
#             # if not config.is_n_sequentialization_reduction_enabled
#             # else config.number_of_sequentializations_per_variant
#         )
#         all_variants = generate_sequentializations(
#             Group.deserialize(c_variant_indexed),
#             n_sequentializations=n_sequentializations,
#         )
#         if len(all_variants) > 0:
#             all_variants_length = len(all_variants)
#
#         for variant in all_variants:
#             alignment = calculate_alignment(
#                 variant, process_tree, InfixType.NOT_AN_INFIX
#             )
#             total_cost += alignment["cost"]
#             deviations += alignment["deviation"]
#
#     except Exception as e:
#         total_cost = 0
#         deviations = 0
#         all_variants_length = 1
#         print("Error calculating variant conformance: " + str(e))
#
#     return {
#         "cost": total_cost / all_variants_length,
#         "deviations": deviations / all_variants_length,
#     }


def calculate_alignment(variant, pt: ProcessTree, infix_type: InfixType):
    # this function uses the specific tree alignment calculation
    trace = Trace()
    for a in variant:
        e = Event()
        e["concept:name"] = a
        trace.append(e)

    align = tree_alignment.apply_from_variants_list([tuple(variant)], pt)
    align = align[0]
    align["deviation"] = align["cost"] > 0

    if infix_type != InfixType.NOT_AN_INFIX:
        align["deviation"] = align["cost"] >= STD_MODEL_LOG_MOVE_COST

    # remove non essential information
    res = {k: align[k] for k in ["alignment", "cost", "deviation"]}

    return res


def index_leafs(variant, indices=None):
    if indices is None:
        indices = Counter()
    if "follows" in variant:
        res = {"follows": []}
        for v in variant["follows"]:
            childs = index_leafs(v, indices)
            res["follows"].append(childs)
        return res
    elif "parallel" in variant:
        res = {"parallel": []}
        for v in variant["parallel"]:
            childs = index_leafs(v, indices)
            res["parallel"].append(childs)
        return res
    else:
        leafs = []
        for activity in variant["leaf"]:
            leafs.append(LabelWithIndex(activity, indices[activity]))
            indices[activity] += 1
        return {"leaf": leafs}


def is_tau(tree_node):
    if (
        (tree_node.children == [] or tree_node.children is None)
        and tree_node.operator is None
        and tree_node.label is None
    ):
        return True

    return False


def calculate_process_tree_edit_distance(tree1: ProcessTree, tree2: ProcessTree) -> int:
    return zss.simple_distance(tree1, tree2, get_children, get_label, label_dist)


def get_children(n):
    return n.children


def get_label(n):
    if n.operator is not None:
        return n.operator.value

    return n.label


def label_dist(a, b):
    if a == b:
        return 0
    return 1


def remove_subtree_from_parents_children(subtree: ProcessTree):
    if subtree.id > 0:
        index = -1
        for i in range(len(subtree.parent.children)):
            child = subtree.parent.children[i]
            if child.id == subtree.id:
                index = i
                break
        if index != -1:
            del subtree.parent.children[index]

    return subtree



def find_tree_node_by_id(tree_node_id: int, pt: ProcessTree):
    if pt.id == tree_node_id:
        return pt
    elif len(pt.children) > 1:
        i = 0
        for child in pt.children:
            if tree_node_id > child.id:
                i += 1
                if len(pt.children) == i:
                    return find_tree_node_by_id(tree_node_id, child)
                else:
                    continue
            elif tree_node_id < child.id:
                return find_tree_node_by_id(tree_node_id, pt.children[i - 1])
            elif child.id == tree_node_id:
                return child
            else:
                return None


def is_activity_in_children(activity_name: str, pt: ProcessTree):
    for child in pt.children:
        if (
            (child.children == [] or child.children is None)
            and child.operator is None
            and child.label is None
            and activity_name == "tau"
        ):
            return True
        elif (
            (child.children == [] or child.children is None)
            and child.operator is None
            and child.label is not None
            and child.label == activity_name
        ):
            return True
    return False


def get_activity_reference_in_children(activity_name: str, pt: ProcessTree):
    for child in pt.children:
        if (
            (child.children == [] or child.children is None)
            and child.operator is None
            and child.label is None
            and activity_name == "tau"
        ):
            return child
        elif (
            (child.children == [] or child.children is None)
            and child.operator is None
            and child.label is not None
            and child.label == activity_name
        ):
            return child
    return None


def remove_and_return_child_activity(
    parent_node: ProcessTree, label: str
) -> ProcessTree:
    index = -1
    child: ProcessTree = None
    for i in range(len(parent_node.children)):
        if parent_node.children[i].label == label:
            index = i
            break

    if index != -1:
        child = copy.deepcopy(parent_node.children[index])
        del parent_node.children[index]

    return child


def remove_and_return_child_subtree(
    parent_node: ProcessTree, parent_id_of_removal_node: int
) -> ProcessTree:
    index = -1
    child: ProcessTree = None

    for i in range(len(parent_node.children)):
        if parent_node.children[i].operator is not None:
            if (
                parent_node.children[i].id == parent_id_of_removal_node
                or find_tree_node_by_id(parent_id_of_removal_node, parent_node)
                is not None
            ):
                index = i
                break

    if index != -1:
        child = copy.deepcopy(parent_node.children[index])
        del parent_node.children[index]

    return child

def remove_and_return_child_by_id(
    parent_node: ProcessTree, id: int
) -> ProcessTree:
    index = -1
    child: ProcessTree = None
    for i in range(len(parent_node.children)):
        if parent_node.children[i].id == id:
            index = i
            break

    if index != -1:
        child = copy.deepcopy(parent_node.children[index])
        del parent_node.children[index]

    return child
