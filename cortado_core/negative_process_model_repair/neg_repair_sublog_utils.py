import multiprocessing
from pm4py import ProcessTree
from pm4py.objects.log.obj import Trace, EventLog, Event
from cortado_core.negative_process_model_repair.frequency_trace import FrequencyTypedTrace
from cortado_core.process_tree_utils.miscellaneous import is_leaf_node, is_subtree
from cortado_core.process_tree_utils.to_petri_net_transition_bordered import (
    apply as pt_to_petri_net,
)
from cortado_core.utils.alignment_utils import (
    alignment_contains_deviation,
)
from cortado_core.utils.parallel_alignments import calculate_alignments_parallel, \
    calculate_alignments_parallel_typed_trace, calculate_alignments_typed_trace


def calculate_sublog_for_negative_process_model_repair_with_frequency(
    pt: ProcessTree,
    fitting_traces: list[FrequencyTypedTrace],
    pool,
):
    # not_infix_log, infix_traces = __split_log_by_infix_type_with_frequency(fitting_traces)
    sublogs = __calculate_sub_log_for_each_node_regular_traces_with_frequency(
        pt, fitting_traces, pool=pool
    )
    return sublogs


def add_alignment_to_sublogs_with_frequency(alignment, sublogs, index, allow_deviations=False):
    if not allow_deviations:
        assert not alignment_contains_deviation(alignment)
    currently_active_pt_nodes = {}

    for step in alignment["alignment"]:
        # executed transition always corresponds to a node in the process tree
        current_pt = step[0][1][0]
        if (current_pt, current_pt.id) in currently_active_pt_nodes:
            if current_pt.id not in sublogs:
                sublogs[current_pt.id] = EventLog()

            sublogs[current_pt.id].append(
                currently_active_pt_nodes[(current_pt, current_pt.id)]
            )

            if len(currently_active_pt_nodes[(current_pt, current_pt.id)]) == 0:
                sublogs[current_pt.id][len(sublogs[current_pt.id])-1].attributes['frequency'] = alignment["frequency"]
                sublogs[current_pt.id][len(sublogs[current_pt.id]) - 1].attributes['index'] = index

            # every pt node occurs at least twice in an alignment, i.e., start and end. Hence when we observe a pt
            # node for the second time, we know it is closed
            assert step[0][1][1] == "closed"
            del currently_active_pt_nodes[(current_pt, current_pt.id)]
        elif not is_leaf_node(current_pt):
            currently_active_pt_nodes[(current_pt, current_pt.id)] = Trace()

        if is_leaf_node(current_pt):
            activity_name = step[1][1]
            if activity_name:
                for active_node, active_node_obj_id in currently_active_pt_nodes:
                    if is_subtree(active_node, current_pt):
                        event = Event()
                        event["concept:name"] = activity_name
                        event["concept:id"] = current_pt.id

                        currently_active_pt_nodes[(active_node, active_node.id)].append(
                            event
                        )
                        currently_active_pt_nodes[(active_node, active_node.id)].attributes['frequency'] = alignment["frequency"]
                        currently_active_pt_nodes[(active_node, active_node.id)].attributes['index'] = index


    return sublogs


def __calculate_sub_log_for_each_node_regular_traces_with_frequency(
    pt: ProcessTree, log: EventLog, pool: multiprocessing.pool.Pool
) -> dict[int, EventLog]:
    """
    Calculates the sublog for each full, already added trace by first computing the alignment and then adding the relevant
    parts to the sublog of the lca.
    Parameters
    ----------
    pt
    log
    pool

    Returns
    -------

    """
    sublogs: dict[int, EventLog] = {}
    # assumption: log is replayable on process tree without deviations
    net, im, fm = pt_to_petri_net(pt)

    if pool is not None:
        alignments = calculate_alignments_parallel_typed_trace(
            log, net, im, fm, parameters={"ret_tuple_as_trans_desc": True}, pool=pool
        )
    else:
        alignments = calculate_alignments_typed_trace(
            log, net, im, fm, parameters={"ret_tuple_as_trans_desc": True}
        )
    index = 0
    for alignment in alignments:
        sublogs = add_alignment_to_sublogs_with_frequency(alignment, sublogs, index)
        index = index + 1

    return sublogs

