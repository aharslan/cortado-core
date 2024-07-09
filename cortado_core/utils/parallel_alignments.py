from typing import List

from pm4py.objects.log.obj import EventLog, Trace
from pm4py.objects.petri_net.obj import PetriNet, Marking
from pm4py.util.typing import AlignmentResult
from pm4py.algo.conformance.alignments.petri_net.algorithm import (
    apply as calculate_alignment,
)
from pm4py.algo.conformance.alignments.petri_net.algorithm import (
    variants as variants_calculate_alignments,
)
from cortado_core.negative_process_model_repair.frequency_trace import FrequencyTypedTrace


def calculate_alignments_typed_trace(
    log: list[FrequencyTypedTrace], net: PetriNet, im: Marking, fm: Marking, parameters
) -> List[AlignmentResult]:
    results = []
    for trace in log:
        result = calculate_alignment_a_star_frequency(trace.trace, net, im, fm, trace.frequency, parameters)
        results.append(result)

    return results


def calculate_alignments_parallel_typed_trace(
    log: list[FrequencyTypedTrace], net: PetriNet, im: Marking, fm: Marking, parameters, pool
) -> List[AlignmentResult]:
    results = []
    for trace in log:
        result = pool.apply_async(
            calculate_alignment_a_star_frequency,
            args=[trace.trace, net, im, fm, trace.frequency],
            kwds={"parameters": parameters},
        )
        results.append(result)

    return [r.get() for r in results]


def calculate_alignment_a_star_frequency(
    trace: Trace, net: PetriNet, im: Marking, fm: Marking, frequency: int, parameters
) -> AlignmentResult:
    """
    Calculates an alignment using the a star search algorithm. This function is necessary as the variants accepted
    by PM4PY are python modules. So, they cannot be pickled and therefore not used directly in a pool.apply_async()-call
    """
    alignment =  calculate_alignment(
        trace,
        net,
        im,
        fm,
        parameters=parameters,
        variant=variants_calculate_alignments.state_equation_a_star,
    )
    alignment['frequency'] = frequency
    return alignment


def calculate_alignments_parallel(
    log: EventLog, net: PetriNet, im: Marking, fm: Marking, parameters, pool
) -> List[AlignmentResult]:
    results = []
    for trace in log:
        result = pool.apply_async(
            calculate_alignment_a_star,
            args=[trace, net, im, fm],
            kwds={"parameters": parameters},
        )
        results.append(result)

    return [r.get() for r in results]


def calculate_alignment_a_star(
    trace: Trace, net: PetriNet, im: Marking, fm: Marking, parameters
) -> AlignmentResult:
    """
    Calculates an alignment using the a star search algorithm. This function is necessary as the variants accepted
    by PM4PY are python modules. So, they cannot be pickled and therefore not used directly in a pool.apply_async()-call
    """
    return calculate_alignment(
        trace,
        net,
        im,
        fm,
        parameters=parameters,
        variant=variants_calculate_alignments.state_equation_a_star,
    )
