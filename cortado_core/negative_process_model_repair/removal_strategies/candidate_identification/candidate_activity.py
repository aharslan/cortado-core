from dataclasses import dataclass
from pm4py import ProcessTree
from pm4py.objects.process_tree.obj import Operator


@dataclass
class CandidateActivity:
    parent_id: int
    operator: Operator
    activity_name: str
    frequency_in_positive_traces: int
    frequency_in_negative_trace: int
    number_of_positive_traces_with_activity: int
    presence_measure: float
    activity_reference: ProcessTree

    def __init__(
        self,
        parent_id,
        operator,
        activity_name,
        frequency_in_positive_traces,
        frequency_in_negative_trace,
        number_of_positive_traces_with_activity,
        activity_reference,
        number_of_positive_traces,
    ):
        self.parent_id: int = parent_id
        self.operator: Operator = operator
        self.activity_name: str = activity_name
        self.frequency_in_positive_traces: int = frequency_in_positive_traces

        self.frequency_in_negative_trace = frequency_in_negative_trace

        self.number_of_positive_traces_with_activity = (
            number_of_positive_traces_with_activity
        )
        self.activity_reference: ProcessTree = activity_reference

        self.calculate_presence_measure(number_of_positive_traces)

    def calculate_presence_measure(self, number_of_positive_traces: int):
        frequency_in_negative_trace = 1

        if number_of_positive_traces == 0:
            number_of_positive_traces = 1
        if self.frequency_in_negative_trace != 0:
            frequency_in_negative_trace = self.frequency_in_negative_trace

        self.presence_measure = (
            (self.number_of_positive_traces_with_activity / number_of_positive_traces)
        )

