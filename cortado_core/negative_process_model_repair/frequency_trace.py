from dataclasses import dataclass
from cortado_core.utils.trace import TypedTrace
from pm4py.objects.log.obj import Trace
from cortado_core.models.infix_type import InfixType


@dataclass
class FrequencyTypedTrace(TypedTrace):
    frequency: Trace

    def __init__(self, trace, infix_type: InfixType, frequency: int):
        self.frequency = frequency
        super().__init__(trace, infix_type)


    def __hash__(self):
        return hash((hash(self.trace), self.infix_type.value, self.frequency))

    def __eq__(self, other):
        return self.trace == other.trace and self.infix_type == other.infix_type and self.frequency == other.frequency
