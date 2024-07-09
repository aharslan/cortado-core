import copy

from cortado_core.process_tree_utils.miscellaneous import get_root
from pm4py import ProcessTree
from pm4py.objects.process_tree.obj import Operator
from cortado_core.negative_process_model_repair.removal_strategies.candidate_identification.candidate_subtree import \
    CandidateSubtree
from cortado_core.negative_process_model_repair.temp_utils import find_tree_node_by_id, \
    remove_subtree_from_parents_children, is_tau, remove_and_return_child_activity, remove_and_return_child_subtree, \
    remove_and_return_child_by_id


class UpdateRule:
    removal_candidate_subtree: CandidateSubtree

    def __init__(self, removal_candidate_subtree: CandidateSubtree):
        assert isinstance(removal_candidate_subtree, CandidateSubtree)
        assert isinstance(removal_candidate_subtree.node_id, int)
        assert removal_candidate_subtree.node_id >= 0

        self.removal_candidate_subtree = removal_candidate_subtree


class SequenceUpdateRule(UpdateRule):

    def __init__(self, removal_candidate_subtree: CandidateSubtree):
        super().__init__(removal_candidate_subtree)

    def apply_rule(self, tree_to_update: ProcessTree) -> ProcessTree:
        assert isinstance(tree_to_update, ProcessTree)
        assert tree_to_update is not None

        tree_node_to_update_reference = find_tree_node_by_id(
            self.removal_candidate_subtree.node_id, tree_to_update
        )
        remove_subtree_from_parents_children(tree_node_to_update_reference)

        return tree_to_update


class ChoiceUpdateRule(UpdateRule):

    def __init__(self, removal_candidate_subtree: CandidateSubtree):
        super().__init__(removal_candidate_subtree)

    def apply_rule(self, tree_to_update: ProcessTree) -> ProcessTree:
        assert isinstance(tree_to_update, ProcessTree)
        assert tree_to_update is not None
        assert self.removal_candidate_subtree.get_child_candidate_id() is not None

        child_id_to_remove = self.removal_candidate_subtree.get_child_candidate_id()

        tree_node_to_update = find_tree_node_by_id(
            child_id_to_remove, tree_to_update
        )
        remove_subtree_from_parents_children(tree_node_to_update)

        return tree_to_update


class LoopUpdateRule(UpdateRule):

    def __init__(self, removal_candidate_subtree: CandidateSubtree):
        super().__init__(removal_candidate_subtree)

    def apply_remove_redundant_redo(self, tree_to_update: ProcessTree) -> ProcessTree:
        assert isinstance(tree_to_update, ProcessTree)
        assert tree_to_update is not None
        assert self.removal_candidate_subtree.reference is not None
        assert len(self.removal_candidate_subtree.reference.children) > 0

        tree_node_to_update = find_tree_node_by_id(
            self.removal_candidate_subtree.node_id, tree_to_update
        )
        temp_node: ProcessTree = self.removal_candidate_subtree.reference.children[0]

        for i in range(len(tree_node_to_update.parent.children)):
            if tree_node_to_update.parent.children[i].id == tree_node_to_update.id:
                tree_node_to_update.parent.children[i] = temp_node
                temp_node.parent = tree_node_to_update.parent
                break

        return tree_to_update

    def apply_optional_redo_mandatory(self, tree_to_update: ProcessTree) -> ProcessTree:
        assert isinstance(tree_to_update, ProcessTree)
        assert tree_to_update is not None

        tree_node_to_update = find_tree_node_by_id(
            self.removal_candidate_subtree.node_id, tree_to_update
        )

        assert is_tau(tree_node_to_update.children[0])

        if is_tau(tree_node_to_update.children[0]):
            left_child = copy.deepcopy(tree_node_to_update.children[0])

            tree_node_to_update.children[0] = copy.deepcopy(tree_node_to_update.children[1])
            tree_node_to_update.children[0].parent = tree_node_to_update

            tree_node_to_update.children[1] = left_child
            tree_node_to_update.children[1].parent = tree_node_to_update

        return tree_to_update

    def apply_repeat_exactly_n(self, tree_to_update: ProcessTree, repetitions_to_encode: list[int]) -> ProcessTree:
        assert isinstance(tree_to_update, ProcessTree)
        assert tree_to_update is not None
        assert repetitions_to_encode is not None
        assert len(repetitions_to_encode) > 0

        tree_node_to_update = find_tree_node_by_id(
            self.removal_candidate_subtree.node_id, tree_to_update
        )

        choice_node = ProcessTree()
        choice_node.operator = Operator.XOR
        choice_node.children = []

        index_in_parents_children = 0
        for i in range(len(tree_node_to_update.parent.children)):
            if tree_node_to_update.parent.children[i].id == tree_node_to_update.id:
                index_in_parents_children = i
                break

        for key in repetitions_to_encode:
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

    def apply_repeat_at_least_n(self, tree_to_update: ProcessTree, do_frequency_neg_variant: int) -> ProcessTree:
        assert isinstance(tree_to_update, ProcessTree)
        assert tree_to_update is not None
        assert do_frequency_neg_variant is not None
        assert do_frequency_neg_variant > 0

        tree_node_to_update = find_tree_node_by_id(
            self.removal_candidate_subtree.node_id, tree_to_update
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

        for i in range(do_frequency_neg_variant):
            temp_node.children.append(copy.deepcopy(left_child_loop))
            temp_node.children.append(copy.deepcopy(right_child_loop))

        loop_node = copy.deepcopy(tree_node_to_update)
        loop_node.parent = temp_node
        temp_node.children.append(loop_node)

        temp_node.parent = tree_node_to_update.parent
        tree_node_to_update.parent.children[index_in_parents_children] = temp_node

        return tree_to_update


class ParallelUpdateRule(UpdateRule):

    def __init__(self, removal_candidate_subtree: CandidateSubtree):
        super().__init__(removal_candidate_subtree)

    def apply_pre_sequeltialization_rule(self, tree_to_update: ProcessTree,
                                         child_ids_to_sequentialize: list[int]) -> ProcessTree:
        assert isinstance(tree_to_update, ProcessTree)
        assert tree_to_update is not None
        assert child_ids_to_sequentialize is not None
        assert len(child_ids_to_sequentialize) > 0

        tree_update_reference_node = find_tree_node_by_id(
            self.removal_candidate_subtree.node_id, tree_to_update
        )

        assert len(tree_update_reference_node.children) >= 3

        sequence: ProcessTree = ProcessTree()
        sequence.operator = Operator.SEQUENCE
        sequence.parent = tree_update_reference_node.parent
        sequence.label = None

        if tree_update_reference_node.parent is not None:
            for i in range(0, len(tree_update_reference_node.parent.children)):
                if tree_update_reference_node.parent.children[i].id == tree_update_reference_node.id:
                    tree_update_reference_node.parent.children[i] = sequence
                    break

        for id in child_ids_to_sequentialize:
            removed_child = remove_and_return_child_by_id(tree_update_reference_node, id)
            if removed_child is not None:
                removed_child.parent = sequence
                sequence.children.append(removed_child)

        sequence.children.append(tree_update_reference_node)
        tree_update_reference_node.parent = sequence

        return get_root(tree_to_update)

    def apply_post_sequeltialization_rule(self, tree_to_update: ProcessTree,
                                          child_ids_to_sequentialize: list[int]) -> ProcessTree:
        assert isinstance(tree_to_update, ProcessTree)
        assert tree_to_update is not None
        assert child_ids_to_sequentialize is not None
        assert len(child_ids_to_sequentialize) > 0

        tree_update_reference_node = find_tree_node_by_id(
            self.removal_candidate_subtree.node_id, tree_to_update
        )

        assert len(tree_update_reference_node.children) >= 3

        sequence: ProcessTree = ProcessTree()
        sequence.operator = Operator.SEQUENCE
        sequence.parent = tree_update_reference_node.parent
        sequence.label = None

        if tree_update_reference_node.parent is not None:
            for i in range(0, len(tree_update_reference_node.parent.children)):
                if tree_update_reference_node.parent.children[i].id == tree_update_reference_node.id:
                    tree_update_reference_node.parent.children[i] = sequence
                    break

        sequence.children.append(tree_update_reference_node)
        tree_update_reference_node.parent = sequence

        for id in child_ids_to_sequentialize:
            removed_child = remove_and_return_child_by_id(tree_update_reference_node, id)
            if removed_child is not None:
                removed_child.parent = sequence
                sequence.children.append(removed_child)

        return get_root(tree_to_update)

    def apply_mid_sequeltialization_rule(self, tree_to_update: ProcessTree,
                                         child_ids_to_sequentialize: dict[list[int], list[int], list[int]],
                                         execution_sequence_of_child_subtrees: list[int]) -> ProcessTree:
        assert isinstance(tree_to_update, ProcessTree)
        assert tree_to_update is not None
        assert child_ids_to_sequentialize is not None
        assert execution_sequence_of_child_subtrees is not None
        assert len(execution_sequence_of_child_subtrees) > 0

        tree_update_reference_node = find_tree_node_by_id(
            self.removal_candidate_subtree.node_id, tree_to_update
        )

        sequence: ProcessTree = ProcessTree()
        sequence.operator = Operator.SEQUENCE
        sequence.parent = tree_update_reference_node.parent
        sequence.label = None

        left_parallel: ProcessTree = ProcessTree()
        left_parallel.operator = Operator.PARALLEL
        left_parallel.parent = sequence
        left_parallel.label = None

        right_parallel: ProcessTree = ProcessTree()
        right_parallel.operator = Operator.PARALLEL
        right_parallel.parent = sequence
        right_parallel.label = None

        sequence.children.append(left_parallel)

        position = 'left'
        # look at the execution_sequence_of_child_subtrees to decide which subtrees to place on the left and right subtrees
        # removing duplicates to make sure repeating subtree is positioned (left or right) where it is observed first..
        execution_sequence_of_child_subtrees_no_duplicates = list(dict.fromkeys(execution_sequence_of_child_subtrees))
        for id in execution_sequence_of_child_subtrees_no_duplicates:

            # if position == 'left' and id in child_ids_to_sequentialize:
            #     position = 'middle'
            # if position == 'middle' and id not in child_ids_to_sequentialize:
            #     position = 'right'

            removed_child = remove_and_return_child_by_id(tree_update_reference_node, id)

            if id in child_ids_to_sequentialize['left'] and removed_child is not None:
                removed_child.parent = left_parallel
                left_parallel.children.append(removed_child)

            if id in child_ids_to_sequentialize['middle'] and removed_child is not None:
                removed_child.parent = sequence
                sequence.children.append(removed_child)

            if id in child_ids_to_sequentialize['right'] and removed_child is not None:
                removed_child.parent = right_parallel
                right_parallel.children.append(removed_child)

        sequence.children.append(right_parallel)

        if len(tree_update_reference_node.children) > 0:
            for child in tree_update_reference_node.children:
                removed_child = remove_and_return_child_by_id(tree_update_reference_node, child.id)
                removed_child.parent = left_parallel
                left_parallel.children.append(removed_child)

        if tree_update_reference_node.parent is None:
            tree_to_update = sequence
        elif tree_update_reference_node.parent is not None:
            for i in range(0, len(tree_update_reference_node.parent.children)):
                if tree_update_reference_node.parent.children[i].id == tree_update_reference_node.id:
                    tree_update_reference_node.parent.children[i] = sequence
                    break

        return get_root(tree_to_update)

    def apply_dynamic_sequeltialization_rule(self, tree_to_update: ProcessTree,
                                             child_ids_to_sequentialize: list[int]) -> ProcessTree:
        assert isinstance(tree_to_update, ProcessTree)
        assert tree_to_update is not None
        assert child_ids_to_sequentialize is not None
        assert len(child_ids_to_sequentialize) > 0

        tree_update_reference_node = find_tree_node_by_id(
            self.removal_candidate_subtree.node_id, tree_to_update
        )

        sequence: ProcessTree = ProcessTree()
        sequence.operator = Operator.SEQUENCE
        sequence.parent = tree_update_reference_node
        sequence.label = None

        for id in child_ids_to_sequentialize:
            removed_child = remove_and_return_child_by_id(tree_update_reference_node, id)
            if removed_child is not None:
                removed_child.parent = sequence
                sequence.children.append(removed_child)

        tree_update_reference_node.children.append(sequence)
        return tree_to_update


