from pm4py.objects.process_tree.utils.generic import parse

from cortado_core.negative_process_model_repair.temp_utils import calculate_process_tree_edit_distance

if __name__ == "__main__":
    # process_tree_1 = parse("+(+('a',+('b','c')),->(+('c',X('f','g'))))")
    # process_tree_2 = parse("+(+('a',+('b','c')),->(+('c',X('f','g')),'d'))")

    process_tree_1 = parse("->(X('d','c','f'),'e')")
    process_tree_2 = parse("->(X('d','g'),'e','h')")

    process_tree_3 = parse("->('d','c','f','e')")

    dist = calculate_process_tree_edit_distance(process_tree_1, process_tree_3)

    dist = dist
