from typing import Optional
import numpy as np
from tabulate import tabulate

DEBUG = 0


class RTriple:
    index = -1
    order = 0  # 0 is an invalid order, they are 1-indexed
    value = []

    def __init__(self, index, value):
        self.index = index
        self.value = value

    def __repr__(self):
        return f"index {self.index}, value {self.value}, order {self.order}"


def countingSort(array: list[RTriple], index_to_sort_on, place):
    size = len(array)
    output = [None] * size
    count = [0] * 10

    for i in range(0, size):
        index = array[i].value[index_to_sort_on] // place
        count[index % 10] += 1

    for i in range(1, 10):
        count[i] += count[i - 1]

    i = size - 1
    while i >= 0:
        index = array[i].value[index_to_sort_on] // place
        output[count[index % 10] - 1] = array[i]
        count[index % 10] -= 1
        i -= 1

    for i in range(0, size):
        array[i] = output[i]


# Radix sort algorithm lifted from https://www.programiz.com/dsa/radix-sort
def radixSort(triples: list[RTriple], index_to_sort_on):
    # Gets the portion of the triple being sorted on based on the index
    # Radix sort is called a fixed 6 times during the algorithm, so this is still linear
    relevant_values = [val.value[index_to_sort_on] for val in triples]
    max_element = max(relevant_values)

    place = 1
    while max_element // place > 0:
        countingSort(triples, index_to_sort_on, place)
        place *= 10

# Calculates the R12 indexes
def get_mods(arr):
    mods = []
    two_mods = []
    for i in range(len(arr) - 2):
        if i % 3 == 1:
            mods.append(i)
        elif i % 3 == 2:
            two_mods.append(i)
    mods.extend(two_mods)
    return mods


def get_r12_order(arr, mods):
    # Triples need to be made based on the R12 indexes first
    orders = np.array([RTriple(mod, arr[mod:mod + 3]) for mod in mods])

    # First pass of sorting - first number in the triple
    radixSort(orders, 0)
    if DEBUG: print("First pass:", orders)

    # Second pass of sorting - second number in the triple
    buckets = [[0, 0]]
    curr = orders[0].value[0]
    bucket_index = 0
    for i, val in enumerate(orders):
        if val.value[0] == curr:
            buckets[bucket_index][1] += 1
        else:
            buckets.append([i, i + 1])
            bucket_index += 1
        curr = val.value[0]
    # orders = []
    for bucket in buckets:
        radixSort(orders[bucket[0]:bucket[1]], 1)
        # orders.extend(bucket)
    if DEBUG: print("Second pass:", orders)

    # Sorting third digit
    bucket_index = 0
    curr = [orders[0].value[0], orders[0].value[1]]
    buckets = [[0, 0]]
    for i, val in enumerate(orders):
        if val.value[0] == curr[0] and val.value[1] == curr[1]:
            buckets[bucket_index][1] += 1
        else:
            buckets.append([i, i + 1])
            bucket_index += 1
        curr = [val.value[0], val.value[1]]
    for bucket in buckets:
        radixSort(orders[bucket[0]:bucket[1]], 2)
    if DEBUG: print("Third pass:", orders)

    # Assign R12 orders to the triples
    bucket_index = 0
    curr = [-1, -1, -1]
    for val in orders:
        if curr != val.value:
            bucket_index += 1
        val.order = bucket_index
        curr = val.value
    return orders

# Convert orders to suffix array
def triple_to_index(arr, mods):
    index_order_map = {triple.index: triple for triple in arr}
    result = [index_order_map[mod].order for mod in mods]
    result.extend([0, 0, 0])
    return result


# R0 equivalent of get_mods for the R12 orders
def get_r0_mods(arr):
    return np.array([index for index in range(len(arr) - 2) if index % 3 == 0])


def create_r0_pairs(arr, mods):
    return np.array([RTriple(mod, [arr[mod], mod]) for mod in mods])


def get_r0_order(arr, r12):
    # First pass
    arr = np.array(arr)
    radixSort(arr, 0)

    # Second pass
    buckets = [[]]
    curr = arr[0].value[0]
    bucket_index = 0
    for val in arr:
        if val.value[0] != curr:
            buckets.append([val])
            bucket_index += 1
        else:
            buckets[bucket_index].append(val)
        curr = val.value[0]
    orders = []

    # Second pass for R0 a little different - use R12 order to break ties
    r12_index_mapping = {r12_index: index for index, r12_index in enumerate(r12)}
    for bucket in buckets:
        if len(bucket) > 1:
            converted = np.array([RTriple(pair.index, [r12_index_mapping[pair.value[1] + 1] + 1, pair.value[1]]) for pair in
                         bucket])
            radixSort(converted, 0)
            orders.extend(converted)
        else:
            orders.extend(bucket)

    # Finalize
    orders = [order.value[1] for order in orders]

    return orders


def compare_orders(t, r12_val, r0_val, inverse_suffix_array) -> bool:
    # Pairs and triples can be compared to find out which is the next correct value in the sequence
    t_r12_val = t[r12_val] if r12_val <= len(t) else 0
    t_r0_val = t[r0_val] if r0_val <= len(t) else 0
    if t_r12_val < t_r0_val:
        return True
    if t_r12_val > t_r0_val:
        return False
    if r12_val % 3 != 0 and r0_val % 3 != 0:
        return inverse_suffix_array[r12_val] < inverse_suffix_array[r0_val]
    return compare_orders(t, r12_val + 1, r0_val + 1, inverse_suffix_array)


def merge_orders(t, r12, r0, remove_sentinel=False):
    inverse_suffix_array = {r12[i]: i for i in range(len(r12))}
    output = []
    r12_p, r0_p = 0, 0
    if remove_sentinel:
        if t[r0[0]] == 0:
            r0_p += 1
        else:
            r12_p += 1
    while r12_p < len(r12) and r0_p < len(r0):
        r12_val = r12[r12_p]
        r0_val = r0[r0_p]

        # R12 is next
        if compare_orders(t, r12_val, r0_val, inverse_suffix_array):
            output.append(r12_val)
            r12_p += 1
        #R0 is next
        else:
            output.append(r0_val)
            r0_p += 1
    output.extend(r12[r12_p:])
    output.extend(r0[r0_p:])
    return output

# Transform orders back into the original string indexes
def remap_suffixes(r12_indexes, t_suffixes, alg_input):
    final_r12 = [r12_indexes[i] for i in t_suffixes]

    final_r0_mods = get_r0_mods(alg_input)
    final_r0_pairs = create_r0_pairs(alg_input, final_r0_mods)
    final_r0_order = get_r0_order(final_r0_pairs, final_r12)
    if DEBUG: print(final_r0_order)

    result = merge_orders(alg_input, final_r12, final_r0_order, remove_sentinel=True)
    return result


def generate_suffix_array(input_string):
    # convert input to digits
    alg_input = [ord(c) for c in input_string]

    # append three zeroes
    alg_input.extend([0, 0, 0])

    original_mods = get_mods(alg_input)
    original_orders = get_r12_order(alg_input, original_mods)
    t_prime = triple_to_index(original_orders, original_mods)
    if DEBUG: print(t_prime)

    t_mods = get_mods(t_prime)
    t_orders = get_r12_order(t_prime, t_mods)
    r12 = [val.index for val in t_orders]
    if DEBUG: print(r12)

    r0_mods = get_r0_mods(t_prime)
    r0_pairs = create_r0_pairs(t_prime, r0_mods)
    r0_order = get_r0_order(r0_pairs, r12)

    merged = merge_orders(t_prime, r12, r0_order)[1:]  # remove sentinel
    final_suffix_array = remap_suffixes(original_mods, merged, alg_input)

    return final_suffix_array


def tabulate_suffix_array(suffix_array, original_input, lcps):
    table = [[index, original_input[index:], lcps[index]] for index in suffix_array]
    print(tabulate(table, headers=["Index", "Suffix", "LCP"], tablefmt='double_outline'))


def generate_lcps(suffix_array, original_input):
    # LCP array algorihtm adapted from https://www.cs.helsinki.fi/u/tpkarkka/teach/15-16/SPA/lecture10-2x4.pdf
    sa_length = len(suffix_array)
    lex = [None] * (sa_length - 1)
    plcp = [None] * (sa_length - 1)
    lcps = [-1] * sa_length
    for i in range(1, sa_length):
        lex[suffix_array[i]] = suffix_array[i - 1]
    l = 0
    for j in range(sa_length - 1):
        # Linear due to l being set to 0 at most n times, and l cannot be set to more than n
        while original_input[j + l] == original_input[lex[j] + l]:
            l += 1
        plcp[j] = l
        l = 0
    if DEBUG: print(plcp)
    for i in range(1, sa_length):
        lcps[suffix_array[i - 1]] = plcp[suffix_array[i]]

    return lcps


def translate_lcps(lcps, suffix_array):
    return [lcps[item] for index, item in enumerate(suffix_array)]


class Edge:
    def __init__(self, parent, child, label):
        self.parent = parent
        self.child = child
        self.edge_label = label


class Node:
    edges: list[Edge]
    parent_edge: Optional[Edge]
    node_label: int
    depth: int
    is_leaf: bool

    def __init__(self, label=-1, is_leaf=False):
        self.edges = []
        self.parent_edge = None
        self.node_label = label
        self.has_sifted = False
        self.is_leaf = is_leaf

    def add_child(self, child, edge_label=""):
        edge = Edge(self, child, edge_label)
        self.edges.append(edge)
        child.parent_edge = edge
        return child

    def get_parent(self):
        return self.parent_edge.parent if self.parent_edge else None


    def sift_up(self):
        parent = self.get_parent()
        if parent and parent.node_label > self.node_label:
            tmp = self.node_label
            tmp2 = parent.parent_edge.edge_label
            self.parent_edge = Edge(parent.parent_edge.parent, self, tmp2[:tmp])
            parent.parent_edge = Edge(self, parent, tmp2[tmp:])

            # Sift_up designed with the knowledge that we will only ever call this function if we are the
            # rightmost child of our parent and have no children ourselves. This keeps everything constant-time
            # as we can always pop the rightmost element rather than locating the item and then removing it.
            parent.edges.pop(-1)
            self.edges.append(parent.parent_edge)
            self.parent_edge.parent.edges.pop(-1)
            self.parent_edge.parent.edges.append(self.parent_edge)
            return self.sift_up()
        return self


def build_suffix_tree(original_input, suffix_array, lcp_array):
    root = Node(label=-1)
    curr: Optional[Node] = root
    n = len(original_input)
    prev_lcp = 0
    for index in range(n):
        lcp = lcp_array[index]
        suffix = suffix_array[index]
        if lcp <= 0 and prev_lcp == 0:
            root.add_child(Node(suffix, is_leaf=True), original_input[suffix:])
        elif lcp > 0 and prev_lcp == 0:
            curr = root.add_child(Node(lcp), original_input[suffix:suffix+lcp]).sift_up()
            curr.add_child(Node(suffix, is_leaf=True), original_input[suffix + lcp:])
        elif lcp == 0 and prev_lcp > 0:
            curr.add_child(Node(suffix, is_leaf=True), original_input[suffix + prev_lcp:])
        elif lcp == prev_lcp:
            curr.add_child(Node(suffix, is_leaf=True), original_input[suffix + prev_lcp:])
        elif lcp > 0:
            if prev_lcp < lcp:
                islice = original_input[-lcp:-prev_lcp]
            else:
                islice = original_input[prev_lcp:prev_lcp+lcp]
            curr = curr.add_child(Node(lcp), islice).sift_up()
            if len(curr.edges) > 0:
                # Current only has children if an LCP node has been sifted up - in which case this suffix
                # belongs to the *original* larger parent LCP instead of the current node
                curr.edges[-1].child.add_child(Node(suffix), original_input[curr.edges[-1].child.node_label + suffix:])
            else:
                curr.add_child(Node(suffix, is_leaf=True), original_input[suffix + lcp:])
        else:
            curr.add_child(Node(suffix, is_leaf=True), original_input[suffix + prev_lcp:])
        prev_lcp = lcp
    return root


def print_suffix_tree(node: Node, leaves_only=False):
    if len(node.edges) == 0:
        print(str(node.node_label) + "*", end='\n')  # signify leaf with *
    elif not leaves_only:
        print(str(node.node_label), end=' ')
    for edge in node.edges:
        if not leaves_only: print(edge.edge_label, end=' ')
        print_suffix_tree(edge.child, leaves_only)
        if not leaves_only: print(node.node_label, end=' ')


def main(orig_input, print_output=True):
    final = generate_suffix_array(orig_input)
    lcps = generate_lcps(final, orig_input)
    tree = build_suffix_tree(orig_input, final, translate_lcps(lcps, final))
    if print_output:
        print("\nOriginal Input:", orig_input)
        tabulate_suffix_array(final, orig_input, lcps)
        print("\nTraversal of generated suffix tree (Eulerian):")
        print_suffix_tree(tree)
        print("\nTraversal of generated suffix tree w/leaves only (Eulerian):")
        print_suffix_tree(tree, True)


if __name__ == "__main__":
    main(input("Enter in a string: "))

