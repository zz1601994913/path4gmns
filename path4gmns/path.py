import ctypes
import collections
import heapq
import os.path
from sys import platform

from .consts import MAX_LABEL_COST

# for precheck on connectivity of each OD pair
# 0: isolated, has neither outgoing links nor incoming links
# 1: has at least one outgoing link
# 2: has at least one incoming link
# 3: has both outgoing and incoming links
_zone_degrees = {}


def _single_source_shortest_path_fifo(G, origin_node_no):
    """ FIFO implementation of MLC using built-in list and indicator array

    The caller is responsible for initializing node_label_cost,
    node_predecessor, and link_predecessor.
    """
    G.node_label_cost[origin_node_no] = 0  # 将起点的cost设置为0
    # node status array
    status = [0] * G.node_size  # 所有点的status设置为0,status表示点是否在SElist中,能够一定程度上减少点被重复使用的次数

    # scan eligible list
    SEList = []
    SEList.append(origin_node_no)
    status[origin_node_no] = 1
    # label correcting
    while SEList:
        from_node = SEList.pop(0)
        status[from_node] = 0
        for link in G.node_list[from_node].outgoing_link_list:
            to_node = link.to_node_seq_no
            new_to_node_cost = (G.node_label_cost[from_node]
                                + link.cost)
            # we only compare cost at the downstream node ToID
            # at the new arrival time t
            if new_to_node_cost < G.node_label_cost[to_node]:
                # update cost label and node/time predecessor
                G.node_label_cost[to_node] = new_to_node_cost
                # pointer to previous physical node index
                # from the current label at current node and time
                G.node_predecessor[to_node] = from_node
                # pointer to previous physical node index
                # from the current label at current node and time
                G.link_predecessor[to_node] = link.link_seq_no
                if not status[to_node]:
                    SEList.append(to_node)
                    status[to_node] = 1


def _single_source_shortest_path_deque(G, origin_node_no):
    """ Deque implementation of MLC using deque list and indicator array

    The caller is responsible for initializing node_label_cost,
    node_predecessor, and link_predecessor.

    Adopted and modified from
    https://github.com/jdlph/shortest-path-algorithms
    """
    G.node_label_cost[origin_node_no] = 0
    # node status array
    status = [0] * G.node_size  # 一开始所有点的status均为0
    # scan eligible list
    SEList = collections.deque()
    SEList.append(origin_node_no)

    # label correcting
    while SEList:
        from_node = SEList.popleft()  # 每次都从左边取
        status[from_node] = 2  # 将更新源点的status改为2表示该点被作为标签更新源点使用过了
        for link in G.node_list[from_node].outgoing_link_list:
            to_node = link.to_node_seq_no
            new_to_node_cost = (G.node_label_cost[from_node]
                                + link.cost)
            # we only compare cost at the downstream node ToID
            # at the new arrival time t
            if new_to_node_cost < G.node_label_cost[to_node]:
                # update cost label and node/time predecessor
                G.node_label_cost[to_node] = new_to_node_cost
                # pointer to previous physical node index
                # from the current label at current node and time
                G.node_predecessor[to_node] = from_node
                # pointer to previous physical node index
                # from the current label at current node and time
                G.link_predecessor[to_node] = link.link_seq_no
                if status[to_node] != 1:  # 0或2
                    if status[to_node] == 2:  # 被作为标签更新源点使用过
                        SEList.appendleft(to_node)  # 添加到最左边首先使用
                    else:
                        SEList.append(to_node)
                    status[to_node] = 1


def _single_source_shortest_path_dijkstra(G, origin_node_no):
    """ Simplified heap-Dijkstra's Algorithm using heapq

    The caller is responsible for initializing node_label_cost,
    node_predecessor, and link_predecessor.

    Adopted and modified from
    https://github.com/jdlph/shortest-path-algorithms
    """
    G.node_label_cost[origin_node_no] = 0
    # node status array
    status = [0] * G.node_size
    # scan eligible list
    SEList = []
    heapq.heapify(SEList)
    heapq.heappush(SEList, (G.node_label_cost[origin_node_no], origin_node_no))

    # label setting
    while SEList:
        (label_cost, from_node) = heapq.heappop(SEList)
        # already scanned, pass it
        if status[from_node] == 1:
            continue  # 因为heap堆结构具有排序的功能，此处是默认没有负数弧，因此已经被作为的源节点的点可以直接略过
        status[from_node] = 1
        for link in G.node_list[from_node].outgoing_link_list:
            to_node = link.to_node_seq_no
            new_to_node_cost = label_cost + link.cost
            # we only compare cost at the downstream node ToID
            # at the new arrival time t
            if new_to_node_cost < G.node_label_cost[to_node]:
                # update cost label and node/time predecessor
                G.node_label_cost[to_node] = new_to_node_cost
                # pointer to previous physical node index
                # from the current label at current node and time
                G.node_predecessor[to_node] = from_node
                # pointer to previous physical node index
                # from the current label at current node and time
                G.link_predecessor[to_node] = link.link_seq_no
                heapq.heappush(SEList, (G.node_label_cost[to_node], to_node))


def single_source_shortest_path(G, origin_node_id,
                                engine_type='p', sp_algm='deque'):
    origin_node_no = G.get_node_no(origin_node_id)  # node_id→Node

    if engine_type.lower() == 'c':  # lower：大写转小写
        G.allocate_for_CAPI()  # 作用不详

    else:
        # just in case user uses C++ and Python path engines in a mixed way
        G.has_capi_allocated = False

        # Initialization for all nodes
        G.node_label_cost = [MAX_LABEL_COST] * G.node_size
        # pointer to previous node index from the current label at current node
        G.node_predecessor = [-1] * G.node_size
        # pointer to previous node index from the current label at current node
        G.link_predecessor = [-1] * G.node_size

        # make sure node_label_cost, node_predecessor, and link_predecessor
        # are initialized even the source node has no outgoing links
        if not G.node_list[origin_node_no].outgoing_link_list:  # 即使origin_node没有出弧也不会报错
            return

        if sp_algm.lower() == 'fifo':
            _single_source_shortest_path_fifo(G, origin_node_no)
        elif sp_algm.lower() == 'deque':
            _single_source_shortest_path_deque(G, origin_node_no)
        elif sp_algm.lower() == 'dijkstra':
            _single_source_shortest_path_dijkstra(G, origin_node_no)
        else:
            raise Exception('Please choose correct shortest path algorithm: '
                            + 'fifo or deque or dijkstra')


def _get_path_cost(G, to_node_id):
    to_node_no = G.node_id_to_no_dict[to_node_id]
    return G.node_label_cost[to_node_no]


def output_path_sequence(G, to_node_id, type='node'):
    """ output shortest path in terms of node sequence or link sequence

    Note that this function returns GENERATOR rather than list.
    """
    path = []
    current_node_seq_no = G.node_id_to_no_dict[to_node_id]

    if type.startswith('node'):
        # retrieve the sequence backwards
        while current_node_seq_no >= 0:
            path.append(current_node_seq_no)
            current_node_seq_no = G.node_predecessor[current_node_seq_no]
        # reverse the sequence
        for node_seq_no in reversed(path):
            yield G.node_no_to_id_dict[node_seq_no]
    else:
        # retrieve the sequence backwards
        current_link_seq_no = G.link_predecessor[current_node_seq_no]
        while current_link_seq_no >= 0:
            path.append(current_link_seq_no)
            current_node_seq_no = G.node_predecessor[current_node_seq_no]
            current_link_seq_no = G.link_predecessor[current_node_seq_no]
        # reverse the sequence
        for link_seq_no in reversed(path):
            yield G.link_list[link_seq_no].get_link_id()


def find_shortest_path(G, from_node_id, to_node_id, seq_type='node'):
    if from_node_id not in G.node_id_to_no_dict.keys():
        raise Exception(f"Node ID: {from_node_id} not in the network")
    if to_node_id not in G.node_id_to_no_dict.keys():
        raise Exception(f"Node ID: {to_node_id} not in the network")

    single_source_shortest_path(G, from_node_id, engine_type='p')

    path_cost = _get_path_cost(G, to_node_id)

    if path_cost == MAX_LABEL_COST:
        return f'distance: infinitity | path: '

    path = ';'.join(
        str(x) for x in output_path_sequence(G, to_node_id, seq_type)
    )

    return f'distance: {path_cost:.2f} | path: {path}'


def find_path_for_agents(G, column_pool, engine_type='c'):
    """ find and set up shortest path for each agent

    the internal node and links will be used to set up the node sequence and
    link sequence respectively

    Note that we do not cache the predecessors and label cost even some agents
    may share the same origin and each call of the single-source path algorithm
    will calculate the shortest path tree from the source node.
    """
    if G.get_agent_count() == 0:
        print('setting up individual agents')
        G.setup_agents(column_pool)

    from_node_id_prev = -1
    for agent in G.agent_list:
        from_node_id = agent.o_node_id
        to_node_id = agent.d_node_id

        # just in case agent has the same origin and destination
        if from_node_id == to_node_id:
            continue

        if from_node_id not in G.node_id_to_no_dict.keys():
            raise Exception(f'Node ID: {from_node_id} not in the network')
        if to_node_id not in G.node_id_to_no_dict.keys():
            raise Exception(f'Node ID: {to_node_id} not in the network')

        # simple caching strategy
        # if the current from_node_id is the same as from_node_id_prev,
        # then there is no need to redo shortest path calculation.
        if from_node_id != from_node_id_prev:
            from_node_id_prev = from_node_id
            single_source_shortest_path(G, from_node_id, engine_type)

        node_path = []
        link_path = []

        current_node_seq_no = G.node_id_to_no_dict[to_node_id]
        # set up the cost
        agent.path_cost = G.node_label_cost[current_node_seq_no]

        # retrieve the sequence backwards
        while current_node_seq_no >= 0:
            node_path.append(current_node_seq_no)
            current_link_seq_no = G.link_predecessor[current_node_seq_no]
            if current_link_seq_no >= 0:
                link_path.append(current_link_seq_no)
            current_node_seq_no = G.node_predecessor[current_node_seq_no]

        # make sure it is a valid path
        if not link_path:
            continue

        agent.node_path = [x for x in node_path]
        agent.link_path = [x for x in link_path]

# ----------------------------------------------------------------------------------
