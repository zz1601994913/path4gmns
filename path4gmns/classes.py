import ctypes
from copy import deepcopy

from .path import find_path_for_agents, find_shortest_path, \
    single_source_shortest_path
from .consts import MAX_LABEL_COST, SMALL_DIVISOR


class Assignment:

    def __init__(self):
        self.agent_types = []
        self.demand_periods = []
        self.demands = []
        # 4-d array
        self.column_pool = {}
        self.network = None
        self.spnetworks = []
        self.accessnetwork = None
        self.memory_blocks = 4
        self.map_atstr_id = {}
        self.map_dpstr_id = {}
        self.map_name_atstr = {}

    def update_agent_types(self, at):
        if at.get_type_str() not in self.map_atstr_id:
            self.map_atstr_id[at.get_type_str()] = at.get_id()
        else:
            raise Exception('agent type is not unique:' + at.get_type_str())

        if at.get_name() not in self.map_name_atstr:
            self.map_name_atstr[at.get_name()] = at.get_type_str()
        else:
            raise Exception('agent type name is not unique:' + at.get_name())

        self.agent_types.append(at)

    def update_demand_periods(self, dp):
        if dp.get_period() not in self.map_dpstr_id:
            self.map_dpstr_id[dp.get_period()] = dp.get_id()
        else:
            raise Exception('demand period is not unique:' + dp.get_period())

        self.demand_periods.append(dp)

    def update_demands(self, d):
        self.demands.append(d)

    def get_agent_type_count(self):
        return len(self.agent_types)

    def get_demand_period_count(self):
        return len(self.demand_periods)

    def get_agent_type_id(self, at_str):
        try:
            return self.map_atstr_id[at_str]
        except KeyError:
            raise Exception('NO agent type: ' + at_str)

    def get_demand_period_id(self, dp_str):
        try:
            return self.map_dpstr_id[dp_str]
        except KeyError:
            raise Exception('NO demand period: ' + dp_str)

    def get_agent_type(self, at_str):
        return self.agent_types[self.get_agent_type_id(at_str)]

    def get_demand_period(self, dp_str):
        return self.demand_periods[self.get_demand_period_id(dp_str)]

    def setup_spnetwork(self):
        spvec = {}

        # z is zone id starting from 1
        for z in self.network.zones:
            if z == -1:
                continue

            for d in self.demands:
                at = self.get_agent_type(d.get_agent_type_str())
                dp = self.get_demand_period(d.get_period())
                if z - 1 < self.memory_blocks:
                    sp = SPNetwork(self.network, at, dp)
                    spvec[(at.get_id(), dp.get_id(), z - 1)] = sp
                    sp.orig_zones.append(z)
                    sp.add_orig_nodes(self.network.get_nodes_from_zone(z))
                    for node_id in self.network.get_nodes_from_zone(z):
                        sp.node_id_to_no[node_id] = self.network.get_node_no(node_id)
                    self.spnetworks.append(sp)
                else:
                    m = (z - 1) % self.memory_blocks
                    sp = spvec[(at.get_id(), dp.get_id(), m)]
                    sp.orig_zones.append(z)
                    sp.add_orig_nodes(self.network.get_nodes_from_zone(z))
                    for node_id in self.network.get_nodes_from_zone(z):
                        sp.node_id_to_no[node_id] = self.network.get_node_no(node_id)

    def _convert_mode(self, mode):
        """convert mode to the corresponding agent type name and string"""
        if mode in self.map_atstr_id:
            at = self.get_agent_type(mode)
            return at.get_name(), mode

        if mode in self.map_name_atstr:
            return mode, self.map_name_atstr[mode]

        # for distance-based shortest path calculation only
        # it shall not be used with any accessibility evaluations
        if mode == 'all':
            return mode, mode

        raise Exception('Please provide a valid mode!')

    def find_shortest_path(self, from_node_id, to_node_id, mode, seq_type='node'):
        """ call find_shortest_path() from path.py

        exceptions will be handled in find_shortest_path()
        """
        # reset agent type str or mode according to user's input
        at_name, _ = self._convert_mode(mode)
        self.network.set_agent_type_name(at_name)

        return find_shortest_path(self.network, from_node_id,
                                  to_node_id, seq_type)

    def find_path_for_agents(self, mode):
        """ find and set up shortest path for each agent """
        # reset agent type str or mode according to user's input
        at_name, _ = self._convert_mode(mode)
        self.network.set_agent_type_name(at_name)

        find_path_for_agents(self.network, self.column_pool)

    def get_accessible_nodes(self, source_node_id, time_budget,
                             mode, time_dependent, tau):
        if source_node_id not in self.network.node_id_to_no_dict.keys():
            raise Exception(f'Node ID: {source_node_id} not in the network')

        assert(time_budget>=0)

        if time_budget == 0:
            return []

        if not self.accessnetwork:
            self.accessnetwork = AccessNetwork(self.network, False)

        # simple caching to avoid duplicate shortest path calculation
        run_sp = False
        if self.accessnetwork.pre_source_node_id != source_node_id:
            self.accessnetwork.set_source_node_id(source_node_id)
            run_sp = True

        at_name, at_str = self._convert_mode(mode)
        if self.accessnetwork.agent_type_name != at_name:
            self.accessnetwork.set_target_mode(at_name)
            at = self.get_agent_type(at_str)
            self.accessnetwork.update_generalized_link_cost(at,
                                                            time_dependent,
                                                            tau)#tau就是demand_period
            run_sp = True

        if run_sp:
            single_source_shortest_path(self.accessnetwork, source_node_id)

        # if max min travel time is less than or equal to time_budget,
        # output the entire node set directly without the following check?
        nodes = []
        for node in self.accessnetwork.get_nodes():
            node_no = node.get_node_no()
            # do not include the source node itself
            if node.get_node_id() == source_node_id:
                continue
            if self.accessnetwork.get_node_label_cost(node_no) <= time_budget:
                nodes.append(node.get_node_id())

        return nodes

    def get_accessible_links(self, source_node_id, time_budget,
                             mode, time_dependent, tau):
        # node id's
        nodes = self.get_accessible_nodes(source_node_id, time_budget,
                                          mode, time_dependent, tau)
        # convert to link id's
        return [
            self.accessnetwork.get_pred_link_id(x) for x in nodes
        ]

    def get_agent_types(self):
        return self.agent_types



class Agent:
    """ individual agent derived from aggragted demand between an OD pair

    agent_id: integer starts from 1
    agent_seq_no: internal agent index starting from 0 used for calculation
    """

    def __init__(self, agent_id, agent_seq_no, agent_type,
                 o_zone_id, d_zone_id):
        """ the attribute of agent """
        self.agent_id = agent_id
        self.agent_seq_no = agent_seq_no
        # vehicle
        self.agent_type = agent_type
        self.o_zone_id = o_zone_id
        self.d_zone_id = d_zone_id
        self.o_node_id = 0
        self.d_node_id = 0
        self.node_path = None
        self.link_path = None
        self.current_link_seq_no_in_path = 0
        self.departure_time_in_min = 0
        # Passenger Car Equivalent (PCE) of the agent
        self.PCE_factor = 1
        self.path_cost = 0
        # self.departure_time_in_simu_interval = int(
        #     self.departure_time_in_min
        #     * 60 /_NUM_OF_SECS_PER_SIMU_INTERVAL
        #     + 0.5)
        self.b_generated = False
        self.b_complete_trip = False
        self.feasible_path_exist_flag = False


class Network:

    def __init__(self):
        self.node_list = []
        self.link_list = []
        self.agent_list = []
        self.node_size = 0
        self.link_size = 0
        self.agent_size = 0
        # key: node id, value: node seq no
        self.node_id_to_no_dict = {}
        # key: node seq no, value: node id
        self.node_no_to_id_dict = {}
        # map link id to link seq no
        self.link_id_dict = {}
        # td:time-dependent, key:simulation time interval,
        # value:agents(list) need to be activated
        self.agent_td_list_dict = {}
        # key: zone id, value: node id list
        self.zone_to_nodes_dict = {}
        self.node_label_cost = None
        self.node_predecessor = None
        self.link_predecessor = None
        # added for CG
        self.zones = None
        self.has_capi_allocated = False
        # the following two are IDs rather than objects
        self._agent_type_size = 1
        self._demand_period_size = 1
        self.agent_type_name = 'all'

    def update(self, agent_type_size, demand_period_size):
        self.node_size = len(self.node_list)
        self.link_size = len(self.link_list)
        self.agent_size = len(self.agent_list)
        # it is needed for setup_spnetwork() and setup_spnetwork_a()
        self.zones = sorted(self.zone_to_nodes_dict.keys())
        self._agent_type_size = agent_type_size
        self._demand_period_size = demand_period_size

    def get_node_no(self, node_id):
        return self.node_id_to_no_dict[node_id]

    def get_nodes_from_zone(self, zone_id):
        return self.zone_to_nodes_dict[zone_id]

    def allocate_for_CAPI(self):
        # execute only on the first call
        if self.has_capi_allocated:
            return

        node_size = self.node_size
        link_size = self.link_size

        # initialization for predecessors and label costs
        node_predecessor = [-1] * node_size
        link_predecessor = [-1] * node_size
        node_label_cost = [MAX_LABEL_COST] * node_size

        # initialize from_node_no_array, to_node_no_array, and link_cost_array
        from_node_no_array = [link.from_node_seq_no for link in self.link_list]
        to_node_no_array = [link.to_node_seq_no for link in self.link_list]
        link_cost_array = [link.cost for link in self.link_list]

        # initialize others as numpy arrays directly
        queue_next = [0] * node_size
        first_link_from = [-1] * node_size
        last_link_from = [-1] * node_size
        sorted_link_no_array = [-1] * link_size

        # internal link index used for shortest path calculation only
        j = 0
        for i, node in enumerate(self.node_list):
            if not node.outgoing_link_list:
                continue
            first_link_from[i] = j
            for link in node.outgoing_link_list:
                # set up the mapping from j to the true link seq no
                sorted_link_no_array[j] = link.link_seq_no
                j += 1
            last_link_from[i] = j

        # setup allowed uses
        # allowed_uses = [''] * link_size
        # self._setup_allowed_use(allowed_uses)
        allowed_uses = [link.allowed_uses for link in self.link_list]

        # set up arrays using ctypes
        int_arr_node = ctypes.c_int * node_size
        int_arr_link = ctypes.c_int * link_size
        double_arr_node = ctypes.c_double * node_size
        double_arr_link = ctypes.c_double * link_size
        # for allowed_uses
        char_arr_link = ctypes.c_wchar_p * link_size

        self.from_node_no_array = int_arr_link(*from_node_no_array)
        self.to_node_no_array = int_arr_link(*to_node_no_array)
        self.first_link_from = int_arr_node(*first_link_from)
        self.last_link_from = int_arr_node(*last_link_from)
        self.sorted_link_no_array = int_arr_link(*sorted_link_no_array)
        self.link_cost_array = double_arr_link(*link_cost_array)
        self.node_label_cost = double_arr_node(*node_label_cost)
        self.node_predecessor = int_arr_node(*node_predecessor)
        self.link_predecessor = int_arr_node(*link_predecessor)
        self.queue_next = int_arr_node(*queue_next)
        self.allowed_uses = char_arr_link(*allowed_uses)

        self.has_capi_allocated = True

    def set_agent_type_name(self, at_name):
        self.agent_type_name = at_name

    def get_agent_count(self):
        return self.agent_size

    def get_zones(self):
        return self.zones

    def get_link(self, seq_no):
        return self.link_list[seq_no]

    def get_links(self):
        return self.link_list

    def get_nodes(self):
        return self.node_list

    def get_node_size(self):
        return self.node_size

    def get_link_size(self):
        return self.link_size



class AgentType:

    def __init__(self, id=0, type='p', name='passenger',
                 vot=10, flow_type=0, pce=1, ffs=60):
        """ default constructor """
        self.id = id
        self.type = type
        self.name = name
        self.vot = vot
        self.flow_type = flow_type
        self.pce = pce
        self.ffs = ffs

    def get_id(self):
        return self.id

    def get_name(self):
        return self.name

    def get_vot(self):
        return self.vot

    def get_type_str(self):
        return self.type

    def get_pce(self):
        return self.pce

    def get_free_flow_speed(self):
        return self.ffs


class DemandPeriod:

    def __init__(self, id=0, period='AM', time_period='0700_0800'):
        self.id = id
        self.period = period
        self.time_period = time_period

    def get_id(self):
        return self.id

    def get_period(self):
        return self.period


class Demand:

    def __init__(self, id=0, period='AM', agent_type='p', file='demand.csv'):
        self.id = id
        self.period = period
        self.agent_type_str = agent_type
        self.file = file

    def get_id(self):
        return self.id

    def get_file_name(self):
        return self.file

    def get_period(self):
        return self.period

    def get_agent_type_str(self):
        return self.agent_type_str


class Node:

    def __init__(self, node_seq_no, node_id, zone_id, x='', y=''):
        """ the attributes of node  """
        # node_seq_no: internal node index used for calculation
        self.node_seq_no = node_seq_no
        # node_id: user defined node id from input
        self.node_id = node_id
        self.outgoing_link_list = []
        self.incoming_link_list = []
        self.zone_id = zone_id
        self.coord_x = x
        self.coord_y = y

    def add_outgoing_link(self, link):
        self.outgoing_link_list.append(link)

    def add_incoming_link(self, link):
        self.incoming_link_list.append(link)

    def get_zone_id(self):
        return self.zone_id

    def update_coordinate(self, args):
        self.coord_x = args[0]
        self.coord_y = args[1]

    def get_coordinate(self):
        return self.coord_x + ' ' + self.coord_y

    def get_node_id(self):
        return self.node_id

    def get_node_no(self):
        return self.node_seq_no



class Link:

    def __init__(self,
                 id,
                 link_seq_no,
                 from_node_no,
                 to_node_no,
                 from_node_id,
                 to_node_id,
                 length,
                 lanes=1,
                 link_type=1,
                 free_speed=60,
                 capacity=49500,
                 allowed_uses='all',
                 geometry='',
                 agent_type_size=1,
                 demand_period_size=1):
        """ the attributes of link """
        self.id = id
        self.link_seq_no = link_seq_no
        self.from_node_seq_no = from_node_no
        self.to_node_seq_no = to_node_no
        self.from_node_id = from_node_id
        self.to_node_id = to_node_id
        # length is mile or km
        self.length = length
        self.lanes = lanes
        # 1:one direction, 2:two way
        self.type = link_type
        # length:km, free_speed: km/h
        self.free_flow_travel_time_in_min = (
                length / max(SMALL_DIVISOR, free_speed) * 60
        )
        # capacity is lane capacity per hour
        self.link_capacity = capacity * lanes
        self.allowed_uses = allowed_uses
        self.geometry = geometry
        self.cost = self.free_flow_travel_time_in_min
        self.flow_volume = 0
        # add for CG
        self.agent_type_size = agent_type_size
        self.demand_period_size = demand_period_size
        self.toll = 0
        self.route_choice_cost = 0
        self.travel_time_by_period = [0] * demand_period_size
        self.flow_vol_by_period = [0] * demand_period_size
        self.vdfperiods = []
        # Peiheng, 04/05/21, not needed for the current implementation
        # self.queue_length_by_slot = [0] * MAX_TIME_PERIODS
        # self.vol_by_period_by_at = [
        #     [0] * demand_period_size for i in range(agent_type_size)
        # ]
        # self.travel_marginal_cost_by_period = [
        #     [0] * demand_period_size for i in range(agent_type_size)
        # ]

    def get_length(self):
        return self.length

    def get_link_id(self):
        return self.id

    def get_seq_no(self):
        return self.link_seq_no

    def get_period_fftt(self, tau):
        try:
            return self.vdfperiods[tau].get_fftt()
        except IndexError:
            raise Exception(f'NO such demand period id: {tau}!'
                            ' Check your input demand_period_id')

    def get_route_choice_cost(self):
        return self.route_choice_cost

    def get_toll(self):
        return self.toll

    def get_free_flow_travel_time(self):
        return self.free_flow_travel_time_in_min





class VDFPeriod:

    def __init__(self, id, alpha=0.15, beta=4, mu=1000,
                 fftt=0.0, cap=99999, phf=-1):
        """ default constructor """
        self.id = id
        # the following four have been defined in class Link
        # they should be exactly the same with those in the corresponding link
        self.alpha = alpha
        self.beta = beta
        self.mu = mu
        # free flow travel time
        self.fftt = fftt
        self.capacity = cap
        self.phf = phf
        self.marginal_base = 1
        self.avg_travel_time = 0
        self.voc = 0

    def get_fftt(self):
        return self.fftt


class SPNetwork(Network):  # 作用不详
    """ attributes related to outputs from shortest path calculations """

    def __init__(self, base, at, dp):
        self.base = base
        # AgentType object
        self.agent_type = at
        # DemandPeriod object
        self.demand_period = dp

        # this is necessary for each instance of SPNetwork
        # to retrieve network topoloy
        if not base.has_capi_allocated:
            base.allocate_for_CAPI()

        # set up attributes unique to each instance
        node_preds = [-1] * base.node_size
        link_preds = [-1] * base.node_size
        node_lables = [MAX_LABEL_COST] * base.node_size
        queue_next = [0] * base.node_size
        link_cost_array = [link.cost for link in base.link_list]

        int_arr_node = ctypes.c_int * base.node_size
        double_arr_node = ctypes.c_double * base.node_size
        double_arr_link = ctypes.c_double * base.link_size

        self.node_predecessor = int_arr_node(*node_preds)
        self.link_predecessor = int_arr_node(*link_preds)
        self.node_label_cost = double_arr_node(*node_lables)
        self.link_cost_array = double_arr_link(*link_cost_array)
        self.queue_next = int_arr_node(*queue_next)

        # node id
        self.orig_nodes = []
        # zone sequence no
        self.orig_zones = []
        self.node_id_to_no = {}
        self.has_capi_allocated = True

    def add_orig_nodes(self, nodes):
        self.orig_nodes.extend(nodes)


class AccessNetwork(Network):
    """ network for accessibility evaluation """

    def __init__(self, base, add_cc=True):
        self.base = base
        self.node_list = self.base.get_nodes()
        self.link_list = self.base.get_links()
        self.map_id_to_no = self.base.node_id_to_no_dict
        self.map_no_to_id = self.base.node_no_to_id_dict
        self.node_size = base.get_node_size()
        self.link_size = base.get_link_size()
        self.centroids = []
        self.agent_type_name = 'all'
        self.pre_source_node_id = -1
        if add_cc:
            self._add_centroids_connectors()
        self.has_capi_allocated = False
        super().allocate_for_CAPI()

    def _get_zone_coord(self, zone_id):
        """ coordinate of each zone is from its first node
        该区域的第一个点的坐标
        """
        node_id = self.get_nodes_from_zone(zone_id)[0]
        node_seq_no = self.base.get_node_no(node_id)
        node = self.base.get_nodes()[node_seq_no]
        return node.coord_x, node.coord_y

    def _add_centroids_connectors(self):  # 连接区域质心与区域内其他点
        # deep copy
        self.node_list = deepcopy(self.node_list)
        self.link_list = deepcopy(self.link_list)

        node_seq_no = self.node_size
        link_seq_no = self.link_size
        # get zones
        for z in self.get_zones():
            if z == -1:
                continue

            # create a centroid
            node_id = 'c_' + str(z)
            centroid = Node(node_seq_no, node_id, z)  # node_seq_no=node_size+1, node_id, zone_id
            centroid.update_coordinate(self._get_zone_coord(z))

            self.node_list.append(centroid)  # 把质点也作为点加入node_list中，该点id为c_区域第一个点的id
            self.centroids.append(centroid)

            self.map_id_to_no[node_id] = node_seq_no
            self.map_no_to_id[node_seq_no] = node_id

            # build connectors

            for i in self.get_nodes_from_zone(z):  # 遍历区域的所有点，与质点连接，length=0
                link_id_f = 'conn_' + str(link_seq_no)
                from_node_no = node_seq_no
                to_node_id = i

                try:
                    to_node_no = self.map_id_to_no[i]
                except KeyError:
                    continue

                # connector from centroid to activity nodes in this zone
                c_forward = Link(link_id_f,
                                 link_seq_no,
                                 from_node_no,
                                 to_node_no,
                                 node_id,
                                 to_node_id,
                                 0)

                # connector from activity nodes in this zone to centroid
                link_id_b = 'conn_' + str(link_seq_no + 1)
                c_backward = Link(link_id_b,
                                  link_seq_no + 1,
                                  to_node_no,
                                  from_node_no,
                                  to_node_id,
                                  node_id,
                                  0)

                self.node_list[from_node_no].add_outgoing_link(c_forward)
                self.node_list[to_node_no].add_outgoing_link(c_backward)

                self.link_list.append(c_forward)
                self.link_list.append(c_backward)

                link_seq_no += 2

            node_seq_no += 1

        self.node_size = len(self.node_list)  # 跟新尺寸
        self.link_size = len(self.link_list)

    def set_target_mode(self, mode):
        """ set up the target mode for accessibility evaluation

        Parameters
        ----------
        mode : agent name which is in settings.yml.
        """
        self.agent_type_name = mode

    def get_zones(self):
        return self.base.get_zones()

    def get_centroids(self):
        return self.centroids

    def get_node_label_cost(self, node_no):
        return self.node_label_cost[node_no]

    def get_sp_distance(self, node_no):
        """ get the shortest path distance """
        dist = 0
        while node_no >= 0:
            link_seq_no = self.link_predecessor[node_no]
            if link_seq_no >= 0:  # link_seq_no!=-1
                dist += self.get_link(link_seq_no).get_length()

            node_no = self.node_predecessor[node_no]

        return dist

    def set_source_node_id(self, node_id):
        self.pre_source_node_id = node_id


    def update_generalized_link_cost(self, at, time_dependent, demand_period_id):
        """ update generalized link costs to calculate accessibility """
        vot = at.get_vot()

        if time_dependent:
            for link in self.get_links():
                # do not update connectors
                #不更新质心连接线
                if link.get_link_id().startswith('conn_'):
                    continue

                self.link_cost_array[link.get_seq_no()] = (
                    link.get_period_fftt(demand_period_id)#获得link在demand_period的自由流旅行时间
                    + link.get_route_choice_cost()
                    + link.get_toll() / min(SMALL_DIVISOR, vot) * 60
                )
        else:
            if at.get_type_str().startswith('p'):
                for link in self.get_links():
                    self.link_cost_array[link.get_seq_no()] = (
                        link.get_free_flow_travel_time()
                        + link.get_route_choice_cost()
                        + link.get_toll() / min(SMALL_DIVISOR, vot) * 60
                    )
            else:
                ffs = at.get_free_flow_speed()

                for link in self.get_links():
                    self.link_cost_array[link.get_seq_no()] = (
                        (link.get_length() / max(SMALL_DIVISOR, ffs) * 60)
                        + link.get_route_choice_cost()
                        + link.get_toll() / min(SMALL_DIVISOR, vot) * 60
                    )

    def get_nodes_from_zone(self, zone_id):
        return self.base.zone_to_nodes_dict[zone_id]
    def get_node_no(self, node_id):
        return self.map_id_to_no[node_id]

    def get_pred_link_id(self, node_id):
        """ return id of the predecessor link to node_id """
        link_no = self.link_predecessor[self.get_node_no(node_id)]
        return self.link_list[link_no].get_link_id()



class UI:
    """ an abstract class only with user interfaces """

    def __init__(self, assignment):
        self._base_assignment = assignment

    def find_shortest_path(self, from_node_id, to_node_id,
                           mode='all', seq_type='node'):
        """ return shortest path between from_node_id and to_node_id

        Parameters
        ----------
        from_node_id
            the starting node id

        to_node_id
            the ending node id

        seq_type
            'node' or 'link'. You will get the shortest path in sequence of
            either node IDs or link IDs. The default is 'node'.

        mode
            the target transportation mode which is defined in settings.yml. It
            can be either agent type or its name. For example, 'w' and 'walk'
            are equivalent inputs.

            The default is 'all', which means that links are open to all modes.
            不太能理解

        Outputs
        -------
        the shortest path between from_node_id and to_node_id

        Exceptions will be thrown if either of them or both are not valid node
        IDs.
        """
        return self._base_assignment.find_shortest_path(
            from_node_id,
            to_node_id,
            mode,
            seq_type
        )

    def find_path_for_agents(self, mode='all'):
        """ find and set up shortest path for each agent """
        return self._base_assignment.find_path_for_agents(mode)

    def get_accessible_nodes(self,
                             source_node_id,
                             time_budget,
                             mode='p',
                             time_dependent=False,
                             demand_period_id=0):
        """ get the accessible nodes from a node given mode and time budget

        Parameters
        ----------
        source_node_id
            the starting node id for evaluation

        time_budget
            the amount of time to travel in minutes

        mode
            the target transportation mode which is defined in settings.yml. It
            can be either agent type or its name. For example, 'w' and 'walk'
            are equivalent inputs. Its default value is 'p' (i.e., mode auto).

        time_dependent
            True or False. Its default value is False.

            If True, the accessibility will be evaluated using the period link
            free-flow travel time (i.e., VDF_fftt). In other words, the
            accessibility is time-dependent.

            If False, the accessibility will be evaluated using the link length
            and the free flow travel speed of each mode.

        demand_period_id
            The sequence number of demand period listed in demand_periods in
            settings.yml. demand_period_id of the first demand_period is 0.

            Use it with time_dependent when there are multiple demand periods.
            Its default value is 0.

        Outputs
        -------
        print out the number of nodes that can be accessible from \
        source_node_id given time_budget and mode, and the node list
        """
        nodes = self._base_assignment.get_accessible_nodes(source_node_id,
                                                           time_budget,
                                                           mode,
                                                           time_dependent,
                                                           demand_period_id)

        node_strs = ';'.join(str(x) for x in nodes)

        print(f'number of accessible nodes is {len(nodes)}')
        print(f'accessible nodes are: {node_strs}')

    def get_accessible_links(self,
                             source_node_id,
                             time_budget,
                             mode='p',
                             time_dependent=False,
                             demand_period_id=0):
        """ get the accessible links from a node given mode and time budget

        Parameters
        ----------
        source_node_id
            the starting node id for evaluation

        time_budget
            the amount of time to travel in minutes

        mode
            the target transportation mode which is defined in settings.yml. It
            can be either agent type or its name. For example, 'w' and 'walk'
            are equivalent inputs.


        Outputs
        -------
        print out the number of links that can be accessible from \
        source_node_id given time_budget and mode, and the link list
        """
        links = self._base_assignment.get_accessible_links(source_node_id,
                                                           time_budget,
                                                           mode,
                                                           time_dependent,
                                                           demand_period_id)

        link_strs = ';'.join(str(x) for x in links)

        print(f'number of accessible links is {len(links)}')
        print(f'accessible links are: {link_strs}')