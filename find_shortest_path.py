import csv
import ctypes
import collections
import heapq

SMALL_DIVISOR = 0.001
MAX_LABEL_COST = 999999

# for precheck on connectivity of each OD pair
# 0: isolated, has neither outgoing links nor incoming links
# 1: has at least one outgoing link
# 2: has at least one incoming link
# 3: has both outgoing and incoming links
_zone_degrees = {}


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

    def allocate_for_CAPI(self):  # 作用不详
        # execute only on the first call
        if self.has_capi_allocated:
            return

    def set_agent_type_name(self, at_name):
        self.agent_type_name = at_name


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


def _auto_setup(assignment):
    """ automatically set up one demand period and one agent type

    The two objects will be set up using the default construnctors using the
    default values. See class DemandPeriod and class AgentType for details
    """
    at = AgentType()
    dp = DemandPeriod()
    d = Demand()

    assignment.update_agent_types(at)
    assignment.update_demand_periods(dp)
    assignment.update_demands(d)


def _convert_str_to_int(str):
    """
    TypeError will take care the case that str is None
    ValueError will take care the case that str is empty
    """
    if not str:
        return None

    try:
        return int(str)
    except ValueError:
        return int(float(str))
    except TypeError:
        return None


def _convert_str_to_float(str):
    """
    TypeError will take care the case that str is None
    ValueError will take care the case that str is empty
    """
    if not str:
        return None

    try:
        return float(str)
    except (TypeError, ValueError):
        return None


def _update_orig_zone(oz_id):
    if oz_id not in _zone_degrees:
        _zone_degrees[oz_id] = 1
    elif _zone_degrees[oz_id] == 2:
        _zone_degrees[oz_id] = 3


def _update_dest_zone(dz_id):
    if dz_id not in _zone_degrees:
        _zone_degrees[dz_id] = 2
    elif _zone_degrees[dz_id] == 1:
        _zone_degrees[dz_id] = 3


# 将settings文件中的agents、demand periods、demands files添加到assignment的agent_types、demand_periods、demands中，读不了用默认值添加
def read_settings(input_dir, assignment):
    try:
        import yaml as ym

        with open(input_dir + '/settings.yml') as file:
            settings = ym.full_load(file)

            # agent types
            agents = settings['agents']
            for i, a in enumerate(agents):
                agent_type = a['type']
                agent_name = a['name']
                agent_vot = a['vot']
                agent_flow_type = a['flow_type']
                agent_pce = a['pce']
                agent_ffs = a['free_speed']

                at = AgentType(i,
                               agent_type,
                               agent_name,
                               agent_vot,
                               agent_flow_type,
                               agent_pce,
                               agent_ffs)

                assignment.update_agent_types(at)

            # demand periods
            demand_periods = settings['demand_periods']
            for i, d in enumerate(demand_periods):
                period = d['period']
                time_period = d['time_period']

                dp = DemandPeriod(i, period, time_period)
                assignment.update_demand_periods(dp)

            # demand files
            demands = settings['demand_files']
            for i, d in enumerate(demands):
                demand_file = d['file_name']
                # demand_format_tpye = d['format_type']
                demand_period = d['period']
                demand_type = d['agent_type']

                demand = Demand(i, demand_period, demand_type, demand_file)
                assignment.update_demands(demand)

    except ImportError:
        # just in case user does not have pyyaml installed
        print('Please intall pyyaml next time!')
        print('Engine will set up one demand period and one agent type using '
              'default values for you, which might NOT reflect your case!\n')
        _auto_setup(assignment)
    except FileNotFoundError:
        # just in case user does not provide settings.yml
        print('Please provide settings.yml next time!')
        print('Engine will set up one demand period and one agent type using '
              'default values for you, which might NOT reflect your case!\n')
        _auto_setup(assignment)
    except Exception as e:
        raise e


def read_nodes(input_dir,
               nodes,
               id_to_no_dict,
               no_to_id_dict,
               zone_to_node_dict):
    """ step 1: read input_node """
    with open(input_dir + '/node.csv', 'r', encoding='utf-8') as fp:
        print('read node.csv')

        reader = csv.DictReader(fp)
        node_seq_no = 0
        for line in reader:
            # set up node_id, which should be an integer
            node_id = _convert_str_to_int(line['node_id'])
            if node_id is None:
                continue

            # set up zone_id, which should be an integer
            zone_id = _convert_str_to_int(line['zone_id'])
            if zone_id is None:
                zone_id = -1

            # treat them as string
            coord_x = line['x_coord']
            coord_y = line['y_coord']

            # construct node object
            node = Node(node_seq_no, node_id, zone_id, coord_x, coord_y)
            nodes.append(node)

            # set up mapping between node_seq_no and node_id
            id_to_no_dict[node_id] = node_seq_no
            no_to_id_dict[node_seq_no] = node_id

            # associate node_id to corresponding zone
            if zone_id not in zone_to_node_dict.keys():
                zone_to_node_dict[zone_id] = []
            zone_to_node_dict[zone_id].append(node_id)

            node_seq_no += 1

        print(f"the number of nodes is {node_seq_no}")

        zone_size = len(zone_to_node_dict)
        # do not count virtual zone with id as -1
        if -1 in zone_to_node_dict.keys():
            zone_size -= 1

        print(f"the number of zones is {zone_size}")


def read_links(input_dir,
               links,
               nodes,
               id_to_no_dict,
               link_id_dict,
               agent_type_size,
               demand_period_size,
               load_demand):
    """ step 2: read input_link """
    with open(input_dir + '/link.csv', 'r', encoding='utf-8') as fp:
        print('read link.csv')

        reader = csv.DictReader(fp)
        link_seq_no = 0
        for line in reader:
            # it can be an empty string
            link_id = line['link_id']

            # check the validility
            from_node_id = _convert_str_to_int(line['from_node_id'])
            if from_node_id is None:
                continue

            to_node_id = _convert_str_to_int(line['to_node_id'])
            if to_node_id is None:
                continue

            length = _convert_str_to_float(line['length'])
            if length is None:
                continue

            # pass validility check

            try:
                from_node_no = id_to_no_dict[from_node_id]
                to_node_no = id_to_no_dict[to_node_id]
            except KeyError:
                print(f"EXCEPTION: Node ID {from_node_no} "
                      f"or/and Node ID {to_node_id} NOT IN THE NETWORK!!")
                continue

            # for the following attributes,
            # if they are not None, convert them to the corresponding types
            # if they are None's, set them using the default values
            lanes = _convert_str_to_int(line['lanes'])
            if lanes is None:
                lanes = 1

            link_type = _convert_str_to_int(line['link_type'])
            if link_type is None:
                link_type = 1

            free_speed = _convert_str_to_int(line['free_speed'])
            if free_speed is None:
                free_speed = 60

            # issue: int??
            capacity = _convert_str_to_int(line['capacity'])
            if capacity is None:
                capacity = 49500

            # if link.csv does not have no column 'allowed_uses',
            # set allowed_uses to 'auto'
            # developer's note:
            # we may need to change this implemenation as we cannot deal with
            # cases a link which is not open to any modes
            try:
                allowed_uses = line['allowed_uses']
                if not allowed_uses:
                    allowed_uses = 'all'
            except KeyError:
                allowed_uses = 'all'

            # if link.csv does not have no column 'geometry',
            # set geometry to ''
            try:
                geometry = line['geometry']
            except KeyError:
                geometry = ''

            link_id_dict[link_id] = link_seq_no

            # construct link ojbect
            link = Link(link_id,
                        link_seq_no,
                        from_node_no,
                        to_node_no,
                        from_node_id,
                        to_node_id,
                        length,
                        lanes,
                        link_type,
                        free_speed,
                        capacity,
                        allowed_uses,
                        geometry,
                        agent_type_size,
                        demand_period_size)

            # VDF Attributes
            for i in range(demand_period_size):
                header_vdf_alpha = 'VDF_alpha' + str(i + 1)
                header_vdf_beta = 'VDF_beta' + str(i + 1)
                header_vdf_mu = 'VDF_mu' + str(i + 1)
                header_vdf_fftt = 'VDF_fftt' + str(i + 1)
                header_vdf_cap = 'VDF_cap' + str(i + 1)
                header_vdf_phf = 'VDF_phf' + str(i + 1)

                # case i: link.csv does not VDF attributes at all
                # case ii: link.csv only has partial VDF attributes
                # under case i, we will set up only one VDFPeriod ojbect using
                # default values
                # under case ii, we will set up some VDFPeriod ojbects up to
                # the number of complete set of VDF_alpha, VDF_beta, and VDF_mu
                try:
                    VDF_alpha = line[header_vdf_alpha]
                    if VDF_alpha:
                        VDF_alpha = float(VDF_alpha)
                except (KeyError, TypeError):
                    if i == 0:
                        # default value will be applied in the constructor
                        VDF_alpha = 0.15
                    else:
                        break

                try:
                    VDF_beta = line[header_vdf_beta]
                    if VDF_beta:
                        VDF_beta = float(VDF_beta)
                except (KeyError, TypeError):
                    if i == 0:
                        # default value will be applied in the constructor
                        VDF_beta = 4
                    else:
                        break

                try:
                    VDF_mu = line[header_vdf_mu]
                    if VDF_mu:
                        VDF_mu = float(VDF_mu)
                except (KeyError, TypeError):
                    if i == 0:
                        # default value will be applied in the constructor
                        VDF_mu = 1000
                    else:
                        break

                try:
                    VDF_fftt = line[header_vdf_fftt]
                    if VDF_fftt:
                        VDF_fftt = float(VDF_fftt)
                except (KeyError, TypeError):
                    # set it up using length and free_speed from link
                    VDF_fftt = length / max(SMALL_DIVISOR, free_speed) * 60

                try:
                    VDF_cap = line[header_vdf_cap]
                    if VDF_cap:
                        VDF_cap = float(VDF_cap)
                except (KeyError, TypeError):
                    # set it up using capacity from link
                    VDF_cap = capacity

                # not a mandatory column
                try:
                    VDF_phf = line[header_vdf_phf]
                    if VDF_phf:
                        VDF_phf = float(VDF_phf)
                except (KeyError, TypeError):
                    # default value will be applied in the constructor
                    VDF_phf = -1

                # construct VDFPeriod object
                vdf = VDFPeriod(i, VDF_alpha, VDF_beta, VDF_mu,
                                VDF_fftt, VDF_cap, VDF_phf)

                link.vdfperiods.append(vdf)

            # set up outgoing links and incoming links
            from_node = nodes[from_node_no]
            to_node = nodes[to_node_no]
            from_node.add_outgoing_link(link)
            to_node.add_incoming_link(link)
            links.append(link)

            # set up zone degrees
            if load_demand:
                oz_id = from_node.get_zone_id()
                dz_id = to_node.get_zone_id()
                _update_orig_zone(oz_id)
                _update_dest_zone(dz_id)

            link_seq_no += 1

        print(f"the number of links is {link_seq_no}")


def read_network(load_demand='true', input_dir='.'):
    assignm = Assignment()
    network = Network()

    read_settings(input_dir, assignm)

    read_nodes(input_dir,
               network.node_list,
               network.node_id_to_no_dict,
               network.node_no_to_id_dict,
               network.zone_to_nodes_dict)

    read_links(input_dir,
               network.link_list,
               network.node_list,
               network.node_id_to_no_dict,
               network.link_id_dict,
               assignm.get_agent_type_count(),
               assignm.get_demand_period_count(),
               load_demand)

    network.update(assignm.get_agent_type_count(),
                   assignm.get_demand_period_count())

    assignm.network = network
    assignm.setup_spnetwork()  # 作用不详

    ui = UI(assignm)

    return ui


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
                                engine_type='c', sp_algm='deque'):
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


network = read_network(load_demand=False)

print('\nshortest path (node id) from node 1 to node 2, '
      + network.find_shortest_path(1, 4))
