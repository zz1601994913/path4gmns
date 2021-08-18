import path4gmns as pg
from time import time

def test_find_shortest_path():
    network = pg.read_network(load_demand=False)

    print('\nshortest path (node id) from node 1 to node 2, '
          + network.find_shortest_path(1, 2))
    print('\nshortest path (link id) from node 1 to node 2, '
          + network.find_shortest_path(1, 2, seq_type='link'))

    # retrieve the shortest path under a specific mode (which must be defined
    # in settings.yaml)
    print('\nshortest path (node id) from node 1 to node 2, '
          + network.find_shortest_path(1, 2, mode='w'))
    print('\nshortest path (link id) from node 1 to node 2, '
          + network.find_shortest_path(1, 2, mode='w', seq_type='link'))


def test_accessibility():
    network = pg.read_network(load_demand=False)

    print('\nstart accessibility evaluation\n')
    st = time()

    # multimodal accessibility evaluation
    pg.evaluate_accessibility(network)
    # accessibility evalutation for a target mode
    # pg.evaluate_accessibility(network, multimodal=False, mode='p')

    print('complete accessibility evaluation.\n')
    print(f'processing time of accessibility evaluation: {time()-st:.2f} s')

    # get accessible nodes and links starting from node 1 with a 5-minitue
    # time window for the default mode auto (i.e., 'p')
    network.get_accessible_nodes(1, 5)
    network.get_accessible_links(1, 5)

    # get accessible nodes and links starting from node 1 with a 15-minitue
    # time window for mode walk (i.e., 'w')
    network.get_accessible_nodes(1, 15, 'w')
    network.get_accessible_links(1, 15, 'w')

    # time-dependent accessibility under the default mode auto (i.e., p)
    # for demand period 0 (i.e., VDF_fftt1 in link.csv will be used in the
    # evaluation)
    # pg.evaluate_accessibility(network, multimodal=False, time_dependent=True)

    # it is equivalent to
    # pg.evaluate_accessibility(network, multimodal=False,
    #                           time_dependent=True, demand_period_id=0)

    # get accessible nodes and links starting from node 1 with a 5-minitue
    # time window for the default mode auto (i.e., 'p') for demand period 0
    # network.get_accessible_nodes(1, 5, time_dependent=True)

    # get accessible nodes and links starting from node 1 with a 15-minitue
    # time window for mode walk (i.e., 'w') for demand period 0
    # network.get_accessible_nodes(1, 15, 'w', time_dependent=True)
def demo_mode(mode):
    print(f'the selected mode is {mode}\n')


    if mode == 1:
        # option 1: find shortest path between O and D on Chicago network
        test_find_shortest_path()
    else:
        # option 6: evaluate multimodal accessibility on Chicago network
        test_accessibility()


if __name__=="__main__":

    demo_mode(6)