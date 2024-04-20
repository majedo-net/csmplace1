import numpy as np
from Node import Cell
from Multilevel import LevelGraph
from verilog2hgr import parse_func


def read_hgr(hgr_file):
    '''
    Reads the text file produced by verilog2hgr.py and returns the header and list of edges
    Each edge is also a list of cell indexes
    '''
    with open(hgr_file) as file:
        edges = [line.split() for line in file]
    header=edges.pop(0)

    for idx, edge in enumerate(edges):
        edges[idx] = [int(vertex)-1 for vertex in edge]
        
    return header, edges

if (__name__ == '__main__'):
    verilog_netlist = 'csmplace1/KSA16_yosys.vg'
    hgr_file = 'csmplace1/KSA16_yosys.hgr'
    # convert verilog netlist to hypergraph and write to file  using verilog2hgr
    parse_func(verilog_netlist)
    header,master_hg = read_hgr(hgr_file)
    num_nets = int(header[0])
    master_num_cells = int(header[1])
    master_cell_array = []
    for cidx in range(master_num_cells):
        width = np.random.rand(1)*5
        height = np.random.rand(1)*5
        master_cell_array.append(Cell(1.0+5*cidx,1, width,height,width*height))
    level0 = LevelGraph(master_hg,master_cell_array)
    del master_cell_array
    print(f'Level: {level0.current_level}')
    print(f'Cells: {level0.Nverts}')
    print(f'Edges: {len(level0.edges[level0.current_level])}')
    print(f'Avg Cell Area: {level0.average_cluster_area}')
    print('===================')
    while level0.Nverts > master_num_cells/2:
        level0.doFCCluster()
        print(f'Level: {level0.current_level}')
        print(f'Cells: {level0.Nverts}')
        print(f'Edges: {len(level0.edges[level0.current_level])}')
        print(f'Avg Cell Area: {level0.average_cluster_area}')
        print('===================')

    print('breakpoint')


