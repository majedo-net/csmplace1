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
    verilog_netlist = 'csmplace1/test.vg'
    hgr_file = 'csmplace1/test.hgr'
    # convert verilog netlist to hypergraph using verilog2hgr
    parse_func(verilog_netlist)
    header,edges = read_hgr(hgr_file)
    num_nets = int(header[0])
    num_cells = int(header[1])
    cell_array = []
    for cidx in range(num_cells):
        cell_array.append(Cell(1.0+5*cidx,1, 4, 4, 16))
    level0 = LevelGraph(cell_array,edges)
    level0.findNeighbors()
    level0.calcAffinity()
    level0.sortDegree()
    level0.cluster()

    print('breakpoint')


