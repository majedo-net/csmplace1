import numpy as np
from Node import Cell
from Multilevel import LevelGraph
from verilog2hgr import parse_func
from spef_extractor.lef_def_parser import LefParser, DefParser


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
        width = 0
        height = 0
        master_cell_array.append(Cell(1.0+5*cidx,1, width,height,width*height)) 

    # parse def file for cell information
    def_file = 'csmplace1/KSA16_yosys.def'
    defparser = DefParser(def_file)
    defparser.parse()
    # build the indexed list of cell types (macros) used in the design
    ctypes = []
    for comp in defparser.components:
        if 'tap' in comp.macro: # tap cells not placed
            continue
        else:
            ctypes.append(comp.macro.replace('hd','hs'))
    if len(ctypes) != len(master_cell_array):
        print('Something is wrong with macro list')

    for cidx in range(master_num_cells):
        lef_file = f'csmplace1/cell_lefs/{ctypes[cidx]}.lef'
        lefparser = LefParser(lef_file)
        lefparser.parse()
        master_cell_array[cidx].w = lefparser.macro_dict[ctypes[cidx]].info['SIZE'][0]
        master_cell_array[cidx].h = lefparser.macro_dict[ctypes[cidx]].info['SIZE'][1]
        master_cell_array[cidx].area = master_cell_array[cidx].w * master_cell_array[cidx].h

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


