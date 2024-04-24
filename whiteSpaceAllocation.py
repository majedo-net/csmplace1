from Node import Node, Cell, Partition, getEdges
import numpy as np

"""
def partitionAux(top, minArea, cut_type):
    Recursive auxiliary function for partition()

    Parameters:
        top: Top node in binary tree
        minArea: Minimum allowable area of partition
        cut_type: 'V' for vertical or 'H' for horizontal

    Return:
        top node with child node pointers set
    if top.area >= minArea:
        left = Node()
        right = Node()

        #Alternate vertical and horizontal cuts
        if cut_type == 'V':
            if top.w % 2 == 0:#If width of partition being divided is odd
                new_w = top.w // 2
                left.tlx = top.tlx
                left.tly = top.tly
                left.w = new_w
                left.h = top.h
                left.area = left.w*left.h
                right.tlx = top.tlx + new_w
                right.tly = top.tly
                right.w = new_w
                right.h = top.h
                right.area = right.w*right.h
            else:#If width of partition being divided is odd
                int_new_w = int(top.w / 2)
                left.tlx = top.tlx
                left.tly = top.tly
                left.w = int_new_w
                left.h = top.h
                left.area = left.w*left.h
                right.tlx = top.tlx + int_new_w
                right.tly = top.tly
                right.w = int_new_w + 1
                right.h = top.h
                right.area = right.w*right.h

            top.left = partitionAux(left, minArea, 'H')#Alternate vertical and horizontal cuts
            top.right = partitionAux(right, minArea, 'H')

        if cut_type == 'H':
            if top.h % 2 == 0:#If height of partition being divided is odd
                new_h = top.h // 2
                left.tlx = top.tlx
                left.tly = top.tly
                left.w = top.w
                left.h = new_h
                left.area = left.w*left.h
                right.tlx = top.tlx 
                right.tly = top.tly + new_h
                right.w = top.w
                right.h = new_h
                right.area = right.w*right.h
            else:#If height of partition being divided is odd
                int_new_h = int(top.h / 2)
                left.tlx = top.tlx
                left.tly = top.tly
                left.w = top.w
                left.h = int_new_h
                left.area = left.w*left.h
                right.tlx = top.tlx 
                right.tly = top.tly + int_new_h
                right.w = top.w
                right.h = int_new_h + 1
                right.area = right.w*right.h

            top.left = partitionAux(left, minArea, 'V')#Alternate horizontal and vertical cuts
            top.right = partitionAux(right, minArea, 'V')

        return top
    else:
        return None
"""

def main():
    #tn = partition(0, 0, 24, 20, 30)
    """
    cells_arr = [Cell(2.0,2.0,np.sqrt(14.0),np.sqrt(14.0),14.0), Cell(5.0,4.0,np.sqrt(13.0),np.sqrt(13.0),13.0), Cell(8.0,2.0,np.sqrt(20.0),np.sqrt(20.0),20.0), Cell(11.0,4.5,np.sqrt(11.0),np.sqrt(11.0),11.0), 
                 Cell(5.0, 6.0, np.sqrt(24.0), np.sqrt(24.0), 24.0), Cell(8.0, 9.0, 4.0, 4.0, 16.0), Cell(11.0, 7.0, 4.0, 4.0, 16.0), Cell(12.5,0.5,np.sqrt(5.0), np.sqrt(5.0), 5.0), Cell(16.0, 4.0, np.sqrt(20.0), np.sqrt(20.0), 20.0),
                 Cell(16.0, 1.0, np.sqrt(10.0), np.sqrt(10.0), 10.0), Cell(20.0, 3.0, 4.0, 4.0, 16.0), Cell(23.0, 1.0, np.sqrt(13.0), np.sqrt(13.0), 13.0), Cell(13.5, 8.8, np.sqrt(15.0), np.sqrt(15.0), 15.0),
                 Cell(17.0, 7.0, np.sqrt(15.0), np.sqrt(15.0), 15.0), Cell(21.0, 6.0, np.sqrt(17.0), np.sqrt(17.0), 17.0), Cell(21.0, 9.0, np.sqrt(12.0), np.sqrt(12.0), 12.0)]
    part = Partition(30.0, cells_arr)
    part.construct(0.0, 0.0, 24.0, 10.0)
    ws = part.calcWhiteSpace(part.topNode)
    part.WSA(part.topNode)
    part.printLeaves(part.topNode)
    part.updateCells(part.topNode)
    for cell in part.cells:
        cell.printSelf()
    """

    cells_arr = [Cell(2.0,2.0,np.sqrt(14.0),np.sqrt(14.0),14.0,[1,2,3]), Cell(5.0,4.0,np.sqrt(13.0),np.sqrt(13.0),13.0,[0,3,4]), Cell(8.0,2.0,np.sqrt(20.0),np.sqrt(20.0),20.0,[0,3,5]), 
                 Cell(11.0,4.5,np.sqrt(11.0),np.sqrt(11.0),11.0,[0,1,2,4,5,6]), Cell(5.0, 6.0, np.sqrt(24.0), np.sqrt(24.0), 24.0,[1,3,6]), Cell(8.0, 9.0, 4.0, 4.0, 16.0,[2,3,6]), 
                 Cell(11.0, 7.0, 4.0, 4.0, 16.0,[3,4,5])]
    edges =  getEdges(cells_arr)
    print(edges)
                 

if __name__ == "__main__":
    main()