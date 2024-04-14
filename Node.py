class Node:
    """
    Class declaration for node in slicing tree for white space
    allocation algorithm in CSMPlace1

    Nodes are rectangles in R^2
    """

    #top-left x and y coordinates of rectangle
    tlx = 0
    tly = 0

    #width, height, and area
    w = 0
    h = 0
    area = 0

    #White space
    whiteSpace = 0

    #Child Nodes
    right = None
    left = None

    #List of cell numbers that lie inside this Node
    cell_nums = None

    def __init__(self, tl_x_=0, tl_y_=0, w_=0, h_=0, area_=0, cells_=None):
        self.tlx = tl_x_
        self.tly = tl_y_
        self.w = w_
        self.h = h_
        self.area = area_
        self.whiteSpace = 0
        self.right = None
        self.left = None
        self.cell_nums = cells_


    def printSelf(self):
        print("Top-left coordinates: x: " + str(self.tlx) + "," + str(self.tly))
        print("Width: " + str(self.w) + "   Height: " + str(self.h))
        print("Area: " + str(self.area) + "   White Space Available:  " + str(self.whiteSpace))

    def doesCellIntersect(self, cell):
        """
        Function to determine if the specified cell intersects
        this partition
        """

        ctlx = cell.cx - (cell.w/2)
        ctly = cell.cy - (cell.h/2)
        cbrx = cell.cx + (cell.w/2)
        cbry = cell.cy + (cell.h/2)

        pbrx = self.tlx + self.w
        pbry = self.tly + self.h

        if (self.tlx < cbrx) and (pbrx > ctlx) and (self.tly < cbry) and (pbry > ctly):
            return True
        
        return False
    
    def isCellInside(self, cell):
        if (cell.cx >= self.tlx) and (cell.cx < (self.tlx + self.w)) and (cell.cy >= self.tly) and (cell.cy < (self.tly + self.h)):
            return True
        return False

class Cell:
    cx = 0.0
    cy = 0.0
    w = 0.0
    h = 0.0
    area = 0.0

    def __init__(self, cx_=0.0, cy_=0, w_=0, h_=0, area_=0.0):
        self.cx = cx_
        self.cy = cy_
        self.w = w_
        self.h = h_
        self.area = area_

    def printSelf(self):
        print("Top-left coordinates: x: " + str(self.cx) + "," + str(self.cy))
        print("Width: " + str(self.w) + "   Height: " + str(self.h))
        print("Area: " + str(self.area))

class Partition:
    """
    Class to represent the partition of a rectangular layout region into 
    smaller rectangles that have a minimum area
    """
    cells = None
    topNode = None
    minArea = 0.0

    def __init__(self, minA, cells_):
        """
        Constructor

        Parameters:
            minA: Minimum allowable area of each partition
            cells_: List of all cells that are iinside the overall rectangular region
        """
        self.cells = cells_
        self.minArea = minA

    def construct(self, ovr_tlx, ovr_tly, ovr_w, ovr_h):
        """
        Function to construct the slicing tree of the partition, including the calculation
        of available white space in each rectangle of the slicing tree

        Parameters:
            ovr_tlx: Top-left x-coordinate of overall rectangular region
            ovr_tly: Top-left y-coordinate of overall rectangular region
            ovr_w: Width of overall rectangular region
            ovr_h: Height of overall rectangular region
        """
        cell_nos = []
        for i in range(len(self.cells)):
            cell_nos.append(i)
        top = Node(ovr_tlx, ovr_tly, ovr_w, ovr_h, ovr_w*ovr_h, cell_nos)
        self.topNode = self.partitionAux(top, 'V')

    def partitionAux(self, top, cut_type):
        """
        Recursive auxiliary function for constructor that builds
        the slicing tree of the partition

        Parameters:
            top: Top node in binary tree
            cut_type: 'V' for vertical or 'H' for horizontal

        Return:
            top node of slicing tree with child node pointers set
        """
        if top.area >= self.minArea:
            left = Node()
            right = Node()

            #Alternate vertical and horizontal cuts
            if cut_type == 'V':
                if top.w % 2 == 0:#If width of partition being divided is even
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

            if cut_type == 'H':
                if top.h % 2 == 0:#If height of partition being divided is even
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

            #Determine which cells from the top are inside the left and right partitions
            left_cell_nums = []
            right_cell_nums = []
            for cell_num in top.cell_nums:
                if left.isCellInside(self.cells[cell_num]):
                    left_cell_nums.append(cell_num)
                if right.isCellInside(self.cells[cell_num]):
                    right_cell_nums.append(cell_num)
            left.cell_nums = left_cell_nums
            right.cell_nums = right_cell_nums

            new_cut_type = 'H'
            if cut_type == 'H':
                new_cut_type = 'V'
            top.left = self.partitionAux(left, new_cut_type)#Alternate vertical and horizontal cuts
            top.right = self.partitionAux(right, new_cut_type)
            return top
        
        else:
            return None
        
    def calcWhiteSpace(self, top):
        """
        Function to calculate the white space available in each node in the
        slicing tree in a bottom-up manner (post-order traversal)

        Parameters:
            top: topNode data member in Partition class (so we don't need
                 a recursive auxiliary function)

        Return:
            None, update to slicing tree is performed in-place
        """
        left_ws = 0.0
        right_ws = 0.0
        if top.left is not None:
            left_ws = self.calcWhiteSpace(top.left)
        if top.right is not None:
            right_ws = self.calcWhiteSpace(top.right)
        
        #If we are at a leaf node
        if (top.left is None) and (top.right is None):
            cell_area = 0.0
            for cell_n in top.cell_nums:
                cell_area += self.cells[cell_n].area
            top.whiteSpace = top.area - cell_area
        else:#We are not at a child node
            top.whiteSpace = left_ws + right_ws

        return top.whiteSpace
    
    def WSA(self, top):
        """
        Function to allocate available white space to overflow partitons
        in a top-down manner and update cut lines accordingly

        Parameters:
            top: topNode data member in Partition class (so we don't need
                 a recursive auxiliary function)

        Return:
            None, update to slicing tree is performed in-place
        """
        if top.left is not None and top.right is not None:
            cut_type = 'H'#Horizontal or vertical cut type for cut line adjustment
            if top.left.tlx != top.right.tlx:#Cut is vertical, not horizontal
                cut_type = 'V'
            #Allocate all white space available from top to left child
            if (top.left.whiteSpace < 0) and (top.right.whiteSpace >= 0):
                olws = top.left.whiteSpace
                orws = top.right.whiteSpace
                top.left.whiteSpace = 0.0
                top.right.whiteSpace = top.whiteSpace
                top.left.area = top.left.area + top.left.whiteSpace - olws
                top.right.area = top.right.area + top.right.whiteSpace - orws
            #Allocate all white space available from top to right child
            elif (top.left.whiteSpace >= 0) and (top.right.whiteSpace < 0):
                olws = top.left.whiteSpace
                orws = top.right.whiteSpace
                top.left.whiteSpace = top.whiteSpace
                top.right.whiteSpace = 0.0
                top.left.area = top.left.area + top.left.whiteSpace - olws
                top.right.area = top.right.area + top.right.whiteSpace - orws
            #Both child nodes have available white space --> distribute
            #available white space from top proportionally
            elif (top.left.whiteSpace >= 0) and (top.right.whiteSpace >= 0):
                olws = top.left.whiteSpace
                orws = top.right.whiteSpace
                if top.right.whiteSpace > 0:
                    r = top.left.whiteSpace/top.right.whiteSpace#ratio of white space available in child nodes
                    top.left.whiteSpace = (r*top.whiteSpace) / (r + 1.0)
                    top.right.whiteSpace = top.whiteSpace / (r + 1.0)
                else:
                    top.left.whiteSpace = 0
                top.left.area = top.left.area + top.left.whiteSpace - olws
                top.right.area = top.right.area + top.right.whiteSpace - orws
            
            #Update cut lines now that we know area adjustment
            if cut_type == 'V':
                top.left.tlx = top.tlx
                top.left.tly = top.tly
                top.left.h = top.h
                top.left.w = top.left.area / top.left.h
                top.right.h = top.h
                top.right.w = top.right.area / top.right.h
                top.right.tlx = top.left.tlx + top.left.w
                top.right.tly = top.left.tly 
            else:
                top.left.tlx = top.tlx
                top.left.tly = top.tly
                top.left.w = top.w
                top.left.h = top.left.area / top.left.w
                top.right.w = top.w
                top.right.h = top.right.area / top.right.w
                top.right.tly = top.left.tly + top.left.h
                top.right.tlx = top.left.tlx 

            #Continue pre-order traversal
            self.WSA(top.left)
            self.WSA(top.right)

        else:
            return

    def printLeaves(self, top):
        """
        Function to print leaf nodes in slicing tree
        """
        if top is None:
            return
    
        if (top.left is None) and (top.right is None):
            top.printSelf()

        self.printLeaves(top.left)
        self.printLeaves(top.right)

    def inOrderTraversal(self, top):
        """
        Function to print nodes in slicing tree in order
        """
        if top is None:
            return
        
        self.inOrderTraversal(top.left)
        top.printSelf()
        self.inOrderTraversal(top.right)




        