#!/usr/bin/python

"""
*   Copyright (C) 2019 Keegan Yao, Ashim Ranjeet and Mukul S. Bansal (mukul.bansal@uconn.edu).
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""



import argparse
from ete3 import Tree
import math
import itertools
import csv
import time



class LCAMapping:
    def __init__(self, Tree):
        self.tree = Tree
        self.array = []
        self.lcaToRMQ()
        self.preprocessST()


    # RMQ preprocessing using dynamic programming - O(n^2)
    def preprocessDP(self):
        n = len(self.levelLst)
        self.array = []
        for i in range(n):
            L1 = []
            for j in range(n):
                if(j < i):
                    L1.append('-')
                elif(j == i):
                    L1.append(i)
                else:
                    if(self.levelLst[L1[j-1]] > self.levelLst[j]):
                        L1.append(j)
                    else:
                        L1.append(L1[j-1])
            self.array.append(L1)
        #print(array)


    # RMQ preprocessing using Stable Table - O(nlogn)
    def preprocessST(self):
        self.array = []
        x = int(math.log2(len(self.levelLst)))
        for i in range(len(self.levelLst)):
            self.array.append([])
        for j in range(x+1):
            for i in range(len(self.levelLst)):
                if(j == 0):
                    self.array[i].append(i)
                else:
                    if(i+ 2**(j-1) >= len(self.levelLst)):
                        self.array[i].append('-')
                    elif(self.array[i+ 2**(j-1)][j-1] == '-'):
                        self.array[i].append('-')
                    elif(self.levelLst[self.array[i][j-1]] <= self.levelLst [ self.array[i+ 2**(j-1)][j-1] ]):
                        self.array[i].append(self.array[i][j-1])
                    else:
                        self.array[i].append(self.array[i+ 2**(j-1)][j-1])

        #return array





    # query for stable table prepreprossing
    def queryST(self, i, j):
        if(i == j ):
            if((i >= 0) and (i < len(self.levelLst)) ):
                return i
            else:
                print('out of bounds')
                return

        k = int(math.log2(j - i))
        x = self.array[i][k]
        y = self.array[j-2**(k)+1][k]
        if(self.levelLst[x] <= self.levelLst[y]):
            value = x
        else:
            value = y
        return value


    def eulerTour(self, node,lv):
        self.eTourLst.append(node)
        self.repLst[node.name] = len(self.eTourLst)-1
        self.levelLst.append(lv)
        if(node.is_leaf()):
            return True
        for nd in node.get_children():
            self.eulerTour(nd,lv+1)
            self.eTourLst.append(node)
            self.levelLst.append(lv)


    # convert LCA into RMQ
    def lcaToRMQ(self):
        self.eTourLst = []
        self.levelLst = []
        self.repLst = {}
        self.eulerTour(self.tree, 0)

    # query for LCA
    def queryLCA(self, node1, node2):
        x = self.repLst.get(node1)
        y = self.repLst.get(node2)
        if x is None:
            print("NONETYPE ERROR")
            print(node1)
        if(x == y):
            indexL = x
        elif(x > y):
            indexL = self.queryST(y, x)
        elif(x < y):
            indexL = self.queryST(x, y)
        return self.eTourLst[indexL]





class compareTree:
    def __init__(self, T, S):
        self.colors = ["Green", "Red", "Blue", "Yellow"]

        # Store the two completed trees
        self.t = T.copy()
        self.t2 = S.copy()

        # Store the EF-RF(+) completions, in the order that they are inputted to the initializer
        self.EF1 = None
        self.EF2 = None

        # Store the RF(+) completions after they have been computed, as a check that the function has been run
        self.RF1 = None
        self.RF2 = None

        # Keep track of the color of each vertex
        self.tcolors = {}
        self.t2colors = {}

        # Store which leaves are colored yellow and which leaves are colored red
        self.yellowLeaves, self.redLeaves = set(), set()

        # Store half of the best possible change to the RF distance given that we push up any number of maximal
        # red and yellow subtrees, pair up whenever possible, and leave N maximal subtrees of color c unpaired
        # at the input vertex v
        self.cost = {}

        # Store the maximum number of red (c=0) and yellow (c=1) subtrees contained within the subtree
        # rooted at each vertex in self.t (self.cMax1) and self.t2 (self.cMax2)
        self.cMax1, self.cMax2 = {}, {}

        # Keep track of parameter values at each of the children which produce the optimal cost value
        # for the parent vertex given fixed parameter values N and c
        self.leftN1, self.rightN1 = {}, {}
        self.leftc1, self.rightc1 = {}, {}

        # Keep track of optimal parameters
        self.optN1, self.optc1 = {}, {}
        self.order1, self.order2 = {}, {}

        self.iNode1, self.iNode2 = 0, 0     # if internal nodes are unnammed, name them based on numbers
        self.swapped = 0                    # to keep track of if t and t2 is T and S or if they are swapped

        self.init_swap, self.is_init = 0, None      # to keep track of the original order the two trees were inputted
        self.resetT()

        self.start_time = time.time()
        self.EF_time, self.RF_time = None, None



    def resetT(self):
        # swapping to set the smaller tree to be t, and larger tree to be t2
        # a map of the name to node based on dictionary is made.

        self.map = {}
        self.t1Leafset = {}
        self.t2Leafset = {}
        t1LeafCount = 0
        t2LeafCount = 0

        self.iNode1 = 0
        for node in self.t.traverse("postorder"):
            # if internal node doesn't have a name, name it inode, increment inode
            if(not node.is_leaf()):
                self.iNode1 = self.iNode1 + 1
                node.name = str(self.iNode1)
            else:
                t1LeafCount +=1
                self.t1Leafset[str(node.name)] = t1LeafCount

        self.iNode2 = 0
        for node in self.t2.traverse("postorder"):
            if(not node.is_leaf()):
                self.iNode2 = self.iNode2 + 1
                node.name = str(self.iNode2)
            else:
                t2LeafCount +=1
                self.t2Leafset[str(node.name)] = t2LeafCount

        if(t1LeafCount  > t2LeafCount):
            temp = self.t
            self.t = self.t2
            self.t2 = temp
            self.tLeafset = self.t2Leafset
            self.swapped = self.swapped + 1

            self.tcolors, self.t2colors = self.t2colors, self.tcolors
            self.iNode1, self.iNode2 = self.iNode2, self.iNode1

            self.cMax1, self.cMax2 = self.cMax2, self.cMax1

            if self.is_init is None:
                self.init_swap = 1
        else:
            self.tLeafset = self.t1Leafset

        self.is_init = 0
        self.t1mapping = LCAMapping(self.t)



    def recolor(self):
        # Method used to make sure that all colors are aligned properly before running the DP recurrence relation
        # In the process of the EF-RF completions, it is necessary to treat some nodes which should be red / yellow
        #     as green nodes to prevent confusion over which subtrees should be duplicated and added at each stage
        # Note that the yellow and red leaves are stored precisely to keep track of which subtrees should be which colors
        #     after the EF-RF(+) completions have been computed

        for node in self.t.traverse("postorder"):
            if node.is_leaf():
                if node.name in self.redLeaves:
                    self.tcolors[node.name] = self.colors[1]
                elif node.name in self.yellowLeaves:
                    self.tcolors[node.name] = self.colors[3]
                else:
                    self.tcolors[node.name] = self.colors[0]
            else:
                numGRB = [0,0,0,0]    # to count the number of green, red, blue and yellow children of the node
                # count the number of child depending on color
                for child in node.children:
                    if(self.tcolors[child.name] == "Green"):
                        numGRB[0] = numGRB[0] + 1
                    elif(self.tcolors[child.name] == "Red"):
                        numGRB[1] = numGRB[1] + 1
                    elif(self.tcolors[child.name] == "Blue"):
                        numGRB[2] = numGRB[2] + 1
                    elif(self.tcolors[child.name] == "Yellow"):
                        numGRB[3] = numGRB[3] + 1
                # if both child are green, make color of node green.
                if(numGRB[0] == 2):
                    self.tcolors[node.name] = self.colors[0]
                    node.add_features(mark = False)
                # if both child are red, then make color of node red
                elif(numGRB[1] == 2):
                    self.tcolors[node.name] = self.colors[1]
                    node.add_features(mark = False)
                # if both child are yellow, then make color of node yellow
                elif(numGRB[3] == 2):
                    self.tcolors[node.name] = self.colors[3]
                    node.add_features(mark = False)
                # if one child is red (or yellow), and the other is green or blue, color node blue and mark node
                elif((numGRB[1] == 1 or numGRB[3] == 1) and (numGRB[0] == 1 or numGRB[2] == 1)):
                    self.tcolors[node.name] = self.colors[2]
                    node.add_features(mark = True)

                # otherwise make node blue, unmarked
                else:
                    self.tcolors[node.name] = self.colors[2]
                    node.add_features(mark = False)

        for node in self.t2.traverse("postorder"):
            if node.is_leaf():
                if node.name in self.redLeaves:
                    self.t2colors[node.name] = self.colors[1]
                elif node.name in self.yellowLeaves:
                    self.t2colors[node.name] = self.colors[3]
                else:
                    self.t2colors[node.name] = self.colors[0]
            else:
                numGRB = [0,0,0,0]    # to count the number of green, red, blue and yellow children of the node
                # count the number of child depending on color
                for child in node.children:
                    if(self.t2colors[child.name] == "Green"):
                        numGRB[0] = numGRB[0] + 1
                    elif(self.t2colors[child.name] == "Red"):
                        numGRB[1] = numGRB[1] + 1
                    elif(self.t2colors[child.name] == "Blue"):
                        numGRB[2] = numGRB[2] + 1
                    elif(self.t2colors[child.name] == "Yellow"):
                        numGRB[3] = numGRB[3] + 1
                # if both child are green, make color of node green.
                if(numGRB[0] == 2):
                    self.t2colors[node.name] = self.colors[0]
                    node.add_features(mark = False)
                # if both child are red, then make color of node red
                elif(numGRB[1] == 2):
                    self.t2colors[node.name] = self.colors[1]
                    node.add_features(mark = False)
                # if both child are yellow, then make color of node yellow
                elif(numGRB[3] == 2):
                    self.t2colors[node.name] = self.colors[3]
                    node.add_features(mark = False)
                # if one child is red (or yellow), and the other is green or blue, color node blue and mark node
                elif((numGRB[1] == 1 or numGRB[3] == 1) and (numGRB[0] == 1 or numGRB[2] == 1)):
                    self.t2colors[node.name] = self.colors[2]
                    node.add_features(mark = True)
                # otherwise make node blue, unmarked
                else:
                    self.t2colors[node.name] = self.colors[2]
                    node.add_features(mark = False)



    def ROT_RF_Plus(self):
        #tree.traverse("postorder") is ETE function which tranverse the tree in postorder

        for node in self.t2.traverse("postorder"):
            if node.is_leaf():     # coloring each leaf red or green
                if(self.tLeafset.get(str(node.name)) != None):
                    node.add_features(mark = False)
                    self.t2colors[node.name] = self.colors[0]
                else:
                    if self.swapped % 2 == 1:
                        self.t2colors[node.name] = self.colors[3]
                        self.yellowLeaves.add(node.name)
                    else:
                        self.t2colors[node.name] = self.colors[1]
                        self.redLeaves.add(node.name)
                    node.add_features(mark = False)
            else:   # coloring the internal nodes
                numGRB = [0,0,0,0]    # to count the number of green, red, blue and yellow children of the node
                # count the number of child depending on color
                for child in node.children:
                    if(self.t2colors[child.name] == "Green"):
                        numGRB[0] = numGRB[0] + 1
                    elif(self.t2colors[child.name] == "Red"):
                        numGRB[1] = numGRB[1] + 1
                    elif(self.t2colors[child.name] == "Blue"):
                        numGRB[2] = numGRB[2] + 1
                    elif(self.t2colors[child.name] == "Yellow"):
                        numGRB[3] = numGRB[3] + 1
                # if both child are green, make color of node green.
                if(numGRB[0] == 2):
                    self.t2colors[node.name] = self.colors[0]
                    node.add_features(mark = False)
                # if both child are red, then make color of node red
                elif(numGRB[1] == 2):
                    self.t2colors[node.name] = self.colors[1]
                    node.add_features(mark = False)
                # if both child are yellow, then make color of node yellow
                elif(numGRB[3] == 2):
                    self.t2colors[node.name] = self.colors[3]
                    node.add_features(mark = False)
                # if one child is red (or yellow), and the other is green or blue, color node blue and mark node
                elif((numGRB[1] == 1 or numGRB[3] == 1) and (numGRB[0] == 1 or numGRB[2] == 1)):
                    self.t2colors[node.name] = self.colors[2]
                    node.add_features(mark = True)
                # otherwise make node blue, unmarked
                else:
                    self.t2colors[node.name] = self.colors[2]
                    node.add_features(mark = False)

        # If there are NO green leaves, then there are also no blue leaves, and hence we cannot just follow the algorithm.
        #     (In this case, there is also no possible EF-RF completion)
        # A pair of trees in this form is technically not an instance of the ROT-RF(+) problem. However, following through
        #     with this case in the ROT_RF(+) method will be convenient for computing the optimal R-RF(+) completions
        if self.t2colors[self.t2.name] == "Red" or self.t2colors[self.t2.name] == "Yellow":
            newT = Tree()
            newT.add_child(self.t)
            z = newT.add_child(self.t2.copy())
            if not z.is_leaf():
                self.iNode1 = self.iNode1 + 1
                z.name = str(self.iNode1)
            self.iNode1 = self.iNode1 + 1
            newT.name = str(self.iNode1)
            self.t = newT
            self.tcolors[z.name] = self.t2colors[self.t2.name]
            self.tcolors[newT.name] = self.colors[2]
            return

            newT2 = Tree()
            newT2.add_child(self.t2)
            z = newT2.add_child(self.t.copy())
            if not z.is_leaf():
                self.iNode2 = self.iNode2 + 1
                z.name = str(self.iNode2)
            self.iNode2 = self.iNode2 + 1
            newT2.name = str(self.iNode2)
            self.t2 = newT2
            self.t2colors[z.name] = self.tcolors[self.t2.name]
            self.t2colors[newT2.name] = self.colors[2]
            self.EF_exists = False
            return


        # If at least one leaf is green, then the optimal RF completion can be formed by following the algorithm
        self.EF_exists = True

        #ETE function to tranverse the tree in postorder (the second tree in this case)
        for node in self.t2.traverse("postorder"):
            # mapping
            # when the node is green
            #print(self.t2colors[node.name], node, "\n\n")
            if(self.t2colors[node.name] == self.colors[0]):
                if(node.is_leaf()):
                    # in the dictionary, map the name of the node to the node in the first tree

                    self.map[node.name] = self.t1mapping.queryLCA(node.name, node.name)
                else:
                    # get a list of all the children, and find the common ancestor of the children
                    lst = []
                    for child in node.children:
                        #print(child)
                        lst.append(self.map.get(child.name))
                    #print("Here is the node and list:\t", node, lst, "\n\n")

                    # ETE function to get common ancestor of a list of leaves or nodes
                    #nd = self.t.get_common_ancestor(lst)
					# now geting the common ancestor using the LCA mapping
                    nd = self.t1mapping.queryLCA(lst[0].name, lst[1].name)

                    # map the name of the node to the ancestor node in the first tree
                    self.map[node.name] = nd

            # when the node is blue
            elif(self.t2colors[node.name] == self.colors[2]):
                lst = []
                # and the children that are green or blue to the list
                for child in node.children:
                    if(self.t2colors[child.name] == self.colors[0] or self.t2colors[child.name] == self.colors[2]):
                        lst.append(self.map.get(child.name))

                # when there is only one node in the list, use that node is the common ancestor
                # otherwise get the greatest ancestor of the list of nodes in the first

                if(len(lst) == 1):
                    if(lst[0].is_leaf()):
                        nd = lst[0]
                    else:
                        #nd = self.t.get_common_ancestor(lst[0].children)
                        nd = self.t1mapping.queryLCA(lst[0].name, lst[0].name)
                    self.map[node.name] = nd
                elif(len(lst) == 2):
                    nd = self.t1mapping.queryLCA(lst[0].name,lst[1].name)
                    self.map[node.name] = nd
                else:
                    print('non binary tree')
        #ETE function to tranverse tree in preorder
        for node in self.t2.traverse("preorder"):
            #tree-add
            # if the node is marked, meaning that one child is red
            if(node.mark == True):
                for child in node.children:
                    # when the child is red, add it that node to the first tree.
                    # if the node is the root, then the node is a sibling, to the current tree,
                    # and you have to create a new root
                    # otherwise, you have to detach the maping of the node in the first tree, and
                    # create a new tree with the detach node and the node you want to add as sibling.
                    # this new tree is then added to location where the node was detached in the first tree.
                    if(self.t2colors[child.name] == self.colors[1] or self.t2colors[child.name] == self.colors[3]):
                        nd = self.map.get(node.name)
                        # Case: red/yellow subtree attached at root
                        if(nd.up == None):
                            newT = Tree()
                            newT.add_child(self.t)
                            z = newT.add_child(child.copy())
                            if not z.is_leaf():
                                self.iNode1 = self.iNode1 + 1
                                z.name = str(self.iNode1)
                            self.iNode1 = self.iNode1 + 1
                            newT.name = str(self.iNode1)
                            self.t = newT
                            self.tcolors[z.name] = self.t2colors[child.name]
                            self.tcolors[newT.name] = self.colors[2]
                        #Case: red/yellow subtree attached at different internal node
                        else:
                            up = nd.up
                            newT = Tree()
                            removed = nd.detach()
                            z = newT.add_child(child.copy())
                            if not z.is_leaf():
                                self.iNode1 = self.iNode1 + 1
                                z.name = str(self.iNode1)
                            self.iNode1 = self.iNode1 + 1
                            newT.name = str(self.iNode1)

                            newT.add_child(removed)
                            up.add_child(newT)

                            self.tcolors[z.name] = self.t2colors[child.name]
                            self.tcolors[newT.name] = self.colors[2]



    def EF_R_RF(self):

        # If the EF-RF completions have been computed before, then they are stored in the compareTree class
        if self.EF1 and self.EF2:
            return self.EF1, self.EF2
        # If not, then continue with computing the EF-RF completions

        self.start_time = time.time()

        # calling the one tree completion function.
        # if there are missing leaves in both tree, then it is called twice
        # otherwise, it is only called once
        # it returns the two completed trees
        self.ROT_RF_Plus()


        # The EF_exists variable stores whether such a completion can be formed at all
        # If no leaves are shared between the two trees, then there cannot be an EF-RF completion.
        #     Instead, we will return the optimal completion, which is one large extraneous clade.

        self.resetT()

        self.ROT_RF_Plus()


        if(self.swapped == 0 or self.swapped == 2):
            self.EF1, self.EF2 = self.t.copy(), self.t2.copy()
        else:
            self.EF1, self.EF2 = self.t2.copy(), self.t.copy()

        self.EF_time = time.time() - self.start_time
        return self.EF1, self.EF2



    def setRedMax(self):
        # Simple postorder traversal counting how many maximal red subtrees are contained
        #     in the subtree rooted at each vertex
        for v in self.t.traverse("postorder"):
            if v.is_leaf():
                if self.tcolors[v.name] == "Red" and self.tcolors[v.up.name] != "Red":
                    self.cMax1[0, v.name] = 1
                else:
                    self.cMax1[0, v.name] = 0
            else:
                tempMax, redChildren = 0, 0
                for child in v.children:
                    if not self.tcolors[child.name] == "Red":
                        tempMax = tempMax + self.cMax1[0, child.name]
                    else:
                        redChildren = redChildren + 1
                if redChildren == 1:
                    tempMax = tempMax + 1
                elif redChildren == 2:
                    if self.tcolors[v.name] == "Red" and self.tcolors[v.up.name] != "Red":
                        tempMax = 1
                    elif self.tcolors[v.name] != "Red":
                        tempMax = 2
                self.cMax1[0, v.name] = tempMax



    def setYellowMax(self):
        # Simple postorder traversal counting how many maximal yellow subtrees are contained
        #     in the subtree rooted at each vertex
        for v in self.t.traverse("postorder"):
            if v.is_leaf():
                if self.tcolors[v.name] == "Yellow" and self.tcolors[v.up.name] != "Yellow":
                    self.cMax1[1, v.name] = 1
                else:
                    self.cMax1[1, v.name] = 0
            else:
                tempMax, yellowChildren = 0, 0
                for child in v.children:
                    if not self.tcolors[child.name] == "Yellow":
                        tempMax = tempMax + self.cMax1[1, child.name]
                    else:
                        yellowChildren = yellowChildren + 1
                if yellowChildren == 1:
                    tempMax = tempMax + 1
                elif yellowChildren == 2:
                    if self.tcolors[v.name] == "Yellow" and self.tcolors[v.up.name] != "Yellow":
                        tempMax = 1
                    elif self.tcolors[v.name] != "Yellow":
                        tempMax = 2
                self.cMax1[1, v.name] = tempMax



    def pairExt(self):
        redYellow = ["Red", "Yellow"]

        # Determine if red (and yellow) subtrees came originally from the first or second input tree
        if self.init_swap == 0:
            original_colors = [0, 1]
        else:
            original_colors = [1, 0]


        # Keep track of which nodes will be deleted (since we are moving a grafted subtree) and
        # in which cases, if any, we need to remove the root node and assign a new root node
        to_delete1, to_delete2 = [], []
        is_old_root1, is_old_root2 = [], []


        # Pair up corresponding extraneous clades
        for v in self.t.traverse("postorder"):
            if v.is_leaf():
                self.order1[0, v.name], self.order1[1, v.name] = [], []
                if self.tcolors[v.name] == "Red" and self.tcolors[v.up.name] != "Red":
                    self.order1[0, v.name].append(v)
                elif self.tcolors[v.name] == "Yellow" and self.tcolors[v.up.name] != "Yellow":
                    self.order1[1, v.name].append(v)
            else:
                self.order1[0, v.name], self.order1[1, v.name] = [], []
                for c in [0,1]:
                    if self.tcolors[v.name] == redYellow[c] and self.tcolors[v.up.name] != redYellow[c]:
                        self.order1[c, v.name].append(v)
                    else:
                        for child in v.children:
                            self.order1[c, v.name].extend(self.order1[c, child.name])

            # Keep track of the local optimal color and number of unpaired subtrees
            C, N = self.optc1[v.name], self.optN1[v.name]
            m = len(self.order1[C, v.name])

            # If both red and yellow subtrees have been pushed up to v, then pair up based on the order they appear.
            # After pairing, let the new order of color C be the order of the remaining unpaired subtrees,
            # and let the order of color 1-C be empty

            if self.order1[0, v.name] and self.order1[1, v.name]:
                for n in range(m - N):

                    newT1, newT2 = Tree(), Tree()

                    ext_left1, ext_right1 = self.order1[0, v.name][n], self.order1[1, v.name][n]
                    ext_left2, ext_right2 = self.map2[ext_left1.name], self.map2[ext_right1.name]

                    if original_colors[0] == 0:
                        up1, up2 = ext_left1.up, ext_right2.up

                        to_delete1.append(ext_right1.up)
                        to_delete2.append(ext_left2.up)

                        if ext_right1.up.up is not None:
                            is_old_root1.append(None)
                        else:
                            is_old_root1.append(ext_right1.up)

                        if ext_left2.up.up is not None:
                            is_old_root2.append(None)
                        else:
                            is_old_root2.append(ext_left2.up)
                    else:
                        up1, up2 = ext_right1.up, ext_left2.up

                        to_delete1.append(ext_left1.up)
                        to_delete2.append(ext_right2.up)

                        if ext_right2.up.up is not None:
                            is_old_root2.append(None)
                        else:
                            is_old_root2.append(ext_right2.up)

                        if ext_left1.up.up is not None:
                            is_old_root1.append(None)
                        else:
                            is_old_root1.append(ext_left1.up)

                    ext_left1, ext_right1 = ext_left1.detach(), ext_right1.detach()
                    ext_left2, ext_right2 = ext_left2.detach(), ext_right2.detach()

                    newT1.add_child(ext_left1)
                    newT1.add_child(ext_right1)
                    up1.add_child(newT1)
                    newT2.add_child(ext_left2)
                    newT2.add_child(ext_right2)
                    up2.add_child(newT2)

                self.order1[C, v.name], self.order1[1-C, v.name] = self.order1[C, v.name][m-N:], []

            # If there are only maximal red subtrees, keep the last N of them
            elif self.order1[0, v.name]:
                if C == 0:
                    self.order1[0, v.name] = self.order1[0, v.name][m-N:]
                else:
                    self.order1[0, v.name] = []

            # If there are only maximal yellow subtrees, keep the last N of them
            elif self.order1[1, v.name]:
                if C == 1:
                    self.order1[1, v.name] = self.order1[1, v.name][m-N:]
                else:
                    self.order1[1, v.name] = []

        # Once we have already merged ALL extraneous clades, then delete the redundant nodes
        # If deleting is done too early, then the subtree order that we keep track of will not properly collapse
        for n in range(len(to_delete1)):
            if is_old_root1[n]:
                self.t = is_old_root1[n].children[0]                    # We have already pruned the second child
                self.t.up = None
            else:
                to_delete1[n].delete()

            if is_old_root2[n]:
                self.t2 = is_old_root2[n].children[0]                    # We have already pruned the second child
                self.t2.up = None
            else:
                to_delete2[n].delete()



    def EF_U_RF(self):
        if self.EF1 and self.EF2:
            return self.EF1, self.EF2

        # Outgroup at some arbitrary green leaf and split accordingly to root
        #     leftover rooted tree in the proper place
        outgroup1, outgroup2 = None, None
        for leaf in self.t2:
            if(self.tLeafset.get(str(leaf.name)) != None):
                outgroup1 = self.t.get_leaves_by_name(str(leaf.name))[0]
                outgroup2 = leaf
                break

        self.t.set_outgroup(outgroup1)
        self.t2.set_outgroup(outgroup2)
        T1, T2 = self.t, self.t2

        if len(self.t.children[0]) == 1:
            self.t = self.t.children[1].detach()
            T1 = T1.children[0].detach()
        else:
            T1 = T1.children[1].detach()
            self.t = self.t.children[0].detach()

        if len(self.t2.children[0]) == 1:
            self.t2 = self.t2.children[1].detach()
            T2 = T2.children[0].detach()
        else:
            T2 = T2.children[1].detach()
            self.t = self.t.children[0].detach()

        T1.up, T2.up, self.t.up, self.t2.up = None, None, None, None

        # Compute EF-RF+ distance, add the green leaf back at the root and return the result
        rooted1, rooted2 = self.EF_R_RF()
        newT1, newT2 = Tree(), Tree()
        newT1.add_child(T1)
        newT1.add_child(rooted1)
        newT2.add_child(T2)
        newT2.add_child(rooted2)
        self.EF1, self.EF2 = newT1, newT2
        self.EF1.unroot()
        self.EF2.unroot()
        return self.EF1, self.EF2



    def U_RF_Plus(self):
        if self.RF1 and self.RF2:
            return self.RF1, self.RF2

        # Outgroup at some arbitrary green leaf and split accordingly to root
        #     leftover rooted tree in the proper place
        outgroup1, outgroup2 = None, None
        for leaf in self.t2:
            if(self.tLeafset.get(str(leaf.name)) != None):
                print(leaf.name)
                outgroup1 = self.t.get_leaves_by_name(str(leaf.name))[0]
                outgroup2 = leaf
                print(outgroup1, outgroup2)
                break

        self.t.set_outgroup(outgroup1)
        self.t2.set_outgroup(outgroup2)
        T1, T2 = self.t, self.t2

        if len(self.t.children[0]) == 1:
            self.t = self.t.children[1].detach()
            T1 = T1.children[0].detach()
        else:
            T1 = T1.children[1].detach()
            self.t = self.t.children[0].detach()

        if len(self.t2.children[0]) == 1:
            self.t2 = self.t2.children[1].detach()
            T2 = T2.children[0].detach()
        else:
            T2 = T2.children[1].detach()
            self.t = self.t.children[0].detach()

        T1.up, T2.up, self.t.up, self.t2.up = None, None, None, None

        # Compute RF+ distance, add the green leaf back at the root and return the result
        rooted1, rooted2 = self.R_RF_Plus()
        newT1, newT2 = Tree(), Tree()
        newT1.add_child(T1)
        newT1.add_child(rooted1)
        newT2.add_child(T2)
        newT2.add_child(rooted2)
        self.RF1, self.RF2 = newT1, newT2
        self.RF1.unroot()
        self.RF2.unroot()
        return self.RF1, self.RF2



    def R_RF_Plus(self):

        # If the RF(+) completions have already been computed, then they have been stored in the compareTree class
        if self.RF1 and self.RF2:
            return self.RF1, self.RF2
        # Otherwise, continue with computing the RF(+) completions

        self.start_time = time.time()

        # First, call the EF-R-RF(+) completion function and preprocess any additional necessary information
        # NOTE: If the EF-R-RF(+) completions have been computed already, then the EF_R_RF() method will be skipped
        self.EF_R_RF()


        # We need to frequently relabel the internal nodes in the tree to be processed to prevent confusion between labels of internal nodes
        self.resetT()

        # The variable self.orig_root_name is only used as a check for the sagephy produced test cases
        self.orig_root_name = self.t.name

        self.recolor()
        self.setRedMax()
        self.setYellowMax()


        # Create LCA Mapping from self.t to self.t2 to determine which clades are matches
        self.t2mapping = LCAMapping(self.t2)
        self.map2 = {}
        for node in self.t.traverse("postorder"):
            if node.is_leaf():
                self.map2[node.name] = self.t2mapping.queryLCA(node.name, node.name)
            else:
                lst = []
                for child in node.children:
                    lst.append(self.map2.get(child.name))
                self.map2[node.name] = self.t2mapping.queryLCA(lst[0].name, lst[1].name)


        # Now run the DP recurrence relation based off of self.t
        for node in self.t.traverse("postorder"):
            children, match = [], []
            redYellow = ["Red", "Yellow"]
            if node.is_leaf():
                for c in [0,1]:
                    self.cost[node.name, 0, c] = 0
                    if self.tcolors[node.name] == redYellow[c] and self.tcolors[node.up.name] != redYellow[c]:
                        self.cost[node.name, 1, c] = 0

            else:
                # First, determine the children of the node, and compute if each child is a match between the two EF-RF completions.
                # If a child is a match, then there is a chance of the match being broken after moving red/yellow subtrees up
                # (when the child is not the root of a maximal red/yellow subtree)
                for child in node.children:
                    children.append(child)
                    if len(child) == len(self.map2[child.name]) and (self.tcolors[child.name] != "Red" and self.tcolors[child.name] != "Yellow"):
                        match.append(1)
                    else:
                        match.append(0)

                # Compute the cost values for all possible values of c and N, given the vertex "node"
                for c in [0,1]:
                    if self.cMax1[c, node.name] == 0:
                        self.cost[node.name, 0, c] = 0
                        self.leftN1[node.name, 0, c], self.rightN1[node.name, 0, c] = 0, 0
                        self.leftc1[node.name, 0, c], self.rightc1[node.name, 0, c] = 0, 0
                        continue

                    for N in range(self.cMax1[c, node.name] + 1):

                        if self.cMax1[c, node.name] > 0:
                            self.cost[node.name, N, c] = math.inf
                        else:
                            self.cost[node.name, N, c] = 0
                        self.leftN1[node.name, N, c], self.rightN1[node.name, N, c] = 0, 0
                        self.leftc1[node.name, N, c], self.rightc1[node.name, N, c] = 0, 0

                        # Determine viable ranges for cL, cR, NL and NR, the input cost variables for each of the children
                        # Iterate over all possible 4-tuples and set the cost at v with variables N, c as the minimum possible value
                        #     over all 4-tuples (cL, cR, NL, NR).
                        # Save the optimal values for each of the four variables to construct the completion later

                        # NOTE: The value of N represents the number of left over unpaired maximal red/yellow subtrees.
                        # Therefore, under the assumption that we always pair up opposite colored subtrees that have already been pushed
                        # up to the given vertex, we determine the possible values of NL and NL (quantities pushed from left, right respectively)
                        # DEPENDING on the values of cL and cR (colors of subtrees pushed from left, right respectively)

                        for cL, cR in itertools.product([0, 1], [0, 1]):

                            # There are only three cases. The fourth case never results in the propoper color being pushed up to the specified node
                            # The three meaningful cases are when
                            #   1.) cL == cR == c,
                            #   2.) cL == c != cR,
                            #   3.) cR == c != cL


                            # If the left and right pushed colors BOTH equal the desired pushed color, then the quantities of pushed subtrees from
                            # each of the children must sum to N (since no extraneous clades can be formed)
                            if cL == cR == c:
                                if N <= min(self.cMax1[cL, children[0].name], self.cMax1[cR, children[1].name]):
                                    for n in range(N + 1):
                                        NL = n
                                        NR = N - n

                                        temp = self.cost[children[0].name, NL, cL] + self.cost[children[1].name, NR, cR]
                                        if NL > 0 and self.tcolors[children[0].name] != redYellow[cL]:
                                            temp = temp + match[0]
                                        if NR > 0 and self.tcolors[children[1].name] != redYellow[cR]:
                                            temp = temp + match[1]
                                        if temp < self.cost[node.name, N, c]:
                                            self.cost[node.name, N, c] = temp
                                            self.leftN1[node.name, N, c], self.rightN1[node.name, N, c] = NL, NR
                                            self.leftc1[node.name, N, c], self.rightc1[node.name, N, c] = cL, cR

                                elif N <= self.cMax1[cL, children[0].name]:
                                    for n in range(self.cMax1[cR, children[1].name] + 1):
                                        NL = N - n
                                        NR = n

                                        temp = self.cost[children[0].name, NL, cL] + self.cost[children[1].name, NR, cR]
                                        if NL > 0 and self.tcolors[children[0].name] != redYellow[cL]:
                                            temp = temp + match[0]
                                        if NR > 0 and self.tcolors[children[1].name] != redYellow[cR]:
                                            temp = temp + match[1]
                                        if temp < self.cost[node.name, N, c]:
                                            self.cost[node.name, N, c] = temp
                                            self.leftN1[node.name, N, c], self.rightN1[node.name, N, c] = NL, NR
                                            self.leftc1[node.name, N, c], self.rightc1[node.name, N, c] = cL, cR

                                elif N <= self.cMax1[cR, children[1].name]:
                                    for n in range(self.cMax1[cL, children[0].name] + 1):
                                        NL = n
                                        NR = N - n

                                        temp = self.cost[children[0].name, NL, cL] + self.cost[children[1].name, NR, cR]
                                        if NL > 0 and self.tcolors[children[0].name] != redYellow[cL]:
                                            temp = temp + match[0]
                                        if NR > 0 and self.tcolors[children[1].name] != redYellow[cR]:
                                            temp = temp + match[1]
                                        if temp < self.cost[node.name, N, c]:
                                            self.cost[node.name, N, c] = temp
                                            self.leftN1[node.name, N, c], self.rightN1[node.name, N, c] = NL, NR
                                            self.leftc1[node.name, N, c], self.rightc1[node.name, N, c] = cL, cR

                                else:
                                    for n in range(self.cMax1[c, node.name] - N + 1):
                                        NL = self.cMax1[cL, children[0].name] - n
                                        NR = self.cMax1[cR, children[1].name] - self.cMax1[c, node.name] + N + n

                                        temp = self.cost[children[0].name, NL, cL] + self.cost[children[1].name, NR, cR]
                                        if NL > 0 and self.tcolors[children[0].name] != redYellow[cL]:
                                            temp = temp + match[0]
                                        if NR > 0 and self.tcolors[children[1].name] != redYellow[cR]:
                                            temp = temp + match[1]
                                        if temp < self.cost[node.name, N, c]:
                                            self.cost[node.name, N, c] = temp
                                            self.leftN1[node.name, N, c], self.rightN1[node.name, N, c] = NL, NR
                                            self.leftc1[node.name, N, c], self.rightc1[node.name, N, c] = cL, cR

                            # If cL == c != cR, then we know that there must be N more subtrees pushed from the left than the right
                            # NOTE: The bound depends on whether or not each of the children is a maximal red/blue subtree
                            elif cL == c and self.cMax1[cL, children[0].name] >= N:
                                bound1 = self.cMax1[cL, children[0].name] - N

                                bound2 = self.cMax1[cR, children[1].name]

                                for n in range(min(bound1, bound2) + 1):
                                    NL = N + n
                                    NR = n

                                    temp = self.cost[children[0].name, NL, cL] + self.cost[children[1].name, NR, cR] - n
                                    if NL > 0 and self.tcolors[children[0].name] != redYellow[cL]:
                                        temp = temp + match[0]
                                    if NR > 0 and self.tcolors[children[1].name] != redYellow[cR]:
                                        temp = temp + match[1]
                                    if temp < self.cost[node.name, N, c]:
                                        self.cost[node.name, N, c] = temp
                                        self.leftN1[node.name, N, c], self.rightN1[node.name, N, c] = NL, NR
                                        self.leftc1[node.name, N, c], self.rightc1[node.name, N, c] = cL, cR

                            # If cR == c != cL, then we know that there must be N more subtrees pushed from the right than the left
                            # NOTE: The bound depends on whether or not each of the children is a maximal red/blue subtree
                            elif cR == c and self.cMax1[cR, children[1].name] >= N:
                                bound1 = self.cMax1[cL, children[0].name]

                                bound2 = self.cMax1[cR, children[1].name] - N

                                for n in range(min(bound1, bound2) + 1):
                                    NL = n
                                    NR = N + n

                                    temp = self.cost[children[0].name, NL, cL] + self.cost[children[1].name, NR, cR] - n
                                    if NL > 0 and self.tcolors[children[0].name] != redYellow[cL]:
                                        temp = temp + match[0]
                                    if NR > 0 and self.tcolors[children[1].name] != redYellow[cR]:
                                        temp = temp + match[1]
                                    if temp < self.cost[node.name, N, c]:
                                        self.cost[node.name, N, c] = temp
                                        self.leftN1[node.name, N, c], self.rightN1[node.name, N, c] = NL, NR
                                        self.leftc1[node.name, N, c], self.rightc1[node.name, N, c] = cL, cR


        # Determine top-down how to pair up extraneous clades
        # We can assume that the optimal number N of unpaired maximal colored trees (at the root) is equal to 0
        # Moreover, we can assume that these unpaired trees are red (0)
        self.optN1[self.t.name], self.optc1[self.t.name] = 0, 0

        for v in self.t.traverse("preorder"):
            children = []
            if not v.is_leaf():
                for child in v.children:
                    children.append(child)
                self.optN1[children[0].name] = self.leftN1[v.name, self.optN1[v.name], self.optc1[v.name]]
                self.optN1[children[1].name] = self.rightN1[v.name, self.optN1[v.name], self.optc1[v.name]]
                self.optc1[children[0].name] = self.leftc1[v.name, self.optN1[v.name], self.optc1[v.name]]
                self.optc1[children[1].name] = self.rightc1[v.name, self.optN1[v.name], self.optc1[v.name]]

        self.pairExt()

        if self.swapped % 2 == 0:
            self.RF1, self.RF2 = self.t, self.t2
        else:
            self.RF1, self.RF2 = self.t2, self.t


        self.RF_time = time.time() - self.start_time
        return self.RF1, self.RF2





def convertTreeIntoBinary(T):
    for node in T.traverse("postorder"):
        if ( not node.is_leaf() ):
            if (len(node.children) > 2):
                ndLst = []
                for child in node.children:
                    ndLst.append(child)
                for i in range(len(ndLst)):
                    ndLst[i] = ndLst[i].detach()
                oldTree = Tree()
                oldTree.add_child(ndLst[0])
                oldTree.add_child(ndLst[1])
                for i in range(2, len(ndLst)):
                    if( (i == len(ndLst) - 1) and (ndLst[i].up == None) ):
                        node.add_child(oldTree)
                        node.add_child(ndLst[i])
                    else:
                        newTree = Tree()
                        newTree.add_child(oldTree)
                        newTree.add_child(ndLst[i])
                        oldTree = newTree




def main():
    # takes an input file (required)
    # -o output file to print to is optional. when not present, print to treminal
    # -r is optional. when present, it prints out rf distance instead of modified trees
    parser = argparse.ArgumentParser()
    parser.add_argument("-ext", "--extraneousfree", action="store_true", help = "This flag signifies computation of the EF-RF(+) completions rather than the more general RF(+) completions.")
    parser.add_argument("-u", "--unrooted", action="store_true", help = "This flag signifies that the input trees are unrooted.")
    parser.add_argument("-i", "--inputfile", type = str, help = "The input file contains the trees in newick format. This argument is required.")
    parser.add_argument("-o", "--outputfile", type = str, help = "The output file to which the RF distance and completions in newick format will be printed")
    #parser.add_argument("-r", "--rfdistance", action="store_true", help = "Type this command to print the RF(-), EF-RF(+) and RF(+) distances instead of the completed trees. If this flag is used, then the -ext flag is not necessary.")
    #parser.add_argument("-a", "--analysisfile", type = str, help = "A csv file which will store every pair of trees, labeled by line number, along with the RF(-) and RF(+) distances, sizes of the intersection and union of input tree leaf sets, and runtime for EF-R-RF(+) and R-RF(+) completions. Note the recorded RF(+) runtime is the runtime to compute the RF(+) distance assuming the EF-RF(+) completions have already been computed.")
    args = parser.parse_args()
    msn = []
    index = 0
    with open(args.inputfile) as openbn:
        for line in openbn:
            x = line.strip(' ').strip('\t').strip('\n').strip(' ').strip('\t').strip('\n').strip(' ').strip('\t').strip('\n')
            if(x != ''):
                if(x[-1] != ";"):
                    x = x + ";"
                t = Tree(x)
                # checking to see if give tree is binary or not
                # with the new convertTreeIntoBinary, should no longer reach the error
                convertTreeIntoBinary(t)
                for node in t.traverse("preorder"):
                    if(len(node.children) > 2):
                        print("Error - non binary tree line: {}\n".format(index+1))
                        return
                msn.append(t)
                index += 1

    if(args.outputfile):
        fileNameW = args.outputfile
        fileW = open(fileNameW, "w")


    i = 0
    for j in range(i+1, len(msn)):
        #print("I, J:\t", i, j)
        #print("Length of tree {}:\t".format(j), len(msn[j]))


        a = msn[i].robinson_foulds(msn[j])[0]
        y = compareTree(msn[i], msn[j])
        if args.unrooted:
            y1 = y.EF_U_RF()
            y2 = y.U_RF_Plus()
        else:
            y1 = y.EF_R_RF()
            #y3 = compareTree(msn[1], msn[j])
            y2 = y.R_RF_Plus()
        z1 = y1[0].robinson_foulds(y1[1])[0]
        z2 = y2[0].robinson_foulds(y2[1])[0]


        if args.outputfile:
            fileW.write("Results for Tree {} and Tree {}\n".format(i+1, j+1))
            fileW.write("Number of leaves in Tree {}: {}\n".format(i+1, len(msn[i])))
            fileW.write("Number of leaves in Tree {}: {}\n".format(j+1, len(msn[j])))
            fileW.write("Union of leaf sets: {}\n".format(len(y2[0])))
            fileW.write("Intersection of leaf sets: {}\n".format(len(y2[0]) - len(y.yellowLeaves) - len(y.redLeaves)))
            fileW.write("RF(-) distance: {}\n".format(a))
            fileW.write("RF(+) distance: {}\n".format(z2))
            fileW.write("EF-RF(+) distance: {}\n".format(z1))
            fileW.write("Completed trees:\n")

            if args.extraneousfree:
                fileW.write(y1[0].write(format = 9) + "\n\n")
                #fileW.write("\n\n")
                fileW.write(y1[1].write(format = 9) + "\n\n--------------------------------------------------------------------------------------------------------------------------\n\n")
            else:
                fileW.write(y2[0].write(format = 9) + "\n\n")
                #fileW.write("\n\n")
                fileW.write(y2[1].write(format = 9) + "\n\n--------------------------------------------------------------------------------------------------------------------------\n\n")

        else:
            print("Results for Tree {} and Tree {}".format(i+1, j+1))
            print("Number of leaves in Tree {}: {}".format(i+1, len(msn[i])))
            print("Number of leaves in Tree {}: {}".format(j+1, len(msn[j])))
            print("Union of leaf sets: {}".format(len(y2[0])))
            print("Intersection of leaf sets: {}".format(len(y2[0]) - len(y.yellowLeaves) - len(y.redLeaves)))
            print("RF(-) distance: {}".format(a))
            print("RF(+) distance: {}".format(z2))
            print("EF-RF(+) distance: {}".format(z1))
            print("Completed trees:")

            if args.extraneousfree:
                print(y1[0].write(format = 9) + "\n\n")
                #print("\n\n")
                print(y1[1].write(format = 9) + "\n\n--------------------------------------------------------------------------------------------------------------------------\n\n")
            else:
                print(y2[0].write(format = 9) + "\n\n")
                #print("\n\n")
                print(y2[1].write(format = 9) + "\n\n--------------------------------------------------------------------------------------------------------------------------\n\n")


    if(args.outputfile):
        fileW.close()



if __name__ == "__main__":
    main()
    #pass
