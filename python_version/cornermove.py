from random import choice, uniform
from math import exp

class CornerMove(object):
    def __init__(self, space, beta):
        self.space = space
        self.movetype = "CornerMove"
        self.beta = beta
    def move(self, polyid, polytyp):
        """
        The corner move
        """
        if polytyp == 'sim':
            poly = self.space.Sims[polyid-1]
            sitevalue = 1 # code for free sim
        elif polytyp == 'sumo':
            poly = self.space.Sumos[polyid-1]
            sitevalue = -1 # code for free sumo
        else:
            raise Exception("polymer type undefined")
        
        possible_moves = []
        for pos in range(len(poly.locs) - 2):
            x1, y1 = poly.locs[pos]
            x2, y2 = poly.locs[pos+1]
            x3, y3 = poly.locs[pos+2]
            if x1 != x3 and y1 != y3:
                if x1 == x2:
                    assert(y1 != y2 and y2 == y3 and x2 != x3)
                    if self.space.space[x3][y1] == 0:
                        possible_moves.append((pos+1, (x3, y1)))
                if x1 != x2:
                    assert(y1 == y2 and x2 == x3 and y2 != y3)
                    if self.space.space[x1][y3] == 0:
                        possible_moves.append((pos+1, (x1, y3)))
                    
        if not possible_moves:
            return
        
        (pointid, (xnew, ynew)) = choice(possible_moves)
        xold, yold = poly.locs[pointid]
        nbr_bond_inc = 0
        if abs(self.space.space[xold][yold]) != 1: # (xold, yold) is in a bond
            nbr_bond_inc -= 1 # will lose a bond
            
        bond_choice = [bpoint for bpoint in self.space.neighbor((xnew, ynew)) if self.space.can_build_bond((xnew, ynew), bpoint)]
        if bond_choice: # (x, y) can build a bond
            nbr_bond_inc += 1 # will create a bond
        
        if uniform(0,1) < self.weight(nbr_bond_inc, self.beta): # 
            # removes the monomer on (xold, yold)
            # handles everything except in its polymer.locs and rspace
            self.space.safe_remove(xold, yold) 
            
            # create a monomer on (xnew, ynew)
            # handles everything except in its polymer.locs and rspace
            self.space.safe_create(xnew, ynew, sitevalue)
            
            # make changes on polymer.locs
            poly.locs[pointid] = (xnew, ynew)
                        
            # handle the reverse checking space (not done)
            assert(self.space.rspace[xnew][ynew][0] == 0 and self.space.rspace[xold][yold][0] != 0)
            self.space.rspace[xnew][ynew][0] = self.space.rspace[xold][yold][0]
            self.space.rspace[xnew][ynew][1] = self.space.rspace[xold][yold][1]
            self.space.rspace[xold][yold][0], self.space.rspace[xold][yold][1] = 0, 0
            
            # try building new bonds
            bond_choice = [bpoint for bpoint in self.space.neighbor((xnew, ynew)) if self.space.can_build_bond((xnew, ynew), bpoint)]
            if bond_choice:
                (xb, yb) = choice(bond_choice)
                self.space.create_bond((xnew, ynew), (xb, yb))
            
            
    def weight(self, nbr_bond_inc, beta):
        """
        The weight function for endmove
        """
        return exp(nbr_bond_inc*beta)