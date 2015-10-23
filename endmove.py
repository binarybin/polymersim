from random import choice, uniform


class EndMove(object):
    def __init__(self, space):
        self.space = space
        self.movetype = "EndMove"
    def move(self, polyid, polytyp, endtyp):
        """
        The end move
        """
        if polytyp == 'sim':
            poly = self.space.Sims[polyid-1]
            sitevalue = 1 # code for free sim
        elif polytyp == 'sumo':
            poly = self.space.Sumos[polyid-1]
            sitevalue = -1 # code for free sumo
        else:
            raise Exception("polymer type undefined")
            
        if endtyp == 'head':
            endloc, nextloc = 0, 1
        elif endtyp == 'tail':
            endloc, nextloc = len(poly.locs)-1, len(poly.locs)-2
        else:
            raise Exception("end type undefined")
            
        assert(self.space.space[poly.locs[endloc][0]][poly.locs[endloc][1]] != 0)
            
        newpoint = choice(list(self.space.neighbor(poly.locs[nextloc]) - set([poly.locs[endloc]])))
        
        x, y = newpoint
        if self.space.space[x][y] != 0: # site occupied, not possible to move
            return
        
        if uniform(0,1) < self.weight(poly, endloc, poly.locs[endloc], newpoint, 1): # 
            xold, yold = poly.locs[endloc]
            # removes the monomer on (xold, yold)
            # handles everything except in its polymer.locs and rspace
            self.space.safe_remove(xold, yold) 
            
            # make changes on polymer.locs
            poly.locs[endloc] = newpoint
            
            # create a monomer on (x, y)
            # handles everything except in its polymer.locs and rspace
            self.space.safe_create(x, y, sitevalue)
            
            # handle the reverse checking space
            assert(self.space.rspace[x][y][0] == 0 and self.space.rspace[xold][yold][0] != 0)
            self.space.rspace[x][y][0] = self.space.rspace[xold][yold][0]
            self.space.rspace[x][y][1] = self.space.rspace[xold][yold][1]
            self.space.rspace[xold][yold][0], self.space.rspace[xold][yold][1] = 0, 0
            
            # try building new bonds
            bond_choice = [bpoint for bpoint in self.space.neighbor(newpoint) if self.space.can_build_bond(newpoint, bpoint)]
            if bond_choice:
                (xb, yb) = choice(bond_choice)
                self.space.create_bond((x, y), (xb, yb))
            
            
    def weight(self, polymer, endid, frompoint, topoint, beta):
        """
        The weight function for endmove
        """
        return 2