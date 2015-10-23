from random import uniform, choice


class SnakeMove(object):
    def __init__(self, space):
        self.space = space
        self.movetype = "SnakeMove"
    def move(self, polyid, polytyp, endtyp):
        """
        The snake move
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
            beginloc = len(poly.locs)-1
        elif endtyp == 'tail':
            endloc, nextloc = len(poly.locs)-1, len(poly.locs)-2
            beginloc = 0
        else:
            raise Exception("end type undefined")
        
        x0, y0 = poly.locs[nextloc]
        x1, y1 = poly.locs[endloc]
        x2, y2 = (2*x1 - x0) % self.space.Lx, (2*y1 - y0) % self.space.Ly
        newpoint = (x2, y2)

        if self.space.space[x2][y2] != 0: # site occupied, not possible to move
            return
        
        if uniform(0,1) < self.weight(poly, endloc, poly.locs[endloc], newpoint, 1): # 
            xold, yold = poly.locs[beginloc]
            xnew, ynew = x2, y2
            # removes the monomer on (xold, yold)
            # handles everything except in its polymer.locs and rspace
            self.space.safe_remove(xold, yold) 
            
            # create a monomer on (xnew, ynew)
            # handles everything except in its polymer.locs and rspace
            self.space.safe_create(xnew, ynew, sitevalue)
            
            # make changes on polymer.locs
            if endtyp == 'head':
                del poly.locs[-1]
                poly.locs.insert(0, newpoint)
            elif endtyp == 'tail':
                del poly.locs[0]
                poly.locs.append(newpoint)
            else:
                raise Exception("end type undefined")
            
            # handle the reverse checking space
            assert(self.space.rspace[xnew][ynew][0] == 0 and self.space.rspace[xold][yold][0] != 0)
            self.space.rspace[xold][yold][0], self.space.rspace[xold][yold][1] = 0, 0
            for idx, (xtemp, ytemp) in enumerate(poly.locs):
                self.space.rspace[xtemp][ytemp][0] = polyid
                self.space.rspace[xtemp][ytemp][1] = idx
            
            # try building new bonds
            bond_choice = [bpoint for bpoint in self.space.neighbor((xnew, ynew)) if self.space.can_build_bond(newpoint, bpoint)]
            if bond_choice:
                (xb, yb) = choice(bond_choice)
                self.space.create_bond((xnew, ynew), (xb, yb))
            
    def weight(self, polymer, endid, frompoint, topoint, beta):
        """
        The weight function for endmove
        """
        return 2