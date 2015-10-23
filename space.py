from __future__ import division
import numpy as np
import itertools

NOBOND = -128

class Space(object):
    def __init__(self, NSim, LSim, NSumo, LSumo, Lx, Ly):
        """
        The simplest initialization, just copy numbers.
        """
        self.NSim, self.NSumo = NSim, NSumo
        self.LSim, self.LSumo = LSim, LSumo
        self.Lx, self.Ly = Lx, Ly
        
        self.space = np.zeros((Lx, Ly), dtype=np.int8)
        
        # bond[x][y][0] = the other element in the bond, xvalue
        # bond[x][y][1] = the other element in the bond, yvalue
        # -128 is reserved as "no bond"
        self.bond = np.ones((Lx, Ly, 2), dtype=np.int16) * NOBOND
        # the reverse check space
        # rspace[x][y][0] = the polymer id
        # rspace[x][y][1] = the monomer id in the polymer
        self.rspace = np.zeros((Lx, Ly, 2), dtype=np.int16)
        self.Sims, self.Sumos = [], []
        self.SimId, self.SumoId = 0, 0
        
    def initialize(self):
        """
        The initialization disatcher
        """
        if max(self.LSim, self.LSumo) <= self.Ly and self.NSim+self.NSumo < (self.Lx):
            self.dilute_init('v') # place polymers vertically ("pointing up")
        elif max(self.LSim, self.LSumo) <= self.Lx and self.NSim+self.NSumo < (self.Ly):
            self.dilute_init('h') # place polymers horizontally ("pointing right")
        else:
            self.dense_init()
        
    def place(self, typ, poly_id, locs):
        """
        Place a polymer
        Append this polymer info to self.Sims and self.Sumos
        Mark space as occupied:
            1 for free sim
            2 for sim with bond
           -1 for free sumo
           -2 for sumo with bond
        """
        locs = [x for x in locs]
        if typ == 'sim':
            self.Sims.append(Polymer(poly_id, locs))
            for idx, (x, y) in enumerate(locs):
                self.space[x][y] = 1
                self.rspace[x][y][0] = poly_id
                self.rspace[x][y][1] = idx
        elif typ == 'sumo':
            self.Sumos.append(Polymer(poly_id, locs))
            for idx, (x, y) in enumerate(locs):
                self.space[x][y] = -1
                self.rspace[x][y][0] = poly_id
                self.rspace[x][y][1] = idx
        else:
            raise Exception("Unrecognized polymer type.")
        
    def dilute_init(self, direction):
        """
        The initialization in the dilute (not-so-dense) limit
        """
        if direction == 'v':
            dl = self.Lx // (self.NSim + self.NSumo)
            for x in range(1, self.NSim+1): # id for polymers, start from 1
                self.place('sim', x, itertools.product([(x-1)*dl], range(self.LSim)))
            for x in range(1, self.NSumo+1):
                self.place('sumo', x, itertools.product([self.Lx-x*dl], range(self.LSumo)))
        elif direction == 'h':
            dl = self.Ly // (self.NSim + self.NSumo)
            for y in range(1, self.NSim+1):
                self.palce('sim', y, itertools.product(range(self.LSim), [(y-1) * dl]))
            for y in range(1, self.NSumo+1):
                self.place('sumo', y, itertools.product(range(self.LSumo), [(self.Ly-y)*dl]))
        
    def dense_init(self):
        """
        The initialization in the dense limit
        """
        raise Exception("Polymer density relatively high, requiring better init method.")
        
    def neighbor(self, point):
        """
        Get a set of neighbors of a point
        """
        x, y = point
        four_points = [(x+1, y), (x-1, y), (x, y+1), (x, y-1)]
        points = [(x%self.Lx, y%self.Ly) for x, y in four_points]
        return set(points)
        
    def safe_remove(self, x, y):
        """
        Removes the monomer at (x, y) and
        safely register the change except on polymer.locs and rspace
        """
        if (self.space[x][y]==0):
            print self.rspace[x][y]
        assert(self.space[x][y] != 0)
        self.space[x][y] = 0
        
        xb, yb = self.bond[x][y][0], self.bond[x][y][1]
        if xb != NOBOND and yb != NOBOND: # we do have a bond to remove 
            assert(self.bond[xb][yb][0] == x and self.bond[xb][yb][1] == y) # assert bond reversibility
            assert(abs(self.space[xb][yb]) == 2) # assert the bond consistency
            self.bond[x][y][0], self.bond[x][y][1] = NOBOND, NOBOND # remove the bond from (x,y)
            self.bond[xb][yb][0], self.bond[xb][yb][1] = NOBOND, NOBOND # remove the bond to (x,y)
            self.space[xb][yb] = np.sign(self.space[xb][yb]) # register the space site to be no-bond
            
            
    def safe_create(self, x, y, polyvalue):
        """
        Create the monomer at (x, y) and
        safely register the change except on polymer.locs and rspace
        """
        assert(self.space[x][y] == 0)
        self.space[x][y] = polyvalue
        
    def create_bond(self, (x1, y1), (x2, y2)):
        """
        Create the bond between (x1, y1) and (x2, y2)
        Note that bond information is registered in self.bond
        """
        assert(self.can_build_bond((x1, y1), (x2, y2))) # assert the bond consistency
        assert((self.bond[x1][y1][0] == NOBOND) and (self.bond[x1][y1][0] == NOBOND)) # assert the bond nonexistance
        assert((self.bond[x2][y2][0] == NOBOND) and (self.bond[x2][y2][0] == NOBOND)) # assert the bond nonexistance
        self.bond[x1][y1][0], self.bond[x1][y1][1] = x2, y2
        self.bond[x2][y2][0], self.bond[x2][y2][1] = x1, y1
        self.space[x1][y1] *= 2
        self.space[x2][y2] *= 2
        
    def can_build_bond(self, (x1, y1), (x2, y2)):
        return self.space[x1][y1] * self.space[x2][y2] == -1
    def exist_bond(self, (x1, y1), (x2, y2)):
        return self.space[x1][y1] * self.space[x2][y2] == -4
        
class Polymer(object):
    def __init__(self, poly_id, locs):
        self.id = poly_id
        self.locs = locs