from __future__ import division
import numpy as np
import itertools

class Space(object):
    def __init__(self, NSim, LSim, NSumo, LSumo, Lx, Ly):
        """
        The simplest initialization, just copy numbers.
        """
        self.NSim, self.NSumo = NSim, NSumo
        self.LSim, self.LSumo = LSim, LSumo
        self.Lx, self.Ly = Lx, Ly
        
        self.space = np.zeros((Lx, Ly), dtype=np.int8)
        
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
            for x, y in locs:
                self.space[x][y] = 1
        elif typ == 'sumo':
            self.Sumos.append(Polymer(poly_id, locs))
            for x, y in locs:
                self.space[x][y] = -1
        else:
            raise Exception("Unrecognized polymer type.")
        
    def dilute_init(self, direction):
        """
        The initialization in the dilute (not-so-dense) limit
        """
        if direction == 'v':
            dl = self.Lx // (self.NSim + self.NSumo)
            for x in range(self.NSim):
                self.place('sim', x, itertools.product([x*dl], range(self.LSim)))
            for x in range(self.NSumo):
                self.place('sumo', x, itertools.product([self.Lx-(1+x)*dl], range(self.LSumo)))
        elif direction == 'h':
            dl = self.Ly // (self.NSim + self.NSumo)
            for y in range(self.NSim):
                self.palce('sim', y, itertools.product(range(self.LSim), [y * dl]))
            for y in range(self.NSumo):
                self.place('sumo', y, itertools.product(range(self.LSumo), [(self.Ly-1-y)*dl]))
        
    def dense_init(self):
        """
        The initialization in the dense limit
        """
        raise Exception("Polymer density relatively high, requiring better init method.")
        
    def thermalize(self):
        """
        """
        pass
    
        
    def neighbor(self, point):
        """
        Get a set of neighbors of a point
        """
        x, y = point
        four_points = [(x+1, y), (x-1, y), (x, y+1), (x, y-1)]
        points = [(x%self.Lx, y%self.Ly) for x, y in four_points]
        return set(points)
        
        
        
class Polymer(object):
    def __init__(self, poly_id, locs):
        self.id = poly_id
        self.locs = locs                