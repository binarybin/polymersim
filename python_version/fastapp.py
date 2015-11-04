#import matplotlib.pyplot as plt
from time import time
from random import choice

from space import Space
from endmove import EndMove
from snakemove import SnakeMove
from cornermove import CornerMove


class App(object):
    def __init__(self, NSim, LSim, NSumo, LSumo, Lx, Ly):
        self.space = Space(NSim, LSim, NSumo, LSumo, Lx, Ly)
        self.space.initialize()
        endmove = EndMove(self.space, 0)
        snakemove = SnakeMove(self.space, 0)
        cornermove = CornerMove(self.space, 0)
        self.movelist = [endmove, snakemove, cornermove]
        self.i = 0
    def move(self, one_move):
        if one_move.movetype == "EndMove":
            return one_move.move(choice(range(1, self.space.NSim+1)), choice(['sumo','sim']), choice(['head', 'tail']))
        elif one_move.movetype == "SnakeMove":
            return one_move.move(choice(range(1, self.space.NSim+1)), choice(['sumo','sim']), choice(['head', 'tail']))
        elif one_move.movetype == "CornerMove":
            return one_move.move(choice(range(1, self.space.NSim+1)), choice(['sumo','sim']))
        else:
            raise Exception("Move type undefined")
    def proceed(self):
        self.i += 1
        for one_move in self.movelist:
            self.move(one_move)
            if self.i == 900000:
                self.movelist[0].beta = 10
                self.movelist[1].beta = 10
                self.movelist[2].beta = 10
        
if __name__ == "__main__":
    app = App(30, 20, 30, 20, 120, 100)
   # fig = plt.figure()
    time1 = time()
    for step in range(100000):
        app.proceed()
    time2 = time()
    print "delta t: ", time2-time1, " seconds"
   # plt.imshow(app.space.space)    
   # plt.show()
