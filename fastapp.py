import matplotlib.pyplot as plt

from random import choice

from space import Space
from endmove import EndMove
from snakemove import SnakeMove
from cornermove import CornerMove


class App(object):
    def __init__(self, NSim, LSim, NSumo, LSumo, Lx, Ly):
        self.space = Space(NSim, LSim, NSumo, LSumo, Lx, Ly)
        self.space.initialize()
        endmove = EndMove(self.space)
        snakemove = SnakeMove(self.space)
        cornermove = CornerMove(self.space)
        self.movelist = [endmove, snakemove, cornermove]
        
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
        for one_move in self.movelist:
            self.move(one_move)
        
if __name__=='__main__':
    app = App(10, 10, 10, 20, 40, 30)
    fig = plt.figure()
    for step in range(30000):
        app.proceed()
    plt.imshow(app.space.space)    
    plt.show()