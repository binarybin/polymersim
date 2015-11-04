import matplotlib
matplotlib.use("TkAgg")
from matplotlib import style
style.use("ggplot")
import matplotlib.animation as animation
from copy import copy

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
        endmove = EndMove(self.space, 0)
        snakemove = SnakeMove(self.space, 0)
        cornermove = CornerMove(self.space, 0)
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
    def animate(self, i):
        tempspace = copy(self.space.space)
        plt.imshow(tempspace)
        print i
        self.proceed()
        if i == 30:
            self.movelist[0].beta = 10
            self.movelist[1].beta = 10
            self.movelist[2].beta = 10
        
if __name__=='__main__':
    app = App(10, 10, 10, 20, 40, 30)
    fig = plt.figure()
    ani = animation.FuncAnimation(fig, app.animate, interval=0, frames=5000, repeat=False) 
    plt.show()
    #ani.save('movie5000.mp4', writer='ffmpeg', fps = 10)