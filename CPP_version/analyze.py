
# coding: utf-8

# In[1]:



# In[7]:
from __future__ import division
from pylab import *
import pylab
from sys import argv


# In[23]:

def get_layers(lines):
    for idx, line in enumerate(lines):
        if line == 'END\n':
            print idx
            size = idx - 3
            break
    matlines = lines[3:size+3]
    matstring = []
    for line in matlines:
        matstring.append(line.split('\t'))
    simlayer = pylab.zeros((size,size), int)
    sumolayer = pylab.zeros((size,size), int)
    for i in range(size):
        for j in range(size):
            xy = matstring[i][j].split(',')
            x = int(xy[0])
            y = int(xy[1])
            simlayer[i][j] = x
            sumolayer[i][j] = y
    return simlayer, sumolayer

def get_polys(lines):
    for idx, line in enumerate(lines):
        if line == 'Two dimensional single layer space\n':
            start = idx+1
        if line == 'No bond info for this geometry\n':
            end = idx - 3
            break

    polylines = lines[start:end]

    polyodd  = [pline for idx, pline in enumerate(polylines) if idx %2 == 0]
    polyeven = [pline for idx, pline in enumerate(polylines) if idx %2 != 0]

    polyinfo = zip(polyodd, polyeven)

    sims = []
    sumos = []

    for element in polyinfo:
        identifier = element[0].split()[0]
        points = element[1].split()
        pointinfo = [info.split('(')[1].split(')')[0].split(',') for info in points]
        pointfull = [(int(info[0]), int(info[1])) for info in pointinfo]
        if identifier == 'Sim':
            sims.append(pointfull)
        elif identifier == 'Sumo':
            sumos.append(pointfull)
        else:
            raise Exception("Something wrong with the type")
    return sims, sumos



def drawline(axs, pt1, pt2):
    lw = 3
    cl = 'r'
    y1, x1 = pt1
    y2, x2 = pt2
    y1, y2 = min(y1, y2), max(y1, y2)
    x1, x2 = min(x1, x2), max(x1, x2)
    if x1 == x2:
        if y2-y1 == 1:
            axs.vlines(x1, y1, y2, linewidth=lw, color = cl)
        else:
            axs.vlines(x1, y1, y1-0.5, linewidth=lw, color = cl)
            axs.vlines(x1, y2, y2+0.5, linewidth=lw, color = cl)
    elif y1 == y2:
        if x2-x1 == 1:
            axs.hlines(y1, x1, x2, linewidth=lw, color = cl)
        else:
            axs.hlines(y1, x1, x1-0.5, linewidth=lw, color = cl)
            axs.hlines(y1, x2, x2+0.5, linewidth=lw, color = cl)

def drawframe(axs, pts):
    lw = 2
    cl = 'w'
    xs = pylab.array([y for x,y in pts])
    ys = pylab.array([x for x,y in pts])
    if max(ys) - min(ys) in [1,3]:
        if max(xs)-min(xs) in [1,3]:
            axs.hlines(min(ys)-0.5, min(xs)-0.5, max(xs)+0.5, linewidth=lw, color = cl)
            axs.hlines(max(ys)+0.5, min(xs)-0.5, max(xs)+0.5, linewidth=lw, color = cl)
            axs.vlines(min(xs)-0.5, min(ys)-0.5, max(ys)+0.5, linewidth=lw, color = cl)
            axs.vlines(max(xs)+0.5, min(ys)-0.5, max(ys)+0.5, linewidth=lw, color = cl)
        else:
            maxx1 = max(xs[xs<10])
            minx1 = min(xs)
            maxx2 = max(xs)
            minx2 = min(xs[xs>10])
            axs.vlines(minx2-0.5, min(ys)-0.5, max(ys)+0.5, linewidth=lw, color = cl)
            axs.vlines(maxx1+0.5, min(ys)-0.5, max(ys)+0.5, linewidth=lw, color = cl)
            axs.hlines(min(ys)-0.5, minx1-0.5, maxx1+0.5, linewidth=lw, color = cl)
            axs.hlines(max(ys)+0.5, minx1-0.5, maxx1+0.5, linewidth=lw, color = cl)
            axs.hlines(min(ys)-0.5, minx2-0.5, maxx2+0.5, linewidth=lw, color = cl)
            axs.hlines(max(ys)+0.5, minx2-0.5, maxx2+0.5, linewidth=lw, color = cl)
    elif max(xs)-min(xs) in [1,3]:
        maxy1 = max(ys[ys<10])
        miny1 = min(ys)
        maxy2 = max(ys)
        miny2 = min(ys[ys>10])
        axs.hlines(miny2-0.5, min(xs)-0.5, max(xs)+0.5, linewidth=lw, color = cl)
        axs.hlines(maxy1+0.5, min(xs)-0.5, max(xs)+0.5, linewidth=lw, color = cl)
        axs.vlines(min(xs)-0.5, miny1-0.5, maxy1+0.5, linewidth=lw, color = cl)
        axs.vlines(max(xs)+0.5, miny1-0.5, maxy1+0.5, linewidth=lw, color = cl)
        axs.vlines(min(xs)-0.5, miny2-0.5, maxy2+0.5, linewidth=lw, color = cl)
        axs.vlines(max(xs)+0.5, miny2-0.5, maxy2+0.5, linewidth=lw, color = cl)
    else:
        maxy1 = max(ys[ys<10])
        miny1 = min(ys)
        maxy2 = max(ys)
        miny2 = min(ys[ys>10])
        maxx1 = max(xs[xs<10])
        minx1 = min(xs)
        maxx2 = max(xs)
        minx2 = min(xs[xs>10])
        axs.vlines(minx2-0.5, miny1-0.5, maxy1+0.5, linewidth=lw, color = cl)
        axs.vlines(maxx1+0.5, miny1-0.5, maxy1+0.5, linewidth=lw, color = cl)
        axs.vlines(minx2-0.5, miny2-0.5, maxy2+0.5, linewidth=lw, color = cl)
        axs.vlines(maxx1+0.5, miny2-0.5, maxy2+0.5, linewidth=lw, color = cl)
        axs.hlines(miny2-0.5, minx1-0.5, maxx1+0.5, linewidth=lw, color = cl)
        axs.hlines(maxy1+0.5, minx1-0.5, maxx1+0.5, linewidth=lw, color = cl)
        axs.hlines(miny2-0.5, minx2-0.5, maxx2+0.5, linewidth=lw, color = cl)
        axs.hlines(maxy1+0.5, minx2-0.5, maxx2+0.5, linewidth=lw, color = cl)
    

def plotsim(axs, sims):
    for pts in sims:
        for i in range(len(pts)-1):
            drawline(axs, pts[i], pts[i+1])

def plotsumo(axs, sumos):
    for pts in sumos:
        drawframe(axs, pts)

def make_plot(sumolayer, sims, sumos, prop):
    fig, axs = plt.subplots(1,1,figsize=(10,10))
    axs.imshow(sumolayer, clim=(-2, 2), interpolation='none')
    fig.canvas.draw()
    plotsumo(axs, sumos)
    plotsim(axs, sims)
    title(prop)
#    fig.show()
    savefig(prop + ".png")
    close()
def analyze(lines, prop):
    simlayer, sumolayer = get_layers(lines)
    sims, sumos = get_polys(lines)
    make_plot(sumolayer, sims, sumos, prop)

    
    
    


# In[27]:

fff = open("filelist.txt", 'r')
for filename in fff.readlines():
    try: 
        f = open("dragres/"+filename[:-1], 'r')
        lines = f.readlines()
        propraw = filename
        nepyc = propraw.split("nsim")[1].split('_')[0]
        nrubi = propraw.split("nsumo")[1].split('_')[0]
        lepyc = propraw.split("lsim")[1].split('_')[0]
        lrubi = propraw.split("lsumo")[1].split('_')[0]
        beta = propraw.split("beta_")[1].split('_')[0]
        gamma = propraw.split("gamma_")[1].split(".txt")[0]
        prop = "NEPYC_"+nepyc+"_NRUBI_"+nrubi+"_LEPYC_"+lepyc+"_LRUBY_"+lrubi+"_beta_"+beta+"_gamma_"+gamma
        analyze(lines, prop)
    except e:
        print filename+"failed"






