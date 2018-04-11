
import sys

import numpy
import matplotlib.pyplot as P
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, action = 'append', help = 'Input file(s)')
    parser.add_argument('-p', '--pdf', type = str, default = None, help = 'Output PDF')

    parser.add_argument('-s', '--show', default = False, action = 'store_true', help = 'Show plot')
    parser.add_argument('--min', default = 0.0, type = float, help = 'Minimum likelihood')
    parser.add_argument('--max', default = 3000.0, type = float, help = 'Maximum likelihood')
    parser.add_argument('--bins', default = 300, type = int, help = 'Histogram bins')
    parser.add_argument('--burnin', default = 0, type = int, help = 'Iterations burnin')
    parser.add_argument('--thin', default = 10, type = int, help = 'Iterations thinning')

    parser.add_argument('--chi', default = None, type = float, help = 'Draw Chi squared line')
    parser.add_argument('--chains', default = 0, type = int, help = 'No. chains')
    
    parser.add_argument('--width', type = int, default = 8, help = 'Figure width (inches)')
    parser.add_argument('--height', type = int, default = 6, help = 'Figure height (inches)')
    parser.add_argument('-I', '--show-iterations', action = 'store_true', default = False, help = 'Show iteration counts')
    
    
    args = parser.parse_args()

    if args.input is None:
        print 'Require 1 or more likelihood files'
        sys.exit(-1)

    fig, ax = P.subplots()

    fig.set_size_inches(args.width, args.height)
    fig.set_tight_layout(True)
    
    divider = make_axes_locatable(ax)
    erry = divider.append_axes('right', 1.2, pad = 0.25, sharey = ax)
    erry.yaxis.tick_right()
    
    totalhist, totaledges = numpy.histogram([], bins = args.bins, range = (args.min, args.max))

    totalcentres = totaledges[:-1] + (totaledges[1] - totaledges[0])/2.0

    likelihood = None
    for i in args.input:
        like = numpy.loadtxt(i)

        rows = len(like)
        if likelihood is None:
            likelihood = like.reshape((rows,1))
        else:
            likelihood = numpy.hstack((likelihood, like.reshape((rows,1))))
            
    rows, cols = likelihood.shape

    for c in range(cols):

        likehist = likelihood[:,c]
        
        ax.plot(likehist[args.burnin::args.thin], linestyle = 'solid', alpha = 0.5)

        hist, edges = numpy.histogram(likehist, bins = args.bins, range = (args.min, args.max))

        normhist = numpy.array(hist, dtype = 'float', copy = True)/numpy.sum(hist)
        
        erry.plot(normhist, totalcentres, linestyle = 'solid', alpha = 0.5)

        totalhist = totalhist + hist

    normtotalhist = numpy.array(totalhist, dtype = 'float', copy = True)/numpy.sum(totalhist)

    #erry.fill_betweenx(totalcentres, normtotalhist, edgecolor = 'none', color = 'grey')
    erry.plot(normtotalhist, totalcentres, 'k-')
    
    ax.set_ylim(args.min, args.max)

    ax.set_ylabel('Negative Log Likelihood')
    ax.set_xlabel('Iterations')

    if not args.show_iterations:
        ax.set_xticks([])
    erry.set_xticks([])

    ax.set_xlim(0, len(likelihood[args.burnin::args.thin,0]))

    if args.chi:
        ax.axhline(y = args.chi, linestyle = 'dashed', color = 'red')
        erry.axhline(y = args.chi, linestyle = 'dashed', color = 'red')
        
    if args.show:
        P.show()

    if args.pdf:
        fig.savefig(args.pdf, format = 'PDF')

        
    

    

    
