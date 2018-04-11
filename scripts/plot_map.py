
import sys
import numpy

import argparse

import matplotlib.pyplot as P
from mpl_toolkits.basemap import Basemap, shiftgrid

def load_points(filename):
    f = open(filename, 'r')
    lines = f.readlines()[1:]
    f.close()

    return zip(*map(lambda x: map(float, x.split())[:3], lines))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Input file')

    parser.add_argument('--flip', action = 'store_true', default = False, help = 'Vertically flip image')
    
    parser.add_argument('--clon', type = float, default = 0.0, help = 'Central Longitude')
    parser.add_argument('--cmap', type = str, default = 'hsv', help = 'Colour map')

    parser.add_argument('--vmin', type = float, default = 0.0, help = 'V Min')
    parser.add_argument('--vmax', type = float, default = 1.0, help = 'V Max ')

    parser.add_argument('--colorbar', action = 'store_true', default = False, help = 'Show colour bar')
    parser.add_argument('--clabel', type = str, default = None, help = 'Colour bar label')

    parser.add_argument('--pdf', type = str, default = None, help = 'PDF output')
    parser.add_argument('--png', type = str, default = None, help = 'PNG output')

    parser.add_argument('--noshow', action = 'store_true', default = False, help = 'No output')

    parser.add_argument('--bordercolor', type = str, default = 'black', help = 'Border colour')

    parser.add_argument('--robinson', action = 'store_true', default = False, help = 'Robinson projection')
    
    args = parser.parse_args()
    

    image = numpy.loadtxt(args.input)

    rows, cols = image.shape

    if args.flip:
        image = numpy.flipud(image)

    print numpy.min(image), numpy.max(image)
    
    lons = numpy.linspace(-180.0, 180.0, cols)
    lats = numpy.linspace(-90.0, 90.0, rows)

    fig, ax = P.subplots()
    fig.set_tight_layout(True)

    if args.robinson:
        m = Basemap(resolution='c', projection='robin', lon_0 = args.clon)
    else:
        m = Basemap(resolution='c', projection='moll', lon_0 = args.clon)

    x, y = numpy.meshgrid(lons, lats)
        
    p = m.pcolor(x, y, image, latlon = True, cmap = args.cmap, vmin = args.vmin, vmax = args.vmax)
        
    parallels = numpy.arange(-90.0, 90.0, 45.0)
    m.drawparallels(parallels, linewidth = 0.1, color = args.bordercolor)
    meridians = numpy.arange(-180.0, 180.0, 30.0)
    m.drawmeridians(meridians, linewidth = 0.1, color = args.bordercolor)
    m.drawcoastlines(linewidth = 0.25, color = args.bordercolor)

    if (args.colorbar):
        if args.clabel is None:
            P.colorbar(orientation='horizontal')
        else:
            P.colorbar(orientation='horizontal', label = args.clabel)

    if args.pdf:
        fig.savefig(args.pdf, format = 'PDF')

    if args.png:
        fig.savefig(args.png, format = 'PNG', dpi = 300)

    if not args.noshow:
        P.show()
