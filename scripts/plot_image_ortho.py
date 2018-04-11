
import sys
import numpy
import math

from optparse import OptionParser
import matplotlib.pyplot as P
from mpl_toolkits.basemap import Basemap, shiftgrid

import matplotlib.colors
import matplotlib

def load_points(filename):
    f = open(filename, 'r')
    lines = f.readlines()[1:]
    f.close()

    return zip(*map(lambda x: map(float, x.split())[:3], lines))

def load_data(filename):

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    data = []

    while len(lines) > 0:

        t = lines[0].split()
        tstar = float(t[0])
        noise = float(t[1])
        npoints = int(t[2])

        points = map(lambda x: map(float, x.split()), lines[1:npoints])

        data.append((tstar, noise, points))
        
        lines = lines[npoints + 1:]

    return data

def sphtocart(point):
    lon, lat, r = point
    theta = lon * math.pi/180.0
    phi = lat * math.pi/180.0
    return (r * math.cos(theta) * math.cos(phi),
            r * math.sin(theta) * math.cos(phi),
            r * math.sin(phi))

def distance(a, b):

    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]

    return math.sqrt(dx*dx + dy*dy + dz*dz)

def totaldistance(points):

    cartesian = map(sphtocart, points)
    dist = 0.0
    for i in range(1, len(cartesian)):
        dist = dist + distance(cartesian[i], cartesian[i - 1])

    return dist

def bottomingpoint(points):

    minradius = 6371.0
    minpoint = None
    for (lon, lat, radius) in points:
        if radius < minradius:
            minpoint = (lon, lat, radius)
            minradius = radius

    return minpoint

def approxQ(tstar, points):
    
    return totaldistance(points)/(tstar * 11.0)
    
def draw_bottoming_points(data, m, minlon, maxlon, minlat, maxlat):

    if data:

        lons = []
        lats = []
        approxQs = []
        for tstar, noise, points in data:

            lon, lat, radius = bottomingpoint(points)
            if (lon >= minlon and lon <= maxlon) and (lat >= minlat and lat <= maxlat):
                lons.append(lon)
                lats.append(lat)

                approxQs.append(approxQ(tstar, points))


        print len(lons), ' points'
        x, y = m(lons,lats)
        m.scatter(x, y, s = 40, c = approxQs, cmap = cmap, edgecolor = 'green')

        print 'Q', numpy.min(approxQs), numpy.max(approxQs)

    else:
        print 'draw_bottoming_points: no data'

def draw_paths(data, m, minlon, maxlon, minlat, maxlat):
    
    if data:

        for tstar, noise, points in data:


            lon0, lat0, r0 = points[0]
            lon1, lat1, r1 = points[-1]

            draw = False

            if (maxlon > 180.0):
                tmaxlon = maxlon - 360.0
                draw = (((lon0 >= minlon and lon0 <= 180.0) or (lon0 > -180.0 and lon0 < tmaxlon)) and
                        ((lon1 >= minlon and lon1 <= 180.0) or (lon0 > -180.0 and lon0 < tmaxlon)) and
                        (lat0 >= minlat and lat0 <= maxlat) and
                        (lat1 >= minlat and lat1 <= maxlat))

            else:
                draw = ((lon0 >= minlon and lon0 <= maxlon) and
                        (lon1 >= minlon and lon1 <= maxlon) and
                        (lat0 >= minlat and lat0 <= maxlat) and
                        (lat1 >= minlat and lat1 <= maxlat))
                
            if (draw):
                m.drawgreatcircle(lon0, lat0, lon1, lat1, color = 'green', alpha = 0.2)

    else:
        print 'draw_bottoming_points: no data'


if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option('-f', '--file', dest='file', default=None, help = 'Input file')
    parser.add_option('--fake', dest = 'fake', default = False, action = 'store_true', help = 'Fake image')
    
    parser.add_option('-p', '--pdf', dest='pdf', default=None, help = 'PDF file to write')
    parser.add_option('--png', dest = 'png', default = None, help = 'PNG file to write')
    
    parser.add_option('-n', '--noshow', dest='show', default=True, action='store_false', help = 'Do not show plot')

    parser.add_option('--vmin', dest='vmin', default=0.0, type='float', help = 'Value min')
    parser.add_option('--vmax', dest='vmax', default=1000.0, type='float', help = 'Value max')
    parser.add_option('--log', dest='log', default = False, action='store_true', help = 'Log of field')
    parser.add_option('--flip', dest='flip', default = False, action='store_true', help = 'Flip image vertically')
    parser.add_option('--auto', dest='auto', default=False, action='store_true', help = 'Use auto range')

    parser.add_option('--gray', dest='gray', action='store_true', default=False, help = 'Grayscale colormap')
    parser.add_option('--puor', dest='puor', action='store_true', default=False, help = 'PuOr colormap')
    parser.add_option('--rdbu', dest='rdbu', action='store_true', default=False, help = 'RdBu colormap')
    parser.add_option('--hot', dest='hot', action='store_true', default=False, help = 'Hot colormap')
    parser.add_option('--cool', dest='cool', action='store_true', default=False, help = 'Cool colormap')
    parser.add_option('--viridis', dest='viridis', action='store_true', default=False, help = 'Viridis colormap')
    parser.add_option('--magma', dest='magma', action='store_true', default=False, help = 'Magma colormap')
    parser.add_option('--reverse', dest='reverse', action='store_true', default=False, help = 'Reverse Colormap')
    
    parser.add_option('--lonsamples', dest = 'lonsamples', type = 'int', default = 1024, help = 'Lon samples')
    parser.add_option('--latsamples', dest = 'latsamples', type = 'int', default = 1024, help = 'Lat samples')

    parser.add_option('--points', dest = 'points', default = None, help = 'Points to plot')
    parser.add_option('--triangles', dest = 'triangles', default = None, help = 'Triangles to plot')
    parser.add_option('--data', dest = 'data', default = None, help = 'Data paths to plot as points or paths')

    parser.add_option('--bottom', dest = 'bottom', default = False, action = 'store_true', help = 'Plot bottoming points')

    parser.add_option('--paths', dest = 'paths', default = False, action = 'store_true', help = 'Plot data paths')

    parser.add_option('--shift', dest = 'shift', default = False, action = 'store_true', help = 'Shift image 180 degrees in longitude')

    parser.add_option('--nocolorbar', dest = 'colorbar', default = True, action = 'store_false', help = 'Disable colorbar')

    parser.add_option('--title', dest = 'title', default = None, type = 'str', help = 'Figure title')
    parser.add_option('--label', dest = 'label', default = 'Q', type = 'str', help = 'Colorbar Label')
    
    options, args = parser.parse_args()


    t = None
    
    if options.fake:

        flons = numpy.linspace(-180.0, 180.0, 129)[:-1] + 180.0/128.0
        flats = numpy.linspace(90.0, -90.0, 65)[:-1] + 90.0/64.0

        print flons[0], flons[-1]

        fflon, fflats = numpy.meshgrid(flons, flats)

        t = 500.0*numpy.exp(-((fflon - 140.0)**2 + (fflats + 35.0)**2)/(2.0*10**2))
    else:
        if options.file == None:
            print 'Missing required argument: file'
            sys.exit(-1)

        t = numpy.loadtxt(options.file)

    if options.log:
        t = numpy.log(t)

    if options.flip:
        #
        # Image is north pole to south pole but base map requires reverse in general
        #
        t = numpy.flipud(t)
        
    rows, cols = t.shape

    if options.auto:
        vmin = numpy.min(t)
        vmax = numpy.max(t)
    else:
        vmin = options.vmin
        vmax = options.vmax

    print numpy.min(t), numpy.max(t)
    print vmin, vmax

    #lonhw = 360.0/(float(cols) * 2.0)
    #lons = numpy.linspace(-180.0, 180.0, cols + 1)[:-1] + lonhw

    #lathw = 180.0/(float(rows) * 2.0)
    #lats = numpy.linspace(-90.0, 90.0, rows + 1)[:-1] + lathw

    if options.shift:
        lons = numpy.linspace(0.0, 360.0, cols)
        t, lons = shiftgrid(180.0, t, lons, start = False)
    else:
        lons = numpy.linspace(-180.0, 180.0, cols)
    #lats = numpy.linspace(-90.0, 90.0, rows)
    lats = numpy.linspace(90.0, -90.0, rows)

    fig, ax = P.subplots(3, 4)
    fig.set_size_inches((8,6))
    P.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    
    if options.title:
        fig.suptitle(options.title)
        
    #P.tight_layout()

    print lons[0], lons[-1]
    print lats[0], lats[-1]
    x, y = numpy.meshgrid(lons, lats)
    cmap = None
    if options.gray:
        #cmap = P.cm.get_cmap('gray')
        cmap = P.cm.get_cmap('binary')
    elif options.puor:
        cmap = P.cm.get_cmap('PuOr')
    elif options.rdbu:
        cmap = P.cm.get_cmap('RdBu')
    elif options.hot:
        cmap = P.cm.get_cmap('hot')
    elif options.cool:
        cmap = P.cm.get_cmap('cool')
    elif options.viridis:
        cmap = P.cm.get_cmap('viridis')
    elif options.magma:
        cmap = P.cm.get_cmap('magma')

    if not cmap is None and options.reverse:
        cmap = matplotlib.colors.ListedColormap(cmap.colors[::-1], cmap.name + '_r') 

    maps = []
    mappable = []
    for r in range(3):
        for c in range(4):
            ax[r][c].axis('off')

    data = None
    if options.data:
        data = load_data(options.data)

    #
    # North Pole
    #
    m = Basemap(ax = ax[0][1], resolution='c', projection='ortho', lat_0=90.0, lon_0 = 0.0)
    p = m.pcolor(x, y, t, latlon = True, cmap = cmap, vmin = vmin, vmax = vmax)

    if options.bottom:
        draw_bottoming_points(data, m, -180.0, 180.0, 0.0, 90.0)

    if options.paths:
        draw_paths(data, m, -180.0, 180.0, 0.0, 90.0)
        
    m.drawcoastlines()
    maps.append(m)
    mappable.append(p)

    #
    # Centre
    #
    m = Basemap(ax = ax[1][1], resolution='c', projection='ortho', lat_0=0.0, lon_0 = 0.0)
    p = m.pcolor(x, y, t, latlon = True, cmap = cmap, vmin = vmin, vmax = vmax)

    if options.bottom:
        draw_bottoming_points(data, m, -90.0, 90.0, -90.0, 90.0)

    if options.paths:
        draw_paths(data, m, -90.0, 90.0, -90.0, 90.0)
        
    m.drawcoastlines()
    maps.append(m)
    mappable.append(p)

    #
    # South Pole
    #
    m = Basemap(ax = ax[2][1], resolution='c', projection='ortho', lat_0=-90.0, lon_0 = 0.0)
    p = m.pcolor(x, y, t, latlon = True, cmap = cmap, vmin = vmin, vmax = vmax)

    if options.bottom:
        draw_bottoming_points(data, m, -180.0, 180.0, -90.0, 0.0)

    if options.paths:
        draw_paths(data, m, -180.0, 180.0, -90.0, 0.0)
        
    m.drawcoastlines()
    maps.append(m)
    mappable.append(p)

    #
    # East
    #
    m = Basemap(ax = ax[1][0], resolution='c', projection='ortho', lat_0=0.0, lon_0 = -90.0)
    p = m.pcolor(x, y, t, latlon = True, cmap = cmap, vmin = vmin, vmax = vmax)

    if options.bottom:
        draw_bottoming_points(data, m, -180.0, 0.0, -90.0, 90.0)

    if options.paths:
        draw_paths(data, m, -180.0, 0.0, -90.0, 90.0)
        
    m.drawcoastlines()
    maps.append(m)
    mappable.append(p)

    #
    # West
    #
    m = Basemap(ax = ax[1][2], resolution='c', projection='ortho', lat_0=0.0, lon_0 = 90.0)
    p = m.pcolor(x, y, t, latlon = True, cmap = cmap, vmin = vmin, vmax = vmax)

    if options.bottom:
        draw_bottoming_points(data, m, 0.0, 180.0, -90.0, 90.0)

    if options.paths:
        draw_paths(data, m, 0.0, 180.0, -90.0, 90.0)
        
    m.drawcoastlines()
    maps.append(m)
    mappable.append(p)

    #
    # Far West
    #
    m = Basemap(ax = ax[1][3], resolution='c', projection='ortho', lat_0=0.0, lon_0 = 180.0)
    p = m.pcolor(x, y, t, latlon = True, cmap = cmap, vmin = vmin, vmax = vmax)

    if options.bottom:
        draw_bottoming_points(data, m, 90.0, 270.0, -90.0, 90.0)

    if options.paths:
        draw_paths(data, m, 90.0, 270.0, -90.0, 90.0)
        
    m.drawcoastlines()
    maps.append(m)
    mappable.append(p)
    

    if (options.colorbar):
        norm = matplotlib.colors.Normalize(vmin = vmin, vmax = vmax)
        cbar = fig.colorbar(mappable[0],
                            cax = ax[2][3],
                            orientation = 'horizontal')
        ax[2][3].axis('on')
        bbox = ax[2][2].get_position()
        yc = (bbox.y0 + bbox.y1)/2.0
        yh = (bbox.y1 - bbox.y0)/2.0
        ax[2][2].set_position([bbox.x0, yc - yh/2.0, bbox.x1 - bbox.x0, yh])
        
        ticker = matplotlib.ticker.MaxNLocator(nbins=3)
        cbar.locator = ticker
        cbar.update_ticks()

        cbar.set_label(options.label)
        
    if (options.points):
        lons, lats, heat = load_points(options.points)
        for m in maps:
            x, y = m(lons,lats)
            m.scatter(x, y, s = 25, c = heat, alpha = 1.0, cmap = cmap)

    if options.pdf:
        fig.savefig(options.pdf, format = 'PDF')
        

    if options.png:
        fig.savefig(options.png, format = 'PNG')

    if options.show:
        P.show()
