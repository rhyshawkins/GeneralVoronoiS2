
import argparse

import numpy

def load_data(filename):

    #
    # Format is event lon, event lat, station lon, station lat
    #
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    data = []
    for i in range(n):

        elat, elon, slat, slon = map(float, lines[i].split())

        data.append((elon, elat, slon, slat))

    return data

def gc_angle(lon1, lat1, lon2, lat2):

    phi1 = (90.0 - lat1) * numpy.pi/180.0
    theta1 = lon1 * numpy.pi/180.0

    phi2 = (90.0 - lat2) * numpy.pi/180.0
    theta2 = lon2 * numpy.pi/180.0

    hsinphi = numpy.sin((phi1 - phi2)/2.0);
    hsintheta = numpy.sin((theta1 - theta2)/2.0);
      
    return 2.0 * numpy.arcsin(numpy.sqrt(hsinphi * hsinphi + numpy.sin(phi1)*numpy.sin(phi2)*(hsintheta*hsintheta)));

def gc_distancekm(lon1, lat1, lon2, lat2, radius):

    return radius * gc_angle(lon1, lat1, lon2, lat2)

def e_distance(entry, radius):
    elon, elat, slon, slat, pv, err = entry

    return gc_distancekm(elon, elat, slon, slat, radius)
    
def gc_midpoint(lon1, lat1, lon2, lat2):

    phi1 = (90.0 - lat1) * numpy.pi/180.0
    theta1 = lon1 * numpy.pi/180.0

    phi2 = (90.0 - lat2) * numpy.pi/180.0
    theta2 = lon2 * numpy.pi/180.0

    dtheta = theta2 - theta1;
    bx = numpy.sin(phi2) * numpy.cos(dtheta);
    by = numpy.sin(phi2) * numpy.sin(dtheta);
    delta = numpy.sin(phi1) + bx;
    
    mlat = 180.0/numpy.pi * numpy.arctan2(numpy.cos(phi1) + numpy.cos(phi2), numpy.sqrt(delta*delta + by*by))
    mlon = 180.0/numpy.pi * (theta1 + numpy.arctan2(by, numpy.sin(phi1) + bx))

    return mlon, mlat

def subdividepath(points):

    newpoints = []

    for i in range(len(points) - 1):

        lon1, lat1 = points[i]
        lon2, lat2 = points[i + 1]
        
        mlon, mlat = gc_midpoint(lon1, lat1, lon2, lat2)

        newpoints.extend([(lon1, lat1), (lonwrap(mlon), mlat)])

    newpoints.append(points[-1])

    return newpoints

def mkpath(elon, elat, slon, slat, radius, samplelength):

    path = [(elon, elat), (slon, slat)]

    while gc_distancekm(path[0][0], path[0][1], path[1][0], path[1][1], radius) > samplelength:

        path = subdividepath(path)

    return path

def lonwrap(lon):
    if (lon < -180.0):
        return lonwrap(lon + 360.0)
    elif (lon > 180.0):
        return lonwrap(lon - 360.0)
    else:
        return lon
    
def randomlonlat():

    x = numpy.random.normal()
    y = numpy.random.normal()
    z = numpy.random.normal()

    l = numpy.sqrt(x*x + y*y + z*z)
    
    phi = numpy.arccos(z/l)
    theta = numpy.arctan2(y, x)

    lon = theta * 180.0/numpy.pi
    lat = 90.0 - phi*180.0/numpy.pi

    #
    # Ensure lon in range -180 .. 180
    #
    lon = lonwrap(lon)

    return lon, lat

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-n', '--npaths', type = int, default = 100, help = 'No. random paths')
    
    parser.add_argument('-i', '--input', type = str, default = None, help = 'Input file')
    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')

    parser.add_argument('-r', '--radius', type = float, default = 6371.0, help = 'Earth radius (km)')
    parser.add_argument('-d', '--maxdistance', type = float, default = 100.0, help = 'Max distance between points')
    
    args = parser.parse_args()

    if (args.input is None):

        data = []
        for i in range(args.npaths):
            elon, elat = randomlonlat()
            slon, slat = randomlonlat()
            data.append((elon, elat, slon, slat))

    else:
        data = load_data(args.input)

    f = open(args.output, 'w')

    f.write('%d\n' % len(data))
            
    for d in data:

        elon, elat, slon, slat = d

        path = mkpath(elon, elat, slon, slat, args.radius, args.maxdistance)
        
        f.write('0.0 0.0 %d\n' % len(path))

        for lon, lat in path:

            f.write('%15.9f %15.9f %10.6f\n' % (lon, lat, args.radius))

    f.close()
        
    
    
