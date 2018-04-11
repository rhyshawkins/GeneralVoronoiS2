
import argparse
import numpy

def randomlonlat():

    x = numpy.random.normal()
    y = numpy.random.normal()
    z = numpy.random.normal()

    l = numpy.sqrt(x*x + y*y + z*z)
    
    phi = numpy.arccos(z/l)
    theta = numpy.arctan2(y, x)

    lon = theta * 180.0/numpy.pi
    lat = 90.0 - phi*180.0/numpy.pi

    return lon, lat

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-o', '--output', type = str, required = True, help = 'Output file')
    parser.add_argument('-n', '--npoints', type = int, default = 100, help = 'No. points')

    args = parser.parse_args()

    f = open(args.output, 'w')
    f.write('%d\n' % args.npoints)
    for i in range(args.npoints):

        lon, lat = randomlonlat()

        f.write('%16.9e %16.9e 0.0 0.0\n' % (lon, lat))

    f.close()

                
    

