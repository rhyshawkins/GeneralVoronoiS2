
def parse(lines):

    if len(lines) > 0:
        tstar, sigma, count = map(float, lines[0].split())
        count = int(count)

        points = map(lambda x: map(float, x.split()), lines[1:count + 1])
        if len(points) != count:
            raise Exception("Points mismatch %d %d" % (len(points), count))
        
        rlines = lines[1 + count:]

        return [(tstar, sigma, points)] + parse(rlines)
    else:
        return []

def save(filename, data):

    f = open(filename, 'w')
    f.write('%d\n' % len(data))
    for (tstar, sigma, points) in data:

        f.write('%15.9f %15.9f %d\n' % (tstar, sigma, len(points)))

        for lon, lat, r in points:
            f.write('%15.9f %15.9f %15.9f\n' % (lon, lat, r))

    f.close()
    
    
f = open('coreattenuation.txt', 'r')
lines = f.readlines()
f.close()

data = parse(lines)

save('coreattenuation_converted.txt', data)

