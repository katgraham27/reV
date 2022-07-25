import numpy as np
from shapely.geometry import Polygon, MultiPolygon, Point, LineString, MultiPoint, MultiLineString
import shapely.geometry

if __name__ == "__main__":
    point = Point(1,2)
    if point.contains(Point(1,2)):
        print("asdf")

    m = MultiPoint([(0, 0), (1, 1), (1,2), (2,2)])

    if m.contains(Point(0.3,1)):
        print("poop")
    if m.contains(Point(1,1)):
        print("kat kat")


    #multiline = MultiLineString(LineString([(0, 0), (1, 1), (1,2), (2,2)]))
    multiline = MultiLineString([((0, 0), (1, 1), (1,2), (2,2))])
    print(multiline)
    if multiline.contains(Point(0,0)):
        print("tru")


    l = LineString([(0, 0), (1, 1), (1,2), (2,2)])
    print("l : ", l)
    if l.contains(Point(1,1)):
        print("linestring has da point")
    print(l.geom_type)
    
    l = l.difference(Point(0,0))
    print("l again: ", l)

    l = l.difference(Point(1,1))
    print("l again: ", l)


    list = len(l.coords)
    list1 = len(point.coords)
    print("list", list)
    print(list1)
    for i in range (list1):
        print("i ", i)
        print("point: ", point.coords[i])

    coordspoint = point.coords[:]
    print("cooords point : ", coordspoint)

    linepoints = l.coords[:]
    print("l point : ", linepoints)

    for i in range(len(l.geoms)):
        asdf = l.geoms[i].coords[:]
        print("asdf", asdf)
    # for i in list:
    #     print("point inn line: ", l.coords[i])
    # will not work:
    # for i in len(point):
    #     if point[i].contains(1,2):
    #         print("hi")

    # if (len(point.geoms)==0):
    #     print("hey")

    # multipoly = MultiPolygon(m)
    # print("multipoly.area: ", multipoly.area)


    # multipoly = MultiPolygon(l)