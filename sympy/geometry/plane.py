"""Geometrical Planes.

Contains
========
Plane

"""
from __future__ import print_function, division

from sympy.core import S, C, sympify, Dummy, nan, Eq, symbols
from sympy.functions.elementary.trigonometric import _pi_coeff as pi_coeff, \
    sqrt
from sympy.core.logic import fuzzy_and
from sympy.core.exprtools import factor_terms
from sympy.simplify.simplify import simplify
from sympy.solvers import solve
from sympy.geometry.exceptions import GeometryError
from .entity import GeometryEntity
from .point3d import Point3D
from .line3d import Line3D, Segment3D, Ray3D
from .util import _symbol

class Plane(GeometryEntity):
    def __new__(cls, p1, pt1=None, pt2=None, normal=[], **kwargs):
        p1 = Point3D(p1)
        if pt1 and pt2 is not None and len(normal) == 0:
            a = [p1.args, pt1, pt2]
            if len(set(a)) == 3:           
                try:
                    pt1 = Point3D(pt1)
                    pt2 = Point3D(pt2)
                except NotImplementedError:
                    raise ValueError('Enter 3 valid points or a single point with' 
                    'the normal vector of the desired plane')
            else:
                raise ValueError('Enter 3 distinct points')
        elif len(normal) == 3 and pt1 is not None and pt2 is None:
            pt1 = Point3D(pt1)
            a = p1.direction_ratio(pt1)
            b = sum([i*j for i, j in zip(a, normal)])
            if b == 0:
                pt2 = Point3D(pt1.x + a[0], pt1.y + a[1], pt1.z + a[2])
            else:
                raise ValueError('The plane you are trying to generate is not'
                'valid as the normal vector will not be normal to the plane'
                'specified by those 2 points')
        elif len(normal) == 3 and pt1 is None and pt2 is None:
            x, y, z = symbols('x y z')
            a = Point3D(x, y, z)
            b = p1.direction_ratio(a)
            c = sum(i*j for i, j in zip(b, normal))
            d = solve(Eq(c).subs([(y, 0), (z, 0)]))
            e = solve(Eq(c).subs([(x, 0), (z, 0)]))
            f = solve(Eq(c).subs([(x, 0), (y, 0)]))
            if d != []:
                pt1 = Point3D(d[0], 0, 0)
                pt2 = Point3D(pt1.x + p1.direction_ratio(pt1)[0], pt1.y
                + p1.direction_ratio(pt1)[1], pt1.z + p1.direction_ratio(pt1)[2])
                return GeometryEntity.__new__(cls, p1, pt1, pt2, **kwargs)
            elif e != []:
                pt1 = Point3D(0, d[0], 0)
                pt2 = Point3D(pt1.x + p1.direction_ratio(pt1)[0], pt1.y
                + p1.direction_ratio(pt1)[1], pt1.z + p1.direction_ratio(pt1)[2])
                return GeometryEntity.__new__(cls, p1, pt1, pt2, **kwargs)
            elif f != []:
                pt1 = Point3D(0, 0, d[0])
                pt2 = Point3D(pt1.x + p1.direction_ratio(pt1)[0], pt1.y
                + p1.direction_ratio(pt1)[1], pt1.z + p1.direction_ratio(pt1)[2])
                return GeometryEntity.__new__(cls, p1, pt1, pt2, **kwargs)
        else:
            raise ValueError('A 3rd Point or keyword "normal" must'
            'be used.')
        return GeometryEntity.__new__(cls, p1, pt1, pt2, **kwargs)