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
from .point import Point
from .line3d import Line3D, Segment3D, Ray3D
from sympy.matrices import Matrix

class Plane(GeometryEntity):
    """
    A plane is a flat, two-dimensional surface. A plane is the two-dimensional 
    analogue of a point (zero-dimensions), a line (one-dimension) and a solid 
    (three-dimensions).A plane can generally be constructed by two types of 
    inputs.They are three non collinear points and a point and the plane's
    normal vector.
        
    Attributes
    ==========

    p1
    normal_vector
    
    Examples
    ========
        
    >>> from sympy import Plane, Point3D
    >>> from sympy.abc import x
    >>> Plane(Point3D(1, 1, 1), Point3D(2, 3, 4), Point3D(2, 2, 2))
    Plane(Point3D(1, 1, 1), [-1, 2, -1])
    >>> Plane((1, 1, 1), (2, 3, 4), (2, 2, 2))
    Plane(Point3D(1, 1, 1), [-1, 2, -1])
    >>> Plane(Point3D(1, 1, 1), normal_vector=[1,4,7])
    Plane(Point3D(1, 1, 1), [1, 4, 7])
    
    """    
    def __new__(cls, p1, pt1=None, pt2=None, normal_vector=[], **kwargs):
        p1 = Point3D(p1)
        if pt1 is None and pt2 is None and len(normal_vector) == 3:
            normal_vector = normal_vector
        elif pt1 and pt2 is not None and len(normal_vector) == 0:
            pt1, pt2 = Point3D(pt1), Point3D(pt2)
            if Point3D.is_collinear(p1, pt1, pt2):
                raise NotImplementedError('Enter three non-collinear points')
            a = p1.direction_ratio(pt1)
            b = p1.direction_ratio(pt2)
            normal_vector = list(Matrix(a).cross(Matrix(b)))
        else:
            raise ValueError('Either provide 3 3D points or a point with a'
            'normal vector')
        return GeometryEntity.__new__(cls, p1, normal_vector, **kwargs)
        
    @property
    def p1(self):
        """The only defining point of the plane.Others can be obtained from the
        arbitrary_point method.

        See Also
        ========

        sympy.geometry.point3d.Point3D
        
        Examples
        ========
        
        >>> from sympy import Point3D, Plane        
        >>> a = Plane(Point3D(1, 1, 1), Point3D(2, 3, 4), Point3D(2, 2, 2))
        >>> a.p1
        Point3D(1, 1, 1)
        
        """
        return self.args[0]
        
    @property
    def normal_vector(self):
        """Normal vector of the given plane.
        
        Examples
        ========
        
        >>> from sympy import Point3D, Plane        
        >>> a = Plane(Point3D(1, 1, 1), Point3D(2, 3, 4), Point3D(2, 2, 2))
        >>> a.normal_vector
        [-1, 2, -1]
        >>> a = Plane(Point3D(1, 1, 1), normal_vector=[1, 4, 7])
        >>> a.normal_vector
        [1, 4, 7]
        
        """
        return self.args[1]
        
    def equation(self, x='x', y='y', z='z'):
        """The equation of the Plane.
        
        Examples
        ========

        >>> from sympy import Point3D, Plane
        >>> a = Plane(Point3D(1, 1, 2), Point3D(2, 4, 7), Point3D(3, 5, 1))
        >>> a.equation()
        -23*x + 11*y - 2*z + 16        
        >>> a = Plane(Point3D(1, 4, 2), normal_vector=[6, 6, 6])
        >>> a.equation()
        6*x + 6*y + 6*z - 42
        
        """        
        x, y, z = symbols("x, y, z")
        a = Point3D(x, y, z)
        b = self.p1.direction_ratio(a)
        c = self.normal_vector
        return (sum(i*j for i, j in zip(b, c)))
    
    def projection(self, pt):
        x, y, z = symbols("x, y, z")
        k = self.equation(x, y, z)
        const = [i for i in k.args if i.is_constant() is True]
        a, b, c = k.coeff(x), k.coeff(y), k.coeff(z)
        d = -(a*pt.args[0] + b*pt.args[1] + const.pop())
        t = d/(a**2 + b**2 + c**2)
        if isinstance(pt, Point):
            return Point3D(pt.args[0] + t*a, pt.args[1] + t*b, t*c)
        if isinstance(pt, Point3D):
            return Point3D(pt.args[0] + t*a, pt.args[1] + t*b, pt.args[2] + t*c)
            
    def is_parallel(self, l):
        """Is the given geometric entity parallel to the plane?
        """
        if isinstance(l, Line3D):
            a = l.direction_ratio
            b = self.normal_vector
            c = sum([i*j for i, j in zip(a, b)])
            if c == 0:
                return True
            else:
                return False
        elif isinstance(l, Plane):
            a = Matrix(l.normal_vector)
            b = Matrix(self.normal_vector)
            if sum(list(a.cross(b))) == 0:
                return True
            else:
                return False
                
    def is_perpendicular(self, l):
        """is the given geometric entity perpendicualar to the given plane?
        """
        if isinstance(l, Line3D):        
            a = Matrix(l.direction_ratio)
            b = Matrix(self.normal_vector)
            if sum(list(a.cross(b))) == 0:
                return True
            else:
                return False
        elif isinstance(l, Plane):
           a = Matrix(l.normal_vector)
           b = Matrix(self.normal_vector)
           if a.dot(b) == 0:
               return True
           else:
               return False
        else:
            return False
            
    def perpendicular_line(self, pt):
        """A line perpendicular to the given plane.
        
        Parameters
        ==========
        
        pt: Point3D
        
        Examples
        ========
        >>> a = Plane(Point3D(1,4,6), normal_vector=[2, 4, 6])
        >>> a.perpendicular_line(Point3D(9,8,7))
        Line3D(Point3D(9, 8, 7), Point3D(11, 12, 13))
        
        """
        a = self.normal_vector
        return Line3D(pt, direction_ratio=a)

    def parallel_plane(self, pt):
        """
        Plane parallel to the given plane and passing through the point pt.
        
        Parameters
        ==========
        
        pt: Point3D
        
        Examples
        ========
        
        >>> a = Plane(Point3D(1,4,6), normal_vector=[2, 4, 6])
        >>> a.parallel_plane(Point3D(2, 3, 5))        
        Plane(Point3D(2, 3, 5), [2, 4, 6])        
        
        """        
        a = self.normal_vector
        return Plane(pt, normal_vector=a)
        