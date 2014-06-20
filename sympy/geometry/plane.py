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
from .line import Line, Segment, Ray
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
        '''Projection of any point on the plane
           Will result in a 3D point'''
        x, y, z = symbols("x, y, z")
        k = self.equation(x, y, z)
        const = [i for i in k.args if i.is_constant() is True]
        a, b, c = k.coeff(x), k.coeff(y), k.coeff(z)
        if const != []:
            d = -(a*pt.args[0] + b*pt.args[1] + const.pop())
            t = d/(a**2 + b**2 + c**2)
        else:
            d = -(a*pt.args[0] + b*pt.args[1])
            t = d/(a**2 + b**2 + c**2)
        if isinstance(pt, Point):
            return Point3D(pt.args[0] + t*a, pt.args[1] + t*b, t*c)
        if isinstance(pt, Point3D):
            return Point3D(pt.args[0] + t*a, pt.args[1] + t*b, pt.args[2] + t*c)

    def projection_line(self, l):
        '''Projection of any line on the plane.
           Result will be a 3D line or Ray or segment'''
        if isinstance(l, Line) or isinstance(l, Line3D):
            a, b = self.projection(l.p1), self.projection(l.p2)
            return Line3D(a, b)
        if isinstance(l, Ray) or isinstance(l, Ray3D):
            a, b = self.projection(l.p1), self.projection(l.p2)
            return Ray3D(a, b)
        if isinstance(l, Segment) or isinstance(l, Segment3D):
            a, b = self.projection(l.p1), self.projection(l.p2)
            return Segment3D(a, b)

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

    def random_point(self, seed=None):
        x, y, z = symbols("x, y, z")
        a = self.equation(x, y, z)
        from sympy import Rational
        import random
        if seed is not None:
            rng = random.Random(seed)
        else:
            rng = random
        for i in range(10):
            c = 2*Rational(rng.random()) - 1
            s = sqrt(1 - c**2)
            a = solve(a.subs([(y, c), (z, s)]))
            if a is []:
                d = Point3D(0, c, s)
            else:
                d = Point3D(a[0], c, s)
            if d in self:
                return d
        raise GeometryError(
            'Having problems generating a point in the plane')

    def intersection(self, o):
        from sympy.geometry.line3d import LinearEntity3D
        from sympy.geometry.line import LinearEntity
        if isinstance(o, Point) or isinstance(o, Point3D):
            if o in self:
                return o
            else:
                return []
        if isinstance(o, LinearEntity3D):        
            t = symbols('t')
            x, y, z = symbols("x y z")            
            if o in self:
                return [o]    
            else:
                a = o.arbitrary_point(t)
                b = self.equation(x, y, z)
                c = solve(b.subs([(x, a.x), (y, a.y), (z, a.z)]))
                if c is []:
                    return []
                a = a.subs(t, c[0])
                if a in o:
                    return a
                else:
                    return []
        if isinstance(o, LinearEntity):
            t = symbols('t')
            x, y, z = symbols("x y z")            
            if o in self:
                return [o]
            else:
                a = self.equation(x, y, z).subs(z,0)
                b = o.equation(x, y)
                c = solve((a,b),[x, y])
                if c is {}:
                    return []
                else:
                    return [Point(c[x], c[y])]
        if isinstance(o, Plane):
            if self.is_parallel(o):
                return []
            else:    
                x, y, z = symbols("x y z")
                a, b= Matrix([self.normal_vector]), Matrix([o.normal_vector])
                c = list(a.cross(b))
                d = self.equation(x, y, z)
                e = o.equation(x, y, z)
                f = solve((d.subs(z,0), e.subs(z,0)),[x, y])
                g = solve((d.subs(y,0), e.subs(y,0)),[x, z])
                h = solve((d.subs(x,0), e.subs(x,0)),[y, z])
                if len(f) == 2:                
                    return [Line3D(Point3D(f[x], f[y], 0), c)]
                if len(g) == 2:
                    return [Line3D(Point3D(g[x], 0, g[z]), c)]
                if len(h) == 2:
                    return [Line3D(Point3D(0, h[y], h[z]), c)]

    def __contains__(self, o):
        from sympy.geometry.line3d import LinearEntity3D
        from sympy.geometry.line import LinearEntity
        x, y, z = symbols("x, y, z")
        k = self.equation(x, y, z)
        if isinstance(o, Point3D):
            d = k.subs([(x, o.x), (y, o.y), (z, o.z)])
            return d == 0
        elif isinstance(o, Point):
            if self.projection(o) == Point3D(o.x, o.y, 0):
                return True
            else:
                return False
        elif isinstance(o, LinearEntity3D):
            d = o.arbitrary_point()
            e = k.subs([(x, d[0]), (y, d[1]), (z, d[2])])
            return e == 0
        elif isinstance(o, LinearEntity):
            if self.projection_line(o) == Line3D(Point3D(o.p1.x, o.p1.y, 0), 
                                                Point3D(o.p2.x, o.p2.y, 0)):
                return True
            else:
                return False                    
        else:
            return False