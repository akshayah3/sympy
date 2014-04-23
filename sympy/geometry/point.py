"""Geometrical Points.

Contains
========
Point

"""

from __future__ import print_function, division

from sympy.core import S, sympify
from sympy.core.compatibility import iterable
from sympy.core.containers import Tuple
from sympy.simplify import simplify, nsimplify
from sympy.geometry.exceptions import GeometryError
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.complexes import im
from .entity import GeometryEntity
from sympy.matrices import Matrix
from sympy.core.numbers import Float
from sympy.core.evaluate import global_evaluate


class Point(GeometryEntity):
    """A point in both 2-dimensional and 3-dimensional Euclidean space.

    Parameters
    ==========

    coords : sequence of 2 or 3 coordinate values.

    Attributes
    ==========

    x
    y
    z
    length

    Raises
    ======

    NotImplementedError
        When trying to create a point with more than three dimensions.
        When `intersection` is called with object other than a Point.
    TypeError
        When trying to add or subtract points with different dimensions where 
        any point is having a dimension greater than 3.

    Notes
    =====

    Currently only 2-dimensional and 3-dimensional points are supported.

    See Also
    ========

    sympy.geometry.line.Segment : Connects two Points

    Examples
    ========

    >>> from sympy.geometry import Point
    >>> from sympy.abc import x
    >>> Point(1, 2)
    Point(1, 2)
    
    >>> Point([1, 2])
    Point(1, 2)
    
    >>> Point(0, x)
    Point(0, x)
    
    >>> Point(1,2,3)
    Point(1,2,3)
    
    >>> Point(1,x,2)
    Point(1,x,2)
    
    Floats are automatically converted to Rational unless the
    evaluate flag is False:

    >>> Point(0.5, 0.25)
    Point(1/2, 1/4)
    >>> Point(0.5, 0.25, evaluate=False)
    Point(0.5, 0.25)

    """
    def __new__(cls, *args, **kwargs):
        eval = kwargs.get('evaluate', global_evaluate[0])
        check = True
        if isinstance(args[0], Point):
            if not eval:
                return args[0]
            args = args[0].args
            check = False
        else:
            if iterable(args[0]):
                args = args[0]
            if len(args) > 3 or len(args) == 1:
                raise NotImplementedError(
                    "only 2-D and 3-D points are currently supported")
        coords = Tuple(*args)
        if check:
            if any(a.is_number and im(a) for a in coords):
                raise ValueError('Imaginary args not permitted.')
        if eval:
            coords = coords.xreplace(dict(
                [(f, simplify(nsimplify(f, rational=True)))
                for f in coords.atoms(Float)]))
        return GeometryEntity.__new__(cls, *coords)

    def __hash__(self):
        return super(Point, self).__hash__()

    def __eq__(self, other):
        ts, to = type(self), type(other)
        if ts is not to:
            return False
        self, other = self.convert(other)    
        return self.args == other.args

    def __lt__(self, other):
        return self.args < other.args

    def __contains__(self, item):
        self, item = self.convert(item)
        return item == self

    @property
    def x(self):
        """
        Returns the X coordinate of the Point.

        Examples
        ========

        >>> from sympy import Point
        >>> p = Point(0, 1)
        >>> p.x
        0
        """
        return self.args[0]

    @property
    def y(self):
        """
        Returns the Y coordinate of the Point.

        Examples
        ========

        >>> from sympy import Point
        >>> p = Point(0, 1)
        >>> p.y
        1
        """
        return self.args[1]
        
    @property
    def z(self):
        """
        Returns the Z coordinate of the Point if provided otherwise it is taken
        as 0 by default.

        Examples
        ========

        >>> from sympy import Point
        >>> p = Point(0, 1, 2)
        >>> p.z
        2
        """
        if len(self.args) == 2:
            return S.Zero
        else:    
            return self.args[2]  
            
    @property
    def length(self):
        """
        Treating a Point as a Line, this returns 0 for the length of a Point.

        Examples
        ========

        >>> from sympy import Point
        >>> p = Point(0, 1)
        >>> p.length
        0
        """
        return S.Zero
    
    def point_3d(self):
        """
        Converts any given point into a 3-D point
        
        Examples
        ========
        
        >>> from sympy import Point
        >>> p = Point(0, 1)
        >>> p.point_3d()
        Point(0, 1, 0)
        
        >>> p = Point(1,1,1)
        >>> p.point_3d()
        Point(1, 1, 1)
        """
        if len(self.args) == 3:
            return self
        else:
            a = list(self.args)
            a.append(S.Zero)
            return Point(a)
            
    def convert(self, point):
        """
        This method converts a point accordingly to a different dimension
        
        Examples
        ========
        
        >>>from sympy import Point
        >>> a = Point(1, 2, 3)
        >>> a.convert(Point(1, 2, 4))
        Point(1, 2, 3) , Point(1, 2, 4)
        
        >>> a = Point(1, 2, 3)
        >>> a.convert(Point(1, 2))
        Point(1, 2, 3) , Point(1, 2, 0)
        
        >>> a = Point(1, 2)
        >>> a.convert(Point(2, 1))
        Point(1, 2) , Point(2, 1)
        """
        point = Point(point)
        if len(self.args) == len(point.args):
            return self, point
            
        elif self.args > point.args:
            a = list(point.args)
            a.append(S.Zero)
            point = Point(a)
            return self , point
        else:
            a = list(self.args)
            a.append(S.Zero)
            self = Point(a)
            return self, point
            
    def direction_ratio(self, point):
        """
        Returns the Direction Ratio between two 3-D points
        
        Examples
        ========
        
        >>> from sympy import Point
        >>> p1 = Point(1, 2, 3)
        >>> p1.direction_ratio(Point(2, 3, 5))
        [1, 1, 2]
        
        >>> p1.direction_ratio(Point(1, 2))
        [0, 0, -3]
        """
        point = Point(point)
        return [(point.x - self.x),(point.y - self.y),(point.z - self.z)]
        
    def direction_cosine(self, point):
        """
        Returns the Direction Cosine between two 3-D points
        
        Examples
        ========
        
        >>> from sympy import Point
        >>> p3, p4 = Point(0, 1, -1), Point(-1, 2)
        >>> p3.direction_cosine(p4)
        [-sqrt(3)/3, sqrt(3)/3, sqrt(3)/3]
        """
        a = self.direction_ratio(point)
        b = sqrt(sum([i**2 for i in a]))
        return [i/b for i in a]
        
    def is_collinear(*points):
        """Is a sequence of points collinear?

        Test whether or not a set of points are collinear. Returns True if
        the set of points are collinear, or False otherwise.

        Parameters
        ==========

        points : sequence of Point

        Returns
        =======

        is_collinear : boolean

        Notes
        =====

        Slope is preserved everywhere on a line, so the slope between
        any two points on the line should be the same. Take the first
        two points, p1 and p2, and create a translated point v1
        with p1 as the origin. Now for every other point we create
        a translated point, vi with p1 also as the origin. Note that
        these translations preserve slope since everything is
        consistently translated to a new origin of p1. Since slope
        is preserved then we have the following equality:

              * v1_slope = vi_slope
              * v1.y/v1.x = vi.y/vi.x (due to translation)
              * v1.y*vi.x = vi.y*v1.x
              * v1.y*vi.x - vi.y*v1.x = 0           (*)

        Hence, if we have a vi such that the equality in (*) is False
        then the points are not collinear. We do this test for every
        point in the list, and if all pass then they are collinear.

        See Also
        ========

        sympy.geometry.line.Line

        Examples
        ========

        >>> from sympy import Point
        >>> from sympy.abc import x
        >>> p1, p2 = Point(0, 0), Point(1, 1)
        >>> p3, p4, p5 = Point(2, 2), Point(x, x), Point(1, 2)
        >>> Point.is_collinear(p1, p2, p3, p4)
        True
        >>> Point.is_collinear(p1, p2, p3, p5)
        False

        >>> p1 = Point(2, 2, 2)

        >>> p2 = Point(-3, -3, -3)

        >>> p3 = Point(0, 0, 0)

        >>> p4 = Point(1, 1, 1)

        >>> p5 = Point(1, 2, 3)

        >>> p1.is_collinear( p2, p3, p4)

        True

        >>> p1.is_collinear(p2, p3, p4, p5)

        False

        """
        # Coincident points are irrelevant and can confuse this algorithm.
        # Use only unique points.
        points = list(set(points))
        if not all(isinstance(p, Point) for p in points):
            raise TypeError('Must pass only Point objects')

        if len(points) == 0:
            return False
        if len(points) <= 2:
            return True  # two points always form a line
        points = [Point(a).point_3d() for a in points]

        # XXX Cross product is used now, but that only extends to three
        #     dimensions. If the concept needs to extend to greater
        #     dimensions then another method would have to be used
        """p1 = points[0].point_3d()
        p2 = points[1].point_3d()
        v1 = p2 - p1
        x1, y1, z1= v1.args
        rv = True"""
        for i in xrange(0, len(points) - 3):

            pv1 = [j - k for j, k in zip(points[i].args,

                                         points[i + 1].args)]

            pv2 = [j - k for j, k in zip(points[i + 1].args,

                                         points[i + 2].args)]

            rank = Matrix([pv1, pv2]).rank()

            if(rank != 1):

                return False

        return True

    def is_concyclic(*points):
        """Is a sequence of points concyclic?

        Test whether or not a sequence of points are concyclic (i.e., they lie
        on a circle).

        Parameters
        ==========

        points : sequence of Points

        Returns
        =======

        is_concyclic : boolean
            True if points are concyclic, False otherwise.

        See Also
        ========

        sympy.geometry.ellipse.Circle

        Notes
        =====

        No points are not considered to be concyclic. One or two points
        are definitely concyclic and three points are conyclic iff they
        are not collinear.

        For more than three points, create a circle from the first three
        points. If the circle cannot be created (i.e., they are collinear)
        then all of the points cannot be concyclic. If the circle is created
        successfully then simply check the remaining points for containment
        in the circle.

        Examples
        ========

        >>> from sympy.geometry import Point
        >>> p1, p2 = Point(-1, 0), Point(1, 0)
        >>> p3, p4 = Point(0, 1), Point(-1, 2)
        >>> Point.is_concyclic(p1, p2, p3)
        True
        >>> Point.is_concyclic(p1, p2, p3, p4)
        False

        >>> p1, p2 = Point(-1, 0, 1), Point(1, 0, 0)
        >>> p3, p4 = Point(0, 1, -1), Point(-1, 2)
        >>> Point.is_concyclic(p1, p2, p3)
        False
        
        """
        if len(points) == 0:
            return False
        if len(points) <= 2:
            return True
        points = [Point(p) for p in points]
        if len(points) == 3:
            return (not Point.is_collinear(*points))

        """try:
            from .ellipse import Circle
            c = Circle(points[0], points[1], points[2])
            for point in points[3:]:
                if point not in c:
                    return False
            return True
        except GeometryError:
            # Circle could not be created, because of collinearity of the
            # three points passed in, hence they are not concyclic.
            return False"""
            
    def coplanar_points(*points):

        """

        This function tests whether passed points are coplanar or not.



        It uses the fact that the triple scalar product of three vectors

        vanishes iff the vectors are coplanar. Which means that the volume

        of the solid described by them will have to be zero for coplanarity.



        Parameters

        ==========



        A set of points 3D points



        Returns

        =======



        boolean



        Examples

        ========



        >>> from sympy import Point 

        >>> p1 = Point(1, 2, 2)

        >>> p2 = Point(2, 7, 2)

        >>> p3 = Point(0, 0, 2)

        >>> p4 = Point(1, 1, 2)

        >>> p5 = Point(1, 2, 2)

        >>> coplanar_points(p1, p2, p3, p4, p5)

        True

        """



        if(len(points) == 0):

            raise Exception("No parameters provided")

        points = [Point(point).point_3d() for point in points]

        if(len(points) < 4):

            return True  # These cases are always True

        for i in xrange(0, len(points) - 3):

            pv1 = [j - k for j, k in zip(points[i].args, points[i + 1].args)]

            pv2 = [j - k for j, k in zip(

                points[i + 1].args,

                points[i + 2].args)]

            pv3 = [j - k for j, k in zip(

                points[i + 2].args,

                points[i + 3].args)]

            pv1, pv2, pv3 = Matrix(pv1), Matrix(pv2), Matrix(pv3)

            stp = pv1.dot(pv2.cross(pv3))

            if stp != 0:

                return False

            return True        

    def distance(self, p):
        """The Euclidean distance from self to point p.

        Parameters
        ==========

        p : Point

        Returns
        =======

        distance : number or symbolic expression.

        See Also
        ========

        sympy.geometry.line.Segment.length

        Examples
        ========

        >>> from sympy.geometry import Point
        >>> p1, p2 = Point(1, 1), Point(4, 5)
        >>> p1.distance(p2)
        5
        
        >>> from sympy.geometry import Point
        >>> p1, p2 = Point(1, 1, 2), Point(4, 5, 6)
        >>> p1.distance(p2)
        sqrt(41)
        
        >>> from sympy.geometry import Point
        >>> p1, p2 = Point(1, 1, 1), Point(4, 5)
        >>> p1.distance(p2)
        sqrt(26)

        >>> from sympy.abc import x, y
        >>> p3 = Point(x, y)
        >>> p3.distance(Point(0, 0))
        sqrt(x**2 + y**2)

        """
        p = Point(p)
        a = self.direction_ratio(p)
        return sqrt(sum(i**2 for i in a))

    def midpoint(self, p):
        """The midpoint between self and point p.

        Parameters
        ==========

        p : Point

        Returns
        =======

        midpoint : Point

        See Also
        ========

        sympy.geometry.line.Segment.midpoint

        Examples
        ========

        >>> from sympy.geometry import Point
        >>> p1, p2 = Point(1, 1), Point(13, 5)
        >>> p1.midpoint(p2)
        Point(7, 3)
        
        >>> p1, p2 = Point(1,2,3), Point(1,2,4)
        >>>p1.midpoint(p2)
        Point(1, 2, 3.5)
        
        >>>p1, p2 = Point(1, 2, 3), Point(1, 2)
        >>> p1.midpoint(p2) = Point(1, 2, 3/2)

        """
        self , p = self.convert(p)
        return Point([simplify((a + b)*S.Half) for a, b in zip(self.args, p.args)])

    def evalf(self, prec=None, **options):
        """Evaluate the coordinates of the point.

        This method will, where possible, create and return a new Point
        where the coordinates are evaluated as floating point numbers to
        the precision indicated (default=15).

        Returns
        =======

        point : Point

        Examples
        ========

        >>> from sympy import Point, Rational
        >>> p1 = Point(Rational(1, 2), Rational(3, 2))
        >>> p1
        Point(1/2, 3/2)
        >>> p1.evalf()
        Point(0.5, 1.5)
                
        >>> p2 = Point(1/2, 3/2, 4/5)
        >>> p2.evalf()
        Point(0.5, 1.5, 0.8)
        
        """
        if prec is None:
            coords = [x.evalf(**options) for x in self.args]
        else:
            coords = [x.evalf(prec, **options) for x in self.args]
        return Point(*coords, evaluate=False)

    n = evalf

    def intersection(self, o):
        """The intersection between this point and another point.

        Parameters
        ==========

        other : Point

        Returns
        =======

        intersection : list of Points

        Notes
        =====

        The return value will either be an empty list if there is no
        intersection, otherwise it will contain this point.

        Examples
        ========

        >>> from sympy import Point
        >>> p1, p2, p3 = Point(0, 0), Point(1, 1), Point(0, 0)
        >>> p1.intersection(p2)
        []
        >>> p1.intersection(p3)
        [Point(0, 0)]

        """
        if isinstance(o, Point):
            if self == o:
                return [self]
            return []

        return o.intersection(self)

    def rotate(self, angle, pt=None):
        """Rotate ``angle`` radians counterclockwise about Point ``pt``.

        See Also
        ========

        rotate, scale

        Examples
        ========

        >>> from sympy import Point, pi
        >>> t = Point(1, 0)
        >>> t.rotate(pi/2)
        Point(0, 1)
        >>> t.rotate(pi/2, (2, 0))
        Point(2, -1)

        """
        from sympy import cos, sin, Point

        c = cos(angle)
        s = sin(angle)

        rv = self
        if pt is not None:
            pt = Point(pt)
            rv -= pt
        x, y = rv.args
        rv = Point(c*x - s*y, s*x + c*y)
        if pt is not None:
            rv += pt
        return rv

    def scale(self, x=1, y=1, z=1, pt=None):
        """Scale the coordinates of the Point by multiplying by
        ``x`` , ``y`` and ``z`` after subtracting ``pt`` -- default is (0, 0, 0) --
        and then adding ``pt`` back again (i.e. ``pt`` is the point of
        reference for the scaling).

        See Also
        ========

        rotate, translate

        Examples
        ========

        >>> from sympy import Point
        >>> t = Point(1, 1)
        >>> t.scale(2)
        Point(2, 1)
        >>> t.scale(2, 2)
        Point(2, 2)

        """
        if pt:
            pt = pt.point_3d()
            self = self.point_3d()
            self = self.__sub__(pt)
            a = Point((self.x)*x, (self.y)*y, (self.z)*z)
            return a.__add__(pt)
        return Point(self.x*x, self.y*y, self.z*z)

    def translate(self, x=0, y=0, z=0):
        """Shift the Point by adding x , y and z to the coordinates of the Point.

        See Also
        ========

        rotate, scale

        Examples
        ========

        >>> from sympy import Point
        >>> t = Point(0, 1)
        >>> t.translate(2)
        Point(2, 1)
        >>> t.translate(2, 2)
        Point(2, 3)
        >>> t + Point(2, 2)
        Point(2, 3)
        >>> a = Point(1, 2)
        >>> a.translate(0, 0, 3)
        Point(1, 2, 3)

        """
        return Point(self.x + x, self.y + y, self.z + z)

    def transform(self, matrix):
        """Return the point after applying the transformation described
        by the 3x3 Matrix, ``matrix``.

        See Also
        ========
        geometry.entity.rotate
        geometry.entity.scale
        geometry.entity.translate
        """
        x, y = self.args
        return Point(*(Matrix(1, 3, [x, y, 1])*matrix).tolist()[0][:2])

    def dot(self, p2):
        """Return dot product of self with another Point."""
        p2 = Point(p2)
        if len(self.args) != len(p2.args):
            x1, y1 = self.args
            x2, y2 = p2.args
            return x1*x2 + y1*y2
        else:
            self = self.point_3d()
            p2 = p2.point_3d()
            x1, y1, z1 = self.args
            x2, y2, z2 = p2.args
            return x1*x2 + y1*y2 + z1*z2
            
    def __add__(self, other):
        """Add other to self by incrementing self's coordinates by those of other.

        See Also
        ========

        sympy.geometry.entity.translate

        """

        if isinstance(other, Point):
            if len(other.args) == len(self.args):
                return Point(*[simplify(a + b) for a, b in
                               zip(self.args, other.args)])
            else:
                self = self.point_3d()
                other = other.point_3d()
                return Point(*[simplify(a + b) for a, b in
                               zip(self.args, other.args)])
        else:
            raise ValueError('Cannot add non-Point, %s, to a Point' % other)

    def __sub__(self, other):
        """Subtract two points, or subtract a factor from this point's
        coordinates."""
        if isinstance(other, Point):
            if len(other.args) == len(self.args):
                return Point(*[simplify(a - b) for a, b in
                               zip(self.args, other.args)])
            else:
                self = self.point_3d()
                other = other.point_3d()
                return Point(*[simplify(a - b) for a, b in
                               zip(self.args, other.args)])
        else:
            raise ValueError('Cannot substract non-Point,'
                      '%s, to a Point' % other)

    def __mul__(self, factor):
        """Multiply point's coordinates by a factor."""
        factor = sympify(factor)
        return Point([x*factor for x in self.args])

    def __div__(self, divisor):
        """Divide point's coordinates by a factor."""
        divisor = sympify(divisor)
        return Point([x/divisor for x in self.args])

    __truediv__ = __div__

    def __neg__(self):
        """Negate the point."""
        return Point([-x for x in self.args])

    def __abs__(self):
        """Returns the distance between this point and the origin."""
        origin = Point([0]*len(self.args))
        return Point.distance(origin, self)
