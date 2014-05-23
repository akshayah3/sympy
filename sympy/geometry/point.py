"""Geometrical Points.

Contains
========
Point
Point3D

"""

from __future__ import print_function, division

from sympy.core import S, sympify
from sympy.core.compatibility import iterable
from sympy.core.containers import Tuple
from sympy.simplify import simplify, nsimplify
from sympy.geometry.exceptions import GeometryError
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.complexes import im
from sympy.geometry.entity import GeometryEntity
from sympy.matrices import Matrix
from sympy.core.numbers import Float
from sympy.core.evaluate import global_evaluate

class Point3D(GeometryEntity):
    """Point class for 3D points.
       
    Parameters
    ==========

    coords : sequence of 3 coordinate values.

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
        When trying to add or subtract points with different dimensions
        (other than 2 and 3).
    """
    
    def __new__(cls, *args, **kwargs):
        eval = kwargs.get('evaluate', global_evaluate[0])
        if isinstance(args[0], Point3D):
            if not eval:
                return args[0]
            args = args[0].args            
        else:
            if iterable(args[0]):
                args = args[0]
            if len(args) > 3 or len(args) < 2:
                raise NotImplementedError(
                    "Only two dimensional and three dimensional points currently supported")
        coords = Tuple(*args)
        if eval:
            coords = coords.xreplace(dict(
                [(f, simplify(nsimplify(f, rational=True)))
                for f in coords.atoms(Float)]))
            return GeometryEntity.__new__(cls, *coords)
        
        
        def __hash__(self):
            return super(Point, self).__hash__()
            
        def __lt__(self, other):
            return self.args < other.args    
        
        def __contains__(self, item):
            return item == self
        
    @property
    def x(self):
        """
        Returns the X coordinate of the Point.

        Examples
        ========

        >>> from sympy.geometry.point import Point, Point3D
        >>> p = Point(0, 1)
        >>> p.x
        0
        >>> p = Point3D(1, 2, 3)
        >>> p.x
        1
        """
        return self.args[0]

    @property
    def y(self):
        """
        Returns the Y coordinate of the Point.

        Examples
        ========

        >>> from sympy.geometry.point import Point, Point3D
        >>> p = Point(0, 1)
        >>> p.y
        1
        >>> p = Point3D(0, 1, 1)
        >>>p.y
        1
        """
        return self.args[1]
    
    @property
    def z(self):
        """
        Returns the Z coordinate of the Point.

        Examples
        ========

        >>> from sympy.geometry.point import Point, Point3D
        >>> p = Point3D(0, 1, 1)
        >>> p.z
        1
        >>> p = Point(1, 2)
        >>> p.z
        0
        """
        if len(self.args) == 3:
            return self.args[2]
        else:
            return 0
            
    @property
    def length(self):
        """
        Treating a Point as a Line, this returns 0 for the length of a Point.

        Examples
        ========

        >>> from sympy.geometry.point import Point, Point3D
        >>> p = Point(0, 1)
        >>> p.length
        0
        >>> p = Point3D(1, 1, 1)
        >>> p.length
        0
        """
        return S.Zero
        
    def intersection(self, o):
        """The intersection between this point and another point.
        """
        if isinstance(o, Point) or isinstance(o, Point3D):
            if self == o:
                return [self]
            return []
            
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

         >>> from sympy.geometry.point import Point
         >>> p1, p2 = Point(1, 1), Point(4, 5)
         >>> p1.distance(p2)
         5
         >>> from sympy.geometry.point import Point3D
         >>> p1, p2 = Point3D(1, 1, 2), Point3D(4, 5, 6)
         >>> p1.distance(p2)
         sqrt(41)
         >>> from sympy.geometry.point import *
         >>> p1, p2 = Point3D(1, 1, 1), Point(4, 5)
         >>> p1.distance(p2)
         sqrt(26)

         >>> from sympy.abc import x, y
         >>> p3 = Point(x, y)
         >>> p3.distance(Point(0, 0))
         sqrt(x**2 + y**2)
         """
        a = self.direction_ratio(p)
        return sqrt(sum(i**2 for i in a)) 
            
    def direction_ratio(self, point):
        """Gives the direction ratio between 2 points
        Parameters
        ==========

        p : Point
 
        Returns
        =======

        list
        
        Examples
        ========
        
        >>> from sympy.geomerty.point import Point, Point3D
        >>> p1 = Point3D(1, 2, 3)
        >>> p1.direction_ratio(Point3D(2, 3, 5))
        [1, 1, 2]
        >>> p1.direction_ratio(Point(1, 2))
        [0, 0, -3]
        
        """
        return [(point.x - self.x),(point.y - self.y),(point.z - self.z)]
        
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

        >>> from sympy.geometry.point import Point3D, Point
        >>> p1, p2 = Point3D(1, 1, 1), Point3D(13, 5, 7)
        >>> p1.midpoint(p2)
        Point3D(7, 3, 4)
        >>> p1, p2 = Point3D(1,2,3), Point3D(1,2,4)
        >>>p1.midpoint(p2)
        Point3D(1, 2, 7/2)
        >>>p1, p2 = Point3D(1, 2, 3), Point(1, 2)
        >>> p1.midpoint(p2)
        Point3D(1, 2, 3/2)

        """
        if isinstance(p, Point):
            p = p.convert()
            return Point3D([simplify((a + b)*S.Half) for a, b in zip(self.args,
                            p.args)])
        else:
            return Point3D([simplify((a + b)*S.Half) for a, b in zip(self.args,
                            p.args)])           
            
    def convert(self):
        return self

                
    def evalf(self, prec=None, **options):
        """Evaluate the coordinates of the point.

        This method will, where possible, create and return a new Point
        where the coordinates are evaluated as floating point numbers to
        the precision indicated (default=15).
        """
       
        if prec is None:
            coords = [x.evalf(**options) for x in self.args]
        else:
            coords = [x.evalf(prec, **options) for x in self.args]
        return Point(*coords, evaluate=False)

    n = evalf
    
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
        
        Examples
        ========

        >>> from sympy.geometry.point import Point, Point3D
        >>> p1, p2 = Point(0, 0), Point(1, 1)
        >>> p3, p4, p5 = Point(2, 2), Point(3, 3), Point(1, 2)
        >>> Point.is_collinear(p1, p2, p3, p4)
        True
        >>> Point.is_collinear(p1, p2, p3, p5)
        False

        >>> p1 = Point3D(2, 2, 2)

        >>> p2 = Point3D(-3, -3, -3)
    
        >>> p3 = Point3D(0, 0, 0)

        >>> p4 = Point3D(1, 1, 1)

        >>> p5 = Point3D(1, 2, 3)

        >>> p1.is_collinear(p2, p3, p4)
        True

        >>> p1.is_collinear(p2, p3, p4, p5)
        False
        
        >>> p6 = Point(0,0)
        p1.is_collinear(p2, p4, p6)
        True
        """
        points = list(set(points))
        
        if not all(isinstance(p, Point3D) for p in points):
            raise TypeError('Must pass only Point objects')

        if len(points) == 0:
            return False
        if len(points) <= 2:
            return True
        points = [i.convert() for i in points if isinstance(i, Point)]                            
        
        for i in xrange(0, len(points) - 3):

            pv1 = [j - k for j, k in zip(points[i].args,

                                         points[i + 1].args)]

            pv2 = [j - k for j, k in zip(points[i + 1].args,

                                         points[i + 2].args)]

            rank = Matrix([pv1, pv2]).rank()

            if(rank != 1):

                return False

        return True
        
    def are_coplanar(*points):
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

        >>> from sympy.geometry.point import Point, Point3D
        >>> p1 = Point3D(1, 2, 2)
        >>> p2 = Point3D(2, 7, 2)
        >>> p3 = Point3D(0, 0, 2)
        >>> p4 = Point3D(1, 1, 2)
        >>> p5 = Point3D(1, 2, 2)
        >>> are_coplanar(p1, p2, p3, p4, p5)
        True
        >>> p6 = Point(0,1)
        >>> are_coplanar(p1, p2, p3, p4, p5, p6)
        False

        """
        if(len(points) == 0):

            raise Exception("No parameters provided")

        if(len(points) < 4):

            return True # These cases are always True
        
        points = list(set(points))    
                   
        points = [i.convert() for i in points]
        
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
            
    def translate(self, x=0, y=0, z=0):
        """Shift the Point by adding x , y and z to the coordinates of the 
        Point.

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

        """
        return Point3D(self.x + x, self.y + y, self.z + z)            
            
    def dot(self, p2):
        """Return dot product of self with another Point."""
        p2 = p2.convert()
        self = self.convert()
        x1, y1, z1 = self.args
        x2, y2, z2 = p2.args
        return x1*x2 + y1*y2 + z1*z2      
        
    def __add__(self, other):
        """Add other to self by incrementing self's coordinates by those of other.

        See Also
        ========

        sympy.geometry.entity.translate

        """
        
        if isinstance(other, Point3D):
            if isinstance(other, Point):
                other = other.convert()
                return Point3D(*[simplify(a + b) for a, b in 
                           zip(self.args, other.args)])
        else:                   
            return Point3D(*[simplify(a + b) for a, b in 
                       zip(self.args, other.args)])
        raise ValueError('Cannot add non-Point, %s, to a Point' % other)

    def __sub__(self, other):
        """Subtract two points, or subtract a factor from this point's
        coordinates."""
        if isinstance(other, Point3D):
            if isinstance(other, Point):
                other = other.convert()
                return Point3D(*[simplify(a - b) for a, b in 
                           zip(self.args, other.args)])
        else:                   
            return Point3D(*[simplify(a - b) for a, b in 
                       zip(self.args, other.args)])
        raise ValueError('Cannot subtract non-Point, %s, to a Point' % other)

    def __mul__(self, factor):
        """Multiply point's coordinates by a factor."""
        factor = sympify(factor)
        return Point3D([x*factor for x in self.args])

    def __div__(self, divisor):
        """Divide point's coordinates by a factor."""
        divisor = sympify(divisor)
        return Point3D([x/divisor for x in self.args])

        
class Point(Point3D):
    """
    A 2D Point class which inherits from the Point3D class.
    This is constructed with backward compatibility in mind.
    """
       
    def __init__(self, *args):
        if len(args) != 2:
            raise('Enter a 2 dimensional point')
            
    def convert(self):
        """ Converts the Point to Point3D
        
        
        """
        return Point3D(self.x, self.y, self.z)         
           
    def midpoint(self, p):
        """The midpoint between self and point p.

        Parameters
        ==========

        p : Point

        Returns
        =======

        midpoint : Point


        Examples
        ========

        >>> from sympy.geometry.point import Point, Point3D
        >>> p1, p2 = Point(1, 1), Point(13, 5)
        >>> p1.midpoint(p2)
        Point(7, 3)
        >>>p1, p2 = Point(1, 2), Point3D(1, 2, 4)
        >>> p1.midpoint(p2) = Point3D(1, 2, 2)
        
        """
        if isinstance(p, Point):
            return Point([simplify((a + b)*S.Half) for a, b in zip(self.args,
                            p.args)])
        else:
            self = self.convert()
            return Point3D([simplify((a + b)*S.Half) for a, b in zip(self.args,
                            p.args)])
           
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

        """
        if len(points) == 0:
            return False
        if len(points) <= 2:
            return True
        points = [Point(p) for p in points]
        if len(points) == 3:
            return (not Point.is_collinear(*points))

        try:
            from .ellipse import Circle
            c = Circle(points[0], points[1], points[2])
            for point in points[3:]:
                if point not in c:
                    return False
            return True
        except GeometryError:
            # Circle could not be created, because of collinearity of the
            # three points passed in, hence they are not concyclic.
            return False   
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
        
    def translate(self, x=0, y=0, z=0):
        """Shift the Point by adding x and y to the coordinates of the Point.

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
        >>>t.translate(2,2,2)
        Point3D(2,3,2)

        """
        if z!= 0:
            return self.convert().translate(x,y,z)
        else:    
            return Point(self.x + x, self.y + y)    

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
       
    def __add__(self, other):
        """Add other to self by incrementing self's coordinates by those of other.

        See Also
        ========

        sympy.geometry.entity.translate

        """
        if isinstance(other, Point3D):
            if isinstance(other, Point):
                return Point(*[simplify(a + b) for a, b in 
                           zip(self.args, other.args)])
            else:                   
                return Point3D.__add__(other, self)
        raise ValueError('Cannot add non-Point, %s, to a Point' % other)

    def __sub__(self, other):
        """Subtract two points, or subtract a factor from this point's
        coordinates."""
        if isinstance(other, Point3D):
            if isinstance(other, Point):
                return Point(*[simplify(a - b) for a, b in 
                           zip(self.args, other.args)])
            else:                   
                return Point3D.__sub__(other, self)
        raise ValueError('Cannot subtract non-Point, %s, to a Point' % other)

    def __mul__(self, factor):
        """Multiply point's coordinates by a factor."""
        factor = sympify(factor)
        return Point([x*factor for x in self.args])

    def __div__(self, divisor):
        """Divide point's coordinates by a factor."""
        divisor = sympify(divisor)
        return Point([x/divisor for x in self.args])
