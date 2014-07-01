"""Hyperbolical geometrical entities.

Contains
* Hyperbola
* Rectangular Hyperbola

"""

from sympy.core import S, C, sympify, pi, Dummy
from sympy.core.logic import fuzzy_bool
from sympy.core.numbers import oo, zoo
from sympy.simplify import simplify, trigsimp
from sympy.functions.elementary.miscellaneous import sqrt, Max, Min
from sympy.functions.elementary.complexes import im
from sympy.geometry.exceptions import GeometryError
from sympy.polys import Poly, PolynomialError
from sympy.solvers import solve
from sympy.utilities.lambdify import lambdify
from sympy.utilities.iterables import uniq
from sympy.utilities.misc import filldedent
from .entity import GeometryEntity
from .point import Point
from .line import LinearEntity, Line
from .util import _symbol, idiff
from sympy.mpmath import findroot as nroot

class Hyperbola(GeometryEntity):
    """A Hyoerbolical GeometryEntity.

    Parameters
    ==========

    center : Point, optional
        Default value is Point(0, 0)
    hradius : number or SymPy expression, optional
    vradius : number or SymPy expression, optional
    eccentricity : number or SymPy expression, optional
        Two of `hradius`, `vradius` and `eccentricity` must be supplied to
        create an Hyperbola. The third is derived from the two supplied.

    Attributes
    ==========

    center
    hradius
    vradius
    eccentricity
    periapsis
    apoapsis
    focus_distance
    foci
    
    Raises
    ======

    GeometryError
        When `hradius`, `vradius` and `eccentricity` are incorrectly supplied
        as parameters.
    TypeError
        When `center` is not a Point.

    See Also
    ========

    Ellipse
    """    
    def __new__(
        cls, center=None, hradius=None, vradius=None, eccentricity=None,
            **kwargs):
        hradius = sympify(hradius)
        vradius = sympify(vradius)

        eccentricity = sympify(eccentricity)

        if center is None:
            center = Point(0, 0)
        else:
            center = Point(center)

        if len(list(filter(None, (hradius, vradius, eccentricity)))) != 2:
            raise ValueError('Exactly two arguments of "hradius", '
                '"vradius", and "eccentricity" must not be None."')

        if eccentricity is not None:
            if hradius is None:
                hradius = vradius / sqrt(1 + eccentricity**2)
            elif vradius is None:
                vradius = hradius * sqrt(1 + eccentricity**2)

        return GeometryEntity.__new__(cls, center, hradius, vradius, **kwargs)
        
    @property
    def center(self):
        """The center of the Hyperbola.

        Returns
        =======

        center : number

        See Also
        ========

        sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.center
        Point(0, 0)

        """
        return self.args[0]

    @property
    def hradius(self):
        """The horizontal radius of the Hyperbola.

        Returns
        =======

        hradius : number

        See Also
        ========

        vradius, major, minor

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.hradius
        3

        """
        return self.args[1]

    @property
    def vradius(self):
        """The vertical radius of the Hyperbola.

        Returns
        =======

        vradius : number

        See Also
        ========

        hradius, major, minor

        Examples
        ========

        >>> from sympy import Point, Hyperbola
        >>> p1 = Point(0, 0)
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.vradius
        1

        """
        return self.args[2]

    @property
    def minor(self):
        """Shorter axis of the Hyperbola (if it can be determined) else vradius.

        Returns
        =======

        minor : number or expression

        See Also
        ========

        hradius, vradius, major

        Examples
        ========

        >>> from sympy import Point, Hyperbola, Symbol
        >>> p1 = Point(0, 0)
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.minor
        1

        >>> a = Symbol('a')
        >>> b = Symbol('b')
        >>> Hyperbola(p1, a, b).minor
        b
        >>> Hyperbola(p1, b, a).minor
        a

        >>> m = Symbol('m')
        >>> M = m + 1
        >>> Hyperbola(p1, m, M).minor
        m

        """
        rv = Min(*self.args[1:3])
        if rv.func is Min:
            return self.vradius
        return rv

    @property
    def major(self):
        """Longer axis of the ellipse (if it can be determined) else hradius.

        Returns
        =======

        major : number or expression

        See Also
        ========

        hradius, vradius, minor

        Examples
        ========

        >>> from sympy import Point, Ellipse, Symbol
        >>> p1 = Point(0, 0)
        >>> e1 = Hyperbola(p1, 3, 1)
        >>> e1.major
        3

        >>> a = Symbol('a')
        >>> b = Symbol('b')
        >>> Hyperbola(p1, a, b).major
        a
        >>> Hyperbola(p1, b, a).major
        b

        >>> m = Symbol('m')
        >>> M = m + 1
        >>> Hyperbola(p1, m, M).major
        m + 1

        """
        rv = Max(*self.args[1:3])
        if rv.func is Max:
            return self.hradius
        return rv
