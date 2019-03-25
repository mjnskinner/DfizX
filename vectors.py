from __future__ import division
from functools import reduce
import math


class Point(object):
    """ Point class: Reprepsents a point in the x, y, z space. """

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return '{0}({1}, {2}, {3})'.format(self.__class__.__name__, self.x,
                                           self.y, self.z)

    def substract(self, point):
        """ Return a Point instance as the displacement of two points. """

        return Point(point.x - self.x, point.y - self.y, point.z - self.z)

    @classmethod
    def from_list(cls, l):
        """ Return a Point instance from a given list """

        x, y, z = map(float, l)
        return cls(x, y, z)


class Vector(Point):
    """ Vector class: Represents a vector in the x, y, z space. """

    def __init__(self, x, y, z):
        self.vector = [x, y, z]
        super(Vector, self).__init__(x, y, z)

    def add(self, number):
        """ Return a Vector instance as the product of the vector and a real
            number. """

        return Vector(self.x + number,self.y+number,self.z+number)

    def multiply(self, number):
        """ Return a Vector instance as the product of the vector and a real
            number. """

        return Vector(self.x * number,self.y*number,self.z*number)
		
    def divide(self, number):
        """ Return a Vector instance as the product of the vector and a real
            number. """
        
        return Vector(self.x / number,self.y/number,self.z/number)

    def magnitude(self):
        """ Return magnitude of the vector. """
        return (math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2))

    def sum(self, vector):
        """ Return a Vector instance as the vector sum of two vectors. """
		
        return (Vector(self.x + vector.x,self.y+vector.y,self.z+vector.z))
            

    def substract(self, vector):
        """ Return a Vector instance as the vector difference of two vectors.
        """

        return (Vector(self.x - vector.x,self.y - vector.y,self.z - vector.z))

    def dot(self, vector, theta=None):
        """ Return the dot product of two vectors. If theta is given then the
        dot product is computed as v1*v1 = |v1||v2|cos(theta). Argument theta
        is measured in degrees. """

        if theta is not None:
            return (self.magnitude() * vector.magnitude() *
                    math.degrees(math.cos(theta)))
        return (self.x*vector.x+self.y*vector.y+self.z*vector.z)

    def cross(self, vector):
        """ Return a Vector instance as the cross product of two vectors """

        return Vector((self.y * vector.z - self.z * vector.y),
                      (self.z * vector.x - self.x * vector.z),
                      (self.x * vector.y - self.y * vector.x))

    def angle(self, vector):
        """ Return the angle between two vectors in degrees. """
        return (math.acos((self.dot(vector) / (self.magnitude() * vector.magnitude()))))

    def parallel(self, vector):
        """ Return True if vectors are parallel to each other. """

        if self.cross(vector).magnitude() == 0:
            return True
        return False

    def perpendicular(self, vector):
        """ Return True if vectors are perpendicular to each other. """

        if abs(self.dot(vector)) <= 0.0000000000001:
            return True
        return False

    def non_parallel(self, vector):
        """ Return True if vectors are non-parallel. Non-parallel vectors are
            vectors which are neither parallel nor perpendicular to each other.
        """

        if (self.parallel(vector) is not True and
                self.is_perpendicular(vector) is not True):
            return True
        return False

    @classmethod
    def from_points(cls, point1, point2):
        """ Return a Vector instance from two given points. """

        if isinstance(point1, Point) and isinstance(point2, Point):
            displacement = point1.substract(point2)
            return cls(displacement.x, displacement.y, displacement.z)
        raise TypeError


def vector_from_direction_loft_speed (throw_direction,throw_loft,throw_speed):
	throw_direction /= 57.3
	throw_loft /= 57.3
	z=math.sin(throw_loft)
	x=math.cos(throw_loft)*math.sin(throw_direction)
	y=math.cos(throw_loft)*math.cos(throw_direction)
	return Vector (x*throw_speed,y*throw_speed,z*throw_speed)
	
def make_unit_vector (vector):
	return vector.divide(vector.magnitude())