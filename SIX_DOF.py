'''
@Authors: Harrison Leece, James Hribal, Max Fung, Nils Heidenreich
@Purpose: Explore 6DOF rocket trajectory, esspecially quaternion rotation
Learning resources: https://eater.net/quaternions
'''

import numpy as np
import oyaml as yaml
import math

class Rotator:
    def __init__(self):
        self.re = 0; self.i = 0; self.j = 0; self.k = 1
        self.body_vector = np.array([[0],[1],[0],[0]])

    '''
    https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Using_quaternions_as_rotations
    This function should take inputs: 'Cartesian, unit rotation-axis (Vector),
    Rotation Angle in radians (Theta) and form a quaternion vector
    '''
    def form_quaternion(self, vector, theta):
        assert self.vector_is_unit(vector), 'Class: Rotator, Fxn: form_quaternion, vector is not a unit quaternion'
        r = np.cos(theta/2)
        i = -1*np.sin(theta/2)*vector[0]
        j = -1*np.sin(theta/2)*vector[1]
        k = -1*np.sin(theta/2)*vector[2]
        quaternion = np.array([[r],[i],[j],[k]])
        return quaternion

    '''
    https://math.stackexchange.com/questions/40164/how-do-you-rotate-a-vector-by-a-unit-quaternion
    '''
    def rotate_body(self, quaternion):
        left = quaternion
        right = np.array([[quaternion[0]], -quaternion[1], -quaternion[2], -quaternion[3]])
        h1 = self.hamilton_product(left, self.body_vector)
        print('H1: {}'.format(h1))
        self.body_vector = self.hamilton_product(h1,right)

    '''
    https://en.wikipedia.org/wiki/Quaternion#Hamilton_product
    https://math.stackexchange.com/questions/40164/how-do-you-rotate-a-vector-by-a-unit-quaternion
    '''
    def hamilton_product(self, vec1, vec2):
        a1 = vec1[0]; a2 = vec2[0]
        b1 = vec1[1]; b2 = vec2[1]
        c1 = vec1[2]; c2 = vec2[2]
        d1 = vec1[3]; d2 = vec2[3]

        r = float(a1*a2 - b1*b2 - c1*c2 - d1*d2)
        x = float(a1*b2 + b1*a2 + c1*d2 - d1*c2)
        y = float(a1*c2 - b1*d2 + c1*a2 + d1*b2)
        z = float(d1*d2 + b1*c2 - c1*b2 + d1*a2)
        return np.array([[r],[x],[y],[z]])

    def report_body_vector(self):
        print(self.body_vector)
    '''
    Convert some arbitrary vector to a unit vector (divide components by the magnitude)
    '''
    def unitify_vector(self, vector):
        mag = np.linalg.norm(vector)
        return vector/mag

    '''
    Checker function to verify a vector of arbitrary length is a unit vector
    Tolerance variable to allow 'close enough' cases to succeed
    '''
    def vector_is_unit(vec):
        squares = [x*x for x in vec]
        vec_sum = np.sum(squares)
        norm = np.sqrt(vec_sum)
        tolerance = .01
        if abs(norm-1) < tolerance:
            return true
        return false
    '''
    https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Using_quaternions_as_rotations
    q' = q2q1
    q1 first, then q2
    USE QUATERNION MULTIPLICATION RULE:
    v*w  where v and w are both quaternions with no real part
    v*w = v x w - v * w
    v*w  where v and w are both quaternions with real part s and t (see wikipedia)
    v*w = (s + v)*(t + w) = (st - v * w)+(sw + tv + v x w)
    '''
    def combine_quaternions(self, q1, q2):
        re1 = q1[0]; re2 = q2[0]
        q1 = q1[1:3]; q2 = q2[1:3]
        cross = np.cross(q2, q1); dot = np.dot(q2, q1)
        re_prime = (re1*re2 - dot)
        temp = re1*q2 + re2*q1 + cross
        q_prime = np.array([[re_prime],[temp[0]],[temp[1]],[temp[2]]])
        return q_prime

class Rocket(Rotator):
    def __init__(self, rocket_data, environment_obj, reference=None, units='english'):
        super().__init__()
        '''
        abs_orientation is a (unit) vector representing the orientation of the rocket
        relative to the launch location axis or global axis
        '''
        self.abs_orientation = np.array([[0],[0],[1]])
        '''
        rel_orientation is a (unit) vector representing the orientation of the rocket
        relative to the velocity vector (defined as [0,0,1] if velocity is zero)
        '''
        self.rel_orientation = np.array([[0],[0],[1]])
        self.velocity = np.array([[0],[0],[0]])
        self.abs_position = np.array([[0],[0],[0]])
        '''
        TODO: Get reltaive position and velocities, this object compared to reference
        object
        '''
        self.relative_position = np.array([[None],[None],[None]])
        self.relative_velocity = np.array([[None],[None],[None]])


    '''
    Calculate the angle of attack (alpha) in radians using the rocket's velocity direction
    and rotation state.
    Return alpha
    '''
    def calculate_alpha(self):
        pass
    '''
    Use environmental data regarding gust velocity and rocket geometry to estimate
    the rotation axis and rotation magnitude (radians) of a rocket
    Return a rotation quaternion for this axis and magnitude
    '''
    def calculate_gust_rotation(self):
        pass
    '''
    Use the angle of attack, drag+lift coefficients and rocket geometry to
    estimate the rotation axis and rotation magnitude (radians) of a rocket
    Return a rotation quaternion for this axis and magnitude
    '''
    def calculate_drag_rotation(self):
        pass
    '''
    Place holder - calcultes tvc rotation
    Not needed for final version
    Return a rotation quaternion for this axis and magnitude
    '''
    def calculate_tvc_rotation(self):
        pass
    '''
    lock rotation of the craft despite forces acting on the body
    useful for constraining rocket to a rail at launch for example
    (Prevents integration of accelerations to velocities)
    '''
    def lock_rotation(self):
        pass
    '''
    Unlocks rotation
    '''
    def unlock_rotation(self):
        pass
class Environment():
    '''
    The environment object helps compartmentalize environmental data (atmospheric
    temperature, pressure, gusts etc...).  The object can then be accessed
    to fetch atmospheric or environmental data for the rotator object desired
    '''

    def __init__(self, rocket_position, seed, units='english'):
        self.EARTH_MASS_SLUG = 4.0948607276025*10**23
        self.EARTH_MASS_KG = 5.972 * 10**24
        self.EARTH_RADIUS_MI = 3958.8
        self.EARTH_RADIUS_FT = 20902000
        self.EARTH_RADIUS_M = 6371000
        self.units = units
        '''
        Use the arguments to fetch data from functions
        '''
        pass

    '''
    For these be sure to check which altitude you are working with.  For now I have
    it as altitude relative to center of the earth
    '''
    def calc_gravity(self, altitude):
        if self.units == 'english':
            altitude = altitude * 0.3048
        g = (self.EARTH_MASS_KG*6.67408*10**(-11))/(altitude**2)
        return g

    def fetch_atm_pressure(self, altitude):
        pass
    def fetch_atm_temperature(self, altitude):
        pass
    def fetch_atm_density(self, altitude):
        pass


if __name__ == '__main__':
    with open('rocket_info.yaml') as rocket_info:
        rocket_data = yaml.load(rocket_info, Loader=yaml.FullLoader)
    rot_tester = Rotator()
    rot_tester.report_body_vector()
    rot_quaternion = np.array([[-.707],[0], [.707],[0]])
    rot_tester.rotate_body(rot_quaternion)
    rot_tester.report_body_vector()

    rocenv = Environment(None, None)
    rocket = Rocket(rocket_data, rocenv)
