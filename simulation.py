'''
@Authors: Harrison Leece, James Hribal, Max Fung, Nils Heidenreich
@Purpose: Simulate a rocket using the 6DOF.py classes
'''

from SIX_DOF import Rocket
from SIX_DOF import Environment
import oyaml as yaml
import numpy as np

if __name__ == '__main__':
    with open('rocket_info.yaml') as rocket_info:
        rocket_data = yaml.load(rocket_info, Loader=yaml.FullLoader)
    rocenv = Environment(None, None)
    rocket = Rocket(rocket_data, rocenv)
    #sim_obj1 is a list containing both a rocket and rocket environment object
    sim_obj1 = [rocket, rocenv]
    #object list -- because objects can be added to the list dynamically and simulated
    #(for example a first stage and second stage from one rocket object)
    object_list = [sim_obj1]
'''
    while (condition):
        simulation
'''
