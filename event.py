# Event.py - Event class
#
# A simple collision event class that can store information about the
# type, time and particles involved in a collision event.
#
# Andrew D. McGuire 2017
# a.mcguire227@gmail.com
#----------------------------------------------------------

class Event(object):
    def __init__(self,ct,t,p1_index,p2_index=None,wall=None):

        '''
        Events are either with the wall or another particle
        
        args:
        ct - collision type ("w" is wall collision, anything else is
        binary particle)
        t        - time of collision
        p1_index - location parameter for the first particle
        p2_index - location parameter for the second particle
        (optional)
        wall - name of wall involved in collision
        '''
        
        # Event attributes
        self.wc_log = False # is it wall collision? (otherwise binary particle col)
        self.t = t
        self.p1_index = p1_index
        self.p2_index = p2_index
        self.wall = wall   
        
        # Set event attributes based on the collision type
        if (ct=="w"):
	    self.wc_log = True
        else:
	    # Warn if second particle index was not passed
	    if (p2_index == None):
	        raise RuntimeWarning("Warning: Event: second particle index undefined")

