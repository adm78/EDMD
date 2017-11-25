# Particle.py - Particle class
#
# A simple class to store particle information and apply
# various transformations to simulate particle movement.
#
# Requires:
# - ???
#
# Andrew D. McGuire 2017
# a.mcguire227@gmail.com
#----------------------------------------------------------

import numpy as np
import random

class Particle(object):

    def __init__(self,x,y,r):

        '''
        Initialise the particle with a random velocity and zero
        acceleration.
        
        args:
        x - particle initial x position
        y - particle in ital y position
        r - particle radius     
        '''

        # particle attributes
        self.pos = np.array([x,y])
        self.radius = r
        self.vel = np.array([random.random()*2.0-1.0,
                             random.random()*2.0-1.0])
        self.acc = np.array([0,0])
        self.acc_old = self.acc
        self.mass = 1

    # particle methods
    def update(self,dt):

        ''' Compute the new acceleration, position and velocity
        vectors using the Velocity Verlet algorithm with time-step
        dt. Constant acceleration is assumed for the moment '''
        
        self.update_acc(dt)
	self.update_pos(dt)
	self.update_vel(dt)
         
    
    def show(self):
        
        ''' Draw the particle as an circle on the canvas. Size is
	controlled by the particle radius.  '''
        
        pass
    


    def highlight(self):
        
	# Highlight the particle red.
	
        pass
    

    def update_acc(self,dt):

        ''' Update acceleration term. Constant acceleration for now so
	no change.  '''
        
	self.acc_old = self.acc
	self.acc = self.acc
    

    def update_pos(self,dt):

        ''' Update the particle position according to a time-step of
	size dt.  '''
        
	pos_1 = self.pos
	pos_1 += dt*self.vel
	pos_2 = self.acc*0.5*(dt**2.0)
	self.pos = pos_1 + pos_2

    

    def update_vel(self,dt):

        ''' Update the velocity according to and its accelerate over
	time interval dt.  '''
	
	v_1 = self.vel
	a_sum = self.acc + self.acc
	v_2 = a_sum*0.5*dt
	self.vel = v_1 + v_2
    

    def reflect_side(self):
	#Invert the x velocity component
	self.vel[0] = - self.vel[0]
    

    def reflect_top(self):
        #Invert the y velocity component
        self.vel[1] = - self.vel[1];
        

    def apply_boundary_cond(self):

	'''Update the particle position to simulate periodic boundary
	   conditions with box bounds (0,xmax), (0,ymax)'''
	
	if (self.pos[0] >= xmax):
	    self.pos[0] = self.pos[0] - xmax
            
	if (self.pos[0] < 0):
	    self.pos[0] = self.pos[0] + xmax

	if (self.pos[1] >= ymax):
	    self.pos[1] = self.pos[1] - ymax
	
	if (self.pos[1] < 0):
	    self.pos[1] = self.pos[1] + ymax
	

    def apply_impulse(self,Jx,Jy):

	''' Compute and apply velocity change due to an impulse
	   applied to the particle.
	   
	   args:
	   Jx - scalar x-component of the impulse vector
	   Jy - scalar y-component of the impulse vector 
	'''
	
	self.vel[0] = self.vel[0] + (Jx/self.mass)
	self.vel[1] = self.vel[1] + (Jy/self.mass)
