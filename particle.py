# Particle.py - Particle class
#
# A simple class to store particle information and apply
# various transformations to simulate particle movement.
#
# Requires:
# - pygame
#
# Andrew D. McGuire 2017
# a.mcguire227@gmail.com
#----------------------------------------------------------
import os
import numpy as np
import random
from pygame.locals import *
import pygame.display

class Particle(pygame.sprite.Sprite):

    def __init__(self,x,y,r):

        '''
        Initialise the particle with a random velocity and zero
        acceleration.
        
        args:
        x - particle initial x position
        y - particle in ital y position
        r - particle radius     
        '''
        pygame.sprite.Sprite.__init__(self) #call Sprite initializer

        # load the particle images
        self.image, self.rect = load_image('particle.png', -1)
        self.orig_image = self.image
        self.alt_image, self.alt_rect  = load_image('particle2.png', -1)
        screen = pygame.display.get_surface()
        self.area = screen.get_rect()
        
        # set particle attributes
        self.pos = np.array([x,y])
        self.rect.topleft = x, y
        self.radius = r
        self.vel = np.array([random.random()*2.0-1.0,
                             random.random()*2.0-1.0])
        self.acc = np.array([0,0])
        self.acc_old = self.acc
        self.mass = 1
        self.highlighted = False
                

    # particle methods
    def update(self):
        self.rect.topleft = self.pos       
    
    def updateProperties(self,dt):

        ''' Compute the new acceleration, position and velocity
        vectors using the Velocity Verlet algorithm with time-step
        dt. Constant acceleration is assumed for the moment '''
        
        self.update_acc(dt)
	self.update_pos(dt)
	self.update_vel(dt)
        self.unhighlight()
         

    def highlight(self):
        
	# Highlight the particle red.
	self.image = self.alt_image

    def unhighlight(self):

        # return particle colour to original state
        self.image = self.orig_image

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


# define the image loader
def load_image(name, colorkey=None):
    
    data_dir = os.getcwd()
    fullname = os.path.join(data_dir, name)
    
    try:
        image = pygame.image.load(fullname)
    except pygame.error:
        print 'Cannot load image:', fullname
        raise SystemExit(str(geterror()))
    
    image = image.convert()
    if colorkey is not None:
        if colorkey is -1:
            colorkey = image.get_at((0,0))
        image.set_colorkey(colorkey, RLEACCEL)
        
    return image, image.get_rect()
