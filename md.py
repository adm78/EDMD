#!/usr/bin/python
# md.py - An Event Driven Molecular Dynamics (EDMD) Program
#
# This script performs a simple, event-drive molecular dynamics
# simulation on a pygame canvas
#
# Requires:
# - pygame
# - particle.py
# - event.py
#
# Andrew D. McGuire 2017
# a.mcguire227@gmail.com
#----------------------------------------------------------
import pygame
import numpy as np
from pygame.locals import *
import pygame.display
from particle import Particle
from event import Event


def setup(xmax,ymax,time,particles,r):

    ''' This function is called upon entry to create the
       simulation canvas which we draw onto and the particle array. 
       Canvas size and particle number is dependent on the window 
       size. '''
    
    part_to_init = int(round(xmax*ymax/3000.0))
    particles = initParticles(part_to_init,r,xmax,ymax)
    return particles

def simulate(xmax, ymax, particles, time):

    ''' This function drives the simulation forward in time.'''
    
    dt_step = 10.0
    particles, time =  doStep(particles, time, dt_step, xmax, ymax)
    return particles, time


def doStep(particles, time, dt, xmax, ymax):

    ''' Advances the particle ensemble over the time interval dt, or
       to the next collision time, whichever comes first.  If a
       collision is detected within (time,time+dt) then it's carried
       out and the sim time is updated.
       
       args:
       dt - time to try and advance simulation by
    '''

    # Compute the time to the next collision
    coll_list = getCollisionList(particles,xmax,ymax)
    if (len(coll_list) < 1):
      dt_col = dt + 1 # just needs to exceed dt
    else:
        dt_col = (coll_list[0]).t

    # Check for collisions in the current time
    if (dt < dt_col):
	# No collision in the time step
	particles = advanceParticles(particles, dt)
	time = time + dt
    else:
	# Collision has occured between the step
	# so, carry it out. Highlighting the particles
	# involved. 
	particles = advanceParticles(particles, dt_col)
        firstEvent = coll_list[0]
	particles = highlightEventParticles(firstEvent,particles)
	particles = performCollision(firstEvent,particles)
        time +=  dt_col
        
    return particles, time
   
def highlightEventParticles(CurrentEvent,particles):

    ''' Highlight the particle(s) involved in an event

       args:
       CurrentEvent - a valid Event object
    '''

    p1 = particles[CurrentEvent.p1_index]
    p1.highlight()
    if (CurrentEvent.p2_index != None):
        p2 = particles[CurrentEvent.p2_index]
        p2.highlight()
    return particles    

def getWallCollisionTime(Part,xmax,ymax):
    
    ''' Compute the first collision time with between
       particle Part and any wall 

    returns object attributes
    t  -  first collision time
    wall - wall associated with first collision

     locals vars
    t_side # side wall collision time
    t_ud   # top or bottom wall collision time
    w_side # which side wall ('r' or 'l')
    w_ud   # top or bottom wall first? ('u' ,d')
    '''

    # side walls
    if (Part.vel[0] > 0):
	t_side = (xmax - Part.pos[0] - Part.radius)/Part.vel[0]
	w_side = 'r'
    elif (Part.vel[0] < 0):
	t_side = (0 - Part.pos[0] + Part.radius)/Part.vel[0]
	w_side = 'l'
    else:
	# particle not moving in x direction
	t_side = None
	w_side = None
        

    # top and bottom
    if (Part.vel[1] > 0):
	t_ud = (ymax - Part.pos[1] - Part.radius)/Part.vel[1]
	w_ud = 'd'
    elif (Part.vel[1] < 0):
	t_ud = (0 - Part.pos[1] + Part.radius)/Part.vel[1]
	w_ud = 'u'
    else:
	# particle not moving in y direction
	t_ud = None
	w_ud = None
        
        
    if (t_side == None and t_ud == None):
	# part is stationary
	t = None
	wall= None
    elif (t_side <= t_ud):
	t = t_side
	wall = w_side
    else:
	t = t_ud
	wall = w_ud
        
    return type('obj', (object,),{'t': t, 'wall': wall})



def  getCollisionTime(Part1, Part2):

    ''' Compute the time until collision between particle Part1 and
       Part2.

       return time as None if no collision time solution found '''

    deltaVel = Part1.vel - Part2.vel
    deltaPos = Part1.pos - Part2.pos
    minDist = Part1.radius + Part2.radius
    a = np.dot(deltaVel, deltaVel)
    b = 2.0*np.dot(deltaPos,deltaVel)
    c = np.dot(deltaPos,deltaPos) - minDist*minDist
    discrim = b*b - 4*a*c

    if ((discrim > 0) and (b < 0)):
	t1 = (-b - (discrim**0.5))/(2*a)
	return t1
    
    return None


def getCollisionList(particles,xmax,ymax):
    
    ''' Returns an array of collision Event objects, ordered by their
    time attribute (smallest to largest, Nones at the end)

    args:
    particles - an array of Particle objects '''
    
    # return
    coll_list = []

    # loop through the particle array
    for i in range(len(particles)):

	wall_collision = getWallCollisionTime(particles[i],xmax,ymax)
	firstEvent = Event('w',wall_collision.t,i,None,wall_collision.wall)

	for j in range(i+1, len(particles)):
	    if (i != j):
		col_time = getCollisionTime(particles[i],particles[j])

		# Replace firstEvent if coll time is smaller than current
		# firstEvent.time
	        if col_time != None:
                    if (col_time < firstEvent.t):
			firstEvent = Event('p',col_time,i,j,None)
            
	
	# Add to the collision list if event is valid
	if (firstEvent.t != None):
            coll_list.append(firstEvent)
    

    # Sort the Event array and return it
    coll_list = sorted(coll_list,key=lambda event: event.t)
    return coll_list


def performCollision(event,particles):
    
    ''' Apply collision operator according according to event
       
       args:
       event - a valid Event object
    '''

    if (event.wc_log):
	# Perform wall collision
	if (event.wall == 'r' or event.wall == 'l'):
	    particles[event.p1_index].reflect_side()
	elif (event.wall == 'u' or event.wall == 'd'):
	    particles[event.p1_index].reflect_top()
	else:
	    raise RuntimeError("invalid collison event detected.")
    else:
	# Perform binary particle collision
	J = impulse(particles[event.p1_index],
		    particles[event.p2_index])
	particles[event.p1_index].apply_impulse(J[0],J[1])
	particles[event.p2_index].apply_impulse(-J[0],-J[1])
    return particles
        
def impulse(Part1,Part2):
    
    ''' Compute the impulse associated with a particle-particle collision
    https:#introcs.cs.princeton.edu/java/assignments/collisions.html

       J = 2*m1*m2*(dv*dr)/(sigma*(m1+m2))
    
       args:
       Part1 - valid Particle object
       Part2 - valid Particle object
    '''
    dr = Part2.pos - Part1.pos
    dv = Part2.vel - Part1.vel 
    sigma = Part1.radius + Part2.radius
    hmm = 2*Part1.mass*Part2.mass/(Part1.mass + Part2.mass)
    J = np.dot(dv,dr)*hmm/sigma
    return [J*dr[0]/sigma,J*dr[1]/sigma]


def advanceParticles(particles,dt):
    
    ''' Advance the ensemble forward in time by dt
    in a straight line trajectory (no collisions) '''
    
    for i in range(len(particles)):
      particles[i].updateProperties(dt)
    return particles

def initParticles(n,r,xmax, ymax):
    
    ''' Intialise n particles with radius r in box with dimensions
       (xmax,ymax) such that there are no overlapping particles

       return:
       particle - an array of Particle objects
    '''
       
    parts = []
    dx = initialSpacing(n, xmax, ymax)
    n_init = 0 
    for i in range(int(round(xmax/dx))):
	for j in range(int(round(ymax/dx))):
	    if (n_init < n):
    		parts.append(Particle(dx*(i+0.5),dx*(j+0.5),r))
		n_init += 1
    return parts


def initialSpacing(n, x, y):
    
    ''' Returns the intialise spacing between particles to put n
       particles on a uniform grid with limits x, y '''
    
    num1 = -(x+y)
    num2sqred = (x+y)**2.0 + 4.0*x*y*(n-1)
    num2 = num2sqred**0.5
    den = 2.0*(n-1)
    dx = (num1 + num2) / den
    return dx

def main():

    # define the simualtion parameters
    particles = []    # particle array (to be filled with instances of the Particle class)
    r = 7             # particle radius to use
    time = 0.0        # global simulation time
    xmax=300              # canvas x-width
    ymax=300              # canvas y-width
    paused_log = True # paused indicator bool

    # set-up the screen
    pygame.init()
    screen = pygame.display.set_mode((xmax, ymax))
    pygame.display.set_caption('EDMD')
    pygame.mouse.set_visible(1)
    background = pygame.Surface(screen.get_size())
    background = background.convert()
    background.fill((255, 255, 255))
    screen.blit(background, (0, 0))
    pygame.display.flip()
    clock = pygame.time.Clock()

    # set-up the particles
    particles = setup(xmax, ymax, time, particles,r)
    allsprites = pygame.sprite.RenderPlain(particles)

    print "EDMD simulation initialised with", len(particles), "particles."
    print "Press ESC or click the 'x' to end the simulation."
    print "Click anywhere on the screen to pause."

    #The main run loop
    quit_log = False
    paused_log = False
    while not quit_log:
        clock.tick(60)
        if not paused_log:
            particles, time = simulate(xmax, ymax, particles, time)

        #Handle Input Events
        for event in pygame.event.get():
            if event.type == QUIT:
                quit_log = True
            elif event.type == KEYDOWN and event.key == K_ESCAPE:
                quit_log = True
            elif event.type == MOUSEBUTTONDOWN:
                paused_log = not paused_log
            elif event.type == MOUSEBUTTONUP:
                # fist.unpunch()
                pass


        if not paused_log:
            allsprites.update()
        
            #Draw Everything
            screen.blit(background, (0, 0))
            allsprites.draw(screen)
            pygame.display.flip()
            

    pygame.quit()
    
if __name__ == '__main__':
    main()
