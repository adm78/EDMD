#!/usr/bin/python
# md.py - An Event Driven Molecular Dynamics (EDMD) Program
#
# This script performs a simple, event-drive molecular dynamics
# simulation on a pygame canvas
#
# Requires:
# - ???
# - particle.py
# - event.py
#
# Andrew D. McGuire 2017
# a.mcguire227@gmail.com
#----------------------------------------------------------
import pygame
import numpy as np
from pygame.locals import *
from particle import Particle
from event import Event
import pygame.display


# pygame.display.init()
# info = pygame.display.Info()
# print info

# # pygame.init()
# print "hello world"

# test_Particle = Particle(10,10,2)

# pygame.display.set_mode((640, 480))
particles = []    # particle array (to be filled with instances of the Particle class)
r = 5             # particle radius to use
time = 0.0        # global simulation time
xmax=300              # canvas x-width
ymax=300              # canvas y-width
paused_log = True # paused indicator bool
 
def setup(xmax,ymax,time,particles):

    ''' This function is called upon entry to create the
       simulation canvas which we draw onto and the particle array. 
       Canvas size and particle number is dependent on the window 
       size. '''
    
    #canvas= createCanvas(xmax, ymax)
    part_to_init = int(round(xmax*ymax/5000.0))
    print "part_to_init = ", part_to_init
    particles = initParticles(part_to_init,r,xmax,ymax)
    return particles

def simulate(xmax, ymax, particles, time):

    ''' This function drives the simulation forward in time. 
       It's continuously called for the
       lifetime of the scripts executions after setup()
       has completed. '''

    # set the background-up
    # background(51)
    # stroke(255)
    # strokeWeight(1)

    # # draw the particle to the canvas
    # for (i = 0 i < particles.length i++) {
    #          particles[i].show()
    # }
    
    # # set up stroke for progress box
    # noStroke()
    # fill(0)
    # rect(0.9*xmax,0.9*ymax,0.1*xmax,0.1*ymax)
    # fill(255)

    # # Step through time unless sim is paused,
    # # reporting status in progress box.
    # if (!(paused_log)) {
    #   text("Running",0.91*xmax,0.9*ymax,0.2*xmax,0.1*ymax)

    dt_step = 1.0
    particles, time =  doStep(particles, time, dt_step, xmax, ymax)

    # }
    # else {
    #   text("Paused",0.91*xmax,0.9*ymax,0.2*xmax,0.1*ymax)
    # }
    writeTime(time)
    return particles, time

def writeTime(time):

    # Write the current simulation time to the process box
    # stroke(255)
    # strokeWeight(1)
    # fill(255)
    # text(time.toFixed(0),0.91*xmax,0.95*ymax,20,20)
    print "time =",time


def doStep(particles, time, dt, xmax, ymax):

    ''' Advances the particle ensemble over the
       time interval dt, or to the next collision time,
       whichever comes first.
       If a collision is detected within (time,time+dt)
       then it's carried out and the sim time is updated.
       
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
	advanceParticles(dt)
	time = time + dt
    else:
	# Collision has occured between the step
	# so, carry it out. Highlighting the particles
	# involved. 
	advanceParticles(dt_col)
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

    ''' Compute the time until collision 
       between particle Part1 and Part2.

       return time as None if no collision 
       time solution found '''

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
	    print("Error: performCollision: invalid event")
	    print(event)
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


def advanceParticles(dt):
    
    ''' Advance the ensemble forward in time by dt
    in a straight line trajectory (no collisions) '''
    
    for i in range(len(particles)):
      particles[i].update(dt)
      

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
	        parts[n_init].show()
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


def addParticle():

    # Add a new Particle object to the particles array
    # new_part = Particle(mouseX,mouseY,r)
    # particles.append(new_part)
    pass

def mousePressed():

    # Act on left mouse press
    paused_log = not paused_log

# testing
particles = setup(xmax, ymax, time, particles)
while True:
    particles, time = simulate(xmax, ymax, particles, time)
