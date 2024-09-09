# Bruce A. Maxwell
# CS 151S
# Fall 2015
#
# Project 9
# Collision handler
#
# modified slightly by Eric Aaron for CS 152, Spring '19
# updated by Bruce Maxwell for CS 152, Fall 2019
#   fixed an issue in Particle-block collision that would return False when the Particle was inside the block
#   updated the Particle-Particle collisions to use the proper ray-circle intersections
#

'''
Jordan Wang
CS152 Section B
Professor Harper
05/05/2023
This module contains methods that account for the colisions of objects within the ideal gas law simulator.
Overall, this program consists of a several different components that ultimately make a gas container
simulation. In order to run this code, the user can enter collision.py into the terminal or press the run
button in Visual Studio Code.
'''

import math
import graphicsPlus as gr
import objects as pho
import random

def length(v):
    """utility math function for calculating Euclidean length of a 2D vector"""
    return math.sqrt(v[0]*v[0] + v[1]*v[1])

def unit(v):
    """utility math function for creating a unit 2D vector"""
    l = math.sqrt(v[0]*v[0] + v[1]*v[1])
    if l > 0.0:
        return (v[0]/l, v[1]/l)
    return v


def collisionTest_Particle_wall( Particle, wall ):
    """Tests if there is a collision with the wall along the path of 
    the Particle. Returns the distance to the collision or 1e+6 (a big number)"""

    # get the Particle's velocity and current position
    v = unit( Particle.getVelocity() )
    ParticleP = Particle.getPosition()

    # get the position of the floor
    wallP = wall.getPosition()

    # a variation on Liang-Barsky clipping
    p1 = -v[0]
    p2 = v[0]
    if ParticleP[0] < wallP[0]:
        # Particle is to the left of the wall, so use the left boundary and Particlex + radius
        q1 = (ParticleP[0]+ Particle.getRadius()) - (wallP[0] - wall.getWidth()*0.5)
        q2 = (wallP[0] - wall.getWidth()*0.5) - (ParticleP[0] - Particle.getRadius())
    else:
        # Particle is on the right, so subtract radius and add wall width
        q1 = (ParticleP[0] - (Particle.getRadius())) - (wallP[0] + wall.getWidth()*0.5)
        q2 = (wallP[0] + wall.getWidth()*0.5) - (ParticleP[0] - (Particle.getRadius()))

    # running parallel to the wall, no collision for a stationary wall
    if p1 == 0.0: 
        return 1e+6

    if p1 < 0: # Particle is heading in a +y direction
        if q1 > 0: # Particle is headed away from the wall
            return 1e+6
        else: # Particle is headed towards the wall
            return q1 / p1
    else: # Particle is heading in a -y direction
        if q2 > 0: # Particle is headed away from the wall
            return 1e+6
        else:
            return q2/p2
    

def collision_Particle_wall(Particle, wall, dt):
    """Main collision function for handling Particle/wall collisions
    Updates the Particle's position and returns true if there was a collision. 
    Returns False if there was no collision (Particle still needs to be udpated). """

    # returns the distance between the Particle and the wall
    tk = collisionTest_Particle_wall( Particle, wall )
    
    d = length(Particle.getVelocity())
    if d == 0.0: # special case if the Particle is not moving at all
        return False

    # check if the collision will happen during dt
    delta = tk / (d*dt)
    if delta <= 1.0:

        # update to the collision
        Particle.update( delta*dt )

        # change the velocity
        v = Particle.getVelocity()
        be = Particle.getElasticity()
        fe = wall.getElasticity()
        Particle.setVelocity( -v[0]*be*fe, v[1] )

        # update after the collision
        Particle.update(dt - delta*dt)

        return True

    # if no collision, calling function handles the update
    return False

def collisionTest_Particle_Particle(Particle1, Particle2):
    """Tests if there is a collision with another Particle along the path of
    the Particle.  Returns the distance to the collision or 1e+6 (a big number)"""
    
    # Concept: hold Particle2 still and test if Particle1 will hit it
    # Ray-circle intersection
    v1 = Particle1.getVelocity()
    p1 = Particle1.getPosition()
    p2 = Particle2.getPosition()
    r = Particle1.getRadius() + Particle2.getRadius()

    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]

    # quadratic equation
    a = v1[0]*v1[0] + v1[1]*v1[1]
    b = dx*v1[0] + dy*v1[1]
    c = dx*dx + dy*dy - r*r

    delta = b*b - a*c
    
    if  delta <= 0: # no intersection, imaginary roots
        return 1e+6

    deltaroot = math.sqrt(delta)
    t1 = (-b + deltaroot) / a
    t2 = (-b - deltaroot) / a

    if t1 < 0 and t2 < 0: # intersection is behind the Particle
        return 1e+6

    # one of these could be negative
    tmin = min(t1, t2)

    # Particle's already intersect, so move Particle1 back to the boundary
    if t1 < 0 or t2 < 0:
        newpx = p1[0] + tmin*v1[0]
        newpy = p1[1] + tmin*v1[1]
        Particle1.setPosition( newpx, newpy )
        return 0.0
        
    dx = tmin*v1[0]
    dy = tmin*v1[1]
    distToImpact = math.sqrt(dx*dx + dy*dy)

    return distToImpact
    

def collision_Particle_Particle(Particle1, Particle2, dt):
    """Main collision function for handling Particle/Particle collisions
    Updates the Particle1's position and returns true if there was a collision.
    Returns False if there was no collision (Particle1 still needs to be udpated).
    Particle2's velocity is changed, but it is not updated by this function"""

    # holds Particle2 steady, tests Particle1's trajectory
    distToImpact = collisionTest_Particle_Particle(Particle1, Particle2)

    # get the magnitude of the velocity of Particle1
    v1 = Particle1.getVelocity()
    v2 = Particle2.getVelocity()
    vmag1 = length( Particle1.getVelocity() )

    # no collision if it's too far away
    if distToImpact > vmag1 * dt:
        return False

    # check for a stationary Particle
    if vmag1 < 1e-6:
        # just update and return a collision
        Particle1.update( dt )
        return True

    delta = distToImpact / (vmag1*dt)

    # don't update backwards, which can happen, strangely enough
    if delta > 0.0:
        Particle1.update( delta*dt )

    p1 = Particle1.getPosition()
    p2 = Particle2.getPosition()
    rvec = unit( (p1[0] - p2[0], p1[1] - p2[1]) ) # reflection vector

    # create the reflection matrix R(th)M(X)R(-th)

    # update Particle1's velocity
    # rotate reflection vector to the Y axis
    tvx =  rvec[0] * v1[0] + rvec[1] * v1[1]
    tvy = -rvec[1] * v1[0] + rvec[0] * v1[1]

    # mirror in X
    tvx = - tvx*Particle1.getElasticity()*Particle2.getElasticity() # need to add the loss factor here

    # rotate back
    vfx = rvec[0] * tvx - rvec[1] * tvy
    vfy = rvec[1] * tvx + rvec[0] * tvy

    Particle1.setVelocity( vfx, vfy )

    # update Particle2's velocity
    tvx =  rvec[0] * v2[0] + rvec[1] * v2[1]
    tvy = -rvec[1] * v2[0] + rvec[0] * v2[1]

    # mirror in X
    tvx = - tvx*Particle1.getElasticity()*Particle2.getElasticity() # need to add the loss factor here

    # rotate back
    vfx = rvec[0] * tvx - rvec[1] * tvy
    vfy = rvec[1] * tvx + rvec[0] * tvy

    Particle2.setVelocity( vfx, vfy )

    # finish updating Particle1
    if delta > 0.0:
        Particle1.update( dt - delta*dt )
    else:
        Particle1.update( dt )

    return True


def collisionTest_Particle_block(Particle, block):
    """Test if a Particle is colliding with any side of a block, and indicate
    which side. Sends out a line along the Particle's velocity vector and
    compares it with all four sides of the object."""

    # get the trajectory and position of the Particle
    v = unit( Particle.getVelocity() )
    ParticleP = Particle.getPosition()
    radius = Particle.getRadius()

    # get the position of the block
    blockP = block.getPosition()

    # a variation on Liang-Barsky clipping
    # expands the block by the size of the Particle before testing
    dx = block.getWidth()
    dy = block.getHeight()

    p = ( -v[0], v[0], -v[1], v[1] )
    q = (ParticleP[0] - (blockP[0] - dx*0.5 - radius),
         (blockP[0] + dx*0.5 + radius) - ParticleP[0],
         ParticleP[1] - (blockP[1] - dy*0.5 - radius),
         (blockP[1] + dy*0.5 + radius) - ParticleP[1] )


    # for all four cases
    tmin = -1e+6
    tmax = 1e+6
    side = -1
    sidemax = -1
    for i in range(4):
        if p[i] == 0.0: # no collision for this side of the block, motion is parallel to it
            if q[i] < 0: # outside the boundary of the box, no collision
                return 1e+6,0
            continue

        tk = q[i] / p[i]

        if p[i] < 0: # outside moving in
            if tk > tmin:
                tmin = tk
                side = i
        else:
            if tk < tmax:
                tmax = tk
                sidemax = i

        if tmax <= tmin: # no intersection with the box
            return 1e+6,0

    if tmin < 0 and tmax < 0: # both intersections behind the Particle
        tmin = 1e+6
    elif tmin < 0 and tmax > 0: # Particle is intersecting the block
        #print("Particle is intersecting")

        # move the Particle back along its velocity to the intersection point
        if v[0] == 0.0 and v[1] == 0.0:
            v = (1.0, 0.0)
            tmin = (blockP[0] - 0.5*dx - radius) - ParticleP[0]

        # move it to the closest side and set distance to impact to 0
        if -tmin < tmax:
            #print("setting position using tmin and velocity %.2f %d" % (tmin, side))
            Particle.setPosition( ParticleP[0] + (tmin+1e-3)*v[0], ParticleP[1] + (tmin+1e-3)*v[1])
        else:
            #print("setting position using tmax and velocity %.2f %d" % (tmax, sidemax))
            Particle.setPosition( ParticleP[0] + (tmax+1e-3)*v[0], ParticleP[1] + (tmax+1e-3)*v[1])
        tmin = 0

    # tmin is the closest intersection on side i
    # 0: coming up from below
    # 1: coming down from above
    # 2: coming from the left
    # 3: coming from the right
    return (tmin, side)


def collision_Particle_block(Particle, block, dt):
    """Main collision code for Particle/block interactions.
    Updates the Particle's position and returns true if there was a collision.                                     
    Returns False if there was no collision (Particle still needs to be udpated)."""

    # get distance to impact
    distToImpact, side = collisionTest_Particle_block( Particle, block )

    # check if the impact is farther away than one step
    vmag = length( Particle.getVelocity() )
    if vmag == 0.0 or distToImpact > vmag * dt:
        return False

    # update the Particle prior to the collision
    delta = distToImpact / (vmag * dt)
    Particle.update( delta * dt )

    # modify the velocities
    v = Particle.getVelocity()
    if side == 0 or side == 1: # left or right wall, so adjust x
        Particle.setVelocity( -v[0]*Particle.getElasticity()*block.getElasticity(), v[1]  )
    elif side == 2 or side == 3: # top or bottom wall, so adjust y
        Particle.setVelocity( v[0], -v[1]*Particle.getElasticity()*block.getElasticity()  )

    # update the Particle post-collision
    Particle.update( (1 - delta) * dt )

    return True


collision_router = {}
collision_router[ ('Particle', 'Particle') ] = collision_Particle_Particle
collision_router[ ('Particle', 'block') ] = collision_Particle_block

def collision(Particle, thing, timestep):
    '''Method checks for collision'''

    return collision_router[('Particle',thing.getType())](Particle,thing,timestep)

def buildObstacles(win):
    '''This function creates the border and obstacles'''
    leftborder = pho.Block(win, 2, 40, 25, 40, color = (10, 10, 10))
    rightborder = pho.Block(win, 2, 40, 75, 40, color = (10, 10, 10))
    topborder = pho.Block(win, 50, 2, 50, 59, color = (10, 10 , 10))
    bottomborder = pho.Block(win, 50, 2, 50, 21, color = (10, 10, 10))
    particlegun = pho.Block(win, 5, 2 , 25, 40, color = (72, 61, 139))
    particlestorage = pho.Block(win, 10, 10, 18, 40, color = (72, 61, 139))
    handletop = pho.Block(win, 2, 1, 77, 42, color = (191,239,255))
    handlebottom = pho.Block(win, 2, 1, 77, 38, color = (191,239,255))
    handleside = pho.Block(win, 1, 3, 77.5, 40, color = (191,239,255))
    wall = [leftborder, rightborder, topborder, bottomborder, particlegun, particlestorage, handletop, handlebottom, handleside]
    # Return the list of Things



    def handle_mouse_move(event):
        global handle_dragging, handle_drag_offset
        handle_dragging = False  # Flag to track handle dragging
        handle_drag_offset = 0  # Offset between mouse and handle position during dragging

        if handle_dragging:
            new_y = event.y - handle_drag_offset
            delta_y = new_y - handleside.getCenter().getY()

            # Move the right border
            rightborder.move(0, delta_y)

            # Move the handle top and handle bottom
            handletop.move(0, delta_y)
            handlebottom.move(0, delta_y)

            # Move the handle side
            handleside[8].move(0, delta_y)


    win.bind("<B1-Motion>", handle_mouse_move)


    return wall

IDXleftborder = 0
IDXrightborder = 1
IDXtopborder = 2
IDXbottomborder = 3
IDXparticlegun = 4
IDXparticlestorage = 5
IDXhandletop = 6
IDXhandlebottom = 7
IDXhandleside = 8

