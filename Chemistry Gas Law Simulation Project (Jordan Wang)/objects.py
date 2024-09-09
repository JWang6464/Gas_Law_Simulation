'''
Jordan Wang
CS152 Section B
Professor Harper
05/05/2023
This module contains several different classes that contain attributes specific to those classes. In addition, each class
contains methods that serve to initialize and vary the preset values of the attributes. This module imports the graphicsPlus
and random modules and utilizes them to develop the visualization aspect of each shape/class as well as check for
potential collisions between objects.
'''

import graphicsPlus as gr
import random

class Thing():
    '''parent class for objects'''
    def __init__(self, win, the_type):
        '''Initializer for the objects/obstacles'''
        self.type= the_type
        self.mass = 1
        self.position = [0, 0]
        self.velocity = [0, 0]
        self.acceleration = [0, 0]
        self.elasticity = 1
        self.win = win
        self.scale = 10
        self.vis = []
        self.color = [0, 0, 0]
        self.drawn = False
    
    def getPosition(self):
        '''returns the position value'''
        return tuple(self.position[:])
    
    def setPosition(self, px, py):
        '''sets the Particles position value to the value inputed in px and py'''
        point = self.getPosition()
        #assign to old_x_position the current x position
        old_x_position = point[0]
        #assign to old_y_position the current y position
        old_y_position = point[1]
        #assign to the x coordinate in self.pos the new x position
        #assign to the y coordinate in self.pos the new y position
        self.position = [px, py]
        newpoint = self.getPosition()
        #assign to dx the change in the x position times self.scale
        dx = (newpoint[0] - old_x_position) * self.scale
        #assign to dy the change in the y position times -self.scale
        dy = (newpoint[1] - old_y_position) * -1 * self.scale
        #for each item in the vis field of self
        for item in self.vis:
            item.move(dx, dy)
    
    def getMass(self):
        '''return the mass value'''
        return float(self.mass)
    def setMass(self, m):
        '''method changes the mass of thing'''
        self.mass = m
    def getVelocity(self):
        '''method returns the velocity value'''
        return tuple(self.velocity[:])
    def setVelocity(self, vx, vy):
        '''method changes the velocity by inputing vx and vy(the velocity of x and y in the x and y in their perspective axis)'''
        self.velocity = [vx, vy]
    def getAcceleration(self):
        '''method returns the velocity value'''
        return tuple(self.acceleration[:])
    def setAcceleration(self, ax, ay):
        ''''method sets the acceleration inputing ax and ay'''
        self.acceleration = [ax, ay]
    def getElasticity(self):
        '''method returns the elasticity value'''
        return float(self.elasticity)
    def setElasticity(self, e):
        '''method changes the elasticity value'''
        self.elasticity = e
    def getScale(self):
        '''Method returns the scale of the object'''
        return float(self.scale)
    
    def getType(self):
        "This method returns the type of the object."
        return str(self.type)
    def getColor(self):
        '''Method gets the color of the object'''
        return tuple(self.color[:])
    def setColor(self, c): #takes in an (r,g,b) tuple
        '''Method sets the color in rgb'''
        self.color = c
        if c != None:
            color_name = gr.color_rgb(c[0], c[1], c[2])
            for x in self.vis:
                x.setFill(color_name)
    def draw(self):
        '''Method draws the shapes in the window'''
        vis = self.vis
        win = self.win
        for shape in vis:
            shape.draw(win)
        self.drawn = True
    
    def undraw(self):
        '''Method allows for the deletion of the shapes'''
        vis = self.vis
        for shape in vis:
            shape.undraw()
        self.drawn = False

    def update(self, dt):
        '''method updates the the thing and all sub things'''
        point = self.getPosition()
        vel = self.getVelocity()
        acc = self.getAcceleration()
        #assign to old_x_position the current x position
        old_x_position = point[0]
        #assign to old_y_position the current y position
        old_y_position = point[1]
        x_vel = vel[0]
        y_vel = vel[1]
        x_acc = acc[0]
        y_acc = acc[1]
        #update the x position to be old_x_position + x_vel*dt + 0.5*x_acc * dt*dt
        update_x_pos = old_x_position + x_vel * dt + 0.5 * x_acc * dt * dt
        self.position[0] = update_x_pos
        #update the y position to be old_y_position + y_vel*dt + 0.5*y_acc * dt*dt
        update_y_pos = old_y_position + y_vel * dt + 0.5 * y_acc * dt * dt
        self.position[1] = update_y_pos
        #assign to dx the change in the x position times the scale factor (self.scale)
        dx = (update_x_pos - old_x_position) * self.scale
        #assign to dy the negative of the change in the y position times the scale factor (self.scale)
        dy = (update_y_pos - old_y_position) * -1 * self.scale
        #for each item in self.vis
        for item in self.vis:
        #call the move method of the graphics object with dx and dy as arguments..
            item.move(dx, dy)
        #update the x velocity by adding the acceleration times dt to its old value
        update_x_vel = x_acc * dt + x_vel
        self.velocity[0] = update_x_vel
        #update the y velocity by adding the acceleration times dt to its old value
        update_y_vel = y_acc * dt + y_vel
        self.velocity[1] = update_y_vel


class Particle(Thing):
    '''This is a class for our object block'''

    def __init__(self, win, radius = 1, x0 = 0, y0 = 0, color = None):
        '''This initializer method takes a graphics window given as the win parameter and creates a default Block object with the width given as dx, the heigth given as the dy, positioned at (0, 0), velocity and acceleration at (0, 0), scale at 10, and a Rectangle 
        object'''
        Thing.__init__(self,win,'Particle')
        self.radius = radius
        self.position = [x0, y0]
        self.refresh()
        self.setColor(color)

    def refresh(self):
        drawn = self.drawn
        if drawn:
            self.undraw()

        win = self.win
        self.vis = [gr.Circle(gr.Point(self.position[0] * self.scale, win.getHeight() - (self.position[1] * self.scale)), self.radius * self.scale)]
        if drawn:
            self.draw()


    def getRadius(self):
        '''Method obtains the radius'''
        return float(self.radius)
    
    def setRadius(self, r):
        '''Method sets the radius of the object'''
        self.radius= r
        self.refresh()


class Block(Thing):
    '''This is the class for the block obstacles'''
    def __init__(self, win, width = 2, height = 1, x0=0, y0=0, color = None):
        Thing.__init__(self, win, "block")
        self.dx = width
        self.dy = height
        self.position = [x0, y0]
        self.reshape()
        self.setColor(color)
    
    def reshape(self):
        '''Method reshapes the block'''
        drawn = self.drawn
        if drawn:
            self.undraw()
        win = self.win
        self.vis = [gr.Rectangle(gr.Point((self.position[0] - (self.dx / 2)) * self.scale, win.getHeight() - ((self.position[1] - (self.dy / 2)) * self.scale)), gr.Point((self.position[0] + (self.dx / 2)) * self.scale, win.getHeight() - ((self.position[1] + (self.dy / 2)) * self.scale)))]
        if drawn:
            self.draw()
    
    def getWidth(self):
        '''Method gets the width'''
        return float(self.dx)
    
    def setWidth(self, w):
        '''Method sets the width'''
        self.dx = w
        self.reshape()

    def getHeight(self):
        '''Method gets the height'''
        return float(self.dy)
    
    def setHeight(self, h):
        '''Method sets the height'''
        self.dy = h
        self.reshape()

    def getPosition(self):
        '''Gets the position of the block'''
        return self.position
    
    def setPosition(self, px, py):
        point = self.getPosition()
        #assign to old_x_position the current x position
        old_x_position = point[0]
        #assign to old_y_position the current y position
        old_y_position = point[1]
        #assign to the x coordinate in self.pos the new x position
        #assign to the y coordinate in self.pos the new y position
        self.position = [px, py]
        newpoint = self.getPosition()
        #assign to dx the change in the x position times self.scale
        dx = (newpoint[0] - old_x_position) * self.scale
        #assign to dy the change in the y position times -self.scale
        dy = (newpoint[1] - old_y_position) * -1 * self.scale
        #for each item in the vis field of self
        for item in self.vis:
            item.move(dx, dy)

        