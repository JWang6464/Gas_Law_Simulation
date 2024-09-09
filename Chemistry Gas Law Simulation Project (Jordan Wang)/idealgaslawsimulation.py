'''
Jordan Wang
CS152 Section B
Professor Harper
05/05/2023
This module includes the simulation of the ideal gas law. It takes in the graphicsPlus, objects, random, collision modules
as well as matplotlib and csv to analyze the data produced at the end. In the beginning, several objects are initialized
and the window is constantly updated with a while loop. In the simulation, the user can interact through the particle
generator and temperature buttons as well as the arrow keys to move the moveable wall. The user can run this module/simulation
by typing idealgaslawsimulation.py into the terminal.
'''


import graphicsPlus as gr
import objects as pho
import random
import collision as col
import matplotlib.pyplot as plt 
import csv


def main():
    win = gr.GraphWin("PinParticle Simulation", 1400, 800)
    win.setBackground("light blue")
    walls = col.buildObstacles(win)
    moveable = [walls[col.IDXrightborder], walls[col.IDXhandlebottom], walls[col.IDXhandletop], walls[col.IDXhandleside]]
    for x in walls:
        x.draw()

    ##Temperature controller GUI objects
    Temperature = 273

    ## Objects and texts are initialized to an empty list
    objects = []
    texts = []

    ## Each object and list on the interface is added to the respective lists, which is then drawn on the interface
    objects.append(pho.Block(win, 4, 2, 100, 60, color = (255,48,48)))
    objects.append(pho.Block(win, 4, 2, 108, 60, color = (0,201,87)))
    temptext = gr.Text(gr.Point(1040, 160), 'Temperature = {} K'.format(Temperature))
    texts.append(temptext)
    particlegen = gr.Text(gr.Point(180, 320), 'Press generator to\nproduce a partice')
    texts.append(particlegen)
    barriertext = gr.Text(gr.Point(880, 400 ), 'Press left and right arrows to move\nwall and decrease total volume')
    texts.append(barriertext)
    quitinstructions = gr.Text(gr.Point(500, 650), 'Press q key to stop simulation and analyze temperature vs. collisions')
    texts.append(quitinstructions)
    title = gr.Text(gr.Point(500, 100), "Chemistry Ideal Gas Law Simulation")
    title.setSize(36)  # Set font size to 36
    title.setStyle("bold italic") # Sets the style to bold italic
    title.setFace("times roman") #Sets the Face/Font to times roman
    texts.append(title)


    for obj in objects:
        obj.draw()
    for t in texts:
        t.draw(win)

    total_collisions = 0

    ## writes the temperature and total collisions into a csv file
    with open('data.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Temperature', 'Total Collisions'])


    dt = 0.1
    frame = 0
    particles = []
    print("Press q to quit")

    ###while the simulation is occuring...
    while True:
        if frame % 10 == 0:
            win.update()
        win.update()
        
        ## moves the wall to the left
        keypressed = win.checkKey()
        if keypressed == 'Left':
            for obj in moveable:
                obj_position = obj.getPosition()
                x_old = obj_position[0]
                y_old = obj_position[1]
                x_new = x_old - 1
                y_new = y_old
                obj.setPosition(x_new, y_new)
        
        ## moves the wall to the right
        if keypressed == 'Right':
            for obj in moveable:
                obj_position = obj.getPosition()
                x_old = obj_position[0]
                y_old = obj_position[1]
                x_new = x_old + 1
                y_new = y_old
                obj.setPosition(x_new, y_new)


        ## gets the location of the mouse click and performs actions based on that location
        mouseclick = win.checkMouse()

        if mouseclick:
            mousex = mouseclick.getX()
            mousey = mouseclick.getY()
            print(mousex, mousey)

            if 130 <= mousex <= 230 and 350 <= mousey <= 450:
                newparticle = pho.Particle(win, 1, 29, 40, color=(255, 20, 147))
                particles.append(newparticle)
                walls[col.IDXparticlegun].setPosition(24, 40)
                newparticle.draw()
                newparticle.setVelocity(random.randint(1, 3), random.randint(-3, 3))
                walls[col.IDXparticlegun].setPosition(25, 40)

            ## alter the velocity of the particles if the red temp button is pressed
            if 980 <= mousex <= 1020 and 192 <= mousey <= 210:
                for particle in particles:
                    velocity = particle.getVelocity()
                    if particle.velocity[0] < 0 and particle.velocity[1] < 0:
                        particle.setVelocity(particle.velocity[0] + 0.5, particle.velocity[1] + 0.5)
                    if particle.velocity[0] < 0 and particle.velocity[1] > 0:
                        particle.setVelocity(particle.velocity[0] + 0.5, particle.velocity[1] - 0.5)
                    if particle.velocity[0] > 0 and particle.velocity[1] < 0:
                        particle.setVelocity(particle.velocity[0] - 0.5, particle.velocity[1] + 0.5)
                    else:
                        particle.setVelocity(particle.velocity[0] - 0.5, particle.velocity[1] - 0.5)
                    Temperature -= 1
                    temptext.setText(f'Temperature = {Temperature} K')

            ## alter the velocity of the particles if the green temp button is pressed
            if 1061 <= mousex <= 1100 and 192 <= mousey <= 211:
                for particle in particles:
                    velocity = particle.getVelocity()
                    if particle.velocity[0] < 0 and particle.velocity[1] < 0:
                        particle.setVelocity(particle.velocity[0] - 0.5, particle.velocity[1] - 0.5)
                    if particle.velocity[0] < 0 and particle.velocity[1] > 0:
                        particle.setVelocity(particle.velocity[0] - 0.5, particle.velocity[1] + 0.5)
                    if particle.velocity[0] > 0 and particle.velocity[1] < 0:
                        particle.setVelocity(particle.velocity[0] + 0.5, particle.velocity[1] - 0.5)
                    else:
                        particle.setVelocity(particle.velocity[0] + 0.5, particle.velocity[1] + 0.5)
                    Temperature += 1
                    temptext.setText(f'Temperature = {Temperature} K')


        for i, particle in enumerate(particles):
            # Check for particle-wall collisions
            for wall in walls:
                collided_with_wall = col.collision(particle, wall, dt)
                # print(f'in loop{collided_with_wall}')
                if collided_with_wall:
                    total_collisions += 1
                    # print('in loop: True')

            # Check for particle-particle collisions
            collided_with_particle = False
            for other_particle in particles[i+1:]:
                if col.collisionTest_Particle_Particle(particle, other_particle) < 1e+6:
                    collided_with_particle = True
                    break

            # Only update the particle if no collision occurred
            if not collided_with_wall and not collided_with_particle:
                particle.update(dt)


        with open('data.csv', 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([Temperature, total_collisions])


        frame += 10  # Increment frame
        
        ##pressing q ends the simulation and graphs the data obtained
        if keypressed == 'q':
            win.close()
            temperatures = []
            total_collisions = []

            with open('data.csv', mode='r') as data_file:
                data_reader = csv.reader(data_file)
                next(data_reader)  # Skip header row
                for row in data_reader:
                    temperatures.append(float(row[0]))
                    total_collisions.append(int(row[1]))

            plt.plot(temperatures, total_collisions)
            plt.xlabel('Temperature')
            plt.ylabel('Total Collisions')
            plt.show()
            break



if __name__ == "__main__":
    main()
