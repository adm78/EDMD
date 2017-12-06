# EDMD
## Event Driven Molecular Dynamics with Python

A simple, Python-based event-driven molecular dynamics simulation for hard-spheres. The simulation moves forward in time by calculating and moving to the time of the next 'event'. An event can be a collision between two particles, a particle and a wall or a graphical update checkpoint. The times associated with potential events are solved and placed in a priority queue. Particles move forward in time according to Newton’s Laws of Motion using a Velocity-Verlavet algorithm until the time of the first event in the priority queue is reached. If this is a collision event, then this is carried out and the queue is rebuilt based on the new system state. All collisions are treated as perfectly elastic. 

This is a Pythonic version of the JavaScript EDMD code that I authored for the [Visual Chemical Engineering (VCE) project](http://visualchemeng.com/). The JavaScript version is hosted on GitHub [here](https://github.com/adm78/visualchemeng_js/blob/master/modules/md/md.js). 

![edmd example image](https://user-images.githubusercontent.com/17439476/33520943-7bb99ff8-d7bc-11e7-9cd6-9e676f3a297d.png)

## Dependencies
- Python 2.*
- [Pygame](https://www.pygame.org/wiki/GettingStarted) (Graphics)
- [Numpy](http://www.numpy.org/)

## Usage
Default conditions
```shell
./md.py
```
Simulate 60 particles within a 500 x 500px simulation box and a maximum time step of size 10.0
```shell
./md.py -n 60 -x 500 -y 500 --dt 10.0
```
Display help and exit
```shell
./md.py --help
```

## Future work
- Intelligent reconstruction of the event queue. This is currently rebuilt after every event, which is inefficient but acceptable for small scale simulations.
- Addition of particle interaction potentials.
