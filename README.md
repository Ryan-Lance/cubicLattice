This is a program I made using the Vpython module, it's a mass and spring cubic lattice structure, where an arbitray number of masses may be added.
 
Right now, it uses Euler's method to update the positions of all the masses. My goal is to switch to Verlet's method, but it's been a while since I've worked on it so making changes isn't easy. I might end up keeping the Euler version as is, then basically restarting with Verlet. We'll see how that goes.

Goals for this project:
- Switch to Verlet method of integration (this might fix the flying away problem)
- Store the history of the mass's positions in arrays on which I can then do Fourier Transforms
- Add a GUI (maybe write using p5 library?)
