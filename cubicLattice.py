from __future__ import division
from visual import *
import numpy as np

# http://vpython.org/

setting = display(range=(5,5,5))

# Lattice Constants
a = 1.0
k_ = 0.01
m = 5
q = 1.0

# Box
L = 4
W = 2
H = 3
numL = int(L/a)
numW = int(W/a)
numH = int(H/a)

# Time Constants
T = 10
t = 0 
dt = 0.001
Nt = int(T/dt)


# Lists to keep track of position, velocity and acceleration of all atoms
acc = [[[vector(0,0,0) for k in range(numH)] for j in range(numW)] for i in range(numL)]
vel = [[[vector(0,0,0) for k in range(numH)] for j in range(numW)] for i in range(numL)]
pos = [[[None for k in range(numH)] for j in range(numW)] for i in range(numL)]

# List for drawing all spheres 
atoms = [[[None for k in range(numH)] for j in range(numW)] for i in range(numL)]

# Lists for all springs
xsprings = [[[None for k in xrange(numH)] for j in xrange(numW)] for i in xrange(numL-1)]
ysprings = [[[None for k in xrange(numH)] for j in xrange(numW-1)] for i in xrange(numL)]
zsprings = [[[None for k in xrange(numH-1)] for j in xrange(numW)] for i in xrange(numL)]


# Place Spheres (and fill pos list)
for i in range(numL):
	for j in range(numW):
		for k in range(numH):
			curPos = vector(i*a,j*a,k*a)
			atoms[i][j][k] = sphere(pos = curPos,radius = 0.2, color = (i/numL+0.1,j/numW+0.1, k/numH+0.1))
			pos[i][j][k] = curPos

#____For Fourier Analysis_____#
# History of spring's stretch
xHist = []
yHist = []
zHist = []

# Place Springs
for i in range(numL):
	for j in range(numW):
		for k in range(numH):
			if i < numL-1:
				curPos = pos[i][j][k]
				adjPos = pos[i+1][j][k]
				xsprings[i][j][k] = helix(pos = curPos,radius = 0.05, thickness = 0.05, color = (0,1,1), axis = adjPos-curPos)
				
			if j < numW-1:
				curPos = pos[i][j][k]
				adjPos = pos[i][j+1][k]
				ysprings[i][j][k] = helix(pos = curPos,radius = 0.05, thickness = 0.05, color = (1,0,1), axis = adjPos-curPos)
				
			if k < numH-1:
				curPos = pos[i][j][k]
				adjPos = pos[i][j][k+1]
				zsprings[i][j][k] = helix(pos = curPos,radius = 0.05, thickness = 0.05, color = (1,1,0), axis = adjPos-curPos)

# Perturbation
pos[0][0][0] += vector(0.1,0.13,0.1)
pos[1][0][0] += vector(0.3,-0.1,0.2)

RCM = vector(0,0,0)

# Add some Newtonian Physics
# F = -k*(x_self - x_other)
while True:
	rate(100)

	setting.center = RCM
	RCM = vector(0,0,0)

	for i in range(numL):
		for j in range(numW):
			for k in range(numH):
				
				Fx = vector(0,0,0)
				Fy = vector(0,0,0)
				Fz = vector(0,0,0)
				
				ax = vector(a,0,0)
				ay = vector(0,a,0)
				az = vector(0,0,a)
				
				curPos = pos[i][j][k]
				
				# Calculate forces
				if i < numL-1:
					nextX = pos[i+1][j][k]
					Fx += (curPos-nextX)+ax
					
				if i > 0:
					prevX = pos[i-1][j][k]
					Fx += (curPos-prevX)-ax
					
				if j < numW-1:
					nextY = pos[i][j+1][k]
					Fy += (curPos-nextY)+ay
					
				if j > 0:
					prevY = pos[i][j-1][k]
					Fy += (curPos-prevY)-ay
					
				if k < numH-1:
					nextZ = pos[i][j][k+1]
					Fz += (curPos-nextZ)+az
					
				if k > 0:
					prevZ = pos[i][j][k-1]
					Fz += (curPos-prevZ)-az
				
				Fnet = -k_ * (Fx + Fy + Fz)
				
				# Eularate
				acc[i][j][k] = Fnet/m
				vel[i][j][k] += acc[i][j][k]
				pos[i][j][k] += vel[i][j][k]
				
				curPos = pos[i][j][k]
	
				# Update object positions
				atoms[i][j][k].pos = curPos
				
				if i < numL-1:
					curPos = pos[i][j][k]
					adjPos = pos[i+1][j][k]
					xsprings[i][j][k].pos = curPos
					xsprings[i][j][k].axis = adjPos-curPos
					
				if j < numW-1:
					curPos = pos[i][j][k]
					adjPos = pos[i][j+1][k]
					ysprings[i][j][k].pos = curPos
					ysprings[i][j][k].axis = adjPos-curPos
					
				if k < numH-1:
					curPos = pos[i][j][k]
					adjPos = pos[i][j][k+1]
					zsprings[i][j][k].pos = curPos
					zsprings[i][j][k].axis = adjPos-curPos

				# Rcm = (1/M)*sum(r_i*m_i)
				#RCM += (curPos*m)/(L*W*H)