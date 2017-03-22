## PROGRAM INFO ##
# ERICH MCMILLAN
# 2017_03_16
# This program will find the shortest and longest routes of a spaceship through a blackhole field, i.e.
# a kessel run. The input for the program will be a file containing the information on the masses, and
# the x and y positions of the blackholes on the playing field. The x, y coordinates should lie within
# -10<=x<=10, and -10<=y<=10. A starting position is chosen along -5<x<5 and y = -10 with a uniform
# distrbution. A starting velocity is chosen along a normally distributed angle with a mean of pi/2 and
# a standard deviation of pi/4, the magnitude of the velocity is uniformly chosen between 2 and 5. Next
# the program steps through the run iteratively using a small discrete time step and recalculating the
# force/acceleration/velocity/position at each step. If the force exceeds 4[Newtons] the run is
# terminated as unsuccessful. After a few thousand iterations the best and worst sucessful runs are
# plotted on a graph and saved to a .png in the current directory.

##	Import science modules	##
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import math as math
from collections import namedtuple


## 	Define Classes			##
# used for setting the parameters for the starting configuration in the kessel run
runParameters = namedtuple("runParameters",
	"minX, maxX, startY, finalY, avgVelAngle, stdVelAngle, minVel, maxVel maxForce resolution maxTime shipMass" )
# contains the information about the ship and its path
class kesselRun:
	# kesselRun constructor
	def __init__(self):
		self.pos = [-1,-1]
		self.path = []
		self.vel = [-1,-1]
		self.acc = [-1,-1]
		self.time = 0
		self.success = False
		self.terminated = False
		self.mass = -1

	# update functions
	def setMass(self, mass):
		self.mass = mass
	def updatePosition(self, position):
		self.path.append(position)
		self.pos = position
		# print("position")
		# print(self.pos)
	def updateVelocity(self, velocity):
		self.vel = velocity
		# print("velocity")
		# print(self.vel)
	def updateAcceleration(self, acceleration):
		self.acc = acceleration
		# print("acceleration")
		# print(self.acc)
	def updateTime(self, deltat):
		self.time += deltat
	def stepRungeKuttaIntegration(self, blackholes, runParameters):
		# Definition:

		# declare variables
		forceonship = [0,0]

		# determine force upon the ship
		for currentbh in blackholes:
			forceonship = np.add(forceonship, currentbh.calculateForce(self.pos, self.mass))

		# print(forceonship)
		# print(np.linalg.norm(forceonship))
		# print((np.linalg.norm(forceonship) > runParameters.maxForce))
		# check if force has exceeded specified max if so terminate run
		self.terminated = (np.linalg.norm(forceonship) > runParameters.maxForce) if True else self.terminated
		if(self.terminated):
			 return

		# update acceleration
		self.updateAcceleration(np.divide(forceonship, self.mass))

		# update velocity
		self.updateVelocity(np.add(self.vel, np.multiply(self.acc, runParameters.resolution)))

		# update position
		self.updatePosition(np.add(self.pos, np.multiply(self.vel, runParameters.resolution)))

		# update time
		self.updateTime(runParameters.resolution)

		# check if ship's y is greater than finalY if so then sucessful run
		self.success = (self.pos[1] > runParameters.finalY) if True else self.sucess

		# print(self.pos[0])
		# print((self.pos[0] < -10.0 or self.pos[0] > 10.0))
		# check if ship is outside the bounds
		self.terminated = (self.pos[0] < -10.0 or self.pos[0] > 10.0) if True else self.terminated
		# print(self.terminated)
		# check if run time has exceeded maxtime if so terminate run
		# self.terminated = (self.time > runParameters.maxTime) if True else self.terminated
		# print(self.terminated)
# contains the information about the blackhole each blackhole shall have its own object
class blackHole:
	# blackHole constructor
	def __init__(self, position, mass):
		self.pos = position
		self.mass = mass
		#self.mass = 6E+9
	# calculates the force on the ship due to blackHole
	def calculateForce(self, shiplocation, shipmass):
		# returns the x and y components of the gravitation force due to the blackhole on the ship
		# the shiplocation in the euclidean plane and the shipmass are required as inputs
		# the function will return the x and y components of the gravitational force on the ship [fx,fy]

		# find the vector between the ship and the blackhole pointing toward the blackhole
		rvector = np.subtract(self.pos[0], shiplocation)
		# print(rvector)
		rmagnitude = np.linalg.norm(rvector)

		# determine force of gravity along each axis
		# gravity = 6.67408E-11
		gravity = 1
		numerator = np.multiply(self.mass, shipmass)
		numerator = np.multiply(numerator, gravity)
		denominator = np.power(rmagnitude, 3)
		gravityvector = np.multiply(rvector,  np.divide(numerator, denominator))
		return gravityvector

##	Function declarations	##
def loadBlackHoles( filename ):
	# loads the information about the locations and masses of the blackholes from the specified .mat file
	# the file shall contain 3 arrays 'hX','hY','hM' where hX is the x position, hY is the y position and
	# hM is the mass of the blackhole

	# declare blackhole array
	blackholes = []

	# open file
	bhinfo = sio.loadmat(filename)
	hX = bhinfo['hX']
	hY = bhinfo['hY']
	hM = bhinfo['hM']

	# load file info
	for i in range(0,hX.size):
		blackholes.append(blackHole([hX[i], hY[i]], hM[i]));

	# clean-up/return
	return blackholes
def generateInitialConfiguration( runParameters ):
	# new instance of kesselRun to store starting configuration
	newKRun = kesselRun()

	# generate starting position
	x = np.random.uniform(runParameters.minX, runParameters.maxX)
	y = runParameters.startY

	# generate starting velocity
	vmag = np.random.uniform(runParameters.minVel, runParameters.maxVel)
	vang = np.random.normal(runParameters.avgVelAngle, runParameters.stdVelAngle)

	# add starting info to the new kessel run
	newKRun.updatePosition([x,y])
	newKRun.updateVelocity([vmag*np.cos(vang), vmag*np.sin(vang)])

	# set ship mass
	newKRun.setMass(runParameters.shipMass)

	# clean-up/return
	return newKRun
def simulateKesselRun( blackholes, runParameters ):
	# Definition: Will step through the kessel run specified after generatiing an intial configuration
	# Inputs: blackholes: an array of blackHole() objects containing the positions and masses of the
	# 	blackholes on the field. runParameters: a runParameters() tuple containing information about how
	#	the run should be set up.
	# Outputs: returns the kesselrun object created, which will contain the path and the information
	# 	about whether the run was a success of not.

	# create a random initial configuration
	currKesselRun = generateInitialConfiguration(runParameters)

	# while max time has not elapsed continue simulation
	while(currKesselRun.terminated == False and currKesselRun.success == False):
		currKesselRun.stepRungeKuttaIntegration(blackholes, runParameters)

	# return kessel run object
	return currKesselRun

## 	Main code				##
runParameters = runParameters(minX = -5, maxX = 5, startY = -10, finalY = 10, avgVelAngle = math.pi/2.0,
	stdVelAngle = math.pi/4.0, minVel = 2, maxVel = 5, maxForce = 5, resolution = .005, maxTime = 5, shipMass = 1)
blackholes = loadBlackHoles("cluster1.mat")
currKRun = simulateKesselRun(blackholes, runParameters)
print()
print()
print(currKRun.terminated)
print(currKRun.success)

x = []
y = []
outfile = open("output.txt", "w")

for p in currKRun.path:
	x.append(p[0])
	y.append(p[1])
	outfile.write(str(p[0]))
	outfile.write("\t")
	outfile.write(str(p[1]))
	outfile.write("\n")


# print(x)
plt.plot(x, y, color = 'blue')

for b in blackholes:
	plt.plot(b.pos[0], b.pos[1], color = 'red', marker = "H", ms = 5)
plt.savefig("output.png", dpi=96);
# print(currKRun.path)
