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
import matplotlib.pyplot as plt
import math as math
from collections import namedtuple


## 	Define Classes			##
# used for setting the parameters for the starting configuration in the kessel run
runParameters = namedtuple("runParameters",
	"minX, maxX, startY, finalY, avgVelAngle, stdVelAngle, minVel, maxVel maxForce resolution maxTime shipMass" )
# contains the information about the ship and its path
class kesselRun:
	def __init__(self):
		self.pos = [-1,-1]
		self.path = []
		self.vel = [-1,-1]
		self.acc = [-1,-1]
		self.time = 0
		self.success = False
		self.terminated = False
		self.mass = -1

	def setMass(self, mass):
		self.mass = mass

	def addPosition(self, position):
		#self.path.append(position)
		np.append(self.path, position)
		self.pos = np.array(position)

	def updateVelocity(self, velocity):
		self.vel = np.array(velocity)

	def updateAcceleration(self, acceleration):
		self.acc = np.array(acceleration)

	def incTime(self, deltat):
		self.time += deltat
# contains the information about the blackhole each blackhole shall have its own object
class blackHole:
	def __init__(self, position, mass):
		self.pos = np.array(position)
		self.mass = mass

	def calculateForce(self, shiplocation, shipmass):
		# returns the x and y components of the gravitation force due to the blackhole on the ship
		# the shiplocation in the euclidean plane and the shipmass are required as inputs
		# the function will return the x and y components of the gravitational force on the ship [fx,fy]

		# find the vector between the ship and the blackhole pointing toward the blackhole
		rvector = np.subtract(self.pos, shiplocation)
		rmagnitude = np.linalg.norm(rvector)

		# determine force of gravity along each axis
		gravity = 6.67408E-11
		gravityvector = np.multiply(rvector, self.mass * shipmass * gravity / rmagnitude^3)
		return gravityvector


##	Function declarations	##
def loadBlackHoles( filename ):
	# loads the information about the locations and masses of the blackholes from the specified text file
	# the data shall be in the configuration x\ty\tmass\n

	# open file
	bhfile = open(filename, "r")

	# load file info
	for line in bhfile:
		print(line)

	# clean-up/return
	return
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
	newKRun.addPosition([x,y])
	newKRun.updateVelocity([vmag, vang])

	# set ship mass
	newKRun.setMass(runParameters.shipMass)

	# clean-up/return
	return newKRun

# def stepRungeKuttaIntegration( currKesselRun, blackholes, runParameters):
# 	# Definition
#
# 	# determine force upon the ship
#
#
# 	# check if force has exceeded specified max if so terminate run
# 	currKesselRun.terminated = forceonship[0] > runParameters.maxForce ? True : False
#
# 	# update the ship acceleration
#
# 	# update the ship velocity
#
# 	# update the ship position
#
# 	# check if ship's y is greater than finalY if so then sucessful run
# 	currKesselRun.sucess = currKesselRun.pos[1] > runParameters.finalY ? True, False
#
# 	# increment ship time by resolution
# 	currKesselRun.incTime(runParameters.resolution)
#
# 	# check if run time has exceeded maxtime if so terminate run
# 	currKesselRun.terminated = currKesselRun.time > runParameters.maxtime ? True, False
#
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
	while(!currKesselRun.terminated && !currKesselRun.sucess):
		stepRungeKuttaIntegration(currKesselRun, blackholes, runParameters)

	# return kessel run object
	return currKesselRun

## 	Main code				##
runParameters = runParameters(minX = -5, maxX = 5, startY = -10, finalY = 10, avgVelAngle = math.pi/2.0,
	stdVelAngle = math.pi/4.0, minVel = 2, maxVel = 5, maxForce = 4, resolution = .005, maxTime = 1000, shipMass = 5)
currKRun = generateInitialConfiguration(runParameters)
blk = blackHole([-4,1], 10)
loadBlackHoles("testblackholes.txt")
