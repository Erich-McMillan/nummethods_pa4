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
		#np.append(self.path, np.array(position))
		self.pos = np.array(position)
		#print(self.path)
		#print(self.pos)
	def updateVelocity(self, velocity):
		self.vel = np.array(velocity)
		#print(self.vel)
	def updateAcceleration(self, acceleration):
		self.acc = np.array(acceleration)
		#print(self.acc)
	def updateTime(self, deltat):
		self.time += deltat

	def stepRungeKuttaIntegration(self, blackholes, runParameters):
		#

		# declare variables
		forceonship = np.array([0,0])

		# determine force upon the ship
		for currentbh in blackholes:
			forceonship = np.add(forceonship, currentbh.calculateForce(self.pos, self.mass))

		# check if force has exceeded specified max if so terminate run
		self.terminated = (forceonship[0] > runParameters.maxForce) if True else False

		# update acceleration
		self.updateAcceleration(np.divide(forceonship, self.mass))

		# update velocity
		self.updateVelocity(np.add(self.vel, np.multiply(self.acc, runParameters.resolution)))

		# update position
		self.updatePosition(np.add(self.pos, np.multiply(self.vel, runParameters.resolution)))

		# update time
		self.updateTime(runParameters.resolution)

		# check if ship's y is greater than finalY if so then sucessful run
		self.success = (self.pos[1] > runParameters.finalY) if True else False

		# check if run time has exceeded maxtime if so terminate run
		self.terminated = (self.time > runParameters.maxTime) if True else False


# contains the information about the blackhole each blackhole shall have its own object
class blackHole:
	# blackHole constructor
	def __init__(self, position, mass):
		self.pos = np.array(position)
		#self.mass = mass
		self.mass = 5.2E+15

	# calculates the force on the ship due to blackHole
	def calculateForce(self, shiplocation, shipmass):
		# returns the x and y components of the gravitation force due to the blackhole on the ship
		# the shiplocation in the euclidean plane and the shipmass are required as inputs
		# the function will return the x and y components of the gravitational force on the ship [fx,fy]

		# find the vector between the ship and the blackhole pointing toward the blackhole
		#print(np.array(map(float, self.pos))
		#print(self.pos)
		#print(shiplocation)
		rvector = np.subtract(self.pos, shiplocation)
		rmagnitude = np.linalg.norm(rvector)

		# determine force of gravity along each axis
		print(self.mass)
		gravity = 6.67408E-11
		numerator = np.multiply(self.mass, shipmass)
		numerator = np.multiply(numerator, gravity)
		denomiator = np.power(rmagnitude, 3)
		gravityvector = np.multiply(rvector,  np.divide(numerator, denomiator))
		print(gravityvector)
		return gravityvector


##	Function declarations	##
def loadBlackHoles( filename ):
	# loads the information about the locations and masses of the blackholes from the specified text file
	# the data shall be in the configuration x\ty\tmass\n

	# declare blackhole array
	blackholes = []

	# open file
	bhfile = open(filename, "r")

	# load file info
	for line in bhfile:
		#values = line.split()
		values = np.array(line.split(), dtype='|S4')
		values = values.astype(np.float)
		blackholes.append(blackHole([values[0],values[1]], values[2]))

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
	newKRun.updateVelocity([vmag, vang])

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
	stdVelAngle = math.pi/4.0, minVel = 2, maxVel = 5, maxForce = 4, resolution = .00005, maxTime = 1000, shipMass = 5)
# currKRun = generateInitialConfiguration(runParameters)
# blk = blackHole([-4,1], 10)
blackholes = loadBlackHoles("testblackholes.txt")
currKRun = simulateKesselRun(blackholes, runParameters)
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

plt.plot(x, y)
plt.savefig("output.png", dpi=96);
#print(currKRun.path)
