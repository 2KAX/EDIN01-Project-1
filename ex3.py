import argparse
import numpy as np
from tabulate import tabulate



def getPrimes(F):  #Function to obtain the elements in F, get the first F prime numbers
	primes = []
	pr = 2  #pr: candidate prime number, starting with 2
	isPrime = True
	while len(primes) < F: #while we have less than F numbers
		for i in range(2,int(np.floor(np.sqrt(pr)))+1):  #check numbers from 2 to √pr
			if pr % i == 0:
				isPrime = False  #if remaining is 0 pr is not a prime number
				break
		if isPrime:
			primes.append(pr) #adding prime to the list
		else:
			isPrime = True
		pr += 1
	return primes, pr #returns array of prime numbers and the smooth boundary B

def findFactors(rmodN, primes):
	factors = []
	fSmooth = False
	binaryRow = np.zeros(len(primes),dtype=int)  #we will also return the binary row with 1s when a factor is present
	curr = rmodN #number we want to factorize
	i = 0
	pr = primes[i]
	while i < len(primes):
		pr = primes[i] #get prime
		if curr % pr == 0: #check if factor
			factors.append(pr) #if so add as factor and remove from number, pr is not changed as a number can be divided more than one time by same prime
			curr = curr/pr
			binaryRow[i] = (binaryRow[i] + 1) % 2 #factor is present and even-> 1 in binary row
			if curr == 1:
				fSmooth = True #if we got to 1 it means we succesfully factorized the number with the primes inside F which means the number is F-smooth
				break
		else:
			i += 1 #if not a factor we continue on prime list
	return factors, fSmooth, binaryRow

def getNextPair(pair) :
	# This function goes through all pairs in NxN : (0,0) -> (1,0) -> (0,1) -> (2,0) -> (1,1) -> (0,2) -> (3,0) -> ...
	(k,j) = pair
	if k == 0 :
		return (k+j+1,0)
	else :
		return (k-1,j+1)

def createBinaryMatrix(numberToFactorize,factorbase) :
	# This functions create a binary matrix with  len(factorbase)+2  rows.

	binaryMatrix = np.array([])
	parameters = []
	numbers = []
	remainders = []
	remainderDecompositions = []
	
	numberOfBinaryRowsNeeded = len(factorbase) + 2

	parameterPair = (0,0)

	while binaryMatrix.shape[0] < numberOfBinaryRowsNeeded : # Until enough rows are included in the binary matrix

		# Test if the next number is smooth
		parameterPair = getNextPair(parameterPair)
		candidate = int(np.floor(np.sqrt(parameterPair[0]*numberToFactorize)) + parameterPair[1])
		remainder = candidate**2 % numberToFactorize
		factors,isRemainderSmooth,binaryRow = findFactors(remainder,factorbase)
		if not isRemainderSmooth :
			continue

		# Test if the binary row is already in the binary matrix
		if len(binaryMatrix) > 0 and (np.any(np.all(binaryMatrix == binaryRow, axis=1))):  #add row if not already there 
			continue

		if len(binaryMatrix) == 0:
			binaryMatrix = np.array([binaryRow])
		else :
			binaryMatrix = np.vstack((binaryMatrix,binaryRow))

		parameters.append(parameterPair)
		numbers.append(candidate)
		remainders.append(remainder)
		remainderDecompositions.append(factors)

		
	return binaryMatrix, parameters, numbers, remainders, remainderDecompositions

def solveEquation(M):
	# This function return several different solutions to the equation M^T * x = 0 where sz(M) = (F+2,F)
	# null_space = np.linalg.lstsq(M, np.zeros(M.shape[0]))[0]
	#print(np.transpose(M))
	# print(str(np.linalg.lstsq(np.transpose(M), np.zeros(M.shape[1]))))
	pass

def factorize(numberToFactorize, cardinalOfFactorbase, showSteps = False):

	factorbase,factorbaseBoudary = getPrimes(cardinalOfFactorbase)
	if showSteps :
		print("Factorbase : " + str(factorbase))

	binaryMatrix, parameters, numbers, remainders, remainderDecompositions = createBinaryMatrix(numberToFactorize,factorbase)
	if showSteps :
		print("Data coresponding to the binary matrix :")
		data = zip(parameters, numbers, remainders, remainderDecompositions)
		headers = ["(k,j)", "r", "r**2 mod N", "decomposition"]
		print(tabulate(data, headers = headers, tablefmt="fancy_grid"))
		print("Binary matrix :")
		print(str(binaryMatrix))

	# TBC solve equation x*binaryMatrix = 0
	pass

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Script that finds 2 prime factors of N")
	parser.add_argument("--N", required=True, type=int)
	args = parser.parse_args()
	N = args.N

	N = 323 # tmp test
	F = 5

	print("Running the program on N = " + str(N))

	factorize(N, F, showSteps = True)

	#p,q = factorize(N, 10)
	#print('N = ',p,' x ',q)