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

	if rmodN == 0 :
		return factors, fSmooth, binaryRow

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
	# This functions create a binary matrix with  len(factorbase)+4  rows.

	binaryMatrix = np.array([])
	parameters = []
	numbers = []
	remainders = []
	remainderDecompositions = []
	
	numberOfBinaryRowsNeeded = len(factorbase) + 4

	parameterPair = (0,0)

	trialsSinceSmoothRemainderFound = 0
	while binaryMatrix.shape[0] < numberOfBinaryRowsNeeded : # Until enough rows are included in the binary matrix
		
		trialsSinceSmoothRemainderFound += 1
		if trialsSinceSmoothRemainderFound > 10000 :
			print("Timeout Error : more than 10000 trials since the last smooth remainder was found, consider increasing the size of the factorbase")
			break

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

		trialsSinceSmoothRemainderFound = 0

		if len(binaryMatrix) == 0:
			binaryMatrix = np.array([binaryRow])
		else :
			binaryMatrix = np.vstack((binaryMatrix,binaryRow))

		parameters.append(parameterPair)
		numbers.append(candidate)
		remainders.append(remainder)
		remainderDecompositions.append(factors)

		
	return binaryMatrix, parameters, numbers, remainders, remainderDecompositions

def solveEquation(binaryMatrix):
	# This function return several different solutions to the equation M^T * x = 0 where sz(M) = (F+2,F)
	
	linearCombinations = np.eye(binaryMatrix.shape[0])
	indexLastFrozenRow = -1
	indexCurrentPrime = 0

	while indexCurrentPrime < binaryMatrix.shape[1] :
		# Find a row that is not frozen and has 1 as value for the current prime and place it under the last frozen row and freeze it
		# (if no row coresponds then skip the next step and increase the index of the current prime)
		for rowIndex in range(indexLastFrozenRow+1,binaryMatrix.shape[0]) :
			if binaryMatrix[rowIndex][indexCurrentPrime] == 1 :
				break
		if binaryMatrix[rowIndex][indexCurrentPrime] != 1 :
			indexCurrentPrime += 1
			continue

		tmp = binaryMatrix[indexLastFrozenRow+1].copy()
		binaryMatrix[indexLastFrozenRow+1] = binaryMatrix[rowIndex].copy()
		binaryMatrix[rowIndex] = tmp.copy()

		tmp = linearCombinations[indexLastFrozenRow+1].copy()
		linearCombinations[indexLastFrozenRow+1] = linearCombinations[rowIndex].copy()
		linearCombinations[rowIndex] = tmp.copy()

		indexLastFrozenRow += 1

		# For all the the none frozen rows, if they have 1 as value for the current prime then add the last frozen row to it
		for rowIndex in range(indexLastFrozenRow+1,binaryMatrix.shape[0]) :
			if binaryMatrix[rowIndex][indexCurrentPrime] == 1 :
				binaryMatrix[rowIndex] += binaryMatrix[indexLastFrozenRow].copy()
				binaryMatrix[rowIndex] %= 2
				linearCombinations[rowIndex] += linearCombinations[indexLastFrozenRow].copy()

	solutions = linearCombinations[indexLastFrozenRow+1:binaryMatrix.shape[0]]%2
	return solutions

def extractResults(solutions, r, remainders, remainderDecompositions) :
	xList = []
	xDecompositionList = []
	xSquaredList = []
	ySquaredList = []
	ySquaredDecompositionList = []
	yList = []
	for solution in solutions :
		x = 1
		xDecomposition = []
		xSquared = 1
		ySquared = 1
		ySquaredDecomposition = []
		y = 1
		for index in range(solution.shape[0]) :
			if solution[index] != 1 :
				continue
			x *= r[index]
			xDecomposition.append(r[index])
			xSquared *= r[index]**2
			ySquared *= remainders[index]
			ySquaredDecomposition += remainderDecompositions[index].copy()
		ySquaredDecomposition.sort()
		for factorIndex in range(0,len(ySquaredDecomposition),2) :
			y *= ySquaredDecomposition[factorIndex]
		xList.append(x)
		xDecompositionList.append(xDecomposition)
		xSquaredList.append(xSquared)
		ySquaredList.append(ySquared)
		ySquaredDecompositionList.append(ySquaredDecomposition)
		yList.append(y)

	return xList, xDecompositionList, xSquaredList, ySquaredList, ySquaredDecompositionList, yList

def gcd(a,b) :
	if a % b == 0 :
		return b
	return gcd(b, a%b)

def factorize(numberToFactorize, cardinalOfFactorbase, showSteps = False):

	print("Generating factorbase...")
	factorbase,factorbaseBoudary = getPrimes(cardinalOfFactorbase)
	if showSteps :
		print("Factorbase : " + str(factorbase))
	
	print("Generating binary matrix...")
	binaryMatrix, parameters, r, remainders, remainderDecompositions = createBinaryMatrix(numberToFactorize,factorbase)
	if showSteps :
		print("Data coresponding to the binary matrix :")
		data = zip(parameters, r, remainders, remainderDecompositions)
		headers = ["(k,j)", "r", "r**2 mod N", "decomposition"]
		print(tabulate(data, headers = headers, tablefmt="fancy_grid"))
		print("Binary matrix :")
		print(str(binaryMatrix))
	
	print("Solving  Mx=0...")
	solutions = solveEquation(binaryMatrix)
	xList, xDecompositionList, xSquaredList, ySquaredList, ySquaredDecompositionList, yList = extractResults(solutions, r, remainders, remainderDecompositions)
	if showSteps :
		print("Data coresponding to the solutions :")
		data = zip(xList, xDecompositionList, xSquaredList, ySquaredList, ySquaredDecompositionList, yList)
		headers = ["x", "decomposition of x", "x**2", "y**2", "decomposition of y**2", "y"]
		print(tabulate(data, headers = headers, tablefmt="fancy_grid"))
	
	print("Extracting solution...")
	isSolutionFound = False
	for solutionIndex in range(solutions.shape[0]) :
		if xList[solutionIndex] - yList[solutionIndex] == 0 :
			continue
		p = gcd(numberToFactorize, xList[solutionIndex] - yList[solutionIndex])
		if p == 1 or p == numberToFactorize :
			continue
		q = numberToFactorize // p
		isSolutionFound = True
		break
	if not isSolutionFound :
		print("Error : None of the pairs of squares resulted in a factorization")
		return (numberToFactorize, 1)
	if showSteps :
		print("Factorisation :")
		print(str(numberToFactorize) + " = " + str(p) + " x " + str(q))

	return (p,q)

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Script that finds 2 prime factors of N")
	parser.add_argument("--N", required=True, type=int)
	parser.add_argument("--F", required=True, type=int)
	parser.add_argument("--showSteps", required=False, type=bool)
	args = parser.parse_args()
	N = args.N
	F = args.F
	showSteps = args.showSteps

	print("Running the factorization program on N = " + str(N))

	print(factorize(N, F, showSteps))