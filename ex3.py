import sys
import argparse
import numpy as np
import pdb



def getPrimes(F):  #Function to obtain the elements in F, get the first F prime numbers
	primes = []
	pr = 1  #pr: candidate prime number, starting with 1
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
	binaryRow = np.zeros(len(primes))  #we will also return the binary row with 1s when a factor is present
	curr = rmodN #number we want to factorize
	i = 1
	pr = primes[i]
	while i < len(primes):
		pr = primes[i] #get prime
		if curr % pr == 0: #check if factor
			factors.append(pr) #if so add as factor and remove from number, pr is not changed as a number can be divided more than one time by same prime
			curr = curr/pr
			binaryRow[i] = 1 #factor is present -> 1 in binary row
			if curr == 1:
				fSmooth = True #if we got to 1 it means we succesfully factorized the number with the primes inside F which means the number is F-smooth
				break
		else:
			i += 1 #if not a factor we continue on prime list
	return factors, fSmooth, binaryRow


def createBinaryMatrix(N,F,primes):
	M = np.array([]) #Empty binary matrix
	L = F + 2

	k = 0
	j = 0

	while M.shape[0] < L: #goes through various j and k
		k += 1
		while M.shape[0] < L:
			j += 1
			r = int(np.floor(np.sqrt(k*N) + j)) #calculate r
			res = r**2 % N #calculate r^2 mod N
			factors,fSmooth,binaryRow = findFactors(res,primes) #get factors and binary row to add to M
			if not fSmooth: #only add binary row if it is f smooth
				continue
			if len(M) == 0: #first row
				M = binaryRow
			else:
				if not ([binaryRow] in M.tolist()): ## NOT WORKING, KEEPS ADDING DUPLICATED ROWS (trying with N=323)   #add row if not already there
					M = np.vstack((M,binaryRow))
	return M

def factorize(N, F):
	primes,B = getPrimes(F)
	M = createBinaryMatrix(N,F,primes)   ## ENDED TRYING TO MAKE createBinaryMatrix work
	return M,0
	#return p,q

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script that finds 2 prime factors of N")
    parser.add_argument("--N", required=True, type=int)
    args = parser.parse_args()

    N = args.N

    p,q = factorize(N, 10)

    print('N = ',p,' x ',q)

    # if __name__ == "__main__":
#   N = sys.argv[1]
#   p,q = factorize(N)
#   print('N = ',p,' x ',q)