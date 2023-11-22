from ex3 import *


print("Running tests")

print("===== GETPRIMES =====")
print()

(factorbase,boundary) = getPrimes(100)
print("factorbase : " + str(factorbase))
print("boundary : " + str(boundary))
print()

print("===== FINDFACTORS =====")
print()

(factorbase,boundary) = getPrimes(4)
print("factorbase : " + str(factorbase))
print()

print("Number : 60")
factors, fSmooth, binaryRow = findFactors(60,factorbase)
print("factors : " + str(factors))
print("fSmooth : " + str(fSmooth))
print("binaryRow : " + str(binaryRow))
print()

print("Number : 17")
factors, fSmooth, binaryRow = findFactors(17,factorbase)
print("factors : " + str(factors))
print("fSmooth : " + str(fSmooth))
print("binaryRow : " + str(binaryRow))
print()

print("Number : 6")
factors, fSmooth, binaryRow = findFactors(6,factorbase)
print("factors : " + str(factors))
print("fSmooth : " + str(fSmooth))
print("binaryRow : " + str(binaryRow))
print()

print("Number : 1")
factors, fSmooth, binaryRow = findFactors(1,factorbase)
print("factors : " + str(factors))
print("fSmooth : " + str(fSmooth))
print("binaryRow : " + str(binaryRow))
print()

print("Number : 0")
factors, fSmooth, binaryRow = findFactors(0,factorbase)
print("factors : " + str(factors))
print("fSmooth : " + str(fSmooth))
print("binaryRow : " + str(binaryRow))
print()

print("Number : 27")
factors, fSmooth, binaryRow = findFactors(27,factorbase)
print("factors : " + str(factors))
print("fSmooth : " + str(fSmooth))
print("binaryRow : " + str(binaryRow))
print()

print("===== BINARY MATRIX =====")
print()

(factorbase,boundary) = getPrimes(10)
M = createBinaryMatrix(323,factorbase)[0]
print(str(M))
print()