from fpylll import IntegerMatrix, SVP
import sys

def svp(B):
	A = IntegerMatrix.from_matrix(B)
	return SVP.shortest_vector(A)

def first_primes(n):
	p = 1
	P = []
	while len(P) < n:
		p = next_prime(p)
		P += [p]
	return P

def is_smooth(x, P):
	if x==0:
		return False
	y = x
	for p in P:
		while p.divides(y):
			y /= p
	return abs(y) == 1

# Test if a factoring relation was indeed found.
def test_Schnorr(N, n, c, prec=6):
	P = first_primes(n)
	#f = list(range(1, n+1))
	#shuffle(f)
	#print(P)

	# Scale up and round
	def sr(x):
		return round(x * 2^prec)

#	diag = [sr(N*f[i]) for i in range(n)] + [sr(N*ln(N))]       #  this is the Basis from https://eprint.iacr.org/2021/232.pdf, doesn't seem to work
#	diag = [sr( sqrt( ln(P[i]) ) ) for i in range(n)] + [sr(N^c*ln(N))] #  this is the Basis from Schnorr 1991 + taking sqrt
	diag = [sr( ln(P[i]) ) for i in range(n)] + [sr(N^c*ln(N))] #  this is the Basis from Schnorr 1991, at least it works for small n

	B = diagonal_matrix(diag, sparse=False)
	for i in range(n):
		B[i, n] = sr(N^c*ln(P[i]))

	b = svp(B)
	e = [b[i] / sr( ln(P[i]) ) for i in range(n)]
#	e = [b[i] / sr( sqrt( ln(P[i]) ) ) for i in range(n)]

	u = 1
	v = 1
	for i in range(n):
		assert e[i] in ZZ
		if e[i] > 0:
			u *= P[i]^e[i]
		if e[i] < 0:
			v *= P[i]^(-e[i])
	print("u, u - v*N", u, u - v*N)
	success = is_smooth(u - v*N, P) 
	if success:
		print("diag of B:", diag)
		print("e = ", e)

	return success

try:
	bits = int(sys.argv[1])
except:
	bits = 50

try:
	n = int(sys.argv[2])
except:
	n = 50

try:
	trials = int(sys.argv[3])
except:
	trials = 10

print("Testing Schnorr's relation finding algorithm with n=%d on RSA-moduli of %d bits, %d trials"%(n, bits, trials))

successes = 0
for i in range(trials):
	p = random_prime(2^(bits/2), false, 2^(bits/2-1))
	q = random_prime(2^(bits/2), false, 2^(bits/2-1))
	N = p*q
	print("N, p, q", N, p, q)
	success = test_Schnorr(N, n, 1.1)
	successes += success
	print(success)
	sys.stdout.flush()

print("\n%d Factoring Relation found out of %d trials"%(successes, trials))
