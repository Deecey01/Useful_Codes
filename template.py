import sys
input = sys.stdin.readline
from math import ceil,gcd,sqrt,floor,log,log2,comb
from io import BytesIO, IOBase
import os
import sys
from collections import defaultdict
import bisect
import heapq

def SieveOfEratosthenes(n):
  prime = [True for i in range(n+1)]
  p = 2
  while (p * p <= n):
      if (prime[p] == True):
          # Update all multiples of p
          for i in range(p * p, n+1, p):
              prime[i] = False
      p += 1
  for p in range(2, n+1):
      if prime[p]:
          print(p)
                       
'''a+b=(a^b)+2(a&b)'''
                       
def printDivisors(n,arr) : 
    i = 1
    while i <= sqrt(n): 
        if (n % i == 0) : 
            if (n / i == i) : 
                arr.append(i)
            else : 
                arr.extend([i,int(n/i)]) 
        i = i + 1
          
def primeFactors(n,arr):
      while n % 2 == 0:
        arr.append(2)
        n= n // 2
      for i in range(3,int(sqrt(n))+1,2):
          while n % i== 0:
            arr.append(i)
            n = n // i
      if n > 2:
          arr.append(n)
                       
              
def power(base, exponent):
    if exponent == 0:
        return 1
    elif exponent == 1:
        return base
    elif exponent % 2 == 0:
        return power(base*base, exponent//2)
    else:
        return base * power(base*base, (exponent-1)//2)
              
def maxSubArraySum(a):
    size=len(a)
    max_so_far = -10**18
    max_ending_here = 0

    for i in range(0, size):
        max_ending_here = max_ending_here + a[i]
        if (max_so_far < max_ending_here):
            max_so_far = max_ending_here

        if max_ending_here < 0:
            max_ending_here = 0
    return max_so_far

def longestIncreasingSubsequence(arr):
    from bisect import bisect_left
    
    tmp = [arr[0]]
    for i in range(1, len(arr)):
        if arr[i] > tmp[-1]:
            tmp.append(arr[i])
        else:
            ind = bisect_left(tmp, arr[i])
            tmp[ind] = arr[i]
    return len(tmp)

def overlap(v):
    ans = 0
    count = 0
    data = []
    for i in range(len(v)):
        data.append([v[i][0], 'x'])
        data.append([v[i][1], 'y'])

    data = sorted(data)

    for i in range(len(data)):
        if (data[i][1] == 'x'):
            count += 1
        if (data[i][1] == 'y'):
            count -= 1

        ans = max(ans, count)
    return (ans)

def ncr(n, r, p):
    # initialize numerator and denominator
    num = den = 1
    for i in range(r):
        num = (num * (n - i)) % p
        den = (den * (i + 1)) % p
    return (num * pow(den, p - 2, p)) % p


#program for KMP Algorithm
def KMPSearch(pat, txt):
	M = len(pat)
	N = len(txt)
	lps = [0]*M
	j = 0 
	computeLPSArray(pat, M, lps)

	i = 0 # index for txt[]
	while (N - i) >= (M - j):
		if pat[j] == txt[i]:
			i += 1
			j += 1

		if j == M:
			print("Found pattern at index " + str(i-j))
			j = lps[j-1]

		elif i < N and pat[j] != txt[i]:
			if j != 0:
				j = lps[j-1]
			else:
				i += 1

# Function to compute LPS array
def computeLPSArray(pat, M, lps):
	len = 0 
	lps[0] = 0 
	i = 1
	while i < M:
		if pat[i] == pat[len]:
			len += 1
			lps[i] = len
			i += 1
		else:
			if len != 0:
				len = lps[len-1]
			else:
				lps[i] = 0
				i += 1
#Disjoint Set Union
class DSU:
    def __init__(self, n):
        # e[i] < 0 means i is a root and -e[i] is the size of the set
        self.e = [-1] * n

    def get(self, x):
        # Path compression
        if self.e[x] < 0:
            return x
        self.e[x] = self.get(self.e[x])
        return self.e[x]

    def same_set(self, x, y):
        return self.get(x) == self.get(y)

    def size(self, x):
        return -self.e[self.get(x)]

    def unite(self, x, y):
        x, y = self.get(x), self.get(y)
        if x == y:
            return False
        # Union by size
        if self.e[x] > self.e[y]:
            x, y = y, x
        self.e[x] += self.e[y]
        self.e[y] = x
        return True

