import numpy as np

n = 1000
m = 150
l = 3
mat = np.zeros((m,n))
fact = [1, 5, 25]
for i in range(l):
	f = fact[i]
	for j in range(m//3):
		col = np.zeros(m)
		col[3 * f * (j // f) + j%f] = 1
		col[3 * f * (j // f) + j%f + f] = 1
		col[3 * f * (j // f) + j%f + 2*f] = 1
		mat[:,i*(m//3)+j] = col

fact1 = [3,7,15,16,17,18,19,20,21,22,23,24,26,27,28,29,41]
fact2 = [50-x for x in fact1]
parts = [3*x for x in fact1]

for i in range(n//(m//3)-l):
	part = parts[i]
	base = part // 3

	f1 = fact1[i]
	for j in range(base):
		col = np.zeros(m)
		col[3 * f1 * (j // f1) + j%f1] = 1
		col[3 * f1 * (j // f1) + j%f1 + f1] = 1
		col[3 * f1 * (j // f1) + j%f1 + 2*f1] = 1
		mat[:,(i+l)*(m//3)+j] = col

	f2 = fact2[i]
	for j in range(m//3 - base):
		col = np.zeros(m)
		col[part + 3 * f2 * (j // f2) + j%f2] = 1
		col[part + 3 * f2 * (j // f2) + j%f2 + f2] = 1
		col[part + 3 * f2 * (j // f2) + j%f2 + 2*f2] = 1
		mat[:,(i+l)*(m//3)+j+base] = col

np.set_printoptions(threshold=1000000)

print(mat)

# print(np.sum(mat,axis=0))
# print(np.sum(mat,axis=1))

# for i in range(1000):
# 	for j in range(i+1,1000):
# 		d = np.dot(mat[:,i],mat[:,j])
# 		if d > 1:
# 			print(d)