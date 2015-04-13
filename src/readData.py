import numpy

f = open("A.txt", "r")

uN = numpy.array([int(x) for x in f.readline().split()])
u = numpy.array([int(x) for x in f.readline().split()]).reshape(tuple(uN))

vN = numpy.array([int(x) for x in f.readline().split()])
v = numpy.array([int(x) for x in f.readline().split()]).reshape(tuple(vN))

wN = numpy.array([int(x) for x in f.readline().split()])
w = numpy.array([int(x) for x in f.readline().split()]).reshape(tuple(wN))

pN = numpy.array([int(x) for x in f.readline().split()])
p = numpy.array([int(x) for x in f.readline().split()]).reshape(tuple(pN))

f.close()

print(uN)
print(vN)
print(wN)
print(pN)


print('\nu:')
for i in range(uN[2]):
    print(u[:, :, i].T)


print('\nv:')
for i in range(vN[2]):
    print(v[:, :, i].T)


print('\nw:')
for i in range(wN[2]):
    print(w[:, :, i].T)


print('\np:')
for i in range(pN[2]):
    print(p[:, :, i].T)
