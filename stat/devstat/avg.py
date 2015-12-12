import sys
f = open(str(sys.argv[1]),'r+b')
a = [0.0,0.0,0.0]
for line in f:
	b = line.split(' ')
	a[0] += float(b[2])
	a[1] += float(b[4])
	a[2] += float(b[6])
#	print a[0],a[1],a[2]
#	print "\n"
	
a[0] = a[0] /  20
a[1] = a[1] /  20
a[2] = a[2] /  20

c = str(a[0]) + ',' + str(a[1]) + ',' + str(a[2])
print c
