import random as r
import math
n = 0
m = 0
o = 0

n = input("Enter DNA size: ")
m = input("Enter oligonucleotide size: ")
o = input("Enter error %: ")
DNA = ""
parts = []
listlen = 0

n = int(n)
m = int(m)
o = int(o)

print(n/m)

i = 0
tmp = 0
tmpstr = ""

while(i < n):
    tmp = r.randint(0, 3)
    if(tmp == 0):
        DNA = DNA + "A"
    elif(tmp == 1):
        DNA = DNA + "T"
    elif(tmp == 2):
        DNA = DNA + "C"
    else:
        DNA = DNA + "G"

    i += 1

print(DNA)

i = 0
j = 0

while i < math.ceil(n/m):
    while(j < m):
        if(len(DNA[i*m+j:i*m+m+j]) == m):
            parts.append(DNA[i*m+j:i*m+m+j])
        j += 1
    j = 0
    i += 1

i = 0
j = 0
listlen = len(parts)
f = open('data.txt', 'a', encoding='utf-8')
f.write(DNA)
f.write("\n")
f.write(parts[0])
f.write("\n")
parts.pop(0)
print(parts)
while i < o/100*listlen:
    parts.pop(r.randint(0, len(parts)))
    i += 1
i = 0
while i < o/100*listlen:
    while j < m:
        tmp = r.randint(0, 3)
        if tmp == 0:
            tmpstr = tmpstr + "A"
        elif tmp == 1:
            tmpstr = tmpstr + "T"
        elif tmp == 2:
            tmpstr = tmpstr + "C"
        else:
            tmpstr = tmpstr + "G"
        j += 1
    j = 0
    parts.append(tmpstr)
    tmpstr = ""
    i += 1
parts.sort()
print(parts)
while(i < len(parts)):
    f.write(parts[i])
    i += 1
    if i < len(parts):
        f.write("\n")
f.close()
