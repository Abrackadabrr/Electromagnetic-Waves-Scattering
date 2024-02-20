import csv


file1 = open('/media/evgen/SecondLinuxDisk/4_level/INM RAS/Experimental/meshes/In_DM/25_50(2050).rwg', 'r')
count = 0


numNodes = int(next(file1))

csvfile = open(str(numNodes) + '_nodes.csv', 'w', newline='\n')
wNodes = csv.writer(csvfile)
wNodes.writerow(["x", 'y', 'z'])

print("Nodes = ", numNodes)
for i in range(numNodes):
    count += 1
    wNodes.writerow(next(file1).split())

numCells = int(next(file1))
print("Cells = ", numCells)
csvfile = open(str(numCells) + '_cells.csv', 'w', newline='\n')
wCells = csv.writer(csvfile)
wCells.writerow(['f', 's', 't', 'fo'])
for i in range(numCells):
    count += 1
    wCells.writerow(next(file1).split())


# Closing files
file1.close()
# -0.082250 -1.250000 0.651075
