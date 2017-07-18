from utilities.mat import *
output=open("test_files\OUTPUT.txt","r+")
input=MAT("test_files\CLUSTERS.mat")
input_file=open("test_files\CLUSTERS.mat","r+")

output.write("#1.2\n")
output.write(str(len(input.matrix))+"\t"+str(len(input.matrix['1']))+"\n")
output.write("NAME\tDescription\t")
top_line=input_file.readline()
tok = top_line.strip().split('\t')
for i in range(0,len(tok)):
    if i!=0:
        output.write(tok[i]+"\t")
output.write("\n")

for i in input.matrix.keys():
    if len(input.matrix[i])!=0:
        output.write(str(i) + "\t")
        for j in input.matrix[i]:
            output.write(str(j) + "\t")
        output.write("\n")