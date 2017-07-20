from utilities.mat import *
output=open("test_files\OUTPUT.txt","r+")
input=MAT("test_files\CLUSTERS.mat")
input_file=open("test_files\CLUSTERS.mat","r+")

output.write("#1.2\n")
output.write(str(len(input.dict)) + "\t" + str(len(input.dict['1'])) + "\n")
output.write("NAME\tDescription\t")
top_line=input_file.readline()
tok = top_line.strip().split('\t')
for i in range(0,len(tok)):
    if i!=0:
        output.write(tok[i]+"\t")
output.write("\n")

for i in input.dict.keys():
    if len(input.dict[i])!=0:
        output.write(str(i) + "\t")
        for j in input.dict[i]:
            output.write(str(j) + "\t")
        output.write("\n")