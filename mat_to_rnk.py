from utilities.mat import *
output=open("test_files\OUTPUT.txt","r+")
input=MAT("test_files\CLUSTERS.mat")
input_file=open("test_files\CLUSTERS.mat","r+")

output.write("\n")


for i in input.dict.keys():

    output.write(str(i) + "\t" + str(input.dict[i][0]))
    output.write("\n")