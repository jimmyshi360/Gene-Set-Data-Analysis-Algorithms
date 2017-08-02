import matplotlib.pyplot as plt

input=open("C:\Users\Jimmy\Documents\dev\projects\enrichments\\files\\test_files\\fdr_vals.txt", "r+")

fisher = []
binomial = []
hypergeo =[]
chisqd =[]
arr=[]
count=0
for line in input.readlines():

    if line !="" and line != "\n":
        items=line.split(",")
        items=map(float,items)
        arr.append(items)
fisher=arr[0]
binomial=arr[1]
hypergeo=arr[2]
chisqd=arr[3]

data_to_plot = [fisher, binomial, hypergeo, chisqd]


fig = plt.figure(1, figsize=(9, 6))

ax = fig.add_subplot(111)

bp = ax.boxplot(data_to_plot)

fig.savefig('fig1.png', bbox_inches='tight')