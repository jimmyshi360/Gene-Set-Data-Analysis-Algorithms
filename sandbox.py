
rankings= [0.0001, 0.00034, 0.0035,0.004, 0.024, 0.034,0.042]
prev_bh_value=0
bh_arr=[]
for i in range(0, len(rankings)):
    bh_value=rankings[i] * len(rankings) / float(i+1)
    bh_value=min(bh_value,1)

    #to preserve monotonicity
    if prev_bh_value!=0:
        prev_bh_value=min(bh_value,prev_bh_value)
        bh_arr.append(prev_bh_value)
    if i == len(rankings)-1:
        bh_arr.append(bh_value)
    prev_bh_value=bh_value
print(bh_arr)