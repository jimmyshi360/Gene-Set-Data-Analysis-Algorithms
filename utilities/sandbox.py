import multiprocessing
import time
import math

def square_numbers(numbers):
    square_result = []
    for n in (0,numbers):
        square_result.append(math.factorial(int (n)))
    return 2

def n_test():
    square_result = []
    for i in range(0,5500):
            square_result.append(math.factorial(int (i)))
    return 2


if __name__ == "__main__":
    numbers = list(range(5500))

    t1 = time.time()

    p1 = multiprocessing.Pool(processes=8)
    r=p1.map(square_numbers,numbers)
    p1.close()

    t2 = time.time()


    print t2-t1