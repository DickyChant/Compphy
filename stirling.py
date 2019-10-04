import numpy as np
import pandas as pd

def lnstirling(n):
    return np.log(2*np.pi*n)*0.5+n*np.log(n/np.e)

def Hilbert(n):
    sum = 0
    for i in range(1,n-1):
        sum += 4*lnstirling(i+1)
    for i in range(1,2*n-1):
        sum -= lnstirling(i+1)
    return sum

if __name__=="__main__":
    for i in range(10):
        print(i+1,Hilbert(i+1),np.exp(Hilbert(i+1)))
    table = [[i+1,Hilbert(i+1),np.exp(Hilbert(i+1))] for i in range(10)]
    table_df = pd.DataFrame(table,columns=["n","ln(H_n)","H_n"])
    table_df.to_latex("Hilbert.tex")
   