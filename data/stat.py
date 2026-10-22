import pickle
import numpy as np

data=pickle.load(open('NC_000018.10_label.pkl', "rb"))
std=np.std(data)
mean=np.mean(data)
m=np.max(data)
print(std,mean,m)
