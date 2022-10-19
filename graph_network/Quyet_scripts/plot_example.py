import matplotlib.pyplot as plt
import numpy as np

def log(x: iter)->iter:
    out = [math.log10(i) for i in x]
    return out

t1 = np.arange(1.0, 5.0, 0.1)
t2 = np.arange(1.0, 5.0, 0.02)

plt.plot(log(t1),log(t1),'r.',label="train")
plt.plot(log(t2),log(t2),'b.',label="test")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title('Histogram of IQ')
plt.show()