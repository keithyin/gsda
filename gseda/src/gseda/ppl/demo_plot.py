
print("start")
import matplotlib
print("end")

# print("backend before:", matplotlib.get_backend())

matplotlib.use("Agg")

import matplotlib.pyplot as plt

print("backend after:", plt.get_backend())

print("before")

fig, axes = plt.subplots(
    1,
    2,
    figsize=(14, 5),
    squeeze=False,
)

print("after")

fig.savefig("test.png")

print("saved")