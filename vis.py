
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, pts, xs, ys):
        self.pts = pts
        self.xs,self.ys = self.unzip()
        fig = plt.figure(figsize=(7, 7))
        ax = plt.axes(xlim=(-30, 30), ylim=(-30, 30))
        self.scatter = ax.scatter(self.xs,self.ys)
        anim = FuncAnimation(fig, self.update, interval=.05)
        plt.show()

    def update(self, i):
        self.scatter.set_offsets(self.pts[i])
        return self.scatter,

    def unzip(self):
        xs, ys = [], []
        for sample in self.pts:
            x, y = zip(*sample)
            xs.append(x)
            ys.append(y)
        return xs,ys