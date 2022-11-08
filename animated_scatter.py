import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, pts, save=False, filename="rope_vid.mp4"):
        self.pts = pts
        self.xs, self.ys = self.unzip()
        fig = plt.figure(figsize=(7, 7))
        ax = plt.axes(xlim=(-30, 30), ylim=(-30, 30))
        self.scatter = ax.scatter(self.xs, self.ys)
        anim = FuncAnimation(fig, self.update, frames=self.pts,interval=.05)
        if save:
            writervideo = animation.FFMpegWriter(fps=60)
            anim.save(filename, writer=writervideo)
        plt.show()
        plt.close()

    def update(self, frame):
        self.scatter.set_offsets(frame)
        return self.scatter,

    def unzip(self):
        xs, ys = [], []
        for sample in self.pts:
            x, y = zip(*sample)
            xs.append(x)
            ys.append(y)
        return xs, ys