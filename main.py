"""
Rope simulation demos.

Author: Kyle Mathenia
Date: 11/7/2022
"""

from rope import Rope
from animated_scatter import AnimatedScatter
import math


def demo1():
    """Damped/less springy/more rope-like with one end free."""
    rope1 = Rope(k=5, c=100, start_pos=[0, 10], end_pos=[20, 10], num_pts=16, unstretched_len=15, fixed_tail=False,\
                 drag_coef=1000)
    rope1.sim(num_steps=1000,step_size=0.05)
    a = AnimatedScatter(rope1.sim_data, save=False, filename="damped.mp4")


def demo2():
    """Springy with both ends constrained."""
    rope1.fixed_tail = True
    rope1.update_first_pt([-15, 10])
    rope1.update_last_pt([15, 10])
    rope1.update_unstretched_len(45)
    rope1.sim(num_steps=1000, step_size=0.05)
    a = AnimatedScatter(rope1.sim_data, save=False, filename="springy_constrained.mp4")


def demo3():
    """Springy with one end free."""
    rope2 = Rope(k=5, c=10, start_pos=[0, 10], end_pos=[20, 10], num_pts=16, \
                 rope_density=0.06, unstretched_len=15, initialization='linear', \
                 fixed_tail=False, gravity=-9.81, integration_method='verlet', drag_coef=100)
    rope2.sim(num_steps=1000, step_size=0.05)
    a = AnimatedScatter(rope2.sim_data, save=False, filename="springy.mp4")


def demo4():
    """One end free and the start of the rope moving around."""
    rope3 = Rope(k=5, c=40, start_pos=[0, 10], end_pos=[20, 10], num_pts=16, \
                     rope_density=0.06, unstretched_len=10, initialization='linear', \
                     fixed_tail=False, gravity=-9.81, integration_method='verlet', drag_coef=100)
    num_steps = 1000
    step_size = 0.05
    for i in range(num_steps):
        head_pos = [15* math.sin(i/40),10* math.sin(i/15) + 10]
        rope3.update_first_pt(head_pos)
        rope3.step(step_size)

    a = AnimatedScatter(rope3.sim_data, save=False, filename="update_head.mp4")


def demo5():
    """Both ends constrained and moving around."""
    rope4 = Rope(k=5, c=20, start_pos=[0, 10], end_pos=[20, 0], num_pts=16, \
                     rope_density=0.06, unstretched_len=10, initialization='linear', \
                     fixed_tail=True, gravity=-9.81, integration_method='verlet', drag_coef=50)
    num_steps = 1000
    step_size = 0.05
    for i in range(num_steps):
        head_pos = [15* math.sin(i/40),10* math.sin(i/15) + 10]
        rope4.update_first_pt(head_pos)
        tail_pos = [20, 20 * math.sin(i / 60)]
        rope4.update_last_pt(tail_pos)
        rope4.step(step_size)
    a = AnimatedScatter(rope4.sim_data, save=False, filename="update_both.mp4")


def demo6():
    """Show how Euler integration can become numerically unstable, and why it is typically not used for rope simulations.
    https://en.wikipedia.org/wiki/Euler_method#/media/File:Instability_of_Euler's_method.svg"""
    rope5 = Rope(k=2, c=0, start_pos=[0, 20], end_pos=[20, 10], num_pts=8, \
                 rope_density=0.06, unstretched_len=15, initialization='linear', \
                 fixed_tail=False, gravity=-9.81, integration_method='euler', drag_coef=0.005)
    rope5.sim(num_steps=300,step_size=0.05)
    a = AnimatedScatter(rope5.sim_data, save=False, filename="euler_instability.mp4")

def main():
    demo1()
    demo2()
    demo3()
    demo4()
    demo5()
    demo6()

if __name__ == "__main__":
    main()