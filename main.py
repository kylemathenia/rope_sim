
import rope
from animated_scatter import AnimatedScatter

rope = rope.Rope(k=1,c=100,start_pos=[-10,10],end_pos=[10,10],num_pts=6,\
                 rope_density=0.06,unstretched_len=10,initialization='linear',\
                 fixed_tail=False,gravity=-9.81,integration_method='verlet',drag_coef=1000)

def demo1():
    rope.sim(num_steps=1000,step_size=0.05)
    a = AnimatedScatter(rope.sim_data, save=False, filename="rope_vid.mp4")

def demo2():
    pass

def demo3():
    pass

def demo4():
    pass

def main():
    demo1()
    demo2()
    demo3()
    demo4()

if __name__ == "__main__":
    main()