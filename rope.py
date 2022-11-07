
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

class PtMass:
    def __init__(self,pos,mass):
        self.pos = pos
        self.new_pos = pos
        self.vel = np.array([0,0])
        self.new_vel = self.vel
        self.mass = mass

    def update_state(self,pos,dt):
        """Use a new position to update the velocity."""
        prev_pos = self.pos
        self.pos = pos
        self.vel = (pos - prev_pos) * dt


class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, pts, save):
        self.pts = pts
        self.xs, self.ys = self.unzip()
        fig = plt.figure(figsize=(7, 7))
        ax = plt.axes(xlim=(-30, 30), ylim=(-30, 30))
        self.scatter = ax.scatter(self.xs, self.ys)
        anim = FuncAnimation(fig, self.update, frames=self.pts,interval=.05)
        if save:
            writervideo = animation.FFMpegWriter(fps=60)
            anim.save('rope.mp4', writer=writervideo)
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

#TODO. Make sure that the points are being edited correctly.

class Rope:
    def __init__(self,k,c,start_pos,end_pos,num_pts,rope_density,unstretched_len,initialization='linear',
                 fixed_tail=True,gravity=-9.81,integration_method='verlet'):
        self.k = k
        self.c = c
        self.mass_per_pt = (rope_density * unstretched_len) / num_pts
        self.first = PtMass(start_pos, self.mass_per_pt)
        self.last = PtMass(end_pos, self.mass_per_pt)
        self.pts = [None]*num_pts
        self.num_pts = num_pts
        self.unstretched_len = unstretched_len
        self.lo_per_seg = self.__find_lo_per_seg() # unstretched length per segment
        self.fixed_tail = fixed_tail
        self.GRAVITY = np.array([0, gravity])
        self.integration_method = integration_method
        self.rope_positions = None
        self.initialization = initialization
        self.__initialize(initialization)

    def step(self,dt):
        for i, cur_pt in enumerate(self.pts):
            if cur_pt == self.first: continue
            elif cur_pt == self.last: continue
            prev_pt = self.pts[i - 1]
            cur_pt = self.pts[i]
            next_pt = self.pts[i+1]
            # Find what the step will be, then apply all at once at the end.
            if self.integration_method == 'verlet':
                self.__step_verlet(cur_pt, prev_pt, next_pt, dt)
            elif self.integration_method == 'euler':
                self.__step_euler(cur_pt, prev_pt, next_pt, dt)
            else: raise KeyError
        if not self.fixed_tail: self.__step_tail(prev_pt,dt)
        self.__apply_step()

    def sim(self,num_steps,step_size):
        self.rope_positions = []
        for i in range(num_steps):
            self.step(step_size)
            cur_pos = []
            for pt in self.pts:
                cur_pos.append(pt.pos)
            self.rope_positions.append(cur_pos)
            end_pos = 10*math.sin(i/20) + 4*math.sin(i/10)
            start_pos = -10 * math.sin(i / 25) + 4 * math.sin(i / 15)
            # self.update_ending_pt(np.array([10*math.sin(i/50),end_pos]),step_size)
            # self.update_first_pt(np.array([-end_pos, start_pos]), step_size)

    def animate(self,save=False):
        if self.rope_positions is None: return
        a = AnimatedScatter(self.rope_positions,save)

    def update_first_pt(self, pos, dt):
        self.first.update_state(pos, dt)

    def update_last_pt(self, pos, dt):
        self.last.update_state(pos, dt)

    def update_unstretched_len(self,length):
        self.unstretched_len = length
        self.lo_per_seg = self.__find_lo_per_seg()

    def demo(self,demo_num=0):
        if demo_num == 0: self.__demo0()
        elif demo_num == 1: self.__demo1()
        elif demo_num == 2: self.__demo2()
        elif demo_num == 3: self.__demo2()
        else: raise KeyError

    def __find_acc(self,cur_pos,cur_vel,prev_pos,prev_vel,next_pos,next_vel,mass):
        gravitational_force = mass * self.GRAVITY
        force_prev_pt = self.__find_force_bet_points(cur_pos, cur_vel, prev_pos, prev_vel)
        force_next_pt = self.__find_force_bet_points(cur_pos, cur_vel, next_pos, next_vel)
        return (force_prev_pt + force_next_pt + gravitational_force) / mass

    def __find_acc_end(self,cur_pos,cur_vel,other_pos,other_vel,mass):
        gravitational_force = mass * self.GRAVITY
        force = self.__find_force_bet_points(cur_pos, cur_vel, other_pos, other_vel)
        return (force + gravitational_force) / mass

    def __find_half_step(self,pt, acc_euler, dt):
        v_half_dt = pt.vel + (acc_euler * 0.5*dt)
        # Wiki says to exclude 0.5 for p_half_dt
        p_half_dt = pt.pos + (v_half_dt * dt)
        return p_half_dt,v_half_dt

    def __step_verlet(self, cur_pt, prev_pt, next_pt, dt):
        acc_euler = self.__find_acc(cur_pt.pos,cur_pt.vel,prev_pt.pos,prev_pt.vel,next_pt.pos,next_pt.vel,cur_pt.mass)
        prev_half_step_p, prev_half_step_vel = self.__find_half_step(prev_pt, acc_euler, dt)
        cur_half_step_p, cur_half_step_vel = self.__find_half_step(cur_pt,acc_euler,dt)
        next_half_step_p, next_half_step_vel = self.__find_half_step(next_pt, acc_euler, dt)

        acc = self.__find_acc(cur_half_step_p, cur_half_step_vel, prev_half_step_p, prev_half_step_vel, next_half_step_p, next_half_step_vel,cur_pt.mass)

        cur_pt.new_vel = cur_half_step_vel + 0.5*(acc * acc_euler)*dt
        cur_pt.new_pos = cur_pt.pos + 0.5*(acc_euler * (dt**2))

    def __step_euler(self, cur_pt, prev_pt, next_pt, dt):
        acc = self.__find_acc(cur_pt.pos,cur_pt.vel,prev_pt.pos,prev_pt.vel,next_pt.pos,next_pt.vel,cur_pt.mass)
        cur_pt.new_vel = cur_pt.vel + (acc * dt)
        cur_pt.new_pos = cur_pt.pos + (cur_pt.vel * dt)

    def __step_tail(self,prev_pt,dt):
        cur_pt = self.last
        gravitational_force = cur_pt.mass * self.GRAVITY
        if self.integration_method == 'verlet':
            acc_euler = self.__find_acc_end(cur_pt.pos, cur_pt.vel, prev_pt.pos, prev_pt.vel,cur_pt.mass)
            prev_half_step_p, prev_half_step_vel = self.__find_half_step(prev_pt, acc_euler, dt)
            cur_half_step_p, cur_half_step_vel = self.__find_half_step(cur_pt, acc_euler, dt)
            acc = self.__find_acc_end(cur_half_step_p, cur_half_step_vel, prev_half_step_p, prev_half_step_vel, cur_pt.mass)
            cur_pt.new_vel = cur_half_step_vel + 0.5 * (acc * acc_euler) * dt
            cur_pt.new_pos = cur_pt.pos + 0.5 * (acc_euler * (dt ** 2))
        elif self.integration_method == 'euler':
            force = self.__find_acc_end(cur_pt.pos, cur_pt.vel, prev_pt.pos, prev_pt.vel,cur_pt.mass)
            acc = (force + gravitational_force) / cur_pt.mass
            cur_pt.new_vel = cur_pt.vel + (acc * dt)
            cur_pt.new_pos = cur_pt.pos + (cur_pt.vel * dt)
        else: raise KeyError

    def __apply_step(self):
        for i in range(1, self.num_pts):
            # self.pts[i].pos = self.pts[i].new_pos
            # self.pts[i].vel = np.array([0.95,1]) * self.pts[i].new_vel
            self.pts[i].pos = self.pts[i].new_pos
            self.pts[i].vel = self.pts[i].new_vel
        if not self.fixed_tail:
            self.last.pos = self.last.new_pos
            self.last.vel = self.last.new_vel

    def __get_unit_vec(self,to_pt, from_pt):
        vec_to_pt = to_pt - from_pt
        return vec_to_pt / math.sqrt(vec_to_pt[0]**2 + vec_to_pt[1]**2)

    def __find_force_bet_points(self, cur_pt_pos,cur_pt_vel,other_pt_pos,other_pt_vel):
        unit_vec_to_pt = self.__get_unit_vec(other_pt_pos, cur_pt_pos)
        spring_force = self.__find_spring_force(cur_pt_pos, other_pt_pos, unit_vec_to_pt)
        visc_force = self.__find_visc_force(cur_pt_vel, other_pt_vel, unit_vec_to_pt)
        return spring_force + visc_force

    def __find_spring_force(self, cur_pt_pos, other_pt_pos, unit_vec_to_pt):
        dist = math.dist(other_pt_pos,cur_pt_pos)
        deflection = dist - self.lo_per_seg
        # return self.k * deflection * unit_vec_to_pt
        # Ropes cannot be in compression
        if deflection < 0:
            return np.array([0,0])
        else:
            return self.k * deflection * unit_vec_to_pt

    def __find_visc_force(self,cur_pt_vel, other_pt_vel, unit_vec_to_pt):
        relative_vel = other_pt_vel - cur_pt_vel
        relative_speed = np.dot(relative_vel,unit_vec_to_pt)
        # component_vel = relative_vel * unit_vec_to_pt
        # return -component_vel * self.c
        if relative_speed < 0: # Pts moving towards each other.
            return np.array([0,0])
        else: # Pts moving apart
            component_vel = relative_speed * unit_vec_to_pt
            return component_vel * self.c

    def __initialize(self,initialization):
        self.pts[0] = self.first
        self.pts[-1] = self.last
        if initialization == "linear": self.__linear_init()
        else: raise KeyError

    def __linear_init(self):
        diff = self.last.pos - self.first.pos
        for i in range(1,self.num_pts-1):
            dec_percent = i/(self.num_pts-1)
            pos = self.first.pos + (diff * dec_percent)
            self.pts[i] = PtMass(pos,self.mass_per_pt)

    def __find_lo_per_seg(self):
        return self.unstretched_len/(self.num_pts-1)

def main():
    pass

if __name__ == '__main__':
    main()
