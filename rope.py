
import numpy as np
import math

class PtMass:
    def __init__(self,pos,mass):
        self.pos = pos
        self.new_pos = pos
        self.prev_pos = pos
        self.vel = np.array([0,0])
        self.new_vel = self.vel
        self.prev_vel = self.vel
        self.mass = mass

    def update_state(self,pos,dt):
        """Use a new position to update the velocity."""
        self.prev_pos = self.pos
        self.pos = pos
        self.new_vel = (pos - self.prev_pos) * dt

    def step_pos(self):
        self.prev_pos = self.pos
        self.pos = self.new_pos

    def step_vel(self):
        self.prev_vel = self.vel
        self.vel = self.new_vel

class Rope:
    def __init__(self,k=1,c=1,start_pos=[-10,10],end_pos=[10,10],num_pts=6,
                 rope_density=0.06,unstretched_len=10,initialization='linear',
                 fixed_tail=True,gravity=-9.81,integration_method='verlet', drag_coef = 0.9):
        self.mass_per_pt = (rope_density * unstretched_len) / num_pts
        self.first = PtMass(np.array(start_pos), self.mass_per_pt)
        self.last = PtMass(np.array(end_pos), self.mass_per_pt)
        self.k,self.c = k, c
        self.pts = [None]*num_pts
        self.num_pts = num_pts
        self.unstretched_len = unstretched_len
        self.lo_per_seg = self.__find_lo_per_seg() # unstretched length per segment
        self.fixed_tail = fixed_tail
        self.drag_coef = drag_coef
        self.GRAVITY = np.array([0, gravity])
        self.integration_method = integration_method
        self.initialization = initialization
        self.reset()
        self.__initialize(initialization)

    def step(self,dt):
        for i, cur_pt in enumerate(self.pts):
            if cur_pt == self.first or cur_pt == self.last: continue

            # Find what the step will be for all the points first, then apply all at once.
            prev_pt, next_pt = self.pts[i-1], self.pts[i+1]
            if self.integration_method == 'verlet':
                self.__step_verlet(cur_pt, prev_pt, dt, next_pt=next_pt)
            elif self.integration_method == 'euler':
                self.__step_euler(cur_pt, prev_pt, dt, next_pt=next_pt)
            else: raise KeyError
        if not self.fixed_tail: self.__step_tail(dt)
        self.__apply_step()

    def sim(self,num_steps=1000,step_size=0.05):
        self.reset()
        for i in range(num_steps):
            self.step(step_size)
            cur_pos = []
            for pt in self.pts:
                cur_pos.append(pt.pos)
            self.sim_data.append(cur_pos)

    def reset(self):
        self.sim_data = []

    def update_first_pt(self, pos, dt):
        self.first.update_state(pos, dt)

    def update_last_pt(self, pos, dt):
        self.last.update_state(pos, dt)

    def update_unstretched_len(self,length):
        self.unstretched_len = length
        self.lo_per_seg = self.__find_lo_per_seg()

    def __step_tail(self,dt):
        if self.integration_method == 'verlet': self.__step_verlet(self.pts[-1], self.pts[-2], dt)
        elif self.integration_method == 'euler': self.__step_euler(self.pts[-1], self.pts[-2], dt)
        else: raise KeyError

    # def __find_acc(self,cur_pos,cur_vel,prev_pos,prev_vel,next_pos,next_vel,mass):
    #     gravitational_force = mass * self.GRAVITY
    #     force_prev_pt = self.__find_force_bet_points(cur_pos, cur_vel, prev_pos, prev_vel)
    #     force_next_pt = self.__find_force_bet_points(cur_pos, cur_vel, next_pos, next_vel)
    #     drag = self.__find_drag(cur_vel)
    #     return (force_prev_pt + force_next_pt + gravitational_force - drag) / mass

    # def __find_acc_end(self,cur_pos,cur_vel,other_pos,other_vel,mass):
    #     gravitational_force = mass * self.GRAVITY
    #     force = self.__find_force_bet_points(cur_pos, cur_vel, other_pos, other_vel)
    #     drag = self.__find_drag(cur_vel)
    #     return (force + gravitational_force - drag) / mass

    def __find_drag(self,vel):
        return (0.5*(vel**2) * self.drag_coef)

    def __find_acc(self,cur_pt, prev_pt, next_pt=None):
        gravitational_force = cur_pt.mass * self.GRAVITY
        if next_pt is None:
            force_prev_pt = self.__find_force_bet_points(cur_pt.pos, cur_pt.vel, prev_pt.pos, prev_pt.vel)
            drag = self.__find_drag(cur_pt.vel)
            return (force_prev_pt + gravitational_force - drag) / cur_pt.mass
        else:
            force_prev_pt = self.__find_force_bet_points(cur_pt.pos, cur_pt.vel, prev_pt.pos, prev_pt.vel)
            force_next_pt = self.__find_force_bet_points(cur_pt.pos, cur_pt.vel, next_pt.pos, next_pt.vel)
            drag = self.__find_drag(cur_pt.vel)
            return (force_prev_pt + force_next_pt + gravitational_force - drag) / cur_pt.mass

    def __step_verlet(self, cur_pt, prev_pt, dt, next_pt=None):
        acc = self.__find_acc(cur_pt,prev_pt,next_pt)
        cur_pt.new_pos = (2*cur_pt.pos) - cur_pt.prev_pos + (dt**2) * acc
        cur_pt.new_vel = (cur_pt.new_pos-cur_pt.pos) * dt

    def __step_euler(self, cur_pt, prev_pt, dt,next_pt=None):
        acc = self.__find_acc(cur_pt,prev_pt,next_pt)
        cur_pt.new_vel = cur_pt.vel + (acc * dt)
        cur_pt.new_pos = cur_pt.pos + (cur_pt.vel * dt)

    def __apply_step(self):
        for i in range(1, self.num_pts-1):
            self.pts[i].step_pos()
            self.pts[i].step_vel()
        if not self.fixed_tail:
            self.pts[-1].step_pos()
            self.pts[-1].step_vel()

    def __get_unit_vec(self,to_pt, from_pt):
        vec_to_pt = to_pt - from_pt
        return vec_to_pt / math.sqrt(vec_to_pt[0]**2 + vec_to_pt[1]**2)

    def __find_force_bet_points(self, cur_pt_pos,cur_pt_vel,other_pt_pos,other_pt_vel):
        unit_vec_to_pt = self.__get_unit_vec(other_pt_pos, cur_pt_pos)
        spring_force = self.__find_spring_force(cur_pt_pos, other_pt_pos, unit_vec_to_pt)
        visc_force = self.__find_visc_force(cur_pt_vel, other_pt_vel, cur_pt_pos, other_pt_pos, unit_vec_to_pt)
        return spring_force + visc_force

    def __find_spring_force(self, cur_pt_pos, other_pt_pos, unit_vec_to_pt):
        dist = math.dist(other_pt_pos,cur_pt_pos)
        deflection = dist - self.lo_per_seg
        # Ropes cannot be in compression
        if deflection < 0: return np.array([0,0])
        else: return self.k * deflection * unit_vec_to_pt

    def __find_visc_force(self,cur_pt_vel, other_pt_vel, cur_pt_pos, other_pt_pos, unit_vec_to_pt):
        relative_vel = other_pt_vel - cur_pt_vel
        relative_speed = np.dot(relative_vel,unit_vec_to_pt)
        component_vel = relative_speed * unit_vec_to_pt
        dist = math.dist(other_pt_pos,cur_pt_pos)
        deflection = dist - self.lo_per_seg
        # Ropes cannot be in compression
        if deflection < 0: return np.array([0,0])
        else: return component_vel * self.c

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
