"""
A class for simulating ropes or rope-like things.

Author: Kyle Mathenia
Date: 11/7/2022
"""

import numpy as np
import math

class PtMass:
    def __init__(self,pos,mass):
        self.pos,self.new_pos,self.prev_pos = pos,pos,pos
        self.vel, self.new_vel, self.prev_vel = np.array([0,0]), np.array([0,0]), np.array([0,0])
        self.mass = mass

    def step_pos(self):
        self.prev_pos = self.pos
        self.pos = self.new_pos

    def step_vel(self):
        self.prev_vel = self.vel
        self.vel = self.new_vel

    def update_pos_derive_vel(self, pos, dt):
        """Update position and derive the velocity from that. Useful for update the ends of the rope."""
        self.prev_pos = self.pos
        self.pos = pos
        self.new_vel = (pos - self.prev_pos) * dt

class Rope:
    def __init__(self,k=5,c=100,start_pos=[0,10],end_pos=[20,10],num_pts=16,\
                 rope_density=0.06,unstretched_len=15,initialization='linear',\
                 fixed_tail=False,gravity=-9.81,integration_method='verlet',drag_coef=1000):
        self._mass_per_pt = (rope_density * unstretched_len) / num_pts
        self._first_pt = PtMass(np.array(start_pos), self._mass_per_pt)
        self._last_pt = PtMass(np.array(end_pos), self._mass_per_pt)
        self.k,self.c = k, c
        self._pts = [None] * num_pts
        self._num_pts = num_pts
        self._unstretched_len = unstretched_len
        self._lo_per_seg = self.__find_lo_per_seg() # unstretched length per segment
        self.fixed_tail = fixed_tail
        self.DRAG_COEF = drag_coef
        self.GRAVITY = np.array([0, gravity])
        self.integration_method = integration_method
        self.initialization = initialization
        self.reset()


    ###########################################################################
    ################################# Public ##################################

    def step(self,dt):
        for i, cur_pt in enumerate(self._pts):
            if cur_pt == self._first_pt or cur_pt == self._last_pt: continue

            # Find what the step will be for all the points first, then apply all at once.
            prev_pt, next_pt = self._pts[i - 1], self._pts[i + 1]
            if self.integration_method == 'verlet':
                self.__step_verlet(cur_pt, prev_pt, dt, next_pt=next_pt)
            elif self.integration_method == 'euler':
                self.__step_euler(cur_pt, prev_pt, dt, next_pt=next_pt)
            else: raise KeyError
        if not self.fixed_tail: self.__step_tail(dt)
        self.__apply_steps()
        self.__append_data()

    def sim(self,num_steps=1000,step_size=0.05):
        self.reset()
        for i in range(num_steps):
            self.step(step_size)

    def reset(self):
        self.sim_data = []
        self.__initialize(self.initialization)

    def update_first_pt(self, pos, dt=0):
        self._first_pt.update_pos_derive_vel(np.array(pos), dt)

    def update_last_pt(self, pos, dt=0):
        self._last_pt.update_pos_derive_vel(np.array(pos), dt)

    def update_unstretched_len(self,length):
        """You may want to dynamically change the unstretched length if the rope is going over a pulley, for example."""
        self._unstretched_len = length
        self._lo_per_seg = self.__find_lo_per_seg()


    ###########################################################################
    ################################# Private #################################

    ######## step funcs ########

    def __step_verlet(self, cur_pt, prev_pt, dt, next_pt=None):
        acc = self.__find_acc(cur_pt,prev_pt,next_pt)
        cur_pt.new_pos = (2*cur_pt.pos) - cur_pt.prev_pos + (dt**2) * acc
        cur_pt.new_vel = (cur_pt.new_pos-cur_pt.pos) * dt

    def __step_euler(self, cur_pt, prev_pt, dt, next_pt=None):
        acc = self.__find_acc(cur_pt,prev_pt,next_pt)
        cur_pt.new_vel = cur_pt.vel + (acc * dt)
        cur_pt.new_pos = cur_pt.pos + (cur_pt.vel * dt)

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

    def __find_force_bet_points(self, cur_pt_pos,cur_pt_vel,other_pt_pos,other_pt_vel):
        unit_vec_to_pt = self.__unit_vec(other_pt_pos - cur_pt_pos)
        spring_force = self.__find_spring_force(cur_pt_pos, other_pt_pos, unit_vec_to_pt)
        visc_force = self.__find_visc_force(cur_pt_vel, other_pt_vel, cur_pt_pos, other_pt_pos, unit_vec_to_pt)
        return spring_force + visc_force

    def __find_spring_force(self, cur_pt_pos, other_pt_pos, unit_vec_to_pt):
        dist = math.dist(other_pt_pos,cur_pt_pos)
        deflection = dist - self._lo_per_seg
        # Ropes cannot be in compression
        if deflection < 0: return np.array([0,0])
        else: return self.k * deflection * unit_vec_to_pt

    def __find_visc_force(self,cur_pt_vel, other_pt_vel, cur_pt_pos, other_pt_pos, unit_vec_to_pt):
        relative_vel = cur_pt_vel - other_pt_vel
        relative_speed = np.dot(relative_vel,unit_vec_to_pt)
        component_vel = relative_speed * unit_vec_to_pt
        dist = math.dist(other_pt_pos,cur_pt_pos)
        deflection = dist - self._lo_per_seg
        # Ropes cannot be in compression
        if deflection < 0: return np.array([0,0])
        else: return -component_vel * self.c

    def __find_drag(self,vel):
        speed = self.__mag(vel)
        return (0.5 * (speed**2) * self.DRAG_COEF) * self.__unit_vec(vel)

    def __step_tail(self,dt):
        if self.integration_method == 'verlet': self.__step_verlet(self._pts[-1], self._pts[-2], dt)
        elif self.integration_method == 'euler': self.__step_euler(self._pts[-1], self._pts[-2], dt)
        else: raise KeyError

    def __apply_steps(self):
        for i in range(1, self._num_pts - 1):
            self._pts[i].step_pos()
            self._pts[i].step_vel()
        if not self.fixed_tail:
            self._pts[-1].step_pos()
            self._pts[-1].step_vel()

    def __append_data(self):
        cur_pos = []
        for pt in self._pts:
            cur_pos.append(pt.pos)
        self.sim_data.append(cur_pos)


    ######## other funcs ########

    def __initialize(self,initialization):
        self._pts[0] = self._first_pt
        self._pts[-1] = self._last_pt
        if initialization == "linear": self.__linear_init()
        else: raise KeyError

    def __linear_init(self):
        diff = self._last_pt.pos - self._first_pt.pos
        for i in range(1, self._num_pts - 1):
            dec_percent = i/(self._num_pts - 1)
            pos = self._first_pt.pos + (diff * dec_percent)
            self._pts[i] = PtMass(pos, self._mass_per_pt)

    def __find_lo_per_seg(self):
        return self._unstretched_len / (self._num_pts - 1)

    def __unit_vec(self, vec):
        mag = self.__mag(vec)
        if mag == 0: return np.array([0,0])
        return vec / mag

    def __mag(self,vec):
        return np.linalg.norm(vec)

def main():
    pass

if __name__ == '__main__':
    main()
