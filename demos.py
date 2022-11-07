import rope
import numpy as np

k = 10
c = 0
start_pos = np.array([-10, 10])
end_pos = np.array([5, 10])
num_pts = 10
rope_density = 0.06  # kg
unstretched_len = 5  # meters
num_steps = 900
step_size = 0.005

# Demo 1
rope = rope.Rope(k,c,start_pos,end_pos,num_pts,rope_density,unstretched_len,initialization='linear',
                 fixed_tail=False,gravity=-9.81,integration_method='euler')
rope.sim(num_steps,step_size)
rope.animate(save=False)

# Demo 2
rope.fixed_tail=False
rope.sim(num_steps,step_size)
rope.animate(save=False)

# Demo 3
rope.fixed_tail=False
rope.sim(num_steps,step_size)
rope.animate(save=False)

rope = rope.Rope(k,c,start_pos,end_pos,num_pts,rope_density,unstretched_len,initialization='linear',
                 fixed_tail=True,gravity=-9.81,integration_method='verlet')
rope.sim(num_steps,step_size)
rope.animate(save=False)