# rope_sim

A simple rope simulation I made for funzies at some point.

Parameters:
- fixed or constrained tail
- Verlet or Euler integration
- spring constant
- viscous damping constant
- drag coefficient
- number of points
- unstretched length
- gravity

```
# Springy
rope = Rope(k=5, c=10, start_pos=[0, 10], end_pos=[20, 10], num_pts=16, unstretched_len=15, fixed_tail=False, drag_coef=100)
rope.sim(num_steps=1000,step_size=0.05)
```

![springy](https://media.giphy.com/media/hbvvXZNo0W6ZkiKgTc/giphy.gif)

```
# More damping and drag
rope.c = 100
rope.drag_coef = 1000
```

![damped](https://media.giphy.com/media/C0dNc5XVss67AAO4Ka/giphy.gif)


```
# Update start and end points and step manually. Step size can vary in size between steps, which may be useful for real-time loops.
rope = Rope(k=5, c=20, start_pos=[0, 10], end_pos=[20, 0], num_pts=16, unstretched_len=10, fixed_tail=True, drag_coef=50)
num_steps, step_size = 1000, 0.05
for i in range(num_steps):
    rope.update_first_pt([15* math.sin(i/40),10* math.sin(i/15) + 10])
    rope.update_last_pt([20, 20 * math.sin(i / 60)])
    rope.step(step_size)
```

![both_moving](https://media.giphy.com/media/F9ZC74TZJJMYTQtj2y/giphy.gif)

```
# Other demos:
```

![head_moving](https://media.giphy.com/media/DJLv5HVLfMpvqkcDoX/giphy.gif)
![constrained](https://media.giphy.com/media/wl56Ia4c77fSffSHWd/giphy.gif)

Optionally, use Euler integration and observe [numerical instability!](https://en.wikipedia.org/wiki/Euler_method#/media/File:Instability_of_Euler's_method.svg)
```
rope = Rope(k=2, c=0, start_pos=[0, 20], end_pos=[20, 10], num_pts=8, unstretched_len=15, integration_method='euler', drag_coef=0.005)
rope.sim(num_steps=300,step_size=0.05)
```

![euler_instability](https://media.giphy.com/media/Wv3by7uBcN779ZBjJf/giphy.gif)
