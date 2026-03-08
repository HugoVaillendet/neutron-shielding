# neutron-shielding
The goal of this project is to simulate neutron shielding of a physical medium from a punctual neutron source.
We will start from the most basic model and grow in complexity further to match more realistic results.

We will use a 2 dimensional polar system for our entire simulation.

## Monte Carlo methods

Monte Carlo methods are a class of problem solving methods relying on random steps within a selected domain.
This type of problem solving was invented by Stanisław Ulam and Nicholas Metropolis [[1]](#1)

## Initial simualtion

At first we structure the code with a neutron module. We implement a new type neutron containing postions $(x, y)$, speeds $(vx, vy)$ with an SoA structure.
We also implement a `scatter_tally` array as we as `PI` as a parameter evaluated at compile time.

We then implement the following procedures : 

### init()

I this subroutine we deallocate then allocate arrays size and assign postions and `scatter_tally` to $0$.
This prevents errors in case multiple calls of this subroutine

We also initialize a `seed` to allow repeatable results.

### boundary_check()
### sample_direction()
### sample_free_path()
### evaluate_step()



## References
<a id="1">[1]</a> 
Metropolis, N., & Ulam, S. (1949).
The Monte Carlo Method. 
Journal of the American Statistical Association, 44(247), 335–341.
https://doi.org/10.1080/01621459.1949.10483310
