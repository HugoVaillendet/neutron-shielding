# neutron-shielding
The goal of this project is to simulate neutron shielding of a physical medium from a punctual neutron source.
We will start from the most basic model and grow in complexity further to match more realistic results.

Our problem will be in 3 dimensions. Our source will emit in a spherical way and our shield will be a cylindrical shape.

Most of the numerical methods used here are inspired form the *Particle-Transport Simulation with the Monte Carlo Method* [[1]](#1)

## Monte Carlo methods

Monte Carlo methods are a class of problem solving methods relying on random steps within a selected domain.
Although previous work predates this article, the modern Markov Chain Monte Carlo was published by Stanisław Ulam and Nicholas Metropolis [[2]](#2)

## Program structure

The program will have a main execution program and a module in a secondary file for neutron type as well as procedures.

The project is compiled using GFortran and with the following flags : `gfortran  neutron.f90 main.f90 -O3 -march=native -ffast-math -fopenmp -funroll-loops -fopt-info-vec -o main.exe`.

The project takes avantage of Fortran's ability to vectorize loops and use SIMD registers as well as OpenMP parallelism. Regarding my specific architecture
AVX-256 registers are used by the program.

We implement a new type `neutron` containing postions $(x, y, z)$, speeds $(v_x, v_y, v_z)$ with an SoA structure.
We also implement a `scatter_tally` array as well as `PI` as a parameter evaluated at compile time.

We then implement the following procedures : 

### init()

I this subroutine we deallocate then allocate arrays size and assign postions and `scatter_tally` to $0$.
This prevents errors in case multiple calls of this subroutine

We also initialize a `seed` to allow repeatable results.

### sample_energy()

Neutrons radiated by a neutron source are not monoenergetic. They follow a distribution called the Watt spectrum :

$$f(E) = \frac{2 e^{ab/4}}{\sqrt{\pi a^3 b}}\, e^{-E/a} \sinh\bigl(\sqrt{bE}\bigr)$$

ou $a$ et $b$ are intrinsic parameters of the considered source.

To sample energy in for our MC simulation, we use an algorithm often called R12" but in the report it's actually the R11 (for Rejection 11) algorithm by Everett and Cashwell [[3]](#3).

Although this algorithm is standard the original paper proved difficult to find and thus I based my implementation of a 
recent article Miao, J. and Jin, M. [[4]](#4) explaining in details why this algorithm works and produces a Watt spectrum sampler.
The algorithm works as follows :

We define $K = 1 + ab/8$, $L = a(K + \sqrt{K^2 - 1}$ and $M = L/a - 1$.\
Then we sample two random varaibles $(\xi_1, \xi_2) \in [0, 1[$.\
Next we set $x = -\log{\xi_1}$ and $y = -\log{\xi_2}$\
**if** $y - M(x + 1)^2 \leq bLx$ **then**\
  accept and return $Lx$\
**else**\
  reject

We test this sampler against the theoretical spectrum for avery large amount of neutrons and considering a U235 source:

![Watt sampler vs theoretical sepctrum](/fig1.png)

and we observe the sampler matches the theory perfectly and that the overwhelming majority of neutrons fall in the fast neutron category.
### boundary_check()

This function is very basic and checks if a a neutron has escaped the shielding boundary and becoming a transfered neutron.
We test the using the square like follows :

$$x^2 + y^2 < R^2\quad \text{and} \quad z^2 < h^2$$

to reduce calculation costs (the sqrt function is costly for CPUs to execute).

### update_energy()

After a scatter event we have to update the energy after collision. We use an elsatic collision model for this.

### sample_direction()

After a scatter event we sample a new random direction for the neutron to go to. We sample a random angle $\theta \in [0, 2\pi]$  as well as a $\phi \in[0, \pi]$ and then apply the new $vx$, $vy$ and $vz$ components.

### sample_free_path()

When using the Monte Carlo method we want to determine the free path, the characteristic distance before a neutron interacts with the medium it is being transported in.
For this we define in collision problems a quantity called the cross-seciton. It describes the probability of a given particle interacting with it's medium.

From Carter and Chashwell's book we have the probability of first collision between $l$ and $l+dl$ along the flight path of the electron:

$$p(l)dl = e^{-\Sigma_t l}\Sigma_tdl$$

with $\Sigma_t$ the macroscopic total cross-section of the medium.
We set $\xi$ as the Cumulative Distribution Function (CDF) along the fight path $l$ and we find:

$$\xi = 1 - e^{-\Sigma_tl}$$

$$l = -\frac{1}{\Sigma_t}\ln{(1-\xi)}$$

and since $1-\xi$ is also uniformly distributed on $[0,1[$:

$$l = -\frac{1}{\Sigma_t}\ln{(\xi)}$$

We obtain our known sampling formula.

### evaluate_step()

First we sample a direction. Then while the neutron is active, we first asign it a regime depending on it's energy and set the cross secitons accordingly.

Then we sample the free path, advance the neutron from path in sampled direction. Test if the neutron escpaed the box.\
If yes we push the neutrons final positions to the `fx`, `fy` and `fz` vectors.\
If no we test for scattering. If yes we apply elastic energy update as well as increment the `scatter_tally`.\
Otherwise we do nothing.

## Observations

## References

<a id="1">[1]</a> 
Carter, L. L., & Cashwell, E. D. (1975).
Particle-Transport Simulation with the Monte Carlo Method. 
Technical Information Center, U.S. Energy Research and Development Administration (ERDA Critical Review Series), Los Alamos Scientific Laboratory.
https://mcnp.lanl.gov/pdf_files/Book_MonteCarlo_1975_Carter_ParticleTransportSimulationwiththeMonteCarloMethod.pdf

<a id="2">[2]</a> 
Metropolis, N., & Ulam, S. (1949).
The Monte Carlo Method. 
Journal of the American Statistical Association, 44(247), 335–341.
https://doi.org/10.1080/01621459.1949.10483310

<a id="3">[3]</a> 
Everett, C. J. and Cashwell, E. (1972).
Monte Carlo Sampler. 
Tech. rep., Los Alamos Scientific Lab., N. Mex.
https://www.osti.gov/servlets/purl/4589395-XhDPtN/

<a id="4">[4]</a> 
Miao, J. and Jin, M. (2024).
Understanding the Sampling Algorithm for Watt Spectrum. 
Transactions of American Nuclear Society, Vol 130 (2024) 1005-1007
https://doi.org/10.48550/arXiv.2402.09454
