# d2_qpgf_d1
This is the calculation program of quasi-periodic Green's function for the Helmholtz equations. The quasi-periodicity is 1-dimension ( x component only ), Green's function is 2-dimensions. 

## Definitions
- quasi-periodic Green's function  
  <img src="https://latex.codecogs.com/gif.latex?\,^q\!G(\mathbf{r})=\sum_{l=-\infty}^{\infty}G(\mathbf{r}+l\mathbf{d})\exp(il\mathbf{k}\cdot\mathbf{d})">  
  <img src="https://latex.codecogs.com/gif.latex?G(\mathbf{r})=\frac{i}{4}H_0^{(1)}(k|\mathbf{r}|)">  
  <img src="https://latex.codecogs.com/gif.latex?\,^q\!G(\mathbf{r})"> is quasi-peridic Green's function  
  <img src="https://latex.codecogs.com/gif.latex?G(\mathbf{r})"> is Green's function of the 2-dimensional Helmholtz equation  
  <img src="https://latex.codecogs.com/gif.latex?\mathbf{k}"> is wave number vector, <img src="https://latex.codecogs.com/gif.latex?|\mathbf{k}|=k">  
  <img src="https://latex.codecogs.com/gif.latex?\mathbf{d}"> is lattice vector, <img src="https://latex.codecogs.com/gif.latex?\mathbf{d}=(d,0)">
  
- quasi-periodic condition  

  <img src="https://latex.codecogs.com/gif.latex?\,^q\!G(\mathbf{r}+l\mathbf{d})=\exp(-il\mathbf{k}\cdot\mathbf{d})\,^q\!G(\mathbf{r})">,
  <img src="https://latex.codecogs.com/gif.latex?l\in\mathbb{Z}">
## Usage of example code  
1. type 'make' command to compile
2. type './example.out' to run  

Please see src/d2_qpgf_d1.h for detail of functions, src/example.c for detail of function usages.

## References
1. Capolino, Filippo, Donald R. Wilton, and William A. Johnson. "Efficient computation of the 2-D Green's function for 1-D periodic structures using the Ewald method." IEEE Transactions on Antennas and Propagation 53.9 (2005): 2977-2984.  
2. Beylkin, Gregory, Christopher Kurcz, and Lucas Monzón. "Fast algorithms for Helmholtz Green's functions." Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences 464.2100 (2008): 3301-3326.
3. Linton, Chris M. "Lattice sums for the Helmholtz equation." SIAM review 52.4 (2010): 630-674.
4. Abramowitz, Milton, Irene A. Stegun, and Robert H. Romer. "Handbook of mathematical functions with formulas, graphs, and mathematical tables." (1988).
5. Faddeeva Package. http://ab-initio.mit.edu/Faddeeva
