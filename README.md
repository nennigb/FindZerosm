`FindZerom`
============


## Aim
This package is suitable to **compute all the roots of an analytic functions** present inside a closed contour **without initial guess**. There are numerous root finding algorithm like Newton-Raphson method, Muller's method, the Secant method or the Nelder-Mead simplex method. All these techniques have in common that they require initial approximations for the zeros to start the algorithm. The approach used here is to build a polynomial with the same zeros as the original function.

The proposed implementation is based on the method called the Cauchy Integration Method or the Argument Principle Method and allows to compute the number of zeros (including its multiplicity) of a function from contour integral. A [_short_ documentation](https://github.com/nennigb/FindZerosm/blob/master/Documentation/FindZerosm.pdf) explains quickly the theoretical background, shows the calling sequence and presents some example of applications and validations. 

## Origin
This program was initially developed for poroelastic silencer applications [Nennig, 2010](https://github.com/nennigb/FindZerosm/blob/master/Documentation/Nennig_et_al_jasa_2010.pdf) 

>B. Nennig, E. Perrey-Debain, and M. Ben Tahar. A mode matching method for modelling dissipative silencers lined with poroelastic materials and containing mean flow. J. Acoust. Soc. Am. , 12 (6) :3308-3320, 10.1121/1.3693655, 2010.

in order to solved dispersion equation. The algorithm has been already used in other applications fields (see [Delves, 1967;Chen, 2000; Kravanja, 2000] and the references therein) and we proposed  here a basic numerical implementation in matlab language. Examplew of applications in acoustics can be found in my [work](https://cv.archives-ouvertes.fr/benoit-nennig).

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

## Basic Usage
This command compute the roots of fhandle inside a circle of radius `R` in the complex plane,
`K = FindZerosm(R,N,fhandle)`
Where `N` is the number of integration points (500 is a good start) and `fhandle` is the anonymous function of which the root are sought. 

