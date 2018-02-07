`FindZerom`
============


## Aim
This package is suitable to **compute all the roots of an analytic functions** present inside a closed contour **without initial guess**. There are numerous root finding algorithm like Newton-Raphson method, Muller's method, the Secant method or the Nelder-Mead simplex method. All these techniques have in common that they require initial approximations for the zeros to start the algorithm. 

The proposed implementation is based on the method called the Cauchy Integration Method or the Argument Principle Method and allows to compute the number of zeros (including its multiplicity) of a function from contour integral. A _short_ documentation explains quickly the theoretical background, shows the calling sequence and presents some example of applications and validations. [here](https://github.com/nennigb/FindZerosm/blob/master/Documentation/FindZerosm.pdf)

## Origin
This program was initially developed for poroelastic silencer applications [Nennig, 2010](https://github.com/nennigb/FindZerosm/blob/master/Documentation/Nennig_et_al_jasa_2010.pdf) 

>B. Nennig, E. Perrey-Debain, and M. Ben Tahar. A mode matching method for modelling dissipative silencers lined with poroelastic materials and containing mean flow. J. Acoust. Soc. Am. , 12 (6) :3308-3320, 10.1121/1.3693655, 2010.

in order to solved dispersion equation. The algorithm has been already used in other applications fields (see [Delves, 1967;Chen, 2000; Kravanja, 2000] and the references therein) and we proposed  here a basic numerical implementation in matlab language. Example of applications can be found in my [work](https://cv.archives-ouvertes.fr/benoit-nennig).

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.





