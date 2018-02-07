This package is suitable to compute the roots of an analytic function present inside a closed contour. There are numerous available numerical techniques for this purpose like Newton-Raphson method, Muller's method, the Secant method or the Nelder-Mead simplex method. All these techniques have in common that they all require initial approximations for the zeros to start the algorithm. The proposed implementation is based on the method called the Cauchy Integration Method  or the Argument Principle Method and allows to compute the number of zeros (including its multiplicity) of $f$ from contour integral.

This program was initially developed for poroelastic silencer applications\cite{Nennig:2010} in order to solved dispersion equation. The algorithm has been already used in other application (see \cite{Delves:1967,Chen:2000,Kravanja:2000}) and we proposed here a basic numerical implementation in matlab language.

This \emph{short} documentation explain quickly the theoretical background, show the calling sequence and present some example of applications and validations.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

