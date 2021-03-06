/*---------------------------------------------------------------------------*\
 =========                 |
 \\      /  F ield         | Unsupported Contributions for OpenFOAM
  \\    /   O peration     |
   \\  /    A nd           | Copyright (C) 2014 F. Contino, S. Backaert,
    \\/     M anipulation  |                    N. Bourgeois, T. Lucchini
-------------------------------------------------------------------------------
 License
     This file is a derivative work of OpenFOAM.

     OpenFOAM is free software: you can redistribute it and/or modify it
     under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
     for more details.

     You should have received a copy of the GNU General Public License
     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemPointISAT.H"
#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//Defined as static to be able to dynamicly change it during simulations
//(all chemPoints refer to the same object)
template<class CompType, class ThermoType>
Foam::scalar Foam::chemPointISAT<CompType, ThermoType>::tolerance_;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::chemPointISAT<CompType, ThermoType>::qrDecompose
(
    const label nCols,
    scalarRectangularMatrix& R
)
{
    scalarField c(nCols);
    scalarField d(nCols);
    scalar scale, sigma, sum;
    
    for (label k=0; k<nCols-1; k++)
    {
        scale = 0.0;
        for (label i=k; i<nCols; i++)
        {
            scale=max(scale, fabs(R[i][k]));
        }
        if (scale == 0.0)
        {
            c[k]=d[k]=0.0;
        }
        else
        {
            for (label i=k; i<nCols; i++)
            {
                R[i][k] /= scale;
            }
            sum = 0.0;
            for (label i=k; i<nCols; i++)
            {
                sum += sqr(R[i][k]);
            }
            sigma = sign(R[k][k])*sqrt(sum);
            R[k][k] += sigma;
            c[k]=sigma*R[k][k];
            d[k]=-scale*sigma;
            for (label j=k+1; j<nCols; j++)
            {
                sum=0.0;
                for ( label i=k; i<nCols; i++)
                {
                    sum += R[i][k]*R[i][j];
                }
                scalar tau = sum/c[k];
                for ( label i=k; i<nCols; i++)
                {
                    R[i][j] -= tau*R[i][k];
                }
            }
        }
    }
    d[nCols-1] = R[nCols-1][nCols-1];
    //form R
    for (label i=0; i<nCols; i++)
    {
        R[i][i] = d[i];
        for ( label j=0; j<i; j++)
        {
            R[i][j]=0.0;
        }
    }
}


template<class CompType, class ThermoType>
void Foam::chemPointISAT<CompType, ThermoType>::qrUpdate
(
    scalarRectangularMatrix& R,
    const label n,
    const Foam::scalarField &u,
    const Foam::scalarField &v
)
{
    label k,i;
    scalarField w(u);
    for (k=n-1;k>=0;k--)
    {
        if (w[k] != 0.0)
        {
            break;
        }
    }
    if (k < 0)
    {
        k=0;
    }
    for (i=k-1;i>=0;i--)
    {
        rotate(R, i,w[i],-w[i+1], n);
        if (w[i] == 0.0)
        {
            w[i]=fabs(w[i+1]);
        }
        else if (fabs(w[i]) > fabs(w[i+1]))
        {
            w[i]=fabs(w[i])*sqrt(1.0+sqr(w[i+1]/w[i]));
        }
        else
        {
            w[i]=fabs(w[i+1])*sqrt(1.0+sqr(w[i]/w[i+1]));
        }
    }
    for (i=0;i<n;i++)
    {
        R[0][i] += w[0]*v[i];
    }
    for (i=0;i<k;i++)
    {
        rotate(R, i,R[i][i],-R[i+1][i], n);
    }
}


template<class CompType, class ThermoType>
void Foam::chemPointISAT<CompType, ThermoType>::rotate
(
    scalarRectangularMatrix& R,
    const label i,
    const scalar a,
    const scalar b,
    label n
)
{
    label j;
    scalar c, fact, s, w, y;
    if (a == 0.0)
    {
        c=0.0;
        s=(b >= 0.0 ? 1.0 : -1.0);
    }
    else if (fabs(a) > fabs(b))
    {
        fact = b/a;
        c=sign(a)/sqrt(1.0+(fact*fact));
        s=fact*c;
    }
    else
    {
        fact=a/b;
        s=sign(b)/sqrt(1.0+(fact*fact));
        c=fact*s;
    }
    for (j=i;j<n;j++)
    {
        y=R[i][j];
        w=R[i+1][j];
        R[i][j]=c*y-s*w;
        R[i+1][j]=s*y+c*w;
    }
}


template<class CompType, class ThermoType>
void Foam::chemPointISAT<CompType, ThermoType>::svd
(
    scalarRectangularMatrix& A,
    label m,
    label n,
    scalarDiagonalMatrix& d,
    scalarRectangularMatrix& V
)
{
    //UPDATED VERSION NR3
    bool flag;
    label i,its,j,jj,k,l,nm;
    scalar anorm,c,f,g,h,s,scale,x,y,z;
    scalarField rv1(n);
    scalar eps = std::numeric_limits<scalar>::epsilon();
    g = scale = anorm = 0.0;

    //Householder reduction to bidiagonal form 
    for ( i = 0; i<n; i++)
    {
        l=i+2; //change from i+1 to i+2
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i < m)
        {
            for (k=i;k<m;k++)
            {
                scale += fabs(A[k][i]);
            }
            if (scale != 0.0)
            {
                for ( k=i;k<m;k++)
                {
                    A[k][i] /= scale;
                    s += A[k][i]*A[k][i];
                }
                f = A[i][i];
                g = -sign(f)*sqrt(s);
                h = f*g-s;
                A[i][i]=f-g;
                for (j=l-1;j<n;j++)
                {
                    for (s=0.0,k=i;k<m;k++)
                    {
                        s += A[k][i]*A[k][j];
                    }
                    f = s/h;
                    for (k=i; k<m;k++)
                    {
                        A[k][j] += f*A[k][i];
                    }
                }
                for (k=i; k<m;k++)
                {
                    A[k][i] *= scale;
                }
            }
        }
        d[i] = scale * g;
        g=s=scale=0.0;
        
        if (i+1 <= m && i+1 != n)
        {
            for (k=l-1; k<n; k++)
            {
                scale += fabs(A[i][k]);
            }
            if (scale != 0.0)
            {
                for (k=l-1; k<n; k++)
                {
                    A[i][k] /= scale;
                    s += A[i][k]*A[i][k];
                }
                f = A[i][l-1];
                g = -sign(f)*sqrt(s);
                h = f*g-s;
                A[i][l-1] = f-g;
                for (k=l-1; k<n; k++)
                {
                    rv1[k]=A[i][k]/h;
                }
                for (j=l-1; j<m; j++)
                {
                    for (s=0.0,k=l-1; k<n; k++)
                    {
                        s += A[j][k]*A[i][k];
                    }
                    for (k=l-1; k<n; k++)
                    {
                        A[j][k] += s*rv1[k];
                    }
                }
                for (k=l-1; k<n; k++)
                {
                    A[i][k] *= scale;
                }
            }
        }
        anorm = max(anorm, (fabs(d[i])+fabs(rv1[i])));
    }//end Householder reduction
    
    //Accumulation of right-hand transformations
    for (i=n-1; i>=0; i--)
    {
        if (i < n-1)
        {
            if (g != 0.0)
            {
                for (j=l; j<n; j++)
                {
                    V[j][i] = (A[i][j]/A[i][l])/g;
                }
                for (j=l; j<n; j++)
                {
                    for (s=0.0,k=l; k<n; k++)
                    {
                        s += A[i][k]*V[k][j];
                    }
                    for (k=l; k<n; k++)
                    {
                        V[k][j] += s*V[k][i];
                    }
                }
            }
            for (j=l; j<n; j++)
            {
                V[i][j]=V[j][i]=0.0;
            }
        }
        V[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    //Accumulation of left-hand transformations
    for (i = min(m,n)-1; i>=0; i--)
    {
        l=i+1;
        g=d[i];
        for (j=l; j<n; j++)
        {
            A[i][j] = 0.0;
        }
        if (g != 0.0)
        {
            g = 1.0/g;
            for (j=l; j<n; j++)
            {
                for (s=0.0, k=l; k<m; k++)
                {
                    s+= A[k][i]*A[k][j];
                }
                f = (s/A[i][i])*g;
                for (k=i; k<m; k++)
                {
                    A[k][j] += f*A[k][i];
                }
            }
            for (j=i; j<m; j++)
            {
                A[j][i] *= g;
            }
        }
        else
        {
            for (j=i; j<m; j++)
            {
                A[j][i]=0.0;
            }
        }
        ++A[i][i];
    }
    
    //Diagonalization of the bidiagonal form :
    //Loop over singular values, and over allowed iteration
    for (k=n-1; k>=0; k--)
    {
        for (its=0; its<30; its++)
        {
            flag=true;
            // Test for splitting (rv1[1] always zero)
            for (l=k; l>=0; l--)
            {
                nm = l-1;
                if (l == 0 || fabs(rv1[l]) <= eps*anorm)
                {
                    flag = false;
                    break;
                }
                if (fabs(d[nm]) <= eps*anorm)
                {
                    break;
                }
            }
            //Cancellation of rv1[l], if l>1
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i=l; i<k+1; i++)
                {
                    f = s*rv1[i];
                    rv1[i] = c*rv1[i];
                    if (fabs(f) <= eps*anorm)
                    {
                        break;
                    }
                    g = d[i];
                    h = pythag(f,g);
                    d[i] = h;
                    h = 1.0/h;
                    c = g*h;
                    s = -f*h;
                    for (j=0; j<m; j++)
                    {
                        y = A[j][nm];
                        z = A[j][i];
                        A[j][nm] = y*c + z*s;
                        A[j][i] = z*c - y*s;
                    }
                }
            }
            
            z = d[k];
            if (l == k) //Convergence
            {
                if (z < 0.0) //Singular value is made nonnegative
                {
                    d[k] = -z;
                    for (j=0; j<n; j++)
                    {
                        V[j][k] = -V[j][k];
                    }
                }
                break;
            }
            if (its == 34)
            {
                WarningIn
                (
                    "SVD::SVD"
                    "(scalarRectangularMatrix& A, const scalar minCondition)"
                )   << "no convergence in 35 SVD iterations"
                    << endl;
            }

            x = d[l];
            nm = k-1;
            y = d[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g = pythag(f,1.0);
            f = ((x-z)*(x+z)+h*((y/(f+sign(f)*g))-h))/x;
            c=s=1.0; 
            //Next QR transformation
            for (j=l; j<=nm; j++)
            {
                i = j+1;
                g = rv1[i];
                y = d[i];
                h = s*g;
                g = c*g;
                z = pythag(f,h);
                rv1[j] = z;
                c = f/z;
                s = h/z;
                f = x*c + g*s;
                g = g*c - x*s;
                h = y*s;
                y *= c;
                for (jj=0; jj<n; jj++)
                {
                    x = V[jj][j];
                    z = V[jj][i];
                    V[jj][j] = x*c + z*s;
                    V[jj][i] = z*c - x*s;
                }
                z = pythag(f,h);
                d[j] = z;
                if (z)
                {
                    z = 1.0/z;
                    c = f*z;
                    s = h*z;
                }
                f = c*g + s*y;
                x = c*y - s*g;
                for (jj=0; jj<m; jj++)
                {
                    y = A[jj][j];
                    z = A[jj][i];
                    A[jj][j] = y*c + z*s;
                    A[jj][i] = z*c - y*s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            d[k] = x;
        }
    }
}//end svd


// pythag function used in svd
// compute (a^2+b^2)^1/2 without descrutive underflow or overflow
template<class CompType, class ThermoType>
Foam::scalar
Foam::chemPointISAT<CompType, ThermoType>::pythag(scalar a, scalar b)
{
    scalar absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if (absa > absb)
    {
        return absa*sqrt(1.0+sqr(absb/absa));
    }
    else
    {
        return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+sqr(absa/absb)));
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class CompType, class ThermoType>
Foam::chemPointISAT<CompType, ThermoType>::chemPointISAT
(
TDACChemistryModel<CompType, ThermoType>& chemistry,
const scalarField& phi,
const scalarField& Rphi,
const scalarRectangularMatrix& A,
const scalarField& scaleFactor,
const scalar& tolerance,
const label& completeSpaceSize,
const dictionary& coeffsDict,
binaryNode<CompType, ThermoType>* node
)
:
    chemistry_(chemistry),
    phi_(phi),
    Rphi_(Rphi),
    A_(A),
    scaleFactor_(scaleFactor),
    node_(node),
    completeSpaceSize_(completeSpaceSize),
    nGrowth_(0),
    nActiveSpecies_(chemistry.mechRed()->NsSimp()),
    completeToSimplifiedIndex_(completeSpaceSize-2),
    simplifiedToCompleteIndex_(nActiveSpecies_),
    timeTag_(chemistry_.time().timeOutputValue()),
    lastTimeUsed_(chemistry_.time().timeOutputValue()),
    toRemove_(false),
    maxNumNewDim_(coeffsDict.lookupOrDefault("maxNumNewDim",0)),
    printProportion_(coeffsDict.lookupOrDefault("printProportion",false))
{
    tolerance_=tolerance;

    bool isMechRedActive = chemistry_.mechRed()->active();
    if (isMechRedActive)
    {
        for (label i=0; i<completeSpaceSize-2; i++)
        {
            completeToSimplifiedIndex_[i] =
                chemistry.completeToSimplifiedIndex()[i];
        }
        for (label i=0; i<nActiveSpecies_; i++)
        {
            simplifiedToCompleteIndex_[i] =
                chemistry.simplifiedToCompleteIndex()[i];
        }
    }
    
    label reduOrCompDim = completeSpaceSize;
    if (isMechRedActive)
    {
        reduOrCompDim = nActiveSpecies_+2;
    }

    //SVD decomposition A= U*D*V^T 
    scalarRectangularMatrix Atmp(A);//A computed in ISAT.C
    scalarRectangularMatrix B(reduOrCompDim,reduOrCompDim,0.0);
    DiagonalMatrix<scalar> diag(reduOrCompDim,0.0);
    svd(Atmp, reduOrCompDim, reduOrCompDim, diag, B);

    //replace the value of vector diag by max(diag, 1/2), first ISAT paper, Pope
    for (label i=0; i<reduOrCompDim; i++)
    {
        diag[i]=max(diag[i], 0.5);
    }

    //rebuild A with max length, tol and scale factor before QR decomposition
    scalarRectangularMatrix Atilde(reduOrCompDim,reduOrCompDim);
    //result stored in Atilde
    multiply(Atilde, Atmp, diag, B.T());

    for (label i=0; i<reduOrCompDim; i++)//on species loop
    {
        for (label j=0; j<reduOrCompDim; j++)//species, T and p loop
        {
            label compi=i;
            if (isMechRedActive)
            {
                compi = simplifiedToCompleteIndex(i);
            }
            //SF*A/tolerance (where SF is diagonal with inverse of scale factors)
            //SF*A is the same as dividing each line by the scale factor
            //corresponding to the species of this line
            Atilde[i][j] /= (tolerance*scaleFactor[compi]);
        }
    }
  
    //The object LT_ (the transpose of the Q) describe the EOA, since we have
    // A^T B^T B A that should be factorized into L Q^T Q L^T and is set in the
    //qrDecompose function
    LT_ = scalarRectangularMatrix(Atilde);

    qrDecompose(reduOrCompDim,LT_);
}


template<class CompType, class ThermoType>
Foam::chemPointISAT<CompType, ThermoType>::chemPointISAT
(
    Foam::chemPointISAT<CompType, ThermoType>& p
)
:
    chemistry_(p.chemistry()),
    phi_(p.phi()),
    Rphi_(p.Rphi()),
    LT_(p.LT()),
    A_(p.A()),
    scaleFactor_(p.scaleFactor()),
    node_(p.node()),
    completeSpaceSize_(p.completeSpaceSize()),
    nGrowth_(p.nGrowth()),
    nActiveSpecies_(p.nActiveSpecies()),
    completeToSimplifiedIndex_(p.completeToSimplifiedIndex()),
    simplifiedToCompleteIndex_(p.simplifiedToCompleteIndex()),
    timeTag_(p.timeTag()),
    lastTimeUsed_(p.lastTimeUsed()),
    toRemove_(p.toRemove()),
    maxNumNewDim_(p.maxNumNewDim())
{
   tolerance_ = p.tolerance();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
bool Foam::chemPointISAT<CompType, ThermoType>::inEOA(const scalarField& phiq)
{
    scalarField dphi=phiq-phi();
    bool isMechRedActive = chemistry_.mechRed()->active();
    label dim = (isMechRedActive) ? nActiveSpecies_ : completeSpaceSize()-2;
    scalar epsTemp=0.0;
    List<scalar> propEps(completeSpaceSize(),0.0);

    for (label i=0; i<completeSpaceSize()-2; i++)
    {
        scalar temp(0.0);
        //When mechanism reduction is inactive OR on active species
        //multiply L by dphi to get the distance in the active species direction
        //else (for inactive species), just multiply the diagonal
        //element and dphi
        if
        (
            !(isMechRedActive)
          ||(isMechRedActive && completeToSimplifiedIndex_[i]!=-1)
        )
        {
            label si=(isMechRedActive) ? completeToSimplifiedIndex_[i] : i;
            for (label j=si; j<dim; j++)//LT is upper triangular
            {
                label sj=(isMechRedActive) ? simplifiedToCompleteIndex_[j] : j;
                temp += LT_[si][j]*dphi[sj];
            }
            temp += LT_[si][nActiveSpecies_]*dphi[completeSpaceSize()-2];
            temp += LT_[si][nActiveSpecies_+1]*dphi[completeSpaceSize()-1];
        }
        else
        {
            temp = dphi[i]/(tolerance_*scaleFactor_[i]);
        }
        epsTemp += sqr(temp);
        if (printProportion_)
        {
            propEps[i] = temp;
        }
    }
    //Temperature
    epsTemp +=
        sqr
        (
            LT_[dim][dim]*dphi[completeSpaceSize()-2]
           +LT_[dim][dim+1]*dphi[completeSpaceSize()-1]
        );
    //Pressure
    epsTemp += sqr(LT_[dim+1][dim+1]*dphi[completeSpaceSize()-1]);

    if (printProportion_)
    {
        propEps[completeSpaceSize()-2] =
        sqr
        (
            LT_[dim][dim]*dphi[completeSpaceSize()-2]
           +LT_[dim][dim+1]*dphi[completeSpaceSize()-1]
        );

        propEps[completeSpaceSize()-1] = sqr(LT_[dim+1][dim+1]*dphi[completeSpaceSize()-1]);
    }
    if (sqrt(epsTemp) > 1.0+tolerance_)
    {
        if (printProportion_)
        {
            scalar max=-1.0;
            label maxIndex=-1;
            for (label i=0; i<completeSpaceSize(); i++)
            {
                if(max < propEps[i])
                {
                    max = propEps[i];
                    maxIndex = i;
                }
            }
            word propName;
            if (maxIndex >= completeSpaceSize()-2)
            {
                if(maxIndex == completeSpaceSize()-2)
                {
                    propName = "T";
                }
                else if(maxIndex == completeSpaceSize()-1)
                {
                    propName = "p";
                }
            }
            else
            {
                propName = chemistry_.Y()[maxIndex].name();
            }
            Info<< "Direction maximum impact to error in ellipsoid: "
                << propName << endl;
            Info<< "Proportion to the total error on the retrieve: "
                << max / (epsTemp+SMALL) << endl;
        }
        return false;
    }
    else
    {
        return true;
    }
}


template<class CompType, class ThermoType>
bool Foam::chemPointISAT<CompType, ThermoType>::checkSolution
(
    const scalarField& phiq,
    const scalarField& Rphiq
)
{
    scalar eps2 = 0.0;
    scalarField dR = Rphiq - Rphi();
    scalarField dphi = phiq - phi();
    const scalarField& scaleFactorV = scaleFactor();
    const scalarRectangularMatrix& Avar = A();
    bool isMechRedActive = chemistry_.mechRed()->active();
    scalar dRl = 0.0;
    label dim = completeSpaceSize()-2;
    if (isMechRedActive)
    {
        dim = nActiveSpecies_;
    }

    //Since we build only the solution for the species, T and p are not included
    for (label i=0; i<completeSpaceSize()-2; i++)
    {
        dRl = 0.0;
        if (isMechRedActive)
        {
            label si = completeToSimplifiedIndex_[i];
            //If this species is active
            if (si!=-1)
            {
                for (label j=0; j<dim; j++)
                {
                    label sj=simplifiedToCompleteIndex_[j];
                    dRl += Avar[si][j]*dphi[sj];
                }
                dRl += Avar[si][nActiveSpecies_]*dphi[completeSpaceSize()-2];
                dRl += Avar[si][nActiveSpecies_+1]*dphi[completeSpaceSize()-1];
            }
            else
            {
                dRl = dphi[i];
            }
        }
        else
        {
            for (label j=0; j<completeSpaceSize(); j++)
            {
                dRl += Avar[i][j]*dphi[j];
            }
        }
        eps2 += sqr((dR[i]-dRl)/scaleFactorV[i]);
    }

    eps2 = sqrt(eps2);
    if (eps2 > tolerance())
    {
        return false;
    }
    else
    {
        // if the solution is in the ellipsoid of accuracy
        return true;
    }
}


template<class CompType, class ThermoType>
bool Foam::chemPointISAT<CompType, ThermoType>::grow(const scalarField& phiq)
{
    scalarField dphi = phiq - phi();
    label dim = completeSpaceSize();
    label initNActiveSpecies(nActiveSpecies_);
    bool isMechRedActive = chemistry_.mechRed()->active();

    if (isMechRedActive)
    {
        label activeAdded(0);
        DynamicList<label> dimToAdd(0);

        //check if the difference of active species is lower than the maximum
        //number of new dimensions allowed
        for (label i=0; i<completeSpaceSize()-2; i++)
        {
            //first test if the current chemPoint has an inactive species
            //corresponding to an active one in the query point
            if
            (
                completeToSimplifiedIndex_[i]==-1
             && chemistry_.completeToSimplifiedIndex()[i]!=-1
            )
            {
                activeAdded++;
                dimToAdd.append(i);
            }
            //then test if an active species in the current chemPoint
            //corresponds to an inactive on the query side
            if
            (
                completeToSimplifiedIndex_[i]!=-1
             && chemistry_.completeToSimplifiedIndex()[i]==-1
            )
            {
                activeAdded++;
                //we don't need to add a new dimension but we count it to have
                //control on the difference through maxNumNewDim
            }
            //finally test if both points have inactive species but
            //with a dphi!=0
            if
            (
                completeToSimplifiedIndex_[i]==-1
             && chemistry_.completeToSimplifiedIndex()[i]==-1
             && dphi[i] != 0.0
            )
            {
                activeAdded++;
                dimToAdd.append(i);
            }
        }

        //if the number of added dimension is too large, growth fail
        if (activeAdded > maxNumNewDim_)
        {
            return false;
        }

        //the number of added dimension to the current chemPoint
        nActiveSpecies_ += dimToAdd.size();
        simplifiedToCompleteIndex_.setSize(nActiveSpecies_);
        forAll(dimToAdd,i)
        {
            label si = nActiveSpecies_ - dimToAdd.size() + i;
            //add the new active species
            simplifiedToCompleteIndex_[si] = dimToAdd[i];
            completeToSimplifiedIndex_[dimToAdd[i]] = si;
        }

        //update LT and A :
        //-add new column and line for the new active species
        //-transfer last two lines of the previous matrix (p and T) to the end
        //  (change the diagonal position)
        //-set all element of the new lines and columns to zero except diagonal
        //  (=1/(tolerance*scaleFactor))
        if (nActiveSpecies_ > initNActiveSpecies)
        {
            scalarRectangularMatrix LTvar = LT_; //take a copy of LT_
            scalarRectangularMatrix Avar = A_; //take a copy of A_
            LT_ = scalarRectangularMatrix
                (
                    nActiveSpecies_+2,nActiveSpecies_+2,0.0
                );
            A_ = scalarRectangularMatrix
                (
                    nActiveSpecies_+2,nActiveSpecies_+2,0.0
                );

            //write the initial active species
            for (label i=0; i<initNActiveSpecies; i++)
            {
                for (label j=0; j<initNActiveSpecies; j++)
                {
                    LT_[i][j] = LTvar[i][j];
                    A_[i][j] = Avar[i][j];
                }
            }

            //write the columns for temperature and pressure
            for (label i=0; i<initNActiveSpecies; i++)
            {
                for (label j=1; j>=0; j--)
                {
                    LT_[i][nActiveSpecies_+j]=LTvar[i][initNActiveSpecies+j];
                    A_[i][nActiveSpecies_+j]=Avar[i][initNActiveSpecies+j];
                    LT_[nActiveSpecies_+j][i]=LTvar[initNActiveSpecies+j][i];
                    A_[nActiveSpecies_+j][i]=Avar[initNActiveSpecies+j][i];
                }
            }
            //end with the diagonal elements for temperature and pressure
            LT_[nActiveSpecies_][nActiveSpecies_]=
                LTvar[initNActiveSpecies][initNActiveSpecies];
            A_[nActiveSpecies_][nActiveSpecies_]=
                Avar[initNActiveSpecies][initNActiveSpecies];
            LT_[nActiveSpecies_+1][nActiveSpecies_+1]=
                LTvar[initNActiveSpecies+1][initNActiveSpecies+1];
            A_[nActiveSpecies_+1][nActiveSpecies_+1]=
                Avar[initNActiveSpecies+1][initNActiveSpecies+1];

            for (label i=initNActiveSpecies; i<nActiveSpecies_;i++)
            {
                LT_[i][i]=
                    1.0
                  / (tolerance_*scaleFactor_[simplifiedToCompleteIndex_[i]]);
                A_[i][i]=1.0;
            }
        }//end if (nActiveSpecies_>initNActiveSpecies)

        dim = nActiveSpecies_+2;
    }//end if (isMechRedActive)
    //beginning of grow algorithm
    scalarField phiTilde(dim, 0.0);
    scalar normPhiTilde = 0.0;
    //p' = L^T.(p-phi)

    for (label i=0; i<dim; i++)
    {
        for (label j=i; j<dim-2; j++)//LT is upper triangular
        {
            label sj = j;
            if (isMechRedActive)
            {
                sj=simplifiedToCompleteIndex_[j];
            }
            phiTilde[i] += LT_[i][j]*dphi[sj];
        }
        phiTilde[i] += LT_[i][dim-2]*dphi[completeSpaceSize()-2];
        phiTilde[i] += LT_[i][dim-1]*dphi[completeSpaceSize()-1];
        normPhiTilde += sqr(phiTilde[i]);
    }
    scalar invSqrNormPhiTilde = 1.0/normPhiTilde;
    normPhiTilde = sqrt(normPhiTilde);
    //gamma = (1/|p'| - 1)/|p'|^2
    scalar gamma = (1/normPhiTilde - 1)*invSqrNormPhiTilde;
    scalarField u(gamma*phiTilde);
    scalarField v(dim,0.0);
    for ( label i=0; i<dim; i++)
    {
        for (register label j=0; j<=i;j++)
        {
            v[i] += phiTilde[j]*LT_[j][i];
        }
    }

    qrUpdate(LT_,dim, u, v);
    nGrowth_++;
    
    return true;
}
