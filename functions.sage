"""
This file contains the functions for the article *Deterministic computation of
linear structures on general Drinfeld modules over finite fields*, by Antoine
Leudière and Renate Scheidler.
"""

import numpy as np

def anderson_motive_coordinates(phi, f):
    """
    Compute the coordinates of ``f`` relative to the Anderson motive of
    the Drinfeld module ``phi``.

    INPUT:

    - ``phi`` -- a Drinfeld module
    - ``f`` --  an element of the Ore polynomial ring of ``phi``

    OUTPUT: the numpy array of the coordinates
    """
    f = phi.ore_polring()(f)
    r = phi.rank()

    def pad(list_):
        if len(list_) < r:
            return np.array(list_ + [0] * (r - len(list_)))
        return np.array(list_)

    if f.degree() < r:
        return pad(list(f))
    else:
        quo, rem = f.right_quo_rem(phi.gen())
        return (T * anderson_motive_coordinates(phi, quo)) + pad(list(rem))


def order(phi, x):
    """
    Compute the order (or annihilator) of the ``x``, that is, the
    minimal monic polynomial a such that ``x`` is of a-torsion for
    ``phi``.

    INPUT:

    - ``phi`` -- a Drinfeld module
    - ``x`` --  a rational point of ``phi``

    OUTPUT: a polynomial in the function ring Fq[T]
    """
    assert x in K

    d = phi.base_over_constants_field().degree(Fq)
    phi_T = phi.gen()
    # Vector of d+1 elements [1, phi_T(x), ..., phi_T^d(x)]:
    x_phi = [x]
    for i in range(1, d+1):
        x_phi.append(phi_T(x_phi[i-1]))
    # Get the coordinates of these, and create a matrix:
    mat_x = matrix([list(K(xi)) for xi in x_phi]).transpose()
    # Get its kernel:
    ker = mat_x.right_kernel().basis()
    # Get the corresponding polynomials:
    gens = map(lambda vec: sum(coeff * T^i for i, coeff in enumerate(vec)), ker)
    # And return their gcd, which is the order of x:
    return gcd(gens)


def cast_matrix(matrix):
    """
    Cast the input matrix to a matrix with entries in Fq. The use case
    is that the input comes from a Frobenius form that Sage computes
    over the algebraic closure, even though the Frobenius form actually
    has entries with Fq.

    This implementation leaves a lot to be desired, but we believe it is
    good enough for this proof of concept, and we intend to propose a
    better integration for the final integration in SageMath.
    
    INPUT: a matrix with entries in the algebraic closure of Fq

    OUTPUT: a matrix with entries in Fq
    """
    Fqbar = Fq.algebraic_closure()

    def cast_coeff(coeff):
        for x in Fq:
            if coeff == Fqbar(x):
                return x

    return matrix.apply_map(cast_coeff)


def maximal_rational_separable_torsion(phi):
    """
    Return the unique monic polynomial a in Fq[T] that verifies:

      (i) a is away from the characteristic;
      (ii) a is the maximal a' such that the (a')-torsion is rational.

    INPUT: a Drinfeld module over a finite field

    OUTPUT: an element of Fq[T], the function ring of ``phi``
    """
    d = phi.base_over_constants_field().degree(Fq)
    t = phi.ore_variable()
    # Commute the coordinates of tau^d - 1 in the Anderson motives:
    coords = anderson_motive_coordinates(phi, t^d - 1)
    # Compute their gcd, which is an element of K[T]:
    gcd_ = gcd(coords)
    # Now compute the Fq[T]-coordinates of that gcd (recall that K[T]
    # is a free Fq[T]-module with rank d and canonical basis
    # inherited from the fixed Fq-basis of K):
    coeffs_Fq = list(filter(lambda coeff: list(coeff), list(gcd_)))
    gcd_coords = [sum(coeffs_Fq[j][i] * T^j for j in range(gcd_.degree()+1)) 
                  for i in range(d)]
    # Return the gcd of those coefficients:
    return gcd(gcd_coords)


def ore_pol_matrix(f, basis):
    """
    Return the matrix of the Ore polynomial `f` (seen as an
    Fq-endomorphism of the base field) relative to the input `basis`.

    INPUT:

    - ``f`` -- an Ore polynomial
    - ``basis`` -- an Fq-basis of the ambient base field 

    OUTPUT: a matrix with entries in Fq
    """
    return matrix([list(f(x)) for x in basis]).transpose()


def ore_pol_kernel(f, basis):
    """
    Return an Fq-basis of the matrix of the Ore polynomial ``f`` (seen
    as an Fq-linear endomorphism of the base field), relative to the
    given basis.

    INPUT:

    - ``f`` -- an Ore polynomial
    - ``basis`` -- an Fq-basis of the ambient base field 

    OUTPUT: a matrix with entries in Fq
    """
    matrix = ore_pol_matrix(f, basis)
    return matrix.right_kernel().basis()


def morphism_kernel_invariants(u, basis, frobenius=False):
    """
    Return the invariant factors of the kernel of the input morphism
    ``u``, relative to the input ``basis``; optionally, return a
    Frobenius decomposition.

    This function corresponds to Algorithm 4 in the paper.

    INPUT:

    - ``u`` --  a morphism of Drinfeld modules over a finite field
    - ``basis`` -- an Fq- basis of the base field
    - ``frobenius`` (optional) -- if `True`, return a Frobenius decomposition

    OUTPUT:

    - If ``frobenius`` is set to ``False`` (default), then return a list
      of polynomials.
    - Otherwise, return a 2-tuple containing a list of polynomials and a
      list of elements in the ambient base field.
    """
    phi = u.domain()

    # Compute a basis of the kernel of the morphism:
    ker_Mu = ore_pol_kernel(u.ore_polynomial(), basis)

    # Now, we compute the action of phi_T relative to said basis:
    N = matrix(ker_Mu).transpose()
    Y = matrix([list(phi.gen()(K(list((bi))))) for bi in ker_Mu]).transpose()
    X = N.solve_right(Y)  # This is the action of phi_T

    # We can now compute the invariant factors of the morphism kernel:
    invariant_factors = [A(inv) for inv in X.rational_form(format='invariants')]

    # If the user asks for the Frobenius decomposition, there is still
    # work to do. Namely, get the Frobenius decomposition relative to
    # the basis of the kernel, and get the representation of these
    # vectors relative to the original basis.
    #
    # Unfortunately, at the moment, SageMath does not compute the
    # change-of-basis matrix along with the Frobenius normal form.
    # Consequently, we have to use the method ``is_similar`` with the
    # flag ``transformation`` set to ``True``, which is very costly.
    if frobenius:
        # Reconstruct the Frobenius normal form:
        F = block_diagonal_matrix([companion_matrix(A(inv)) for inv in invariant_factors])

        # Compute the change-of-basis matrix:
        _, _S = F.is_similar(X, transformation=True)
        S = cast_matrix(_S)  # Cast the entries of the matrix to Fq

        # Convert the elements of the change-of-basis matrix to elements
        # represented in the original basis:
        frobenius_matrix = N * S

        # Pick the correct columns:
        cols = frobenius_matrix.columns()
        frobenius_decomposition = []
        counter = 0
        for deg in [inv.degree() for inv in invariant_factors]:
            frobenius_decomposition.append(K(list(cols[counter])))
            counter += deg
        return invariant_factors, frobenius_decomposition
    else:
        return invariant_factors


def module_of_points_invariants(phi, basis, frobenius=False):
    """
    Return the invariant factors of module of points of the input
    Drinfeld module ``phi``, relative to the input ``basis``;
    optionally, return a Frobenius decomposition.

    This function corresponds to Algorithm 5 in the paper.

    INPUT:

    - ``phi`` --  a Drinfeld module over a finite field
    - ``basis`` -- an Fq- basis of the base field
    - ``frobenius`` (optional) -- if `True`, return a Frobenius decomposition

    OUTPUT:

    - If ``frobenius`` is set to ``False`` (default), then return a list
      of polynomials.
    - Otherwise, return a 2-tuple containing a list of polynomials and a
      list of elements in the ambient base field.
    """
    return morphism_kernel_invariants(phi.hom(0), basis, frobenius=frobenius)
