// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Header file for Spectral Poisson elements
#ifndef OOMPH_SPECTRAL_POISSON_ELEMENTS_HEADER
#define OOMPH_SPECTRAL_POISSON_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// OOMPH-LIB headers
#include "poisson_elements.h"
#include "../generic/Qspectral_elements.h"

namespace oomph
{
  //======================================================================
  /// QSpectralPoissonElement elements are linear/quadrilateral/brick-shaped
  /// Poisson elements with isoparametric spectral interpolation for the
  /// function. Note that the implementation is PoissonEquations<DIM> does
  /// not use sum factorisation for the evaluation of the residuals and is,
  /// therefore, not optimal for higher dimensions.
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class QSpectralPoissonElement
    : public virtual QSpectralElement<DIM, NNODE_1D>,
      public virtual PoissonEquations<DIM>
  {
  private:
    /// Static array of ints to hold number of variables at
    /// nodes: Initial_Nvalue[n]
    static const unsigned Initial_Nvalue;

  public:
    /// Constructor: Call constructors for QSpectralElement and
    /// Poisson equations
    QSpectralPoissonElement()
      : QSpectralElement<DIM, NNODE_1D>(), PoissonEquations<DIM>()
    {
    }

    /// Broken copy constructor
    QSpectralPoissonElement(
      const QSpectralPoissonElement<DIM, NNODE_1D>& dummy) = delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const QSpectralPoissonElement<DIM,NNODE_1D>&) =
      delete;*/

    ///  Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue;
    }

    /// Output function:
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream& outfile)
    {
      PoissonEquations<DIM>::output(outfile);
    }

    ///  Output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      PoissonEquations<DIM>::output(outfile, n_plot);
    }


    /// C-style output function:
    ///  x,y,u   or    x,y,z,u
    void output(FILE* file_pt)
    {
      PoissonEquations<DIM>::output(file_pt);
    }


    ///  C-style output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      PoissonEquations<DIM>::output(file_pt, n_plot);
    }


    /// Output function for an exact solution:
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      PoissonEquations<DIM>::output_fct(outfile, n_plot, exact_soln_pt);
    }


    /// Output function for a time-dependent exact solution.
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
    /// (Calls the steady version)
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      PoissonEquations<DIM>::output_fct(outfile, n_plot, time, exact_soln_pt);
    }


  protected:
    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_poisson(const Vector<double>& s,
                                                    Shape& psi,
                                                    DShape& dpsidx,
                                                    Shape& test,
                                                    DShape& dtestdx) const;


    /// Shape, test functions & derivs. w.r.t. to global coords. at
    /// integration point ipt. Return Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_poisson(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const;

    /// Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return Jacobian of mapping (J). Also compute
    /// derivatives of dpsidx, dtestdx and J w.r.t. nodal coordinates.
    inline double dshape_and_dtest_eulerian_at_knot_poisson(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      RankFourTensor<double>& d_dpsidx_dX,
      Shape& test,
      DShape& dtestdx,
      RankFourTensor<double>& d_dtestdx_dX,
      DenseMatrix<double>& djacobian_dX) const;
  };


  // Inline functions:


  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QSpectralPoissonElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_poisson(const Vector<double>& s,
                                      Shape& psi,
                                      DShape& dpsidx,
                                      Shape& test,
                                      DShape& dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian(s, psi, dpsidx);

    // Loop over the test functions and derivatives and set them equal to the
    // shape functions
    unsigned nnod = this->nnode();
    for (unsigned i = 0; i < nnod; i++)
    {
      test[i] = psi[i];
      for (unsigned j = 0; j < DIM; j++)
      {
        dtestdx(i, j) = dpsidx(i, j);
      }
    }

    // Return the jacobian
    return J;
  }

  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QSpectralPoissonElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_poisson(const unsigned& ipt,
                                              Shape& psi,
                                              DShape& dpsidx,
                                              Shape& test,
                                              DShape& dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

    // Set the pointers of the test functions
    test = psi;
    dtestdx = dpsidx;

    // Return the jacobian
    return J;
  }

  //======================================================================
  /// Define the shape functions (psi) and test functions (test) and
  /// their derivatives w.r.t. global coordinates (dpsidx and dtestdx)
  /// and return Jacobian of mapping (J). Additionally compute the
  /// derivatives of dpsidx, dtestdx and J w.r.t. nodal coordinates.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QSpectralPoissonElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_poisson(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      RankFourTensor<double>& d_dpsidx_dX,
      Shape& test,
      DShape& dtestdx,
      RankFourTensor<double>& d_dtestdx_dX,
      DenseMatrix<double>& djacobian_dX) const
  {
    // Call the geometrical shape functions and derivatives
    const double J = this->dshape_eulerian_at_knot(
      ipt, psi, dpsidx, djacobian_dX, d_dpsidx_dX);

    // Set the pointers of the test functions
    test = psi;
    dtestdx = dpsidx;
    d_dtestdx_dX = d_dpsidx_dX;

    // Return the jacobian
    return J;
  }

  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the QSpectralPoissonElement elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QSpectralPoissonElement<DIM, NNODE_1D>>
    : public virtual QSpectralElement<DIM - 1, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : QSpectralElement<DIM - 1, NNODE_1D>() {}
  };


  //=======================================================================
  /// Face geometry for the 1D QPoissonElement elements: Point elements
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<QSpectralPoissonElement<1, NNODE_1D>>
    : public virtual PointElement
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : PointElement() {}
  };

} // namespace oomph

#endif
