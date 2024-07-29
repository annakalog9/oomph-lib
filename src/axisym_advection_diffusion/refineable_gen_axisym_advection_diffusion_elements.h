// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
// Header file for elements that solve the advection diffusion equation
// and that can be refined.

#ifndef OOMPH_REFINEABLE_GEN_AXISYM_ADVECTION_DIFFUSION_ELEMENTS_HEADER
#define OOMPH_REFINEABLE_GEN_AXISYM_ADVECTION_DIFFUSION_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// oomph-lib headers
#include "../generic/refineable_quad_element.h"
#include "../generic/refineable_brick_element.h"
#include "../generic/error_estimator.h"
#include "gen_axisym_advection_diffusion_elements.h"

namespace oomph
{
  //======================================================================
  /// A version of the GeneralisedAxisymAdvectionDiffusion
  /// equations that can be
  /// used with non-uniform mesh refinement. In essence, the class overloads
  /// the fill_in_generic_residual_contribution_cons_axisym_adv_diff()
  /// function so that contributions
  /// from hanging nodes (or alternatively in-compatible function values)
  /// are taken into account.
  //======================================================================
  class RefineableGeneralisedAxisymAdvectionDiffusionEquations
    : public virtual GeneralisedAxisymAdvectionDiffusionEquations,
      public virtual RefineableElement,
      public virtual ElementWithZ2ErrorEstimator
  {
  public:
    /// Empty Constructor
    RefineableGeneralisedAxisymAdvectionDiffusionEquations()
      : GeneralisedAxisymAdvectionDiffusionEquations(),
        RefineableElement(),
        ElementWithZ2ErrorEstimator()
    {
    }


    /// Broken copy constructor
    RefineableGeneralisedAxisymAdvectionDiffusionEquations(
      const RefineableGeneralisedAxisymAdvectionDiffusionEquations& dummy) =
      delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const
     RefineableGeneralisedAxisymAdvectionDiffusionEquations&) = delete;*/

    /// Number of 'flux' terms for Z2 error estimation
    unsigned num_Z2_flux_terms()
    {
      return 2;
    }

    /// Get 'flux' for Z2 error recovery:
    /// Standard flux.from GeneralisedAxisymAdvectionDiffusion equations
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
      this->get_flux(s, flux);
    }


    /// Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      // Set size of Vector: u
      values.resize(1);

      // Find number of nodes
      const unsigned n_node = nnode();

      // Find the index at which the unknown is stored
      const unsigned u_nodal_index = this->u_index_cons_axisym_adv_diff();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise value of u
      values[0] = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        values[0] += this->nodal_value(l, u_nodal_index) * psi[l];
      }
    }

    /// Get the function value u in Vector.
    /// Note: Given the generality of the interface (this function
    /// is usually called from black-box documentation or interpolation
    /// routines), the values Vector sets its own size in here.
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      // Set size of Vector:
      values.resize(1);

      // Find out how many nodes there are
      const unsigned n_node = nnode();

      // Find the nodal index at which the unknown is stored
      const unsigned u_nodal_index = this->u_index_cons_axisym_adv_diff();

      // Shape functions
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise the value of u
      values[0] = 0.0;

      // Calculate value
      for (unsigned l = 0; l < n_node; l++)
      {
        values[0] += this->nodal_value(t, l, u_nodal_index) * psi[l];
      }
    }

    /// Fill in the geometric Jacobian, which in this case is r
    double geometric_jacobian(const Vector<double>& x)
    {
      return x[0];
    }


    ///  Further build: Copy source function pointer from father element
    void further_build()
    {
      RefineableGeneralisedAxisymAdvectionDiffusionEquations*
        cast_father_element_pt =
          dynamic_cast<RefineableGeneralisedAxisymAdvectionDiffusionEquations*>(
            this->father_element_pt());

      // Set the values of the pointers from the father
      this->Source_fct_pt = cast_father_element_pt->source_fct_pt();
      this->Wind_fct_pt = cast_father_element_pt->wind_fct_pt();
      this->Conserved_wind_fct_pt =
        cast_father_element_pt->conserved_wind_fct_pt();
      this->Diff_fct_pt = cast_father_element_pt->diff_fct_pt();
      this->Pe_pt = cast_father_element_pt->pe_pt();
      this->PeSt_pt = cast_father_element_pt->pe_st_pt();

      // Set the ALE status
      this->ALE_is_disabled = cast_father_element_pt->ALE_is_disabled;
    }

  protected:
    /// Add the element's contribution to the elemental residual vector
    /// and/or Jacobian matrix
    /// flag=1: compute both
    /// flag=0: compute only residual vector
    void fill_in_generic_residual_contribution_cons_axisym_adv_diff(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag);
  };


  //======================================================================
  /// Refineable version of QGeneralisedAxisymAdvectionDiffusionElement.
  /// Inherit from the standard QGeneralisedAxisymAdvectionDiffusionElement
  /// and the
  /// appropriate refineable geometric element and the refineable equations.
  //======================================================================
  template<unsigned NNODE_1D>
  class RefineableQGeneralisedAxisymAdvectionDiffusionElement
    : public QGeneralisedAxisymAdvectionDiffusionElement<NNODE_1D>,
      public virtual RefineableGeneralisedAxisymAdvectionDiffusionEquations,
      public virtual RefineableQElement<2>
  {
  public:
    /// Empty Constructor:
    RefineableQGeneralisedAxisymAdvectionDiffusionElement()
      : RefineableElement(),
        RefineableGeneralisedAxisymAdvectionDiffusionEquations(),
        RefineableQElement<2>(),
        QGeneralisedAxisymAdvectionDiffusionElement<NNODE_1D>()
    {
    }


    /// Broken copy constructor
    RefineableQGeneralisedAxisymAdvectionDiffusionElement(
      const RefineableQGeneralisedAxisymAdvectionDiffusionElement<NNODE_1D>&
        dummy) = delete;

    /// Broken assignment operator
    /*void operator=(const
     RefineableQGeneralisedAxisymAdvectionDiffusionElement< NNODE_1D>&) =
     delete;*/

    /// Number of continuously interpolated values: 1
    unsigned ncont_interpolated_values() const
    {
      return 1;
    }

    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return QGeneralisedAxisymAdvectionDiffusionElement<
        NNODE_1D>::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return QGeneralisedAxisymAdvectionDiffusionElement<
        NNODE_1D>::vertex_node_pt(j);
    }

    /// Rebuild from sons: empty
    void rebuild_from_sons(Mesh*& mesh_pt) {}

    /// Order of recovery shape functions for Z2 error estimation:
    /// Same order as shape functions.
    unsigned nrecovery_order()
    {
      return (NNODE_1D - 1);
    }

    ///  Perform additional hanging node procedures for variables
    /// that are not interpolated by all nodes. Empty.
    void further_setup_hanging_nodes() {}
  };

  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the
  /// RefineableQuadGeneralisedAxisymAdvectionDiffusionElement elements:
  /// The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<
    RefineableQGeneralisedAxisymAdvectionDiffusionElement<NNODE_1D>>
    : public virtual QElement<1, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : QElement<1, NNODE_1D>() {}
  };

} // namespace oomph

#endif
