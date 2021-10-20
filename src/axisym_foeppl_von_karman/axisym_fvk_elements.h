// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
// Header file for axisymmetric FoepplvonKarman elements
#ifndef OOMPH_AXISYM_FOEPPLVONKARMAN_ELEMENTS_HEADER
#define OOMPH_AXISYM_FOEPPLVONKARMAN_ELEMENTS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include "../generic/nodes.h"
#include "../generic/Qelements.h"
#include "../generic/oomph_utilities.h"
#include "../generic/element_with_external_element.h"

namespace oomph
{
  //=============================================================
  /// A class for all isoparametric elements that solve the
  /// axisum Foeppl von Karman equations.
  ///
  /// This contains the generic maths. Shape functions, geometric
  /// mapping etc. must get implemented in derived class.
  //=============================================================
  class AxisymFoepplvonKarmanEquations : public virtual FiniteElement
  {
  public:
    ///  Function pointer to pressure function fct(r,f(r)) --
    /// r is a Vector!
    typedef void (*AxisymFoepplvonKarmanPressureFctPt)(const double& r,
                                                       double& f);

    ///  Constructor (must initialise the Pressure_fct_pt and
    /// Airy_forcing_fct_pt to null). Also set physical parameters to their
    /// default values.
    AxisymFoepplvonKarmanEquations()
      : Pressure_fct_pt(0), Airy_forcing_fct_pt(0)
    {
      // Set all the physical constants to the default value (zero)
      Eta_pt = &Default_Physical_Constant_Value;
      Linear_bending_model = false;
    }

    /// Broken copy constructor
    AxisymFoepplvonKarmanEquations(
      const AxisymFoepplvonKarmanEquations& dummy) = delete;

    /// Broken assignment operator
    void operator=(const AxisymFoepplvonKarmanEquations&) = delete;

    /// FvK parameter
    const double& eta() const
    {
      return *Eta_pt;
    }

    /// Pointer to FvK parameter
    double*& eta_pt()
    {
      return Eta_pt;
    }

    ///  Return the index at which the i-th unknown value
    /// is stored. The default value, i, is appropriate for single-physics
    /// problems. By default, these are:
    /// 0: w
    /// 1: laplacian w
    /// 2: phi
    /// 3: laplacian phi
    /// 4-5: smooth first derivatives
    /// In derived multi-physics elements, this function should be overloaded
    /// to reflect the chosen storage scheme. Note that these equations require
    /// that the unknown is always stored at the same index at each node.
    virtual inline unsigned nodal_index_fvk(const unsigned& i = 0) const
    {
      return i;
    }

    /// Output with default number of plot points
    void output(std::ostream& outfile)
    {
      const unsigned n_plot = 5;
      output(outfile, n_plot);
    }

    ///  Output FE representation of soln: r,w,sigma_r_r,sigma_phi_phi
    /// at n_plot plot points
    void output(std::ostream& outfile, const unsigned& n_plot);

    /// C_style output with default number of plot points
    void output(FILE* file_pt)
    {
      const unsigned n_plot = 5;
      output(file_pt, n_plot);
    }

    ///  C-style output FE representation of soln: r,w at
    /// n_plot plot points
    void output(FILE* file_pt, const unsigned& n_plot);

    /// Output exact soln: r,w_exact at n_plot plot points
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

    ///  Output exact soln: r,w_exact at
    /// n_plot plot points (dummy time-dependent version to
    /// keep intel compiler happy)
    virtual void output_fct(
      std::ostream& outfile,
      const unsigned& n_plot,
      const double& time,
      FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      throw OomphLibError(
        "There is no time-dependent output_fct() for Foeppl von Karman"
        "elements ",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Get error against and norm of exact solution
    void compute_error(std::ostream& outfile,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       double& error,
                       double& norm);


    /// Dummy, time dependent error checker
    void compute_error(std::ostream& outfile,
                       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                       const double& time,
                       double& error,
                       double& norm)
    {
      throw OomphLibError(
        "There is no time-dependent compute_error() for Foeppl von Karman"
        "elements",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Access function: Pointer to pressure function
    AxisymFoepplvonKarmanPressureFctPt& pressure_fct_pt()
    {
      return Pressure_fct_pt;
    }

    /// Access function: Pointer to pressure function. Const version
    AxisymFoepplvonKarmanPressureFctPt pressure_fct_pt() const
    {
      return Pressure_fct_pt;
    }

    /// Access function: Pointer to Airy forcing function
    AxisymFoepplvonKarmanPressureFctPt& airy_forcing_fct_pt()
    {
      return Airy_forcing_fct_pt;
    }

    /// Access function: Pointer to Airy forcing function. Const version
    AxisymFoepplvonKarmanPressureFctPt airy_forcing_fct_pt() const
    {
      return Airy_forcing_fct_pt;
    }

    ///  Get pressure term at (Eulerian) position r. This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the strength of the pressure function might be determined by
    /// another system of equations.
    inline virtual void get_pressure_fvk(const unsigned& ipt,
                                         const double& r,
                                         double& pressure) const
    {
      // If no pressure function has been set, return zero
      if (Pressure_fct_pt == 0)
      {
        pressure = 0.0;
      }
      else
      {
        // Get pressure strength
        (*Pressure_fct_pt)(r, pressure);
      }
    }

    ///  Get Airy forcing term at (Eulerian) position r. This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the strength of the pressure function might be determined by
    /// another system of equations.
    inline virtual void get_airy_forcing_fvk(const unsigned& ipt,
                                             const double& r,
                                             double& airy_forcing) const
    {
      // If no Airy forcing function has been set, return zero
      if (Airy_forcing_fct_pt == 0)
      {
        airy_forcing = 0.0;
      }
      else
      {
        // Get Airy forcing strength
        (*Airy_forcing_fct_pt)(r, airy_forcing);
      }
    }

    /// Get gradient of deflection: gradient[i] = dw/dr_i */
    void get_gradient_of_deflection(const Vector<double>& s,
                                    Vector<double>& gradient) const
    {
      // Find out how many nodes there are in the element
      const unsigned n_node = nnode();

      // Get the index at which the unknown is stored
      const unsigned w_nodal_index = nodal_index_fvk(0);

      // Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidr(n_node, 1);

      // Call the derivatives of the shape and test functions
      dshape_eulerian(s, psi, dpsidr);

      // Initialise to zero
      gradient[0] = 0.0;

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        gradient[0] += this->nodal_value(l, w_nodal_index) * dpsidr(l, 0);
      }
    }

    /// Fill in the residuals with this element's contribution
    void fill_in_contribution_to_residuals(Vector<double>& residuals);


    // hierher Jacobian not yet implemented
    // void fill_in_contribution_to_jacobian(Vector<double> &residuals,
    //                                      DenseMatrix<double> &jacobian);

    ///  Return FE representation of vertical displacement, w_fvk(s)
    /// at local coordinate s
    inline double interpolated_w_fvk(const Vector<double>& s) const
    {
      // Find number of nodes
      const unsigned n_node = nnode();

      // Get the index at which the unknown is stored
      const unsigned w_nodal_index = nodal_index_fvk(0);

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise value of u
      double interpolated_w = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_w += this->nodal_value(l, w_nodal_index) * psi[l];
      }

      return (interpolated_w);
    }

    ///  Compute in-plane stresses. Return boolean to indicate success
    /// (false if attempt to evaluate stresses at zero radius)
    bool interpolated_stress(const Vector<double>& s,
                             double& sigma_r_r,
                             double& sigma_phi_phi);


    ///  Self-test: Return 0 for OK
    unsigned self_test();

    ///  Sets a flag to signify that we are solving the linear,
    /// pure bending equations, and pin all the nodal values that will
    /// not be used in this case
    void use_linear_bending_model()
    {
      // Set the boolean flag
      Linear_bending_model = true;

      // Get the index of the first FvK nodal value
      unsigned first_fvk_nodal_index = nodal_index_fvk();

      // Get the total number of FvK nodal values (assuming they are stored
      // contiguously) at node 0 (it's the same at all nodes anyway)
      unsigned total_fvk_nodal_indices = 6;

      // Get the number of nodes in this element
      unsigned n_node = nnode();

      // Loop over the appropriate nodal indices
      for (unsigned index = first_fvk_nodal_index + 2;
           index < first_fvk_nodal_index + total_fvk_nodal_indices;
           index++)
      {
        // Loop over the nodes in the element
        for (unsigned inod = 0; inod < n_node; inod++)
        {
          // Pin the nodal value at the current index
          node_pt(inod)->pin(index);
        }
      }
    }


  protected:
    ///  Shape/test functions and derivs w.r.t. to global coords at
    /// local coord. s; return  Jacobian of mapping
    virtual double dshape_and_dtest_eulerian_axisym_fvk(
      const Vector<double>& s,
      Shape& psi,
      DShape& dpsidr,
      Shape& test,
      DShape& dtestdr) const = 0;


    ///  Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return  Jacobian of mapping
    virtual double dshape_and_dtest_eulerian_at_knot_axisym_fvk(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidr,
      Shape& test,
      DShape& dtestdr) const = 0;

    /// Pointer to FvK parameter
    double* Eta_pt;

    /// Pointer to pressure function:
    AxisymFoepplvonKarmanPressureFctPt Pressure_fct_pt;

    /// Pointer to Airy forcing function
    AxisymFoepplvonKarmanPressureFctPt Airy_forcing_fct_pt;

    /// Default value for physical constants
    static double Default_Physical_Constant_Value;

    ///  Flag which stores whether we are using a linear,
    /// pure bending model instead of the full non-linear Foeppl-von Karman
    bool Linear_bending_model;
  };


  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// Axisym FoepplvonKarmanElement elements are 1D
  /// Foeppl von Karman elements with isoparametric interpolation for the
  /// function.
  //======================================================================
  template<unsigned NNODE_1D>
  class AxisymFoepplvonKarmanElement
    : public virtual QElement<1, NNODE_1D>,
      public virtual AxisymFoepplvonKarmanEquations
  {
  private:
    ///  Static int that holds the number of variables at
    /// nodes: always the same
    static const unsigned Initial_Nvalue;

  public:
    ///  Constructor: Call constructors for QElement and
    /// AxisymFoepplvonKarmanEquations
    AxisymFoepplvonKarmanElement()
      : QElement<1, NNODE_1D>(), AxisymFoepplvonKarmanEquations()
    {
    }

    /// Broken copy constructor
    AxisymFoepplvonKarmanElement(
      const AxisymFoepplvonKarmanElement<NNODE_1D>& dummy) = delete;

    /// Broken assignment operator
    void operator=(const AxisymFoepplvonKarmanElement<NNODE_1D>&) = delete;

    ///   Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue;
    }


    ///  Output function:
    ///  r,w,sigma_r_r,sigma_phi_phi
    void output(std::ostream& outfile)
    {
      AxisymFoepplvonKarmanEquations::output(outfile);
    }

    ///   Output function:
    ///   r,w,sigma_r_r,sigma_phi_phi at n_plot plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      AxisymFoepplvonKarmanEquations::output(outfile, n_plot);
    }

    ///  C-style output function:
    ///  r,w
    void output(FILE* file_pt)
    {
      AxisymFoepplvonKarmanEquations::output(file_pt);
    }

    ///   C-style output function:
    ///   r,w at n_plot plot points
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      AxisymFoepplvonKarmanEquations::output(file_pt, n_plot);
    }

    ///  Output function for an exact solution:
    ///  r,w_exact at n_plot plot points
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      AxisymFoepplvonKarmanEquations::output_fct(
        outfile, n_plot, exact_soln_pt);
    }

    ///  Output function for a time-dependent exact solution.
    ///  r,w_exact at n_plot plot points
    /// (Calls the steady version)
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      AxisymFoepplvonKarmanEquations::output_fct(
        outfile, n_plot, time, exact_soln_pt);
    }


  protected:
    ///  Shape, test functions & derivs. w.r.t. to global coords.
    /// Return Jacobian.
    inline double dshape_and_dtest_eulerian_axisym_fvk(const Vector<double>& s,
                                                       Shape& psi,
                                                       DShape& dpsidr,
                                                       Shape& test,
                                                       DShape& dtestdr) const;

    ///  Shape, test functions & derivs. w.r.t. to global coords. at
    /// integration point ipt. Return Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_axisym_fvk(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidr,
      Shape& test,
      DShape& dtestdr) const;
  };


  // Inline functions:

  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned NNODE_1D>
  double AxisymFoepplvonKarmanElement<
    NNODE_1D>::dshape_and_dtest_eulerian_axisym_fvk(const Vector<double>& s,
                                                    Shape& psi,
                                                    DShape& dpsidr,
                                                    Shape& test,
                                                    DShape& dtestdr) const

  {
    // Call the geometrical shape functions and derivatives
    const double J = this->dshape_eulerian(s, psi, dpsidr);

    // Set the test functions equal to the shape functions
    test = psi;
    dtestdr = dpsidr;

    // Return the jacobian
    return J;
  }


  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned NNODE_1D>
  double AxisymFoepplvonKarmanElement<NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_axisym_fvk(const unsigned& ipt,
                                                 Shape& psi,
                                                 DShape& dpsidr,
                                                 Shape& test,
                                                 DShape& dtestdr) const
  {
    // Call the geometrical shape functions and derivatives
    const double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidr);

    // Set the pointers of the test functions
    test = psi;
    dtestdr = dpsidr;

    // Return the jacobian
    return J;
  }


  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// FSI Axisym FoepplvonKarmanElement elements are 1D
  /// Foeppl von Karman elements with isoparametric interpolation for the
  /// function. Gets traction from adjacent fluid element(s) of type
  /// FLUID_ELEMENT.
  //======================================================================
  template<unsigned NNODE_1D, class FLUID_ELEMENT>
  class FSIAxisymFoepplvonKarmanElement
    : public virtual AxisymFoepplvonKarmanElement<NNODE_1D>,
      public virtual ElementWithExternalElement
  {
  public:
    /// Constructor
    FSIAxisymFoepplvonKarmanElement()
      : Q_pt(&AxisymFoepplvonKarmanEquations::Default_Physical_Constant_Value)
    {
      // Set source element storage: one interaction with an external
      // element -- the fluid bulk element that provides the pressure
      this->set_ninteraction(1);
    }

    /// Empty virtual destructor
    virtual ~FSIAxisymFoepplvonKarmanElement() {}

    ///  Return the ratio of the stress scales used to non-dimensionalise
    /// the fluid and elasticity equations.
    const double& q() const
    {
      return *Q_pt;
    }

    ///  Return a pointer the ratio of stress scales used to
    /// non-dimensionalise the fluid and solid equations.
    double*& q_pt()
    {
      return Q_pt;
    }

    ///  How many items of Data does the shape of the object depend on?
    /// All nodal data
    virtual unsigned ngeom_data() const
    {
      return this->nnode();
    }

    ///  Return pointer to the j-th Data item that the object's
    /// shape depends on.
    virtual Data* geom_data_pt(const unsigned& j)
    {
      return this->node_pt(j);
    }

    ///  Overloaded position function: Return 2D position vector:
    /// (r(zeta),z(zeta)) of material point whose "Lagrangian coordinate"
    /// is given by zeta. Here r=zeta!
    void position(const Vector<double>& zeta, Vector<double>& r) const
    {
      const unsigned t = 0;
      this->position(t, zeta, r);
    }

    ///  Overloaded position function: Return 2D position vector:
    /// (r(zeta),z(zeta)) of material point whose "Lagrangian coordinate"
    /// is given by zeta.
    void position(const unsigned& t,
                  const Vector<double>& zeta,
                  Vector<double>& r) const
    {
      // Find number of nodes
      const unsigned n_node = this->nnode();

      // Get the index at which the poisson unknown is stored
      const unsigned w_nodal_index = this->nodal_index_fvk(0);

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      this->shape(zeta, psi);

      // Initialise
      double interpolated_w = 0.0;
      double interpolated_r = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_w += this->nodal_value(t, l, w_nodal_index) * psi[l];
        interpolated_r += this->node_pt(l)->x(t, 0) * psi[l];
      }

      // Assign
      r[0] = interpolated_r;
      r[1] = interpolated_w;
    }

    ///  j-th time-derivative on object at current time:
    /// \f$ \frac{d^{j} r(\zeta)}{dt^j} \f$.
    void dposition_dt(const Vector<double>& zeta,
                      const unsigned& j,
                      Vector<double>& drdt)
    {
      // Find number of nodes
      const unsigned n_node = this->nnode();

      // Get the index at which the poisson unknown is stored
      const unsigned w_nodal_index = this->nodal_index_fvk(0);

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      this->shape(zeta, psi);

      // Initialise
      double interpolated_dwdt = 0.0;
      double interpolated_drdt = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the timestepper
        const TimeStepper* time_stepper_pt = node_pt(l)->time_stepper_pt();

        // If we are doing an unsteady solve then calculate the derivative
        if (!time_stepper_pt->is_steady())
        {
          // Get the number of values required to represent history
          const unsigned n_time = time_stepper_pt->ntstorage();

          // Loop over history values
          for (unsigned t = 0; t < n_time; t++)
          {
            // Add the contribution to the derivative
            interpolated_dwdt += time_stepper_pt->weight(1, t) *
                                 this->nodal_value(t, l, w_nodal_index) *
                                 psi[l];
          }
        }
      }

      // Assign
      drdt[0] = interpolated_drdt;
      drdt[1] = interpolated_dwdt;
    }


    ///  Overload pressure term at (Eulerian) position r.
    /// Adds fluid traction to pressure imposed by "pressure fct pointer"
    /// (which can be regarded as applying an external (i.e.
    /// "on the other side" of the fluid) pressure
    inline virtual void get_pressure_fvk(const unsigned& ipt,
                                         const double& r,
                                         double& pressure) const
    {
      pressure = 0.0;

      // Get underlying version
      AxisymFoepplvonKarmanEquations::get_pressure_fvk(ipt, r, pressure);

      // Get FSI parameter
      const double q_value = q();

      // Get fluid element
      FLUID_ELEMENT* ext_el_pt =
        dynamic_cast<FLUID_ELEMENT*>(external_element_pt(0, ipt));
      Vector<double> s_ext(external_element_local_coord(0, ipt));

      // Outer unit normal is vertically upward (in z direction)
      // (within an essentiall flat) model for the wall)
      Vector<double> normal(2);
      normal[0] = 0.0;
      normal[1] = 1.0;

      // Get traction
      Vector<double> traction(3);
      ext_el_pt->traction(s_ext, normal, traction);

      // Add z-component of traction
      pressure -= q_value * traction[1];
    }


    /// Output integration points (for checking of fsi setup)
    void output_integration_points(std::ostream& outfile)
    {
      // How many nodes do we have?
      unsigned nnod = this->nnode();
      Shape psi(nnod);

      // Get the index at which the unknown is stored
      const unsigned w_nodal_index = this->nodal_index_fvk(0);

      // Loop over the integration points
      const unsigned n_intpt = this->integral_pt()->nweight();
      outfile << "ZONE I=" << n_intpt << std::endl;
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get shape fct
        Vector<double> s(1);
        s[0] = this->integral_pt()->knot(ipt, 0);
        shape(s, psi);

        // Initialise
        double interpolated_w = 0.0;
        double interpolated_r = 0.0;

        // Loop over the local nodes and sum
        for (unsigned l = 0; l < nnod; l++)
        {
          interpolated_w += this->nodal_value(l, w_nodal_index) * psi[l];
          interpolated_r += this->node_pt(l)->x(0) * psi[l];
        }

        // Get fluid element
        FLUID_ELEMENT* ext_el_pt =
          dynamic_cast<FLUID_ELEMENT*>(external_element_pt(0, ipt));
        Vector<double> s_ext(external_element_local_coord(0, ipt));

        // Get veloc
        Vector<double> veloc(3);
        ext_el_pt->interpolated_u_axi_nst(s_ext, veloc);
        Vector<double> x(2);
        ext_el_pt->interpolated_x(s_ext, x);

        outfile << interpolated_r << " " << interpolated_w << " " << veloc[0]
                << " " << veloc[1] << " " << x[0] << " " << x[1] << " "
                << std::endl;
      }
    }


    /// Output adjacent fluid elements (for checking of fsi setup)
    void output_adjacent_fluid_elements(std::ostream& outfile,
                                        const unsigned& nplot)
    {
      // Loop over the integration points
      const unsigned n_intpt = this->integral_pt()->nweight();
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get fluid element
        FLUID_ELEMENT* ext_el_pt =
          dynamic_cast<FLUID_ELEMENT*>(external_element_pt(0, ipt));

        // Dump it
        ext_el_pt->output(outfile, nplot);
      }
    }

    ///  Perform any auxiliary node update fcts of the adjacent
    /// fluid elements
    void update_before_external_interaction_geometric_fd()
    {
      node_update_adjacent_fluid_elements();
    }

    ///  Perform any auxiliary node update fcts of the adjacent
    /// fluid elements
    void reset_after_external_interaction_geometric_fd()
    {
      node_update_adjacent_fluid_elements();
    }


    ///  Perform any auxiliary node update fcts of the adjacent
    /// fluid elements
    void update_in_external_interaction_geometric_fd(const unsigned& i)
    {
      node_update_adjacent_fluid_elements();
    }


    ///  Perform any auxiliary node update fcts of the adjacent
    /// fluid elements
    void reset_in_external_interaction_geometric_fd(const unsigned& i)
    {
      node_update_adjacent_fluid_elements();
    }


    ///  Update the nodal positions in all fluid elements that affect
    /// the traction on this element
    void node_update_adjacent_fluid_elements()
    {
      // Don't update elements repeatedly
      std::map<FLUID_ELEMENT*, bool> done;

      // Number of integration points
      unsigned n_intpt = integral_pt()->nweight();

      // Loop over all integration points in wall element
      for (unsigned iint = 0; iint < n_intpt; iint++)
      {
        // Get fluid element that affects this integration point
        FLUID_ELEMENT* el_f_pt =
          dynamic_cast<FLUID_ELEMENT*>(external_element_pt(0, iint));

        // Is there an adjacent fluid element?
        if (el_f_pt != 0)
        {
          // Have we updated its positions yet?
          if (!done[el_f_pt])
          {
            // Update nodal positions
            el_f_pt->node_update();
            done[el_f_pt] = true;
          }
        }
      }
    }


    ///  Output FE representation of soln:
    /// r,w,dwdt,sigma_r_r,sigma_phi_phi at n_plot plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      // Vector of local coordinates
      Vector<double> s(1);

      // Tecplot header info
      outfile << "ZONE\n";

      // Loop over plot points
      unsigned num_plot_points = nplot_points(n_plot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        get_s_plot(iplot, n_plot, s);

        // Get velocity
        Vector<double> drdt(2);
        dposition_dt(s, 1, drdt);

        // Get stress
        double sigma_r_r = 0.0;
        double sigma_phi_phi = 0.0;
        bool success = this->interpolated_stress(s, sigma_r_r, sigma_phi_phi);
        if (success)
        {
          outfile << this->interpolated_x(s, 0) << " "
                  << this->interpolated_w_fvk(s) << " " << drdt[0] << " "
                  << drdt[1] << " " << sigma_r_r << " " << sigma_phi_phi
                  << std::endl;
        }
      }
    }


  protected:
    ///  Pointer to the ratio, \f$ Q \f$ , of the stress used to
    /// non-dimensionalise the fluid stresses to the stress used to
    /// non-dimensionalise the solid stresses.
    double* Q_pt;
  };


} // namespace oomph

#endif
