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
// Header file for UnsteadyHeat elements
#ifndef OOMPH_UNSTEADY_HEAT_ELEMENTS_HEADER
#define OOMPH_UNSTEADY_HEAT_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// OOMPH-LIB headers
#include "../generic/projection.h"
#include "../generic/nodes.h"
#include "../generic/Qelements.h"
#include "../generic/oomph_utilities.h"


namespace oomph
{
  /// Base class so that we don't need to know the dimension just to set the
  /// source function!
  class UnsteadyHeatEquationsBase : public virtual FiniteElement
  {
  public:
    ///  Function pointer to source function fct(t,x,f(x,t)) --
    /// x is a Vector!
    typedef void (*UnsteadyHeatSourceFctPt)(const double& time,
                                            const Vector<double>& x,
                                            double& u);

    /// Access function: Pointer to source function
    virtual UnsteadyHeatSourceFctPt& source_fct_pt() = 0;
  };

  //=============================================================
  /// A class for all isoparametric elements that solve the
  /// UnsteadyHeat equations.
  /// \f[
  /// \frac{\partial^2 u}{\partial x_i^2}=\frac{\partial u}{\partial t}+f(t,x_j)
  /// \f]
  /// This contains the generic maths. Shape functions, geometric
  /// mapping etc. must get implemented in derived class.
  /// Note that this class assumes an isoparametric formulation, i.e. that
  /// the scalar unknown is interpolated using the same shape funcitons
  /// as the position.
  //=============================================================
  template<unsigned DIM>
  class UnsteadyHeatEquations : public virtual UnsteadyHeatEquationsBase
  {
  public:
    ///  Function pointer to source function fct(t,x,f(x,t)) --
    /// x is a Vector!
    typedef void (*UnsteadyHeatSourceFctPt)(const double& time,
                                            const Vector<double>& x,
                                            double& u);


    ///  Constructor: Initialises the Source_fct_pt to null and
    /// sets flag to use ALE formulation of the equations.
    /// Also set Alpha (thermal inertia) and Beta (thermal conductivity)
    /// parameters to defaults (both one for natural scaling)
    UnsteadyHeatEquations() : Source_fct_pt(0), ALE_is_disabled(false)
    {
      // Set Alpha and Beta parameter to default (one for natural scaling of
      // time)
      Alpha_pt = &Default_alpha_parameter;
      Beta_pt = &Default_beta_parameter;
    }


    /// Broken copy constructor
    UnsteadyHeatEquations(const UnsteadyHeatEquations& dummy) = delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const UnsteadyHeatEquations&) = delete;*/

    ///  Return the index at which the unknown value
    /// is stored. The default value, 0, is appropriate for single-physics
    /// problems, when there is only one variable, the value that satisfies the
    /// unsteady heat equation.
    /// In derived multi-physics elements, this function should be overloaded
    /// to reflect the chosen storage scheme. Note that these equations require
    /// that the unknown is always stored at the same index at each node.
    virtual inline unsigned u_index_ust_heat() const
    {
      return 0;
    }

    ///  du/dt at local node n.
    /// Uses suitably interpolated value for hanging nodes.
    double du_dt_ust_heat(const unsigned& n) const
    {
      // Get the data's timestepper
      TimeStepper* time_stepper_pt = this->node_pt(n)->time_stepper_pt();

      // Initialise dudt
      double dudt = 0.0;

      // Loop over the timesteps, if there is a non Steady timestepper
      if (!time_stepper_pt->is_steady())
      {
        // Find the index at which the variable is stored
        const unsigned u_nodal_index = u_index_ust_heat();

        // Number of timsteps (past & present)
        const unsigned n_time = time_stepper_pt->ntstorage();

        // Add the contributions to the time derivative
        for (unsigned t = 0; t < n_time; t++)
        {
          dudt +=
            time_stepper_pt->weight(1, t) * nodal_value(t, n, u_nodal_index);
        }
      }
      return dudt;
    }

    ///  Disable ALE, i.e. assert the mesh is not moving -- you do this
    /// at your own risk!
    void disable_ALE()
    {
      ALE_is_disabled = true;
    }


    ///  (Re-)enable ALE, i.e. take possible mesh motion into account
    /// when evaluating the time-derivative. Note: By default, ALE is
    /// enabled, at the expense of possibly creating unnecessary work
    /// in problems where the mesh is, in fact, stationary.
    void enable_ALE()
    {
      ALE_is_disabled = false;
    }

    /// Compute norm of fe solution
    void compute_norm(double& norm);

    /// Output with default number of plot points
    void output(std::ostream& outfile)
    {
      unsigned nplot = 5;
      output(outfile, nplot);
    }


    ///  Output FE representation of soln: x,y,u or x,y,z,u at
    /// n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& nplot);

    /// C_style output with default number of plot points
    void output(FILE* file_pt)
    {
      unsigned n_plot = 5;
      output(file_pt, n_plot);
    }


    ///  C-style output FE representation of soln: x,y,u or x,y,z,u at
    /// n_plot^DIM plot points
    void output(FILE* file_pt, const unsigned& n_plot);


    /// Output exact soln: x,y,u_exact or x,y,z,u_exact at nplot^DIM plot points
    void output_fct(std::ostream& outfile,
                    const unsigned& nplot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);


    ///  Output exact soln: x,y,u_exact or x,y,z,u_exact at
    /// nplot^DIM plot points (time-dependent version)
    virtual void output_fct(
      std::ostream& outfile,
      const unsigned& nplot,
      const double& time,
      FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt);


    /// Get error against and norm of exact solution
    void compute_error(std::ostream& outfile,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       double& error,
                       double& norm);


    /// Get error against and norm of exact solution
    void compute_error(std::ostream& outfile,
                       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                       const double& time,
                       double& error,
                       double& norm);


    /// Access function: Pointer to source function
    UnsteadyHeatSourceFctPt& source_fct_pt()
    {
      return Source_fct_pt;
    }


    /// Access function: Pointer to source function. Const version
    UnsteadyHeatSourceFctPt source_fct_pt() const
    {
      return Source_fct_pt;
    }


    ///  Get source term at continous time t and (Eulerian) position x.
    /// Virtual so it can be overloaded in derived multiphysics elements.
    virtual inline void get_source_ust_heat(const double& t,
                                            const unsigned& ipt,
                                            const Vector<double>& x,
                                            double& source) const
    {
      // If no source function has been set, return zero
      if (Source_fct_pt == 0)
      {
        source = 0.0;
      }
      else
      {
        // Get source strength
        (*Source_fct_pt)(t, x, source);
      }
    }

    /// Alpha parameter (thermal inertia)
    const double& alpha() const
    {
      return *Alpha_pt;
    }

    /// Pointer to Alpha parameter (thermal inertia)
    double*& alpha_pt()
    {
      return Alpha_pt;
    }


    /// Beta parameter (thermal conductivity)
    const double& beta() const
    {
      return *Beta_pt;
    }

    /// Pointer to Beta parameter (thermal conductivity)
    double*& beta_pt()
    {
      return Beta_pt;
    }

    /// Get flux: flux[i] = du/dx_i
    void get_flux(const Vector<double>& s, Vector<double>& flux) const
    {
      // Find out how many nodes there are in the element
      unsigned n_node = nnode();

      // Find the index at which the variable is stored
      unsigned u_nodal_index = u_index_ust_heat();

      // Set up memory for the shape and test functions
      Shape psi(n_node);
      DShape dpsidx(n_node, DIM);

      // Call the derivatives of the shape and test functions
      dshape_eulerian(s, psi, dpsidx);

      // Initialise to zero
      for (unsigned j = 0; j < DIM; j++)
      {
        flux[j] = 0.0;
      }

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over derivative directions
        for (unsigned j = 0; j < DIM; j++)
        {
          flux[j] += nodal_value(l, u_nodal_index) * dpsidx(l, j);
        }
      }
    }


    /// Compute element residual Vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      fill_in_generic_residual_contribution_ust_heat(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }


    /// Compute element residual Vector and element Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_ust_heat(residuals, jacobian, 1);
    }


    /// Return FE representation of function value u(s) at local coordinate s
    inline double interpolated_u_ust_heat(const Vector<double>& s) const
    {
      // Find number of nodes
      unsigned n_node = nnode();

      // Find the index at which the variable is stored
      unsigned u_nodal_index = u_index_ust_heat();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise value of u
      double interpolated_u = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += nodal_value(l, u_nodal_index) * psi[l];
      }

      return (interpolated_u);
    }


    ///  Return FE representation of function value u(s) at local
    /// coordinate s at previous time t (t=0: present)
    inline double interpolated_u_ust_heat(const unsigned& t,
                                          const Vector<double>& s) const
    {
      // Find number of nodes
      unsigned n_node = nnode();

      // Find the index at which the variable is stored
      unsigned u_nodal_index = u_index_ust_heat();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise value of u
      double interpolated_u = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += nodal_value(t, l, u_nodal_index) * psi[l];
      }

      return (interpolated_u);
    }


    /// Return FE representation of function value du/dt(s) at local coordinate
    /// s
    inline double interpolated_du_dt_ust_heat(const Vector<double>& s) const
    {
      // Find number of nodes
      unsigned n_node = nnode();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Initialise value of du/dt
      double interpolated_dudt = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_dudt += du_dt_ust_heat(l) * psi[l];
      }

      return (interpolated_dudt);
    }


    ///  Self-test: Return 0 for OK
    unsigned self_test();


  protected:
    ///  Shape/test functions and derivs w.r.t. to global coords at
    /// local coord. s; return  Jacobian of mapping
    virtual double dshape_and_dtest_eulerian_ust_heat(
      const Vector<double>& s,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const = 0;


    ///  Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return  Jacobian of mapping
    virtual double dshape_and_dtest_eulerian_at_knot_ust_heat(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const = 0;

    ///  Compute element residual Vector only (if flag=and/or element
    /// Jacobian matrix
    virtual void fill_in_generic_residual_contribution_ust_heat(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag);

    /// Pointer to source function:
    UnsteadyHeatSourceFctPt Source_fct_pt;

    ///  Boolean flag to indicate if ALE formulation is disabled when
    /// time-derivatives are computed. Only set to true if you're sure
    /// that the mesh is stationary.
    bool ALE_is_disabled;

    /// Pointer to Alpha parameter (thermal inertia)
    double* Alpha_pt;

    /// Pointer to Beta parameter (thermal conductivity)
    double* Beta_pt;

  private:
    ///  Static default value for the Alpha parameter:
    /// (thermal inertia): One for natural scaling
    static double Default_alpha_parameter;

    ///  Static default value for the Beta parameter (thermal
    /// conductivity): One for natural scaling
    static double Default_beta_parameter;
  };


  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  //======================================================================
  /// QUnsteadyHeatElement elements are linear/quadrilateral/brick-shaped
  /// UnsteadyHeat elements with isoparametric interpolation for the function.
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class QUnsteadyHeatElement : public virtual QElement<DIM, NNODE_1D>,
                               public virtual UnsteadyHeatEquations<DIM>
  {
  private:
    ///  Static array of ints to hold number of variables at
    /// nodes: Initial_Nvalue[n]
    static const unsigned Initial_Nvalue;

  public:
    ///  Constructor: Call constructors for QElement and
    /// UnsteadyHeat equations
    QUnsteadyHeatElement()
      : QElement<DIM, NNODE_1D>(), UnsteadyHeatEquations<DIM>()
    {
    }

    /// Broken copy constructor
    QUnsteadyHeatElement(const QUnsteadyHeatElement<DIM, NNODE_1D>& dummy) =
      delete;

    /// Broken assignment operator
    /*void operator=(const QUnsteadyHeatElement<DIM,NNODE_1D>&) = delete;*/

    ///   Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue;
    }

    ///  Output function:
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream& outfile)
    {
      UnsteadyHeatEquations<DIM>::output(outfile);
    }


    ///   Output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      UnsteadyHeatEquations<DIM>::output(outfile, n_plot);
    }


    ///  C-style output function:
    ///  x,y,u   or    x,y,z,u
    void output(FILE* file_pt)
    {
      UnsteadyHeatEquations<DIM>::output(file_pt);
    }


    ///   C-style output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      UnsteadyHeatEquations<DIM>::output(file_pt, n_plot);
    }


    ///  Output function for an exact solution:
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      UnsteadyHeatEquations<DIM>::output_fct(outfile, n_plot, exact_soln_pt);
    }


    ///  Output function for a time-dependent exact solution.
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
    /// (Calls the steady version)
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      UnsteadyHeatEquations<DIM>::output_fct(
        outfile, n_plot, time, exact_soln_pt);
    }


  protected:
    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_ust_heat(const Vector<double>& s,
                                                     Shape& psi,
                                                     DShape& dpsidx,
                                                     Shape& test,
                                                     DShape& dtestdx) const;


    ///  Shape/test functions and derivs w.r.t. to global coords at
    /// integration point ipt; return  Jacobian of mapping
    inline double dshape_and_dtest_eulerian_at_knot_ust_heat(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const;
  };


  // Inline functions:


  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QUnsteadyHeatElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_ust_heat(const Vector<double>& s,
                                       Shape& psi,
                                       DShape& dpsidx,
                                       Shape& test,
                                       DShape& dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian(s, psi, dpsidx);

    // Loop over the test functions and derivatives and set them equal to the
    // shape functions
    for (unsigned i = 0; i < NNODE_1D; i++)
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
  double QUnsteadyHeatElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_ust_heat(const unsigned& ipt,
                                               Shape& psi,
                                               DShape& dpsidx,
                                               Shape& test,
                                               DShape& dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

    // Set the test functions equal to the shape functions
    //(sets internal pointers)
    test = psi;
    dtestdx = dpsidx;

    // Return the jacobian
    return J;
  }


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the QUnsteadyHeatElement elements: The spatial
  /// dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points
  /// along their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QUnsteadyHeatElement<DIM, NNODE_1D>>
    : public virtual QElement<DIM - 1, NNODE_1D>
  {
  public:
    ///  Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : QElement<DIM - 1, NNODE_1D>() {}
  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Face geometry for the 1D QUnsteadyHeatElement elements: Point elements
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<QUnsteadyHeatElement<1, NNODE_1D>>
    : public virtual PointElement
  {
  public:
    ///  Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : PointElement() {}
  };


  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //==========================================================
  /// UnsteadyHeat upgraded to become projectable
  //==========================================================
  template<class UNSTEADY_HEAT_ELEMENT>
  class ProjectableUnsteadyHeatElement
    : public virtual ProjectableElement<UNSTEADY_HEAT_ELEMENT>
  {
  public:
    ///  Constructor [this was only required explicitly
    /// from gcc 4.5.2 onwards...]
    ProjectableUnsteadyHeatElement() {}

    ///  Specify the values associated with field fld.
    /// The information is returned in a vector of pairs which comprise
    /// the Data object and the value within it, that correspond to field fld.
    Vector<std::pair<Data*, unsigned>> data_values_of_field(const unsigned& fld)
    {
#ifdef PARANOID
      if (fld != 0)
      {
        std::stringstream error_stream;
        error_stream << "UnsteadyHeat elements only store a single field so "
                        "fld must be 0 rather"
                     << " than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Create the vector
      unsigned nnod = this->nnode();
      Vector<std::pair<Data*, unsigned>> data_values(nnod);

      // Loop over all nodes
      for (unsigned j = 0; j < nnod; j++)
      {
        // Add the data value associated field: The node itself
        data_values[j] = std::make_pair(this->node_pt(j), fld);
      }

      // Return the vector
      return data_values;
    }

    ///  Number of fields to be projected: Just one
    unsigned nfields_for_projection()
    {
      return 1;
    }

    ///  Number of history values to be stored for fld-th field.
    /// (Note: count includes current value!)
    unsigned nhistory_values_for_projection(const unsigned& fld)
    {
#ifdef PARANOID
      if (fld != 0)
      {
        std::stringstream error_stream;
        error_stream << "UnsteadyHeat elements only store a single field so "
                        "fld must be 0 rather"
                     << " than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return this->node_pt(0)->ntstorage();
    }

    /// Number of positional history values
    /// (Note: count includes current value!)
    unsigned nhistory_values_for_coordinate_projection()
    {
      return this->node_pt(0)->position_time_stepper_pt()->ntstorage();
    }

    ///  Return Jacobian of mapping and shape functions of field fld
    /// at local coordinate s
    double jacobian_and_shape_of_field(const unsigned& fld,
                                       const Vector<double>& s,
                                       Shape& psi)
    {
#ifdef PARANOID
      if (fld != 0)
      {
        std::stringstream error_stream;
        error_stream << "UnsteadyHeat elements only store a single field so "
                        "fld must be 0 rather"
                     << " than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      unsigned n_dim = this->dim();
      unsigned n_node = this->nnode();
      Shape test(n_node);
      DShape dpsidx(n_node, n_dim), dtestdx(n_node, n_dim);
      double J =
        this->dshape_and_dtest_eulerian_ust_heat(s, psi, dpsidx, test, dtestdx);
      return J;
    }


    ///  Return interpolated field fld at local coordinate s, at time
    /// level t (t=0: present; t>0: history values)
    double get_field(const unsigned& t,
                     const unsigned& fld,
                     const Vector<double>& s)
    {
#ifdef PARANOID
      if (fld != 0)
      {
        std::stringstream error_stream;
        error_stream << "UnsteadyHeat elements only store a single field so "
                        "fld must be 0 rather"
                     << " than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // Find the index at which the variable is stored
      unsigned u_nodal_index = this->u_index_ust_heat();

      // Local shape function
      unsigned n_node = this->nnode();
      Shape psi(n_node);

      // Find values of shape function
      this->shape(s, psi);

      // Initialise value of u
      double interpolated_u = 0.0;

      // Sum over the local nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += this->nodal_value(t, l, u_nodal_index) * psi[l];
      }
      return interpolated_u;
    }


    /// Return number of values in field fld: One per node
    unsigned nvalue_of_field(const unsigned& fld)
    {
#ifdef PARANOID
      if (fld != 0)
      {
        std::stringstream error_stream;
        error_stream << "UnsteadyHeat elements only store a single field so "
                        "fld must be 0 rather"
                     << " than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return this->nnode();
    }


    /// Return local equation number of value j in field fld.
    int local_equation(const unsigned& fld, const unsigned& j)
    {
#ifdef PARANOID
      if (fld != 0)
      {
        std::stringstream error_stream;
        error_stream << "UnsteadyHeat elements only store a single field so "
                        "fld must be 0 rather"
                     << " than " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      const unsigned u_nodal_index = this->u_index_ust_heat();
      return this->nodal_local_eqn(j, u_nodal_index);
    }


    ///  Output FE representation of soln: x,y,u or x,y,z,u at
    /// and history values at n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& nplot)
    {
      unsigned el_dim = this->dim();
      // Vector of local coordinates
      Vector<double> s(el_dim);

      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, nplot, s);

        for (unsigned i = 0; i < el_dim; i++)
        {
          outfile << this->interpolated_x(s, i) << " ";
        }
        outfile << this->interpolated_u_ust_heat(s) << " ";
        outfile << this->interpolated_du_dt_ust_heat(s) << " ";


        // History values of coordinates
        unsigned n_prev =
          this->node_pt(0)->position_time_stepper_pt()->ntstorage();
        for (unsigned t = 1; t < n_prev; t++)
        {
          for (unsigned i = 0; i < el_dim; i++)
          {
            outfile << this->interpolated_x(t, s, i) << " ";
          }
        }

        // History values of velocities
        n_prev = this->node_pt(0)->time_stepper_pt()->ntstorage();
        for (unsigned t = 1; t < n_prev; t++)
        {
          outfile << this->interpolated_u_ust_heat(t, s) << " ";
        }
        outfile << std::endl;
      }


      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);
    }
  };


  //=======================================================================
  /// Face geometry for element is the same as that for the underlying
  /// wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<ProjectableUnsteadyHeatElement<ELEMENT>>
    : public virtual FaceGeometry<ELEMENT>
  {
  public:
    FaceGeometry() : FaceGeometry<ELEMENT>() {}
  };


  //=======================================================================
  /// Face geometry of the Face Geometry for element is the same as
  /// that for the underlying wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<FaceGeometry<ProjectableUnsteadyHeatElement<ELEMENT>>>
    : public virtual FaceGeometry<FaceGeometry<ELEMENT>>
  {
  public:
    FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT>>() {}
  };


} // namespace oomph

#endif
