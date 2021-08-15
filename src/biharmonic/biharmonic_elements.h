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
#ifndef OOMPH_BIHARMONIC_ELEMENTS_HEADER
#define OOMPH_BIHARMONIC_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#ifdef OOMPH_HAS_MPI
// mpi headers
#include "mpi.h"
#endif

// Generic C++ headers
#include <iostream>
#include <math.h>

// oomph-lib headers
#include "../generic/matrices.h"
#include "../generic/elements.h"
#include "../generic/hermite_elements.h"


namespace oomph
{
  //=============================================================================
  /// \short Biharmonic Equation Class - contains the equations
  //=============================================================================
  template<unsigned DIM>
  class BiharmonicEquations : public virtual FiniteElement
  {
  public:
    /// source function type definition
    typedef void (*SourceFctPt)(const Vector<double>& x, double& u);


    /// Constructor (must initialise the Source_fct_pt to null)
    BiharmonicEquations() : Source_fct_pt(0) {}


    ~BiharmonicEquations(){};

    /// Access function: Nodal function value at local node n
    /// Uses suitably interpolated value for hanging nodes.
    virtual double u(const unsigned& n, const unsigned& k) const = 0;


    /// gets source strength
    virtual void get_source(const unsigned& ipt,
                            const Vector<double>& x,
                            double& source) const
    {
      // if source function is not provided, i.e. zero, then return zero
      if (Source_fct_pt == 0)
      {
        source = 0.0;
      }

      // else get source strength
      else
      {
        (*Source_fct_pt)(x, source);
      }
    }


    /// wrapper function, adds contribution to residual
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // create a dummy matrix
      DenseDoubleMatrix dummy(1);

      // call the generic residuals functions with flag set to zero
      fill_in_generic_residual_contribution_biharmonic(residuals, dummy, 0);
    }


    /// wrapper function, adds contribution to residual and generic
    void fill_in_contribution_to_jacobian(Vector<double>& residual,
                                          DenseMatrix<double>& jacobian)
    {
      // call generic routine with flag set to 1
      fill_in_generic_residual_contribution_biharmonic(residual, jacobian, 1);
    }


    /// output with nplot points
    void output(std::ostream& outfile, const unsigned& nplot)
    {
      // Vector of local coordinates
      Vector<double> s(DIM);

      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, nplot, s);

        for (unsigned i = 0; i < DIM; i++)
        {
          outfile << this->interpolated_x(s, i) << " ";
        }

        outfile << interpolated_u_biharmonic(s) << std::endl;
      }

      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);
    }


    /// Output at default number of plot points
    void output(std::ostream& outfile)
    {
      FiniteElement::output(outfile);
    }

    /// C-style output
    void output(FILE* file_pt)
    {
      FiniteElement::output(file_pt);
    }

    /// C_style output at n_plot points
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FiniteElement::output(file_pt, n_plot);
    }


    /// output fluid velocity field
    void output_fluid_velocity(std::ostream& outfile, const unsigned& nplot)
    {
      // Vector of local coordinates
      Vector<double> s(DIM);

      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, nplot, s);

        for (unsigned i = 0; i < DIM; i++)
        {
          outfile << this->interpolated_x(s, i) << " ";
        }

        Vector<double> dudx(2, 0.0);
        interpolated_dudx(s, dudx);

        outfile << dudx[1] << "  " << -dudx[0] << std::endl;
      }

      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);
    }


    /// output with nplot points
    void interpolated_dudx(const Vector<double>& s, Vector<double>& dudx)
    {
      // Find out how many nodes there are
      unsigned n_node = this->nnode();

      // Find out how many values there
      unsigned n_value = this->node_pt(0)->nvalue();

      // set up memory for shape functions
      Shape psi(n_node, n_value);
      DShape dpsidx(n_node, n_value, DIM);

      // evaluate dpsidx at local coordinate s
      dshape_eulerian(s, psi, dpsidx);

      // initialise storage for d2u_interpolated
      dudx[0] = 0.0;
      dudx[1] = 0.0;

      // loop over nodes, degrees of freedom, and dimension to calculate
      // interpolated_d2u
      for (unsigned n = 0; n < n_node; n++)
      {
        for (unsigned k = 0; k < n_value; k++)
        {
          for (unsigned d = 0; d < DIM; d++)
          {
            dudx[d] += this->node_pt(n)->value(k) * dpsidx(n, k, d);
          }
        }
      }
    }


    /// output analytic solution
    void output_fct(std::ostream& outfile,
                    const unsigned& nplot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      // Vector of local coordinates
      Vector<double> s(DIM);

      // Vector for coordintes
      Vector<double> x(DIM);

      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);

      // Exact solution Vector (here a scalar)
      Vector<double> exact_soln(1);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, nplot, s);

        // Get x position as Vector
        this->interpolated_x(s, x);

        // Get exact solution at this point
        (*exact_soln_pt)(x, exact_soln);

        // Output x,y,...,u_exact
        for (unsigned i = 0; i < DIM; i++)
        {
          outfile << x[i] << " ";
        }
        outfile << exact_soln[0] << std::endl;
      }

      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);
    }

    /// \short Output exact solution specified via function pointer
    /// at a given time and at a given number of plot points.
    /// Function prints as many components as are returned in solution Vector.
    /// Implement broken FiniteElement base class version
    void output_fct(std::ostream& outfile,
                    const unsigned& nplot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      FiniteElement::output_fct(outfile, nplot, time, exact_soln_pt);
    }


    /// computes the error
    void compute_error(std::ostream& outfile,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       double& error,
                       double& norm)
    {
      // Initialise
      error = 0.0;
      norm = 0.0;

      // Vector of local coordinates
      Vector<double> s(DIM);

      // Vector for coordintes
      Vector<double> x(DIM);

      // Find out how many nodes there are in the element
      unsigned n_node = this->nnode();

      Shape psi(n_node);

      // Set the value of n_intpt
      unsigned n_intpt = this->integral_pt()->nweight();

      // Tecplot header info
      outfile << this->tecplot_zone_string(3);

      // Tecplot
      // outfile << "ZONE" << std::endl;

      // Exact solution Vector (here a scalar)
      Vector<double> exact_soln(1);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Assign values of s
        for (unsigned i = 0; i < DIM; i++)
        {
          s[i] = this->integral_pt()->knot(ipt, i);
        }

        // Get the integral weight
        double w = this->integral_pt()->weight(ipt);

        // Get jacobian of mapping
        double J = this->J_eulerian(s);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Get x position as Vector
        this->interpolated_x(s, x);

        // Get FE function value
        double u_fe = interpolated_u_biharmonic(s);

        // Get exact solution at this point
        (*exact_soln_pt)(x, exact_soln);

        // Output x,y,...,error
        for (unsigned i = 0; i < DIM; i++)
        {
          outfile << x[i] << " ";
        }
        outfile << exact_soln[0] << " " << fabs(exact_soln[0] - u_fe)
                << std::endl;

        // Add to error and norm
        norm += exact_soln[0] * exact_soln[0] * W;
        error += (exact_soln[0] - u_fe) * (exact_soln[0] - u_fe) * W;
      }

      this->write_tecplot_zone_footer(outfile, 3);
    }


    /// \short Plot the error when compared
    /// against a given time-dependent exact solution \f$ {\bf f}(t,{\bf x})
    /// \f$. Also calculates the norm of the error and that of the exact
    /// solution. Call broken base-class version.
    void compute_error(std::ostream& outfile,
                       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                       const double& time,
                       double& error,
                       double& norm)
    {
      FiniteElement::compute_error(outfile, exact_soln_pt, time, error, norm);
    }

    /// calculates interpolated u at s
    double interpolated_u_biharmonic(const Vector<double>& s)
    {
      // initialise storage for u_interpolated
      double uu = 0;

      // number of nodes
      unsigned n_node = this->nnode();

      // number of degrees of freedom per node
      unsigned n_value = this->node_pt(0)->nvalue();

      // set up memory for shape functions
      Shape psi(n_node, n_value);

      // find shape fn at position s
      this->shape(s, psi);

      // calculate interpolated u
      for (unsigned n = 0; n < n_node; n++)
      {
        for (unsigned k = 0; k < n_value; k++)
        {
          uu += u(n, k) * psi(n, k);
        }
      }

      // return interpolated_u
      return uu;
    }


    /// self test wrapper
    unsigned self_test()
    {
      bool passed = true;

      // Check lower-level stuff
      if (FiniteElement::self_test() != 0)
      {
        passed = false;
      }

      // Return verdict
      if (passed)
      {
        return 0;
      }
      else
      {
        return 1;
      }
    }


    /// return number of second derivate degrees of freedom
    unsigned get_d2_dof()
    {
      return d2_dof[DIM - 1];
    }


    /// \short The number of "DOF types" that degrees of freedom in this element
    /// are sub-divided into (for block preconditioning)
    unsigned ndof_types() const
    {
      return this->required_nvalue(1);
    }


    /// \short Create a list of pairs for all unknowns in this element,
    /// so that the first entry in each pair contains the global equation
    /// number of the unknown, while the second one contains the number
    /// of the "DOF types" that this unknown is associated with.
    /// (Function can obviously only be called if the equation numbering
    /// scheme has been set up.) (for block preconditioning)
    void get_dof_numbers_for_unknowns(
      std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
    {
      // number of nodes
      int n_node = this->nnode();

      // number of degrees of freedom at each node
      int n_value = this->node_pt(0)->nvalue();

      // temporary pair (used to store dof lookup prior to being added to list
      std::pair<unsigned long, unsigned> dof_lookup;

      // loop over the nodes
      for (int n = 0; n < n_node; n++)
      {
        // loop over the degree of freedom
        for (int k = 0; k < n_value; k++)
        {
          // determine local eqn number
          int local_eqn_number = this->nodal_local_eqn(n, k);

          // if local equation number is less than zero then nodal dof pinned
          // then ignore
          if (local_eqn_number >= 0)
          {
            // store dof lookup in temporary pair
            dof_lookup.first = this->eqn_number(local_eqn_number);
            dof_lookup.second = k;
            dof_lookup_list.push_front(dof_lookup);
            // add to list
          }
        }
      }
    }


    /// Access functions for the source function pointer
    SourceFctPt& source_fct_pt()
    {
      return Source_fct_pt;
    }

    /// Access functions for the source function pointers (const version)
    SourceFctPt source_fct_pt() const
    {
      return Source_fct_pt;
    }


  protected:
    /// \short Compute element residual Vector only (if JFLAG=and/or element
    /// Jacobian matrix
    virtual void fill_in_generic_residual_contribution_biharmonic(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned JFLAG);


    /// Pointer to source function:
    SourceFctPt Source_fct_pt;


    /// Array to hold local eqn numbers: Local_eqn[n] (=-1 for BC)
    Vector<int> Local_eqn;

  private:
    // number of degrees of freedom of second derivative
    static const unsigned d2_dof[];
  };


  // declares number of degrees of freedom of second derivative
  template<unsigned DIM>
  const unsigned BiharmonicEquations<DIM>::d2_dof[3] = {1, 3, 6};


  //=============================================================================
  ///  biharmonic element class
  //=============================================================================
  template<unsigned DIM>
  class BiharmonicElement : public virtual QHermiteElement<DIM>,
                            public virtual BiharmonicEquations<DIM>
  {
  public:
    /// access function for value, kth dof of node n
    inline double u(const unsigned& n, const unsigned& k) const
    {
      return this->node_pt(n)->value(k);
    }


    ///\short  Constructor: Call constructors for QElement and
    /// Poisson equations
    BiharmonicElement() : QHermiteElement<DIM>(), BiharmonicEquations<DIM>() {}


    ~BiharmonicElement(){};


    /// \short  Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return DIM * 2;
    }


    /// Output
    void output(std::ostream& outfile)
    {
      BiharmonicEquations<DIM>::output(outfile);
    }

    /// output wrapper
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      BiharmonicEquations<DIM>::output(outfile, n_plot);
    }

    /// C-style output
    void output(FILE* file_pt)
    {
      BiharmonicEquations<DIM>::output(file_pt);
    }

    /// C_style output at n_plot points
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      BiharmonicEquations<DIM>::output(file_pt, n_plot);
    }


    /// analytic solution wrapper
    void output_fct(std::ostream& outfile,
                    const unsigned& nplot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      BiharmonicEquations<DIM>::output_fct(outfile, nplot, exact_soln_pt);
    }


    /// Final override
    void output_fct(std::ostream& outfile,
                    const unsigned& nplot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      BiharmonicEquations<DIM>::output_fct(outfile, nplot, time, exact_soln_pt);
    }


    /// computes error
    void compute_error(std::ostream& outfile,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       double& error,
                       double& norm)
    {
      BiharmonicEquations<DIM>::compute_error(
        outfile, exact_soln_pt, error, norm);
    }

    /// Call the equations-class overloaded unsteady error calculation
    void compute_error(std::ostream& outfile,
                       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                       const double& time,
                       double& error,
                       double& norm)
    {
      BiharmonicEquations<DIM>::compute_error(
        outfile, exact_soln_pt, time, error, norm);
    }
  };


} // namespace oomph
#endif
