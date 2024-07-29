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
// Include guard to prevent multiple inclusions of the header
#ifndef OOMPH_NAVIER_STOKES_FLUX_CONTROL_ELEMENTS
#define OOMPH_NAVIER_STOKES_FLUX_CONTROL_ELEMENTS

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// OOMPH-LIB headers
#include "../generic/nodes.h"
#include "../navier_stokes/navier_stokes_surface_power_elements.h"

namespace oomph
{
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////


  //======================================================================
  /// A template free base class for an element to imposes an applied
  /// boundary pressure to the Navier-Stokes equations in order to
  /// control a volume flux when used in conjunction with a
  /// NetFluxControlElement or
  /// NetFluxControlElementForWomersleyPressureControl).
  //======================================================================
  class TemplateFreeNavierStokesFluxControlElementBase
    : public virtual GeneralisedElement
  {
  public:
    /// Empty constructor
    TemplateFreeNavierStokesFluxControlElementBase() {}

    /// Empty virtual destructor
    virtual ~TemplateFreeNavierStokesFluxControlElementBase() {}

    /// Pure virtual function to calculate integral of the volume flux
    virtual double get_volume_flux() = 0;

    /// Function adds to the external data the Data object whose
    /// single value is the pressure applied by the element
    void add_pressure_data(Data* pressure_data_pt)
    {
      Pressure_data_id = add_external_data(pressure_data_pt);
    }

  protected:
    /// Access function gives id of external Data object whose
    /// single value is the pressure applied by the element
    unsigned& pressure_data_id()
    {
      return Pressure_data_id;
    }

  private:
    /// Id of external Data object whose single value is the
    /// pressure applied by the elements
    unsigned Pressure_data_id;
  };


  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////


  //======================================================================
  /// A class for an element that controls the net fluid flux across a
  /// boundary by the imposition of an unknown applied pressure to the
  /// Navier-Stokes equations. This element is used with a mesh of
  /// NavierStokesFluxControlElement elements which are attached
  /// to the boundary.
  /// Note: fill_in_contribution_to_jacobian() does not calculate
  /// Jacobian contributions for this element as they are calculated by
  /// NavierStokesFluxControlElement::fill_in_contribution_to_jacobian(...)
  //======================================================================
  class NetFluxControlElement : public virtual GeneralisedElement
  {
  public:
    /// Constructor takes a mesh of
    /// TemplateFreeNavierStokesFluxControlElementBase elements
    /// that impose the pressure to control the flux, plus a pointer to
    /// the double which contains the desired flux value
    NetFluxControlElement(Mesh* flux_control_mesh_pt,
                          double* prescribed_flux_value_pt)
      : Flux_control_mesh_pt(flux_control_mesh_pt),
        Prescribed_flux_value_pt(prescribed_flux_value_pt)
    {
      // Construct Pressure_data_pt
      Pressure_data_pt = new Data(1);

      // Add the new Data to internal Data for this element
      add_internal_data(Pressure_data_pt);

      // There's no need to add the external data for this element since
      // this elements Jacobian contributions are calculated by the
      // NavierStokesFluxControlElements

      // Loop over elements in the Flux_control_mesh to add this element's
      // Data to the external Data in the elements in the flux control mesh
      unsigned n_el = Flux_control_mesh_pt->nelement();
      for (unsigned e = 0; e < n_el; e++)
      {
        // Get pointer to the element
        GeneralisedElement* el_pt = Flux_control_mesh_pt->element_pt(e);

        // Perform cast to TemplateFreeNavierStokesFluxControlElementBase
        // pointer
        TemplateFreeNavierStokesFluxControlElementBase* flux_el_pt =
          dynamic_cast<TemplateFreeNavierStokesFluxControlElementBase*>(el_pt);

        flux_el_pt->add_pressure_data(Pressure_data_pt);
      }

      // Default value for Dof_number_for_unknown, indiating that it's
      // uninitialised
      Dof_number_for_unknown = UINT_MAX;
    }


    /// Empty Destructor - Data gets deleted automatically
    ~NetFluxControlElement() {}

    /// Broken copy constructor
    NetFluxControlElement(const NetFluxControlElement& dummy) = delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const NetFluxControlElement&) = delete;*/

    /// Spatial dimension of the problem
    unsigned dim() const
    {
      return Dim;
    }

    /// Function to return a pointer to the Data object whose
    /// single value is the pressure applied by the
    /// NavierStokesFluxControlElement elements
    Data* pressure_data_pt() const
    {
      return Pressure_data_pt;
    }


    /// Add the element's contribution to its residual vector:
    /// i.e. the flux constraint.
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic routine
      fill_in_generic_residual_contribution_flux_control(residuals);
    }

    /// This function returns the residuals, but adds nothing to the
    /// Jacobian as this element's Jacobian contributions are calculated by
    /// the NavierStokesFluxControlElements which impose the traction
    /// used to control the flux.
    inline void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                                 DenseMatrix<double>& jacobian)
    {
      // Call the generic routine
      fill_in_generic_residual_contribution_flux_control(residuals);
    }


    /// The number of "DOF types" that degrees of freedom in this element
    /// are sub-divided into - it's set to Dof_number_for_unknown+1
    /// because it's expected this element is added to a fluid mesh
    /// containing navier stokes elements
    unsigned ndof_types() const
    {
#ifdef PARANOID
      if (Dof_number_for_unknown == UINT_MAX)
      {
        std::ostringstream error_message;
        error_message << "Dof_number_for_unknown hasn't been set yet!\n"
                      << "Please do so using the dof_number_for_unknown()\n"
                      << "access function\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return Dof_number_for_unknown + 1;
    }

    /// Function to set / get the nodal value of the "DOF type" to which
    /// the degree of freedom in this element (the pressure that enforces
    /// the required volume flux!) is added to.
    /// This should be set to the Navier-Stokes pressure DOF type
    /// (usually the dimension of the problem, for example, in 3D, the DOF types
    /// for single-physics Navier-Stokes elements are usually
    /// labelled 0, 1, 2, 3 for u, v and w
    /// velocities and pressure respectively. It is important to note that this
    /// is dimension dependent, so should not be hard coded in!! In
    /// particularly, this should not simply be set to the dimension of the
    /// problem if there is further splitting of the velocity DOF types) if this
    /// element is added to a fluid mesh containing Navier-Stokes elements.
    unsigned& dof_number_for_unknown()
    {
      return Dof_number_for_unknown;
    }

    /// Create a list of pairs for all unknowns in this element,
    /// so that the first entry in each pair contains the global equation
    /// number of the unknown, while the second one contains the number
    /// of the "DOF type" that this unknown is associated with.
    /// (Function can obviously only be called if the equation numbering
    /// scheme has been set up.) The single degree of freedom is given the
    /// DOF type number of Dof_number_for_unknown since it's expected this
    /// unknown is added to the Navier-Stokes pressure DOF block (it is also
    /// assumed that the user has set the Dof_number_for_unknown variable to
    /// the velocity DOF type using the function dof_number_for_unknown()).
    void get_dof_numbers_for_unknowns(
      std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
    {
#ifdef PARANOID
      if (Dof_number_for_unknown == UINT_MAX)
      {
        std::ostringstream error_message;
        error_message << "Dof_number_for_unknown hasn't been set yet!\n"
                      << "Please do so using the dof_number_for_unknown()\n"
                      << "access function\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // pair to store dof lookup prior to being added to list
      std::pair<unsigned, unsigned> dof_lookup;

      dof_lookup.first = this->eqn_number(0);
      dof_lookup.second = Dof_number_for_unknown;

      // add to list
      dof_lookup_list.push_front(dof_lookup);
    }

  protected:
    /// This function returns the residuals for the
    /// flux control master element.
    void fill_in_generic_residual_contribution_flux_control(
      Vector<double>& residuals)
    {
      // Initialise volume flux
      double volume_flux = 0.0;

      // Loop over elements in Flux_control_mesh_pt and calculate flux
      unsigned n_el = Flux_control_mesh_pt->nelement();
      for (unsigned e = 0; e < n_el; e++)
      {
        // Get a pointer to the element
        GeneralisedElement* el_pt = Flux_control_mesh_pt->element_pt(e);

        // Cast to NavierStokesFluxControlElement
        TemplateFreeNavierStokesFluxControlElementBase* flux_control_el_pt = 0;
        flux_control_el_pt =
          dynamic_cast<TemplateFreeNavierStokesFluxControlElementBase*>(el_pt);

#ifdef PARANOID
        if (flux_control_el_pt == 0)
        {
          throw OomphLibError("Element must be used with a mesh of "
                              "NavierStokesFluxControlElements",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Add the elemental volume flux
        volume_flux += flux_control_el_pt->get_volume_flux();
      }

      residuals[0] += *Prescribed_flux_value_pt - volume_flux;
    }


  private:
    /// Data object whose single value is the pressure
    /// applied by the elements in the Flux_control_mesh_pt
    Data* Pressure_data_pt;

    /// Mesh of elements which impose the pressure which controls
    /// the net flux
    Mesh* Flux_control_mesh_pt;

    /// Pointer to the value that stores the prescribed flux
    double* Prescribed_flux_value_pt;

    /// The id number of the "DOF type" to which the degree
    /// of freedom in this element is added to. This should be set to the
    /// number id of the Navier-Stokes pressure DOF block (which is dimension
    /// dependent!) if this element is added to a fluid mesh
    /// containing navier stokes elements
    unsigned Dof_number_for_unknown;

    /// spatial dim of NS system
    unsigned Dim;
  };


  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////


  //======================================================================
  /// A class of element to impose an applied boundary pressure to
  /// Navier-Stokes elements to control to control a volume flux. A mesh of
  /// these elements are used in conjunction with a NetFluxControlElement.
  /// The template arguement ELEMENT is a Navier-Stokes "bulk" element.
  ///
  /// Note: This element calculates Jacobian contributions for both itself
  /// and also for the NetFluxControlElement with respect to its unknowns.
  //======================================================================
  template<class ELEMENT>
  class NavierStokesFluxControlElement
    : public virtual TemplateFreeNavierStokesFluxControlElementBase,
      public virtual NavierStokesSurfacePowerElement<ELEMENT>
  {
  public:
    /// Constructor, which takes a "bulk" element and face index
    NavierStokesFluxControlElement(
      FiniteElement* const& element_pt,
      const int& face_index,
      const bool& called_from_refineable_constructor = false)
      : NavierStokesSurfacePowerElement<ELEMENT>(element_pt, face_index)
    {
#ifdef PARANOID
      {
        // Check that the element is not a refineable 3d element
        if (!called_from_refineable_constructor)
        {
          ELEMENT* elem_pt = new ELEMENT;
          // If it's three-d
          if (elem_pt->dim() == 3)
          {
            // Is it refineable
            if (dynamic_cast<RefineableElement*>(elem_pt))
            {
              // Throw Error
              std::ostringstream error_message;
              error_message
                << "This element does not work properly with refineable bulk \n"
                << "elements in 3D. Please use the refineable version\n"
                << "instead.\n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
      }
#endif

      // Set the dimension from the dimension of the first node (since Dim is
      // private in the parent class)
      Dim = this->node_pt(0)->ndim();
    }

    /// Destructor should not delete anything
    ~NavierStokesFluxControlElement() {}

    /// This function returns just the residuals
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function using a dummy matrix argument
      fill_in_generic_residual_contribution_fluid_traction(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// This function returns the residuals and the Jacobian
    /// including the Jacobian contribution from the flux control
    /// master element with respect to dof in this
    /// element
    inline void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                                 DenseMatrix<double>& jacobian)
    {
      // Call the generic routine
      fill_in_generic_residual_contribution_fluid_traction(
        residuals, jacobian, 1);
    }

    /// Function to get the integral of the volume flux
    double get_volume_flux()
    {
      return NavierStokesSurfacePowerElement<ELEMENT>::get_volume_flux();
    }

  protected:
    /// Access function that returns the local equation numbers
    /// for velocity components.
    /// u_local_eqn(n,i) = local equation number or < 0 if pinned.
    /// The default is to asssume that n is the local node number
    /// and the i-th velocity component is the i-th unknown stored at the node.
    virtual inline int u_local_eqn(const unsigned& n, const unsigned& i)
    {
      return this->nodal_local_eqn(n, i);
    }

    /// Function to compute the shape and test functions and to return
    /// the Jacobian of mapping
    inline double shape_and_test_at_knot(const unsigned& ipt,
                                         Shape& psi,
                                         Shape& test) const
    {
      // Find number of nodes
      unsigned n_node = this->nnode();
      // Calculate the shape functions
      this->shape_at_knot(ipt, psi);
      // Set the test functions to be the same as the shape functions
      for (unsigned i = 0; i < n_node; i++)
      {
        test[i] = psi[i];
      }
      // Return the value of the jacobian
      return this->J_eulerian_at_knot(ipt);
    }


    /// This function returns the residuals for the traction function
    /// flag=1(or 0): do (or don't) compute the Jacobian as well.
    /// This function also calculates the Jacobian contribution for the
    /// NetFluxControlElement
    void fill_in_generic_residual_contribution_fluid_traction(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
    {
      // Find out how many nodes there are
      unsigned n_node = this->nnode();

      // Set up memory for the shape and test functions
      Shape psif(n_node), testf(n_node);

      // Set the value of n_intpt
      unsigned n_intpt = this->integral_pt()->nweight();

      // Integers to store local equation numbers
      int local_eqn = 0;

      // Get the pressure at the outflow
      double pressure = this->external_data_pt(pressure_data_id())->value(0);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the integral weight
        double w = this->integral_pt()->weight(ipt);

        // Find the shape and test functions and return the Jacobian
        // of the mapping
        double J = shape_and_test_at_knot(ipt, psif, testf);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Get the outer unit normal
        Vector<double> unit_normal(Dim);
        this->outer_unit_normal(ipt, unit_normal);

        // Calculate the traction
        Vector<double> traction(Dim);
        for (unsigned i = 0; i < Dim; i++)
        {
          traction[i] = -pressure * unit_normal[i];
        }

        // Loop over the test functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the velocity components
          for (unsigned i = 0; i < Dim; i++)
          {
            local_eqn = u_local_eqn(l, i);

            /*IF it's not a boundary condition*/
            if (local_eqn >= 0)
            {
              // Add the user-defined traction terms
              residuals[local_eqn] += traction[i] * testf[l] * W;

              // Calculate the Jacobian if required. It is assumed
              // that traction DOES NOT depend upon velocities
              // or pressures in the Navier Stokes elements, but
              // depend in the Data value which holds the
              // pressure.
              if (flag)
              {
                // Get equation number of the pressure data unknown
                int local_unknown =
                  this->external_local_eqn(pressure_data_id(), 0);

                // IF it's not a boundary condition
                if (local_unknown >= 0)
                {
                  // Add to Jacobian for this element
                  double jac_contribution = -unit_normal[i] * testf[l] * W;
                  jacobian(local_eqn, local_unknown) += jac_contribution;

                  // Add to Jacobian for master element
                  jacobian(local_unknown, local_eqn) += jac_contribution;
                }
              }
            }
          } // End of loop over dimension
        } // End of loop over shape functions
      }
    }

  protected:
    /// The highest dimension of the problem
    unsigned Dim;
  };


  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////


  //======================================================================
  /// A class of element to impose an applied boundary pressure to
  /// Navier-Stokes elements to control to control a volume flux. A mesh of
  /// these elements are used in conjunction with a NetFluxControlElement.
  /// The template arguement ELEMENT is a Navier-Stokes "bulk" element.
  ///
  /// Note: This element calculates Jacobian contributions for both itself
  /// and also for the NetFluxControlElement with respect to its unknowns.
  ///
  /// THIS IS THE REFINEABLE VERSION.
  //======================================================================
  template<class ELEMENT>
  class RefineableNavierStokesFluxControlElement
    : public virtual NavierStokesFluxControlElement<ELEMENT>,
      public virtual NonRefineableElementWithHangingNodes
  {
  public:
    /// Constructor, which takes a "bulk" element and the face index
    RefineableNavierStokesFluxControlElement(FiniteElement* const& element_pt,
                                             const int& face_index)
      : NavierStokesSurfacePowerElement<ELEMENT>(element_pt, face_index),
        // we're calling this from the constructor of the refineable version.
        NavierStokesFluxControlElement<ELEMENT>(element_pt, face_index, true)
    {
    }

    /// Destructor should not delete anything
    ~RefineableNavierStokesFluxControlElement() {}


    /// Number of continuously interpolated values are the
    /// same as those in the bulk element.
    unsigned ncont_interpolated_values() const
    {
      return dynamic_cast<ELEMENT*>(this->bulk_element_pt())
        ->ncont_interpolated_values();
    }

    /// This function returns just the residuals
    inline void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function using a dummy matrix argument
      refineable_fill_in_generic_residual_contribution_fluid_traction(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// This function returns the residuals and the Jacobian
    /// including the Jacobian contribution from the flux control
    /// master element with respect to dof in this
    /// element
    inline void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                                 DenseMatrix<double>& jacobian)
    {
      // Call the generic routine
      refineable_fill_in_generic_residual_contribution_fluid_traction(
        residuals, jacobian, 1);
    }

  protected:
    /// This function returns the residuals for the traction function
    /// flag=1(or 0): do (or don't) compute the Jacobian as well.
    /// This function also calculates the Jacobian contribution for the
    /// NetFluxControlElement
    void refineable_fill_in_generic_residual_contribution_fluid_traction(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
    {
      // Get the indices at which the velocity components are stored
      unsigned u_nodal_index[this->Dim];
      for (unsigned i = 0; i < this->Dim; i++)
      {
        u_nodal_index[i] =
          dynamic_cast<ELEMENT*>(this->bulk_element_pt())->u_index_nst(i);
      }

      // Pointer to hang info object
      HangInfo* hang_info_pt = 0;

      // Find out how many nodes there are
      unsigned n_node = this->nnode();

      // Set up memory for the shape and test functions
      Shape psif(n_node), testf(n_node);

      // Set the value of n_intpt
      unsigned n_intpt = this->integral_pt()->nweight();

      // Integers to store local equation numbers
      int local_eqn = 0;

      // Get the pressure at the outflow
      double pressure =
        this->external_data_pt(this->pressure_data_id())->value(0);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the integral weight
        double w = this->integral_pt()->weight(ipt);

        // Find the shape and test functions and return the Jacobian
        // of the mapping
        double J = this->shape_and_test_at_knot(ipt, psif, testf);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Get the outer unit normal
        Vector<double> unit_normal(this->Dim);
        this->outer_unit_normal(ipt, unit_normal);

        // Calculate the traction
        Vector<double> traction(this->Dim);
        for (unsigned i = 0; i < this->Dim; i++)
        {
          traction[i] = -pressure * unit_normal[i];
        }


        // Number of master nodes and storage for the weight of the shape
        // function
        unsigned n_master = 1;
        double hang_weight = 1.0;

        // Loop over the nodes for the test functions/equations
        //----------------------------------------------------
        for (unsigned l = 0; l < n_node; l++)
        {
          // Local boolean to indicate whether the node is hanging
          bool is_node_hanging = this->node_pt(l)->is_hanging();

          // If the node is hanging
          if (is_node_hanging)
          {
            hang_info_pt = this->node_pt(l)->hanging_pt();

            // Read out number of master nodes from hanging data
            n_master = hang_info_pt->nmaster();
          }
          // Otherwise the node is its own master
          else
          {
            n_master = 1;
          }

          // Loop over the master nodes
          for (unsigned m = 0; m < n_master; m++)
          {
            // Loop over velocity components for equations
            for (unsigned i = 0; i < this->Dim; i++)
            {
              // Get the equation number
              // If the node is hanging
              if (is_node_hanging)
              {
                // Get the equation number from the master node
                local_eqn = this->local_hang_eqn(
                  hang_info_pt->master_node_pt(m), u_nodal_index[i]);
                // Get the hang weight from the master node
                hang_weight = hang_info_pt->master_weight(m);
              }
              // If the node is not hanging
              else
              {
                // Local equation number
                local_eqn = this->nodal_local_eqn(l, u_nodal_index[i]);

                // Node contributes with full weight
                hang_weight = 1.0;
              }

              // If it's not a boundary condition...
              if (local_eqn >= 0)
              {
                // Add the user-defined traction terms
                residuals[local_eqn] +=
                  traction[i] * testf[l] * W * hang_weight;

                // Calculate the Jacobian if required. It is assumed
                // that traction DOES NOT depend upon velocities
                // or pressures in the Navier Stokes elements, but
                // depend in the Data value which holds the
                // pressure.
                if (flag)
                {
                  // Get equation number of the pressure data unknown
                  int local_unknown =
                    this->external_local_eqn(this->pressure_data_id(), 0);

                  // IF it's not a boundary condition
                  if (local_unknown >= 0)
                  {
                    // Add to Jacobian for this element
                    double jac_contribution =
                      -unit_normal[i] * testf[l] * W * hang_weight;
                    jacobian(local_eqn, local_unknown) += jac_contribution;

                    // Add to Jacobian for master element
                    jacobian(local_unknown, local_eqn) += jac_contribution;
                  }
                }
              }
            } // End of loop over dimension
          } // End of loop over master nodes
        } // End of loop over nodes
      }
    }
  };


} // namespace oomph

#endif
