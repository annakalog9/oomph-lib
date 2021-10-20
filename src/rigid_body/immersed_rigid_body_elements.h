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
// Header file for Rigid Body Elements that are immersed in an external fluid
#ifndef OOMPH_IMMERSED_RIGID_BODY_ELEMENTS_HEADER
#define OOMPH_IMMERSED_RIGID_BODY_ELEMENTS_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include "../generic/elements.h"
#include "../generic/triangle_mesh.h"
#include "../generic/fsi.h"

namespace oomph
{
  //=====================================================================
  // DOXYERROR: I removed a \f from the formula \f$\mbox{\boldmath$F$} =
  // \mbox{\boldmath$F$}^{*}/\mu U\f$ below because it broke doxygen. Someone
  // who knows the maths should check that this still makes sense.
  /// Class that solves the equations of motion for a general
  /// two-dimensional rigid body subject to a particular imposed force and
  /// torque distribution and immersed within an external fluid. The body's
  /// position is entirely specified by the location of its centre of mass,
  /// \f$\mbox{\boldmath$X$}\f$, and a single angle, \f$\phi\f$, that represents
  /// a possible rotation. The equations of motion are then simply Newton's
  /// second law for the conservation of linear momentum in two directions and
  /// angular momentum about the single possible axis of rotation.
  ///
  /// The non-dimensionalisation is based on the viscous scales of
  /// the surrounding fluid, in which case, the governing equations are
  /// \f[  Re St^{2} M \frac{\mbox{d}^{2} \mbox{\boldmath$X$}}{\mbox{d} t^{2}}
  /// = \mbox{\boldmath$F$} + \oint \mbox{\boldmath$t$} \mbox{d} s
  /// \quad\mbox{and}\quad Re St^{2} I \frac{\mbox{d}^{2}\Phi}{\mbox{d} t^{2}}
  /// = \mbox{\boldmath$T$} + \oint \mbox{\boldmath$q$} \mbox{d} s,
  /// \f]
  /// where \f$\mbox{\boldmath$F$} = \mbox{\boldmath$F$}^{*}/\mu U\f$ is the
  /// external force per unit length;
  /// \f$\mbox{\boldmath$t$}\f$ is the net force per unit length
  /// applied on the rigid body
  /// by the surrounding fluid;
  /// \f$\mbox{\boldmath$T$} = \mbox{\boldmath$T$}^{*}/(\mu UL)\f$
  /// is the external torque per unit length; and \f$\mbox{\boldmath$q$}\f$ is
  /// the net torque per unit length applied on the body by the fluid.
  /// \f$M\f$ is a scaled mass the density ratio multiplied by a shape parameter
  /// and \f$I\f$ is a scaled moment of inertia. Finally, \f$Re\f$ and \f$St\f$
  /// are the Reynolds and Strouhal numbers of the surrounding fluid.
  /// Note that these equations may be used without the external fluid,
  /// in which case the non-dimensionalisation doesn't make a lot of sense,
  /// but may be re-interpreted accordingly.
  ///
  /// A Data object whose three values
  /// represent the x and y displacements of the body's centre of mass and
  /// its rotation about the centre of mass may be passed in as external
  /// data or, if not, it will be constructed internally.
  ///
  /// For general usage, an underlying geometric object must passed to the
  /// constructor and the position will be determined by applying rigid body
  /// motions to the underlying object based on the initial_centre_of_mass
  /// and initial_phi (angle), which default to zero. If these defaults are
  /// not suitable, the values must be set externally. In addition a mass
  /// and moment of inertia should also be set externally.
  ///
  /// If added to a mesh in the Problem (in its incarnation as a
  /// GeneralisedElement) the displacement/rotation of the body
  /// is computed in response to (i) user-specifiable applied forces
  /// and a torque and (ii) the net drag (and associated torque) from
  /// a mesh of elements that can exert a drag onto the body (typically
  /// Navier-Stokes FaceElements that apply a viscous drag to an
  /// immersed body, represented by the body.)
  //=====================================================================
  class ImmersedRigidBodyElement : public GeneralisedElement, public GeomObject
  {
  public:
    ///  Function pointer to function that specifies
    /// external force
    typedef void (*ExternalForceFctPt)(const double& time,
                                       Vector<double>& external_force);

    ///  Function pointer to function that specifies
    /// external torque
    typedef void (*ExternalTorqueFctPt)(const double& time,
                                        double& external_torque);

  protected:
    ///  Default constructor that intialises everything to
    /// zero. This is expected to be called only from derived clases
    /// such as the ImmersedRigidBodyTriangleMeshPolygon that can provided
    /// their own position() functions.
    ImmersedRigidBodyElement(TimeStepper* const& time_stepper_pt,
                             Data* const& centre_displacement_data_pt = 0)
      : Initial_Phi(0.0),
        Mass(0.0),
        Moment_of_inertia(0.0),
        Centre_displacement_data_pt(centre_displacement_data_pt),
        Geom_object_pt(0),
        External_force_fct_pt(0),
        External_torque_fct_pt(0),
        Drag_mesh_pt(0),
        G_pt(&Default_Gravity_vector),
        Re_pt(&Default_Physical_Constant_Value),
        St_pt(&Default_Physical_Ratio_Value),
        ReInvFr_pt(&Default_Physical_Constant_Value),
        Density_ratio_pt(&Default_Physical_Ratio_Value),
        Include_geometric_rotation(true)
    {
      this->initialise(time_stepper_pt);
    }

  public:
    ///  Constructor that takes an underlying geometric object:
    /// and timestepper.
    ImmersedRigidBodyElement(GeomObject* const& geom_object_pt,
                             TimeStepper* const& time_stepper_pt,
                             Data* const& centre_displacement_data_pt = 0)
      : Initial_Phi(0.0),
        Mass(0.0),
        Moment_of_inertia(0.0),
        Centre_displacement_data_pt(centre_displacement_data_pt),
        Geom_object_pt(geom_object_pt),
        External_force_fct_pt(0),
        External_torque_fct_pt(0),
        Drag_mesh_pt(0),
        G_pt(&Default_Gravity_vector),
        Re_pt(&Default_Physical_Constant_Value),
        St_pt(&Default_Physical_Ratio_Value),
        ReInvFr_pt(&Default_Physical_Constant_Value),
        Density_ratio_pt(&Default_Physical_Ratio_Value),
        Include_geometric_rotation(true)
    {
      this->initialise(time_stepper_pt);
    }

    /// Set the rotation of the object to be included
    void set_geometric_rotation()
    {
      Include_geometric_rotation = true;
    }

    ///  Set the rotation of the object to be ignored (only really
    /// useful if you have a circular shape)
    void unset_geometric_rotation()
    {
      Include_geometric_rotation = false;
    }

    /// Access function for the initial angle
    double& initial_phi()
    {
      return Initial_Phi;
    }

    /// Access function for the initial centre of mass
    double& initial_centre_of_mass(const unsigned& i)
    {
      return Initial_centre_of_mass[i];
    }

    /// Access function for the initial centre of mass (const version)
    const double& initial_centre_of_mass(const unsigned& i) const
    {
      return Initial_centre_of_mass[i];
    }

    /// Overload the position to apply the rotation and translation
    void position(const Vector<double>& xi, Vector<double>& r) const
    {
      Vector<double> initial_x(2);
      Geom_object_pt->position(xi, initial_x);
      this->apply_rigid_body_motion(0, initial_x, r);
    }

    /// Overload to include the time history of the motion of the object
    void position(const unsigned& t,
                  const Vector<double>& xi,
                  Vector<double>& r) const
    {
      Vector<double> initial_x(2);
      Geom_object_pt->position(xi, initial_x);
      this->apply_rigid_body_motion(t, initial_x, r);
    }


    /// Work out the position derivative, including rigid body motion
    void dposition_dt(const Vector<double>& zeta,
                      const unsigned& j,
                      Vector<double>& drdt);

    ///  Destuctor: Cleanup if required
    ~ImmersedRigidBodyElement()
    {
      if (Displacement_data_is_internal)
      {
        delete Centre_displacement_data_pt;
        // Null out the pointer stored as internal data
        this->internal_data_pt(0) = 0;
      }
    }

    /// Access to dimensionless "mass" shape parameter that must be set by hand
    /// for non polygonal shapes
    double& mass_shape()
    {
      return Mass;
    }

    /// Access to dimensionless polar "moment of inertia" shape parameter
    double& moment_of_inertia_shape()
    {
      return Moment_of_inertia;
    }

    ///  Pointer to Data for centre of gravity displacement.
    /// Values: 0: x-displ; 1: y-displ; 2: rotation angle.
    Data*& centre_displacement_data_pt()
    {
      return Centre_displacement_data_pt;
    }

    /// x-displacement of centre of mass
    double& centre_x_displacement()
    {
      return *(Centre_displacement_data_pt->value_pt(0));
    }

    /// y-displacement of centre of mass
    double& centre_y_displacement()
    {
      return *(Centre_displacement_data_pt->value_pt(1));
    }

    /// rotation of centre of mass
    double& centre_rotation_angle()
    {
      return *(Centre_displacement_data_pt->value_pt(2));
    }

    /// Get current centre of gravity
    Vector<double> centre_of_gravity()
    {
      Vector<double> cog(2);
      for (unsigned i = 0; i < 2; i++)
      {
        cog[i] = Initial_centre_of_mass[i] +
                 this->Centre_displacement_data_pt->value(i);
      }
      return cog;
    }

    /// Pin the i-th coordinate of the centre of mass
    void pin_centre_of_mass_coordinate(const unsigned& i)
    {
      Centre_displacement_data_pt->pin(i);
    }

    /// Unpin the i-th coordinate of the centre of mass
    void unpin_centre_of_mass_coordinate(const unsigned& i)
    {
      Centre_displacement_data_pt->unpin(i);
    }

    /// Pin the rotation angle
    void pin_rotation_angle()
    {
      Centre_displacement_data_pt->pin(2);
    }

    /// Unpin the rotation angle
    void unpin_rotation_angle()
    {
      Centre_displacement_data_pt->unpin(2);
    }

    /// Output position velocity and acceleration of centre of gravity
    void output_centre_of_gravity(std::ostream& outfile);

    /// Get the contribution to the residuals
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Get generic function without jacobian terms
      get_residuals_rigid_body_generic(
        residuals, GeneralisedElement::Dummy_matrix, false);
    }


    /// Get residuals including contribution to jacobian
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Get generic function, but don't bother to get jacobian terms
      bool flag = false;
      get_residuals_rigid_body_generic(residuals, jacobian, flag);
      // Get the internal effects by FD, because a change in internal
      // data can change the fluid loading, so this is one of those
      // subtle interactions
      this->fill_in_jacobian_from_internal_by_fd(residuals, jacobian);
      // Get the effect of the fluid loading on the rigid body
      this->fill_in_jacobian_from_external_by_fd(residuals, jacobian);
    }

    ///  Update the positions of the nodes in fluid elements
    /// adjacent to the rigid body, defined as being elements in the
    /// drag mesh.
    inline void node_update_adjacent_fluid_elements()
    {
      // If there is no drag mesh then do nothing
      if (Drag_mesh_pt == 0)
      {
        return;
      }
      // Otherwise update the elements adjacent to the rigid body
      //(all the elements adjacent to those in the drag mesh)
      else
      {
        unsigned nel = Drag_mesh_pt->nelement();
        for (unsigned e = 0; e < nel; e++)
        {
          dynamic_cast<FaceElement*>(Drag_mesh_pt->element_pt(e))
            ->bulk_element_pt()
            ->node_update();
        }
      }
    }


    ///  After an external data change, update the nodal positions
    inline void update_in_external_fd(const unsigned& i)
    {
      node_update_adjacent_fluid_elements();
    }

    /// Do nothing to reset within finite-differencing of  external data
    inline void reset_in_external_fd(const unsigned& i) {}

    /// After all external data finite-differencing, update nodal
    /// positions
    inline void reset_after_external_fd()
    {
      node_update_adjacent_fluid_elements();
    }

    /// After an internal data change, update the nodal positions
    inline void update_in_internal_fd(const unsigned& i)
    {
      node_update_adjacent_fluid_elements();
    }

    /// Do nothing to reset within finite-differencing of internal data
    inline void reset_in_internal_fd(const unsigned& i) {}

    /// After all internal data finite-differencing, update nodal
    /// positions
    inline void reset_after_internal_fd()
    {
      node_update_adjacent_fluid_elements();
    }

    ///  Get force and torque from specified fct pointers and
    /// drag mesh
    void get_force_and_torque(const double& time,
                              Vector<double>& force,
                              double& torque);

    ///  Access to function pointer to function that specifies
    /// external force
    ExternalForceFctPt& external_force_fct_pt()
    {
      return External_force_fct_pt;
    }

    ///  Access to function pointer to function that specifies
    /// external torque
    ExternalTorqueFctPt& external_torque_fct_pt()
    {
      return External_torque_fct_pt;
    }

    ///  Access fct to mesh containing face elements that allow
    /// the computation of the drag on the body
    Mesh* const& drag_mesh_pt()
    {
      return Drag_mesh_pt;
    }

    ///  Function to set the drag mesh and add the appropriate load
    /// and geometric data as external data to the Rigid Body
    void set_drag_mesh(Mesh* const& drag_mesh_pt);

    ///  Function to clear the drag mesh and all associated external data
    void flush_drag_mesh()
    {
      // Delete the hijacked data created as the load
      this->delete_external_hijacked_data();
      // Flush the external data
      this->flush_external_data();
      // Set the Drag_mesh pointer to null
      Drag_mesh_pt = 0;
    }

    /// The position of the object depends on one data item
    unsigned ngeom_data() const
    {
      return 1;
    }

    ///  Return pointer to the j-th (only) Data item that the object's
    /// shape depends on.
    Data* geom_data_pt(const unsigned& j)
    {
      return this->Centre_displacement_data_pt;
    }

    ///  Access function to the direction of gravity
    Vector<double>*& g_pt()
    {
      return G_pt;
    }

    ///  Access function for gravity
    const Vector<double>& g() const
    {
      return *G_pt;
    }

    ///  Access function for the pointer to the fluid Reynolds number
    double*& re_pt()
    {
      return Re_pt;
    }

    /// Access function for the fluid Reynolds number
    const double& re() const
    {
      return *Re_pt;
    }

    ///  Access function for the pointer to the fluid Strouhal number
    double*& st_pt()
    {
      return St_pt;
    }

    /// Access function for the fluid Strouhal number
    const double& st() const
    {
      return *St_pt;
    }

    ///  Access function for pointer to the fluid inverse Froude number
    /// (dimensionless gravitational loading)
    double*& re_invfr_pt()
    {
      return ReInvFr_pt;
    }

    /// Access to the fluid inverse Froude number
    const double& re_invfr()
    {
      return *ReInvFr_pt;
    }

    ///  Access function for the pointer to the density ratio
    double*& density_ratio_pt()
    {
      return Density_ratio_pt;
    }

    /// Access function for the the density ratio
    const double& density_ratio() const
    {
      return *Density_ratio_pt;
    }

  protected:
    ///  Helper function to adjust the position in
    /// response to changes in position and angle of the solid
    /// about the centre of mass
    inline void apply_rigid_body_motion(const unsigned& t,
                                        const Vector<double>& initial_x,
                                        Vector<double>& r) const
    {
      // Scale relative to the centre of mass
      double X = initial_x[0] - Initial_centre_of_mass[0];
      double Y = initial_x[1] - Initial_centre_of_mass[1];

      // Find the original angle and radius
      double phi_orig = atan2(Y, X);
      double r_orig = sqrt(X * X + Y * Y);

      // Updated position vector
      r[0] = Initial_centre_of_mass[0] +
             this->Centre_displacement_data_pt->value(t, 0);

      r[1] = Initial_centre_of_mass[1] +
             this->Centre_displacement_data_pt->value(t, 1);

      // Add in the rotation terms if there are to be included
      if (Include_geometric_rotation)
      {
        r[0] += r_orig *
                cos(phi_orig + this->Centre_displacement_data_pt->value(t, 2));
        r[1] += r_orig *
                sin(phi_orig + this->Centre_displacement_data_pt->value(t, 2));
      }
      else
      {
        r[0] += r_orig * cos(phi_orig);
        r[1] += r_orig * sin(phi_orig);
      }
    }


  private:
    ///  Return the equation number associated with the i-th
    /// centre of gravity displacment
    /// 0: x-displ; 1: y-displ; 2: rotation angle.
    inline int centre_displacement_local_eqn(const unsigned& i)
    {
      if (Displacement_data_is_internal)
      {
        return this->internal_local_eqn(Index_for_centre_displacement, i);
      }
      else
      {
        return this->external_local_eqn(Index_for_centre_displacement, i);
      }
    }

    ///  Initialisation function
    void initialise(TimeStepper* const& time_stepper_pt);

    /// Get residuals and/or Jacobian
    void get_residuals_rigid_body_generic(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian,
                                          const bool& flag);

    /// Storage for the external data that is formed from hijacked data
    /// that must be deleted by this element
    std::list<unsigned> List_of_external_hijacked_data;

    /// Delete the storage for the external data formed from hijacked data
    void delete_external_hijacked_data()
    {
      for (std::list<unsigned>::iterator it =
             List_of_external_hijacked_data.begin();
           it != List_of_external_hijacked_data.end();
           ++it)
      {
        // Delete the associated external data
        delete this->external_data_pt(*it);
        this->external_data_pt(*it) = 0;
      }
      // Now clear the list
      List_of_external_hijacked_data.clear();
    }

  protected:
    /// X-coordinate of initial centre of gravity
    Vector<double> Initial_centre_of_mass;

    /// Original rotation angle
    double Initial_Phi;

    // Mass of body
    double Mass;

    /// Polar moment of inertia of body
    double Moment_of_inertia;

    ///  Data for centre of gravity displacement.
    /// Values: 0: x-displ; 1: y-displ; 2: rotation angle.
    Data* Centre_displacement_data_pt;

  private:
    /// Underlying geometric object
    GeomObject* Geom_object_pt;

    ///  Function pointer to function that specifies
    /// external force
    ExternalForceFctPt External_force_fct_pt;

    ///  Function pointer to function that specifies
    /// external torque
    ExternalTorqueFctPt External_torque_fct_pt;

    ///  Mesh containing face elements that allow the computation of
    /// the drag on the body
    Mesh* Drag_mesh_pt;

    /// The direction of gravity
    Vector<double>* G_pt;

    /// Reynolds number of external fluid
    double* Re_pt;

    /// Strouhal number of external fluid
    double* St_pt;

    /// Reynolds number divided by Froude number of external fluid
    double* ReInvFr_pt;

    /// Density ratio of the solid to the external fluid
    double* Density_ratio_pt;

    /// Static default value for physical constants
    static double Default_Physical_Constant_Value;

    /// Static default value for physical ratios
    static double Default_Physical_Ratio_Value;

    /// Static default value for gravity
    static Vector<double> Default_Gravity_vector;

    ///  Index for the data (internal or external) that contains the
    /// centre-of-gravity displacement
    unsigned Index_for_centre_displacement;

    /// Boolean flag to indicate whether data is internal
    bool Displacement_data_is_internal;

    /// Boolean to indicate that the rotation variable does not affect the
    /// boundary shape
    bool Include_geometric_rotation;
  };


  //=====================================================================
  /// Class upgrading a TriangleMeshPolygon to a "hole" for use during
  /// triangle mesh generation. For mesh generation purposes, the main (and
  /// only) addition to the base class is the provision of the coordinates of a
  /// hole inside the polygon. To faciliate the movement of the "hole" through
  /// the domain we also provide a Data object whose three values represent the
  /// x and y displacements of its centre of gravity and the polygon's rotation
  /// about its centre of gravity. If added to a mesh in the Problem (in its
  /// incarnation as a GeneralisedElement) the displacement/rotation of the
  /// polygon is computed in response to (i) user-specifiable applied forces and
  /// a torque and (ii) the net drag (and associated torque) from a mesh of
  /// elements that can exert a drag onto the polygon (typically Navier-Stokes
  /// FaceElements that apply a viscous drag to an immersed body, represented by
  /// the polygon.)
  //=====================================================================
  class ImmersedRigidBodyTriangleMeshPolygon : public TriangleMeshPolygon,
                                               public ImmersedRigidBodyElement
  {
  public:
    ///  Constructor: Specify coordinates of a point inside the hole
    /// and a vector of pointers to TriangleMeshPolyLines
    /// that define the boundary segments of the polygon.
    /// Each TriangleMeshPolyLine has its own boundary ID and can contain
    /// multiple (straight-line) segments. The optional final argument
    /// is a pointer to a Data object whose three values represent
    /// the two displacements of and the rotation angle about the polygon's
    /// centre of mass.
    ImmersedRigidBodyTriangleMeshPolygon(
      const Vector<double>& hole_center,
      const Vector<TriangleMeshCurveSection*>& boundary_polyline_pt,
      TimeStepper* const& time_stepper_pt,
      Data* const& centre_displacement_data_pt = 0);

    ///  Empty Destuctor
    ~ImmersedRigidBodyTriangleMeshPolygon() {}

    /// Overload (again) the position to apply the rotation and translation
    void position(const Vector<double>& xi, Vector<double>& r) const
    {
      Vector<double> initial_x(2);
      this->get_initial_position(xi, initial_x);
      this->apply_rigid_body_motion(0, initial_x, r);
    }


    /// Overload (again) the position to apply the rotation and translation
    void position(const unsigned& t,
                  const Vector<double>& xi,
                  Vector<double>& r) const
    {
      Vector<double> initial_x(2);
      this->get_initial_position(xi, initial_x);
      this->apply_rigid_body_motion(t, initial_x, r);
    }


    ///  Update the reference configuration by re-setting the original
    /// position of the vertices to their current ones, re-set the
    /// original position of the centre of mass, and the displacements
    /// and rotations relative to it
    void reset_reference_configuration();

  private:
    ///  Get the initial position of the polygon
    void get_initial_position(const Vector<double>& xi, Vector<double>& r) const
    {
      // Find the number of polylines (boundaries)
      unsigned n_poly = this->npolyline();

      // The boundary coordinate will be contiguous from polyline to
      // polyline and each polyline will have the scaled arclength coordinate
      // in the range 0->1.

      // Find the maximum coordinate
      double zeta_max = Zeta_vertex[n_poly - 1].back();

      // If we are above the maximum complain
      if (xi[0] > zeta_max)
      {
        std::ostringstream error_message;
        error_message << "Value of intrinsic coordinate " << xi[0]
                      << "greater than maximum " << zeta_max << "\n";
        throw OomphLibError(error_message.str(),
                            "TriangleMeshPolygon::position()",
                            OOMPH_EXCEPTION_LOCATION);
      }

      // If equal to the maximum, return the final vertex
      if (xi[0] == zeta_max)
      {
        unsigned n_vertex = this->polyline_pt(n_poly - 1)->nvertex();
        r = this->polyline_pt(n_poly - 1)->vertex_coordinate(n_vertex - 1);
        return;
      }

      // Otherwise

      // Find out which polyline we are in
      // If we've got here this should be less than n_poly
      unsigned p = static_cast<unsigned>(floor(xi[0]));

      // If we are "above" the last polyline then throw an error
      // This should have been caught by the trap above
      if (p >= n_poly)
      {
        std::ostringstream error_message;
        error_message
          << "Something has gone wrong.\n"
          << "The integer part of the input intrinsic coordinate is " << p
          << "\nwhich is equal to or greater than the number of polylines, "
          << n_poly << ".\n"
          << "This should have triggered an earlier error\n";


        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Cache the appropriate polyline
      TriangleMeshPolyLine* const line_pt = this->polyline_pt(p);

      // If we are at the first vertex in the polyline, return it
      if (xi[0] == Zeta_vertex[p][0])
      {
        r = line_pt->vertex_coordinate(0);
        return;
      }

      // Otherwise loop over the other points to find the appropriate
      // segment

      // Find the number of vertices in the chosen polyline
      unsigned n_vertex = line_pt->nvertex();
      // Now start from the first node in the appropriate polyline and loop up
      for (unsigned v = 1; v < n_vertex; v++)
      {
        // First time that zeta falls below the vertex coordinate
        // we have something
        if (xi[0] < Zeta_vertex[p][v])
        {
          double fraction = (xi[0] - Zeta_vertex[p][v - 1]) /
                            (Zeta_vertex[p][v] - Zeta_vertex[p][v - 1]);
          Vector<double> first = line_pt->vertex_coordinate(v - 1);
          Vector<double> last = line_pt->vertex_coordinate(v);
          r.resize(2);
          for (unsigned i = 0; i < 2; i++)
          {
            r[i] = first[i] + fraction * (last[i] - first[i]);
          }
          return;
        }
      }
    }


    /// Helper function to assign the values of the (scaled) arc-length
    /// to each node of each polyline. The direction will be the natural
    /// order of the vertices within the polyline.
    void assign_zeta()
    {
      // Find the number of polylines
      unsigned n_poly = this->npolyline();

      // Allocate appropriate storage for the zeta values
      Zeta_vertex.resize(n_poly);

      // Temporary storage for the vertex coordinates
      Vector<double> vertex_coord_first;
      Vector<double> vertex_coord_next;

      // Set the initial value of zeta
      double zeta_offset = 0.0;

      // Loop over the polylines
      for (unsigned p = 0; p < n_poly; ++p)
      {
        // Cache the pointer to the polyline
        TriangleMeshPolyLine* const line_pt = this->polyline_pt(p);

        // Find the number of vertices in the polyline
        unsigned n_vertex = line_pt->nvertex();

        // Allocate storage and set initial value
        Zeta_vertex[p].resize(n_vertex);
        Zeta_vertex[p][0] = 0.0;

        // Loop over the vertices in the polyline and calculate the length
        // between each for use as the intrinsic coordinate
        vertex_coord_first = line_pt->vertex_coordinate(0);
        for (unsigned v = 1; v < n_vertex; v++)
        {
          vertex_coord_next = line_pt->vertex_coordinate(v);
          double length =
            sqrt(pow(vertex_coord_next[0] - vertex_coord_first[0], 2.0) +
                 pow(vertex_coord_next[1] - vertex_coord_first[1], 2.0));
          Zeta_vertex[p][v] = Zeta_vertex[p][v - 1] + length;
          vertex_coord_first = vertex_coord_next;
        }

        // Now scale the length to unity and add the offset
        Zeta_vertex[p][0] += zeta_offset;
        for (unsigned v = 1; v < n_vertex; v++)
        {
          Zeta_vertex[p][v] /= Zeta_vertex[p][n_vertex - 1];
          Zeta_vertex[p][v] += zeta_offset;
        }
        zeta_offset += 1.0;
      } // End of loop over polylines
    }

    ///  Vector of intrisic coordinate values at the nodes
    Vector<Vector<double>> Zeta_vertex;
  };

} // namespace oomph
#endif
