//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
//Header file for triangular/tetrahedaral Navier Stokes elements

#ifndef OOMPH_TNAVIER_STOKES_ELEMENTS_HEADER
#define OOMPH_TNAVIER_STOKES_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


//OOMPH-LIB headers 
#include "../generic/Telements.h"
#include "navier_stokes_elements.h"

namespace oomph
{

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// NOTE: TRI/TET CROZIER RAVIARTS REQUIRE BUBBLE FUNCTIONS! THEY'RE NOT
// STRAIGHTFORWARD GENERALISATIONS OF THE Q-EQUIVALENTS (WHICH ARE
// LBB UNSTABLE!)
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


/* //========================================================================== */
/* ///TCrouzeix_Raviart elements are Navier--Stokes elements with quadratic */
/* ///interpolation for velocities and positions, but a discontinuous linear */
/* ///pressure interpolation */
/* //========================================================================== */
/* template <unsigned DIM> */
/* class TCrouzeixRaviartElement : public virtual TElement<DIM,3>,  */
/*  public virtual NavierStokesEquations<DIM> */
/* { */
/*   private: */

/*  /// Static array of ints to hold required number of variables at nodes */
/*  static const unsigned Initial_Nvalue[]; */
 
/*   protected: */
 
/*  /// \short Velocity shape and test functions and their derivs  */
/*  /// w.r.t. to global coords  at local coordinate s (taken from geometry) */
/*  ///Return Jacobian of mapping between local and global coordinates. */
/*  inline double dshape_and_dtest_eulerian_nst(const Vector<double> &s, 
    Shape &psi,  */
/*                                          DShape &dpsidx, Shape &test,  */
/*                                          DShape &dtestdx) const; */

/*  /// \short Velocity shape and test functions and their derivs  */
/*  /// w.r.t. to global coords at ipt-th integation point (taken from geometry) */
/*  ///Return Jacobian of mapping between local and global coordinates. */
/*  inline double dshape_and_dtest_eulerian_at_knot_nst(const unsigned &ipt,  */
/*                                                  Shape &psi,  */
/*                                                  DShape &dpsidx, Shape &test,  */
/*                                                  DShape &dtestdx) const; */

/*  /// Pressure shape functions at local coordinate s */
/*  inline void pshape_nst(const Vector<double> &s, Shape &psi) const; */

/*  /// Pressure shape and test functions at local coordinte s */
/*  inline void pshape_nst(const Vector<double> &s, Shape &psi, Shape &test) const; */

/*  /// Unpin all internal pressure dofs */
/*  void unpin_all_internal_pressure_dofs(); */

/*  /// Return the local equation numbers for the pressure values. */
/*  inline int p_local_eqn(const unsigned &n) */
/*   {return this->internal_local_eqn(n,0);} */

/* public: */

/*  /// Constructor, there are DIM+1 internal values (for the pressure) */
/*  TCrouzeixRaviartElement() : TElement<DIM,3>(), NavierStokesEquations<DIM>() */
/*   { */
/*    //Allocate and add DIM+1 Internal data objects, each storing one */
/*    //unknown (the pressure) */
/*    for(unsigned i=0;i<DIM+1;i++) {this->add_internal_data(new Data(1));} */
/*   } */
 
/*  /// \short Number of values (pinned or dofs) required at local node n.  */
/*  virtual unsigned required_nvalue(const unsigned &n) const; */

 
/*  /// \short Return the pressure values at internal dof i_internal */
/*  /// (Discontinous pressure interpolation -- no need to cater for hanging  */
/*  /// nodes).  */
/*  double p(const unsigned &i_internal) const */
/*   {return *(this->internal_data_pt(i_internal)->value_pt(0));} */

/*  /// Return number of pressure values */
/*  unsigned npres() const {return DIM+1;}  */


/*  /// \short Set the integers that are used to label velocities and pressures */
/*  /// to the data structures stored at the nodes or in internal data. The */
/*  /// arguments are an Vector of DIM integers representing the velocity  */
/*  /// components and an integer to represent the fluid pressure. */
/*  void pass_fluid_variable_identifiers(const Vector<int>& u, const int &p) */
/*   { */
/*    //Get the number of nodes */
/*    unsigned n_node = FiniteElement::nnode(); */
/*    //Loop over the nodes */
/*    for(unsigned l=0;l<n_node;l++) */
/*     { */
/*      //Loop over the dimensions */
/*      for(unsigned i=0;i<DIM;i++) */
/*       { */
/*        //The type u[i] is stored at data position i */
/*        this->Node_pt[l]->set_identifier(i,u[i]); */
/*       } */
/*     } */
     
/*    //Get the number of pressure dofs */
/*    unsigned n_pres = npres(); */
/*    //Loop over the internal dofs */
/*    for(unsigned l=0;l<n_pres;l++) */
/*     { */
/*      //The pressure is the first entry in each internal dof */
/*      this->internal_data_pt(l)->set_identifier(0,p); */
/*     } */
/*   } */
 
/*  /// Pin p_dof-th pressure dof and set it to value specified by p_value. */
/*  void fix_pressure(const unsigned &p_dof, const double &p_value) */
/*   { */
/*    this->internal_data_pt(p_dof)->pin(0); */
/*    *(this->internal_data_pt(p_dof)->value_pt(0)) = p_value; */
/*   } */

/*  /// \short Add to the set paired_load_data */
/*  /// pairs of pointers to data objects and unsignedegers that */
/*  /// index the values in the data object that affect the load (traction),  */
/*  /// as specified in the get_load() function.  */
/*  void identify_load_data(set<pair<Data*,unsigned> > &paired_load_data); */

/*  /// Redirect output to NavierStokesEquations output */
/*  void output(std::ostream &outfile) {NavierStokesEquations<DIM>::output(outfile);} */

/*  /// Redirect output to NavierStokesEquations output */
/*  void output(std::ostream &outfile, const unsigned &nplot) */
/*   {NavierStokesEquations<DIM>::output(outfile,nplot);} */

/*  /// Redirect output to NavierStokesEquations output */
/*  void output(FILE* file_pt) {NavierStokesEquations<DIM>::output(file_pt);} */

/*  /// Redirect output to NavierStokesEquations output */
/*  void output(FILE* file_pt, const unsigned &n_plot) */
/*   {NavierStokesEquations<DIM>::output(file_pt,n_plot);} */


/*  /// \short Full output function:  */
/*  /// x,y,[z],u,v,[w],p,du/dt,dv/dt,[dw/dt],dissipation */
/*  /// in tecplot format. Default number of plot points */
/*  void full_output(std::ostream &outfile) */
/*   {NavierStokesEquations<DIM>::full_output(outfile);} */

/*  /// \short Full output function:  */
/*  /// x,y,[z],u,v,[w],p,du/dt,dv/dt,[dw/dt],dissipation */
/*  /// in tecplot format. nplot points in each coordinate direction */
/*  void full_output(std::ostream &outfile, const unsigned &nplot) */
/*   {NavierStokesEquations<DIM>::full_output(outfile,nplot);} */

/* }; */

/* //Inline functions */

/* //======================================================================= */
/* /// 2D */
/* /// Derivatives of the shape functions and test functions w.r.t. to global */
/* /// (Eulerian) coordinates. Return Jacobian of mapping between */
/* /// local and global coordinates. */
/* //======================================================================= */
/* template<> */
/* inline double TCrouzeixRaviartElement<2>::dshape_and_dtest_eulerian( */
/*                                   const Vector<double> &s, Shape &psi,  */
/*                                   DShape &dpsidx, Shape &test,  */
/*                                   DShape &dtestdx) const */
/* { */
/*  //Call the geometrical shape functions and derivatives   */
/*  double J = this->dshape_eulerian(s,psi,dpsidx); */
/*  //Loop over the test functions and derivatives and set them equal to the */
/*  //shape functions */
/*  for(unsigned i=0;i<6;i++) */
/*   { */
/*    test[i] = psi[i];  */
/*    dtestdx(i,0) = dpsidx(i,0); */
/*    dtestdx(i,1) = dpsidx(i,1); */
/*   } */
/*  //Return the jacobian */
/*  return J; */
/* } */


/* //======================================================================= */
/* /// 3D */
/* /// Derivatives of the shape functions and test functions w.r.t. to global */
/* /// (Eulerian) coordinates. Return Jacobian of mapping between */
/* /// local and global coordinates. */
/* //======================================================================= */
/* template<> */
/* inline double TCrouzeixRaviartElement<3>::dshape_and_dtest_eulerian( */
/*                                   const Vector<double> &s, Shape &psi,  */
/*                                   DShape &dpsidx, Shape &test,  */
/*                                   DShape &dtestdx) const */
/* { */
/*  //Call the geometrical shape functions and derivatives   */
/*  double J = this->dshape_eulerian(s,psi,dpsidx); */
/*  //Loop over the test functions and derivatives and set them equal to the */
/*  //shape functions */
/*  for(unsigned i=0;i<10;i++) */
/*   { */
/*    test[i] = psi[i];  */
/*    dtestdx(i,0) = dpsidx(i,0); */
/*    dtestdx(i,1) = dpsidx(i,1); */
/*    dtestdx(i,2) = dpsidx(i,2); */
/*   } */
/*  //Return the jacobian */
/*  return J; */
/* } */



/* //======================================================================= */
/* /// 2D */
/* /// Derivatives of the shape functions and test functions w.r.t. to global */
/* /// (Eulerian) coordinates. Return Jacobian of mapping between */
/* /// local and global coordinates. */
/* //======================================================================= */
/* template<> */
/* inline double TCrouzeixRaviartElement<2>::dshape_and_dtest_eulerian_at_knot( */
/*                                               const unsigned &ipt, Shape &psi,  */
/*                                               DShape &dpsidx, Shape &test,  */
/*                                               DShape &dtestdx) const */
/* { */
/*  //Call the geometrical shape functions and derivatives   */
/*  double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx); */
/*  //Loop over the test functions and derivatives and set them equal to the */
/*  //shape functions */
/*  for(unsigned i=0;i<6;i++) */
/*   { */
/*    test[i] = psi[i];  */
/*    dtestdx(i,0) = dpsidx(i,0); */
/*    dtestdx(i,1) = dpsidx(i,1); */
/*   } */
/*  //Return the jacobian */
/*  return J; */
/* } */


/* //======================================================================= */
/* /// 3D */
/* /// Derivatives of the shape functions and test functions w.r.t. to global */
/* /// (Eulerian) coordinates. Return Jacobian of mapping between */
/* /// local and global coordinates. */
/* //======================================================================= */
/* template<> */
/* inline double TCrouzeixRaviartElement<3>::dshape_and_dtest_eulerian_at_knot( */
/*                                               const unsigned &ipt, Shape &psi,  */
/*                                               DShape &dpsidx, Shape &test,  */
/*                                               DShape &dtestdx) const */
/* { */
/*  //Call the geometrical shape functions and derivatives   */
/*  double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx); */
/*  //Loop over the test functions and derivatives and set them equal to the */
/*  //shape functions */
/*  for(unsigned i=0;i<10;i++) */
/*   { */
/*    test[i] = psi[i];  */
/*    dtestdx(i,0) = dpsidx(i,0); */
/*    dtestdx(i,1) = dpsidx(i,1); */
/*    dtestdx(i,2) = dpsidx(i,2); */
/*   } */
/*  //Return the jacobian */
/*  return J; */
/* } */



/* //======================================================================= */
/* /// 2D : */
/* /// Pressure shape functions */
/* //======================================================================= */
/* template<> */
/* inline void TCrouzeixRaviartElement<2>::pshape(const Vector<double> &s,  */
/*                                               Shape &psi) const */
/* { */
/*  psi[0] = 1.0; */
/*  psi[1] = s[0]; */
/*  psi[2] = s[1]; */
/* } */

/* //======================================================================= */
/* /// Pressure shape and test functions */
/* //======================================================================= */
/* template<> */
/* inline void TCrouzeixRaviartElement<2>::pshape(const Vector<double> &s,  */
/*                                               Shape &psi,  */
/*                                               Shape &test) const */
/* { */
/*  //Call the pressure shape functions */
/*  pshape(s,psi); */
/*  //Loop over the test functions and set them equal to the shape functions */
/*  for(unsigned i=0;i<3;i++) test[i] = psi[i]; */
/* } */


/* //======================================================================= */
/* /// 3D : */
/* /// Pressure shape functions */
/* //======================================================================= */
/* template<> */
/* inline void TCrouzeixRaviartElement<3>::pshape(const Vector<double> &s, */
/*                                               Shape &psi) */
/* const */
/* { */
/*  psi[0] = 1.0; */
/*  psi[1] = s[0]; */
/*  psi[2] = s[1]; */
/*  psi[3] = s[2]; */
/* } */


/* //======================================================================= */
/* /// Pressure shape and test functions */
/* //======================================================================= */
/* template<> */
/* inline void TCrouzeixRaviartElement<3>::pshape(const Vector<double> &s, */
/*                                               Shape &psi,  */
/*                                               Shape &test) const */
/* { */
/*  //Call the pressure shape functions */
/*  pshape(s,psi); */
/*  //Loop over the test functions and set them equal to the shape functions */
/*  for(unsigned i=0;i<4;i++) test[i] = psi[i]; */
/* } */


/* //======================================================================= */
/* /// Face geometry of the 2D Crouzeix_Raviart elements */
/* //======================================================================= */
/* template<> */
/* class FaceGeometry<TCrouzeixRaviartElement<2> >: public virtual QElement<1,3> */
/* { */
/*   public: */
/*  FaceGeometry() : QElement<1,3>() */
/*   { */
/*    oomph_info << "Careful: FaceGeometries for TCrouzeixRaviart have not been tested" */
/*         << std::endl; */
/*   } */
/* }; */

/* //======================================================================= */
/* /// Face geometry of the 3D Crouzeix_Raviart elements */
/* //======================================================================= */
/* template<> */
/* class FaceGeometry<TCrouzeixRaviartElement<3> >: public virtual TElement<2,3> */
/* { */
 
/*   public: */
/*  FaceGeometry() : TElement<2,3>() */
/*   { */
/*    oomph_info << "Careful: FaceGeometries for TCrouzeixRaviart have not been tested" */
/*         << std::endl; */
/*   } */
/* }; */


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////



//=======================================================================
/// Taylor--Hood elements are Navier--Stokes elements
/// with quadratic interpolation for velocities and positions and
/// continous linear pressure interpolation
//=======================================================================
template <unsigned DIM>
class TTaylorHoodElement : public virtual TElement<DIM,3>,
 public virtual NavierStokesEquations<DIM>
{
  private:
 
 /// Static array of ints to hold number of variables at node
 static const unsigned Initial_Nvalue[];
 
  protected:

 /// \short Static array of ints to hold conversion from pressure
 /// node numbers to actual node numbers
 static const unsigned Pconv[];
 
 /// \short Velocity shape and test functions and their derivs
 /// w.r.t. to global coords  at local coordinate s (taken from geometry)
 /// Return Jacobian of mapping between local and global coordinates.
 inline double dshape_and_dtest_eulerian_nst(const Vector<double> &s, 
                                             Shape &psi,
                                             DShape &dpsidx, Shape &test,
                                             DShape &dtestdx) const;

 /// \short Velocity shape and test functions and their derivs
 /// w.r.t. to global coords  at local coordinate s (taken from geometry)
 /// Return Jacobian of mapping between local and global coordinates.
 inline double dshape_and_dtest_eulerian_at_knot_nst(const unsigned &ipt,
                                                     Shape &psi,
                                                     DShape &dpsidx, 
                                                     Shape &test,
                                                     DShape &dtestdx) const;
 

 /// \short Compute the pressure shape and test functions and derivatives 
 /// w.r.t. global coords at local coordinate s.
 /// Return Jacobian of mapping between local and global coordinates.
 virtual double dpshape_and_dptest_eulerian_nst(const Vector<double> &s, 
                                                Shape &ppsi, 
                                                DShape &dppsidx, 
                                                Shape &ptest, 
                                                DShape &dptestdx) const;

 /// Unpin all pressure dofs
 void unpin_all_nodal_pressure_dofs();

 ///  Pin all nodal pressure dofs
 void pin_all_nodal_pressure_dofs();
   
 ///  Unpin the proper nodal pressure dofs
 void unpin_proper_nodal_pressure_dofs();


public:

 /// Constructor, no internal data points
 TTaylorHoodElement() : TElement<DIM,3>(),  NavierStokesEquations<DIM>() {}
 

 /// Broken copy constructor
 TTaylorHoodElement(const TTaylorHoodElement<DIM>& dummy) 
  { 
   BrokenCopy::broken_copy("TTaylorHoodElement");
  } 
 
 /// Broken assignment operator
 void operator=(const TTaylorHoodElement<DIM>&) 
  {
   BrokenCopy::broken_assign("TTaylorHoodElement");
  }

 /// \short Number of values (pinned or dofs) required at node n. Can
 /// be overwritten for hanging node version
 inline virtual unsigned required_nvalue(const unsigned &n) const
  {return Initial_Nvalue[n];}

 /// Test whether the pressure dof p_dof hanging or not?
 //bool pressure_dof_is_hanging(const unsigned& p_dof)
 // {return this->node_pt(Pconv[p_dof])->is_hanging(DIM);}


 /// Pressure shape functions at local coordinate s
 inline void pshape_nst(const Vector<double> &s, Shape &psi) const;

 /// Pressure shape and test functions at local coordinte s
 inline void pshape_nst(const Vector<double> &s, Shape &psi, 
                        Shape &test) const;

 /// \short Which nodal value represents the pressure?
 unsigned p_index_nst() {return DIM;}

 /// \short Pointer to n_p-th pressure node
 //Node* pressure_node_pt(const unsigned &n_p)
 //{return this->Node_pt[Pconv[n_p]];}

 /// Return the local equation numbers for the pressure values.
 inline int p_local_eqn(const unsigned &n)
  {return this->nodal_local_eqn(Pconv[n],DIM);}

 /// \short Access function for the pressure values at local pressure
 /// node n_p (const version)
 double p_nst(const unsigned &n_p) const
  {return this->nodal_value(Pconv[n_p],DIM);}
 
 /// \short Set the value at which the pressure is stored in the nodes
 int p_nodal_index_nst() const {return static_cast<int>(DIM);}

 /// Return number of pressure values
 unsigned npres_nst() const;

 /// Pin p_dof-th pressure dof and set it to value specified by p_value.
 void fix_pressure(const unsigned &p_dof, const double &p_value)
  {
   this->node_pt(Pconv[p_dof])->pin(DIM);
   this->node_pt(Pconv[p_dof])->set_value(DIM,p_value);
  }


 /// \short Build FaceElements that apply the Robin boundary condition
 /// to the pressure advection diffusion problem required by 
 /// Fp preconditioner
 void build_fp_press_adv_diff_robin_bc_element(const unsigned& 
                                               face_index)
 {
  this->Pressure_advection_diffusion_robin_element_pt.push_back(
   new FpPressureAdvDiffRobinBCElement<TTaylorHoodElement<DIM> >(
    this, face_index));
 }

 /// \short Add to the set \c paired_load_data pairs containing
 /// - the pointer to a Data object
 /// and
 /// - the index of the value in that Data object
 /// .
 /// for all values (pressures, velocities) that affect the
 /// load computed in the \c get_load(...) function. 
 void identify_load_data(
  std::set<std::pair<Data*,unsigned> > &paired_load_data);

 /// \short  Add to the set \c paired_pressure_data pairs 
 /// containing
 /// - the pointer to a Data object
 /// and
 /// - the index of the value in that Data object
 /// .
 /// for all pressure values that affect the
 /// load computed in the \c get_load(...) function.
 void identify_pressure_data(
  std::set<std::pair<Data*,unsigned> > &paired_pressure_data);

 /// Redirect output to NavierStokesEquations output
 void output(std::ostream &outfile) 
  {NavierStokesEquations<DIM>::output(outfile);}

 /// Redirect output to NavierStokesEquations output
 void output(std::ostream &outfile, const unsigned &nplot)
  {NavierStokesEquations<DIM>::output(outfile,nplot);}

 /// Redirect output to NavierStokesEquations output
 void output(FILE* file_pt) {NavierStokesEquations<DIM>::output(file_pt);}

 /// Redirect output to NavierStokesEquations output
 void output(FILE* file_pt, const unsigned &n_plot)
  {NavierStokesEquations<DIM>::output(file_pt,n_plot);}

 /// \short The number of "blocks" that degrees of freedom in this element
 /// are sub-divided into: Velocity and pressure.
 unsigned ndof_types()
  {
   return DIM+1;
  }
 
 /// \short Create a list of pairs for all unknowns in this element,
 /// so that the first entry in each pair contains the global equation
 /// number of the unknown, while the second one contains the number
 /// of the "block" that this unknown is associated with.
 /// (Function can obviously only be called if the equation numbering
 /// scheme has been set up.) Velocity=0; Pressure=1
 void get_dof_numbers_for_unknowns(
  std::list<std::pair<unsigned long,unsigned> >& block_lookup_list)
  {
   // number of nodes
   unsigned n_node = this->nnode();
   
   // temporary pair (used to store block lookup prior to being added to list)
   std::pair<unsigned,unsigned> block_lookup;
   
   // loop over the nodes
   for (unsigned n = 0; n < n_node; n++)
    {
     // find the number of values at this node
     unsigned nv = this->node_pt(n)->nvalue(); 

     //loop over these values
     for (unsigned v = 0; v < nv; v++)
      {
       // determine local eqn number
       int local_eqn_number = this->nodal_local_eqn(n, v);
       
       // ignore pinned values - far away degrees of freedom resulting 
       // from hanging nodes can be ignored since these are be dealt
       // with by the element containing their master nodes
       if (local_eqn_number >= 0)
        {
         // store block lookup in temporary pair: Global equation number
         // is the first entry in pair
         block_lookup.first = this->eqn_number(local_eqn_number);
         
         // set block numbers: Block number is the second entry in pair
         block_lookup.second = v;
         
         // add to list
         block_lookup_list.push_front(block_lookup);
        }
      }
    }
  }
};




//Inline functions

//==========================================================================
/// 2D :
/// Number of pressure values
//==========================================================================
template<>
inline unsigned TTaylorHoodElement<2>::npres_nst() const
{
 return 3;
}

//==========================================================================
/// 3D :
/// Number of pressure values
//==========================================================================
template<>
inline unsigned TTaylorHoodElement<3>::npres_nst() const
{
 return 4;
}



//==========================================================================
/// 2D :
/// Derivatives of the shape functions and test functions w.r.t to
/// global (Eulerian) coordinates. Return Jacobian of mapping between
/// local and global coordinates.
//==========================================================================
template<>
inline double TTaylorHoodElement<2>::dshape_and_dtest_eulerian_nst(
                                              const Vector<double> &s,
                                              Shape &psi,
                                              DShape &dpsidx, Shape &test,
                                              DShape &dtestdx) const
{
 //Call the geometrical shape functions and derivatives
 double J = this->dshape_eulerian(s,psi,dpsidx);
 //Loop over the test functions and derivatives and set them equal to the
 //shape functions
 for(unsigned i=0;i<6;i++)
  {
   test[i] = psi[i];
   dtestdx(i,0) = dpsidx(i,0);
   dtestdx(i,1) = dpsidx(i,1);
  }
 //Return the jacobian
 return J;
}


//==========================================================================
/// 3D :
/// Derivatives of the shape functions and test functions w.r.t to
/// global (Eulerian) coordinates. Return Jacobian of mapping between
/// local and global coordinates.
//==========================================================================
template<>
inline double TTaylorHoodElement<3>::dshape_and_dtest_eulerian_nst(
                                              const Vector<double> &s,
                                              Shape &psi,
                                              DShape &dpsidx, Shape &test,
                                              DShape &dtestdx) const
{
 //Call the geometrical shape functions and derivatives
 double J = this->dshape_eulerian(s,psi,dpsidx);
 //Loop over the test functions and derivatives and set them equal to the
 //shape functions
 for(unsigned i=0;i<10;i++)
  {
   test[i] = psi[i];
   dtestdx(i,0) = dpsidx(i,0);
   dtestdx(i,1) = dpsidx(i,1);
   dtestdx(i,2) = dpsidx(i,2);
  }
 //Return the jacobian
 return J;
}


//==========================================================================
/// 2D :
/// Derivatives of the shape functions and test functions w.r.t to
/// global (Eulerian) coordinates. Return Jacobian of mapping between
/// local and global coordinates.
//==========================================================================
template<>
inline double TTaylorHoodElement<2>::dshape_and_dtest_eulerian_at_knot_nst(
 const unsigned &ipt,Shape &psi, DShape &dpsidx, Shape &test,
 DShape &dtestdx) const
{
 //Call the geometrical shape functions and derivatives
 double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);
 //Loop over the test functions and derivatives and set them equal to the
 //shape functions
 for(unsigned i=0;i<6;i++)
  {
   test[i] = psi[i];
   dtestdx(i,0) = dpsidx(i,0);
   dtestdx(i,1) = dpsidx(i,1);
  }
 //Return the jacobian
 return J;
}


//==========================================================================
/// 3D :
/// Derivatives of the shape functions and test functions w.r.t to
/// global (Eulerian) coordinates. Return Jacobian of mapping between
/// local and global coordinates.
//==========================================================================
template<>
inline double TTaylorHoodElement<3>::dshape_and_dtest_eulerian_at_knot_nst(
 const unsigned &ipt,Shape &psi, DShape &dpsidx, Shape &test,
 DShape &dtestdx) const
{
 //Call the geometrical shape functions and derivatives
 double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);
 //Loop over the test functions and derivatives and set them equal to the
 //shape functions
 for(unsigned i=0;i<10;i++)
  {
   test[i] = psi[i];
   dtestdx(i,0) = dpsidx(i,0);
   dtestdx(i,1) = dpsidx(i,1);
   dtestdx(i,2) = dpsidx(i,2);
  }
 //Return the jacobian
 return J;
}



//==========================================================================
/// 2D :
/// Pressure shape and test functions and derivs w.r.t. to Eulerian coords.
/// Return Jacobian of mapping between local and global coordinates.
//==========================================================================
template<>
 inline double TTaylorHoodElement<2>::dpshape_and_dptest_eulerian_nst(
  const Vector<double> &s, 
  Shape &ppsi, 
  DShape &dppsidx, 
  Shape &ptest, 
  DShape &dptestdx) const
 {
  
  ppsi[0] = s[0];
  ppsi[1] = s[1];
  ppsi[2] = 1.0-s[0]-s[1];
  
  dppsidx(0,0)=1.0;
  dppsidx(0,1)=0.0;

  dppsidx(1,0)=0.0;
  dppsidx(1,1)=1.0;

  dppsidx(2,0)=-1.0;
  dppsidx(2,1)=-1.0;

    // Allocate memory for the inverse 2x2 jacobian
  DenseMatrix<double> inverse_jacobian(2);


  //Get the values of the shape functions and their local derivatives
  Shape psi(6);
  DShape dpsi(6,2);
  dshape_local(s,psi,dpsi);
    
  // Now calculate the inverse jacobian
  const double det = local_to_eulerian_mapping(dpsi,inverse_jacobian);
  
  // Now set the values of the derivatives to be derivs w.r.t. to the
  // Eulerian coordinates
  transform_derivatives(inverse_jacobian,dppsidx);
  
  // Loop over the test functions and set them equal to the shape functions
  for(unsigned i=0;i<3;i++)
   {
    ptest[i] = ppsi[i];
    dptestdx(i,0) = dppsidx(i,0);
    dptestdx(i,1) = dppsidx(i,1);
   }
  
  // Return the determinant of the jacobian
  return det;
    
 }


//==========================================================================
/// 3D :
/// Pressure shape and test functions and derivs w.r.t. to Eulerian coords.
/// Return Jacobian of mapping between local and global coordinates.
//==========================================================================
template<>
 inline double TTaylorHoodElement<3>::dpshape_and_dptest_eulerian_nst(
  const Vector<double> &s, 
  Shape &ppsi, 
  DShape &dppsidx, 
  Shape &ptest, 
  DShape &dptestdx) const
 {
 
  ppsi[0] = s[0];
  ppsi[1] = s[1];
  ppsi[2] = s[2];
  ppsi[3] = 1.0-s[0]-s[1]-s[2];
  
  dppsidx(0,0)=1.0;
  dppsidx(0,1)=0.0;
  dppsidx(0,2)=0.0;

  dppsidx(1,0)=0.0;
  dppsidx(1,1)=1.0;
  dppsidx(1,2)=0.0;

  dppsidx(2,0)=0.0;
  dppsidx(2,1)=0.0;
  dppsidx(2,2)=1.0;

  dppsidx(3,0)=-1.0;
  dppsidx(3,1)=-1.0;
  dppsidx(3,2)=-1.0;


  //Get the values of the shape functions and their local derivatives
  Shape psi(10);
  DShape dpsi(10,3);
  dshape_local(s,psi,dpsi);

  // Allocate memory for the inverse 3x3 jacobian
  DenseMatrix<double> inverse_jacobian(3);
  
  // Now calculate the inverse jacobian
  const double det = local_to_eulerian_mapping(dpsi,inverse_jacobian);
  
  // Now set the values of the derivatives to be derivs w.r.t. to the
  // Eulerian coordinates
  transform_derivatives(inverse_jacobian,dppsidx);
  
  // Loop over the test functions and set them equal to the shape functions
  for(unsigned i=0;i<4;i++)
   {
    ptest[i] = ppsi[i];
    dptestdx(i,0) = dppsidx(i,0);
    dptestdx(i,1) = dppsidx(i,1);
    dptestdx(i,2) = dppsidx(i,2);
   }
  
  // Return the determinant of the jacobian
  return det;
    
 }



//==========================================================================
/// 2D :
/// Pressure shape functions
//==========================================================================
template<>
inline void TTaylorHoodElement<2>::pshape_nst(const Vector<double> &s, Shape &psi)
const
{
 psi[0] = s[0];
 psi[1] = s[1];
 psi[2] = 1.0-s[0]-s[1];
}

//==========================================================================
/// 3D :
/// Pressure shape functions
//==========================================================================
template<>
inline void TTaylorHoodElement<3>::pshape_nst(const Vector<double> &s, Shape &psi)
const
{
 psi[0] = s[0];
 psi[1] = s[1];
 psi[2] = s[2];
 psi[3] = 1.0-s[0]-s[1]-s[2];
}


//==========================================================================
/// 2D :
/// Pressure shape and test functions
//==========================================================================
template<>
inline void TTaylorHoodElement<2>::pshape_nst(const Vector<double> &s, Shape &psi,
                                           Shape &test) const
{
 //Call the pressure shape functions
 pshape_nst(s,psi);
 //Loop over the test functions and set them equal to the shape functions
 for(unsigned i=0;i<3;i++) test[i] = psi[i];
}


//==========================================================================
/// 3D :
/// Pressure shape and test functions
//==========================================================================
template<>
inline void TTaylorHoodElement<3>::pshape_nst(const Vector<double> &s, Shape &psi,
                                           Shape &test) const
{
 //Call the pressure shape functions
 pshape_nst(s,psi);
 //Loop over the test functions and set them equal to the shape functions
 for(unsigned i=0;i<4;i++) test[i] = psi[i];
}


//=======================================================================
/// Face geometry of the 2D Taylor_Hood elements
//=======================================================================
template<>
class FaceGeometry<TTaylorHoodElement<2> >: public virtual TElement<1,3>
{
  public:

  /// Constructor: Call constructor of base
  FaceGeometry() : TElement<1,3>() {}
/*   { */
/*    throw OomphLibWarning( */
/*     "Careful: FaceGeometries for TTaylorHood have not been tested", */
/*     "FaceGeometry::FaceGeometry()", */
/*     OOMPH_EXCEPTION_LOCATION); */
/*   } */
};

//=======================================================================
/// Face geometry of the 3D Taylor_Hood elements
//=======================================================================
template<>
class FaceGeometry<TTaylorHoodElement<3> >: public virtual TElement<2,3>
{
 
  public:

 /// Constructor: Call constructor of base
  FaceGeometry() : TElement<2,3>() {}
/*   { */
/*    throw OomphLibWarning( */
/*     "Careful: FaceGeometries for TTaylorHood have not been tested", */
/*     "FaceGeometry::FaceGeometry()", */
/*     OOMPH_EXCEPTION_LOCATION); */
/*   } */
};

}

#endif
