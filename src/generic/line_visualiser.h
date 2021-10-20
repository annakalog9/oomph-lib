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
// Header file for line visualiser

// Include guard to prevent multiple inclusions of the header
#ifndef OOMPH_LINE_VISUALISER_HEADER
#define OOMPH_LINE_VISUALISER_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include <fstream>
#include <iostream>


namespace oomph
{
  //====================================================================
  ///  Class to aid visualisation of the values on a set
  /// of points. NOTE: in a distributed problem, output is only done
  /// on processor 0.
  //====================================================================
  class LineVisualiser
  {
  public:
    ///  Constructor: Pass pointer to mesh and coordinates of
    /// desired plot points: coord_vec[j][i] is the i-th coordinate of
    /// the j-th plot point. Optional final parameter specifies the
    /// maximum search radius in bin when locating the plot points.
    /// Defaults to DBL_MAX and will therefore keep searching through
    /// the entire bin structure if a point cannot be found. It's worth
    /// limiting this to the order of the size of a typical element
    /// in the mesh in parallel computations with fine meshes, otherwise
    /// the setup can take forever.
    LineVisualiser(Mesh* mesh_pt,
                   const Vector<Vector<double>>& coord_vec,
                   const double& max_search_radius = DBL_MAX)
      : Max_search_radius(max_search_radius),
        Comm_pt(mesh_pt->communicator_pt())
    {
      // Do the actual work
      setup(mesh_pt, coord_vec);
    }


    ///  Constructor reading centerline file
    /// - Open "file_name" and extract 3 first doubles of each line
    /// - Skip lines which does not begin with a number. Scaling
    /// factor allows points defined in input file to be scaled.
    LineVisualiser(Mesh* mesh_pt,
                   const std::string file_name,
                   const double& scale = 1.0)
      : Max_search_radius(DBL_MAX), Comm_pt(mesh_pt->communicator_pt())
    {
      setup_from_file(mesh_pt, file_name, scale);
    }


    ///  Constructor reading centerline file
    /// - Open "file_name" and extract 3 first doubles of each line
    /// - Skip lines which does not begin with a number. Scaling
    /// factor allows points defined in input file to be scaled.
    /// Second  parameter specifies the
    /// maximum search radius in bin when locating the plot points. It's worth
    /// setting this to the order of the size of a typical element
    /// in the mesh in parallel computations with fine meshes, otherwise
    /// the setup can take forever.
    LineVisualiser(Mesh* mesh_pt,
                   const double& max_search_radius,
                   const std::string file_name,
                   const double& scale = 1.0)
      : Max_search_radius(max_search_radius),
        Comm_pt(mesh_pt->communicator_pt())
    {
      setup_from_file(mesh_pt, file_name, scale);
    }


    ///  Output function: output each plot point.
    /// NOTE: in a distributed problem, output is only done
    /// on processor 0.
    void output(std::ostream& outfile)
    {
      // Get data in array
      Vector<Vector<double>> data(Nplot_points);
      get_output_data(data);


      // Loop over the points
      for (unsigned i = 0; i < Nplot_points; i++)
      {
        // Get the size of the line
        unsigned n = data[i].size();

        // Loop over the values on the line
        for (unsigned j = 0; j < n; j++)
        {
          outfile << data[i][j] << " ";
        }
        if (n > 0)
        {
          outfile << std::endl;
        }
      }
    }

    ///  Output data function: store data associated with each
    /// plot point in data array
    void get_output_data(Vector<Vector<double>>& data)
    {
      // Resize output data array
      data.resize(Nplot_points);

      // Check if mesh is distributed and a communication pointer
      // exists.
      if (Comm_pt != 0)
      {
        int nproc = Comm_pt->nproc();

        if (nproc > 1)
        {
#ifdef OOMPH_HAS_MPI

          // Declaration of MPI variables
          MPI_Status stat;
          int tag = 0;
          int my_rank = Comm_pt->my_rank();


          // Buffer
          unsigned buff_size;

          // Create array which contains data found in every process
          Vector<Vector<double>> vec(Nplot_points);

          // Loop over the points to fill in vec
          for (unsigned i = 0; i < Nplot_points; i++)
          {
            // Check if the point was found in the mesh
            if (Plot_point[i].first != NULL) // success
            {
              // Check if the point is halo
              if (!((*Plot_point[i].first).is_halo()))
              {
                // Get the line of output data from the element
                // (specified by .first), at its local coordinate
                // (specified by .second)
                Plot_point[i].first->point_output_data(Plot_point[i].second,
                                                       vec[i]);
              }
            }
          }


          // Analyse which plot points have been found
          // locally and concatenate the data:

          // This contains the flat-packed doubles to be sent
          // for all located plot points
          Vector<double> local_values;

          // Number of values to be sent for each plot point
          // (almost certainly the same for all plot points, but...)
          // size_values[i] gives the number of doubles to be
          // sent for plot point i.
          Vector<unsigned> size_values;

          // Each processor indicates if it has found a given plot point.
          // Once this is gathered on the root processor we know
          // exactly which data we'll receive from where.
          Vector<unsigned> tmp_proc_point_found_plus_one(Nplot_points, 0);

          // Loop over the plot points
          for (unsigned i = 0; i < Nplot_points; i++)
          {
            unsigned ndata = vec[i].size();
            if (ndata != 0)
            {
              // Store the number of fields
              size_values.push_back(ndata);

              // Update found vector
              tmp_proc_point_found_plus_one[i] = my_rank + 1;

              // Store values
              for (unsigned j = 0; j < ndata; j++)
              {
                local_values.push_back(vec[i][j]);
              }
            }
          }

          // Gather information on root

          // Find out who's found the points
          Vector<unsigned> proc_point_found_plus_one(Nplot_points, 0);
          MPI_Reduce(&tmp_proc_point_found_plus_one[0],
                     &proc_point_found_plus_one[0],
                     Nplot_points,
                     MPI_UNSIGNED,
                     MPI_MAX,
                     0,
                     Comm_pt->mpi_comm());


          // Main process write data
          if (my_rank == 0)
          {
            // Collect all the data
            Vector<Vector<double>> received_data(nproc - 1);
            Vector<Vector<unsigned>> received_size(nproc - 1);
            Vector<unsigned> counter_d(nproc - 1, 0);
            Vector<unsigned> counter_s(nproc - 1, 0);

            // Loop over processors that send their points
            for (int i = 1; i < nproc; i++)
            {
              // Receive sizes of data
              MPI_Recv(&buff_size,
                       1,
                       MPI_UNSIGNED,
                       i,
                       tag,
                       Comm_pt->mpi_comm(),
                       &stat);
              received_size[i - 1].resize(std::max(unsigned(1), buff_size));
              MPI_Recv(&received_size[i - 1][0],
                       buff_size,
                       MPI_UNSIGNED,
                       i,
                       tag,
                       Comm_pt->mpi_comm(),
                       &stat);

              // Receive actual data
              MPI_Recv(&buff_size,
                       1,
                       MPI_UNSIGNED,
                       i,
                       tag,
                       Comm_pt->mpi_comm(),
                       &stat);
              received_data[i - 1].resize(std::max(unsigned(1), buff_size));
              MPI_Recv(&received_data[i - 1][0],
                       buff_size,
                       MPI_DOUBLE,
                       i,
                       tag,
                       Comm_pt->mpi_comm(),
                       &stat);
            }

            // Analyse data for each point
            for (unsigned i = 0; i < Nplot_points; i++)
            {
              // Somebody has found it
              if (proc_point_found_plus_one[i] != 0)
              {
                // Root processor has found it
                if (proc_point_found_plus_one[i] == 1)
                {
                  // Copy directly from vec vector
                  data[i] = vec[i];
                }
                // Another (non-root) processor has found it
                else
                {
                  unsigned line_i = proc_point_found_plus_one[i] - 2;

                  // Resize data line
                  data[i].resize(received_size[line_i][counter_s[line_i]]);

                  // Copy values
                  for (unsigned j = 0;
                       j < received_size[line_i][counter_s[line_i]];
                       j++)
                  {
                    data[i][j] = received_data[line_i][counter_d[line_i] + j];
                  }

                  // Increase counter
                  counter_d[line_i] += received_size[line_i][counter_s[line_i]];
                  counter_s[line_i]++;
                }
              } // end somebody has found it -- no output at all if nobody
              // has found the point (e.g. outside mesh)
            }
          }
          // Send data to root
          else
          {
            // Send the number of fields to the main process
            buff_size = size_values.size();
            MPI_Send(&buff_size, 1, MPI_UNSIGNED, 0, tag, Comm_pt->mpi_comm());

            // Send the sizes of fields to the main process
            if (buff_size == 0) size_values.resize(1);
            MPI_Send(&size_values[0],
                     buff_size,
                     MPI_UNSIGNED,
                     0,
                     tag,
                     Comm_pt->mpi_comm());

            // Send the number of data fields to the main process
            buff_size = local_values.size();
            MPI_Send(&buff_size, 1, MPI_UNSIGNED, 0, tag, Comm_pt->mpi_comm());

            // Send the data to the main process
            if (buff_size == 0) local_values.resize(1);
            MPI_Send(&local_values[0],
                     buff_size,
                     MPI_DOUBLE,
                     0,
                     tag,
                     Comm_pt->mpi_comm());
          }

#endif // Serial version
        }
      }
      else
      {
        // Loop over the points
        for (unsigned i = 0; i < Nplot_points; i++)
        {
          // Check if the point was found in the mesh
          if (Plot_point[i].first != NULL) // success
          {
            // Copy line into data array
            Plot_point[i].first->point_output_data(Plot_point[i].second,
                                                   data[i]);
          }
          else // not found -- keep empty block there for debugging
          {
            // oomph_info << "Point " << i << " not found\n";
          }
        }
      }
    }


    ///  Update plot points coordinates (in preparation of remesh,
    /// say).
    void update_plot_points_coordinates(Vector<Vector<double>>& coord_vec)
    {
      // Resize coord_vec
      coord_vec.resize(Nplot_points);

      // Check that the communication pointer is initialised and the
      // problem is distributed.
      if (Comm_pt != 0)
      {
        int nproc = Comm_pt->nproc();
        if (nproc > 1)
        {
#ifdef OOMPH_HAS_MPI

          // Declaration of MPI variables
          MPI_Status stat;
          int tag; // cgj: tag should be initialised before use
          int my_rank = Comm_pt->my_rank();


          // Buffer
          unsigned buff_size;

          // Create array which contains data found in every process
          Vector<Vector<double>> vec(Nplot_points);

          for (unsigned i = 0; i < Nplot_points; i++)
          {
            if (Plot_point[i].first != NULL)
            {
              if (!((*Plot_point[i].first).is_halo()))
              {
                unsigned dim = Plot_point[i].second.size();

                vec[i].resize(dim);

                for (unsigned j = 0; j < dim; j++)
                {
                  vec[i][j] = Plot_point[i].first->interpolated_x(
                    Plot_point[i].second, j);
                }
              }
            }
          }


          // Analyse which plot points have been found
          // locally and concatenate the data:

          // This contains the flat-packed doubles to be sent
          // for all located plot points
          Vector<double> local_values;

          // Number of values to be sent for each plot point
          // (almost certainly the same for all plot points, but...)
          // size_values[i] gives the number of doubles to be
          // sent for plot point i.
          Vector<unsigned> size_values;

          // Each processor indicates if it has found a given plot point.
          // Once this is gathered on the root processor we know
          // exactly which data we'll receive from where.
          Vector<unsigned> tmp_proc_point_found_plus_one(Nplot_points, 0);

          // Loop over the plot points
          for (unsigned i = 0; i < Nplot_points; i++)
          {
            unsigned ndata = vec[i].size();
            if (ndata != 0)
            {
              // Store the number of fields
              size_values.push_back(ndata);

              // Update found vector
              tmp_proc_point_found_plus_one[i] = my_rank + 1;


              // Store values
              for (unsigned j = 0; j < ndata; j++)
              {
                local_values.push_back(vec[i][j]);
              }
            }
          }

          // Gather information on root

          // Find out who's found the points
          Vector<unsigned> proc_point_found_plus_one(Nplot_points, 0);
          MPI_Reduce(&tmp_proc_point_found_plus_one[0],
                     &proc_point_found_plus_one[0],
                     Nplot_points,
                     MPI_UNSIGNED,
                     MPI_MAX,
                     0,
                     Comm_pt->mpi_comm());

          // Main process write data
          if (my_rank == 0)
          {
            // Collect all the data
            Vector<Vector<double>> received_data(nproc - 1);
            Vector<Vector<unsigned>> received_size(nproc - 1);
            Vector<unsigned> counter_d(nproc - 1, 0);
            Vector<unsigned> counter_s(nproc - 1, 0);

            // Loop over processors that send their points
            for (int i = 1; i < nproc; i++)
            {
              // Receive sizes of data
              MPI_Recv(&buff_size,
                       1,
                       MPI_UNSIGNED,
                       i,
                       tag,
                       Comm_pt->mpi_comm(),
                       &stat);
              received_size[i - 1].resize(std::max(unsigned(1), buff_size));
              MPI_Recv(&received_size[i - 1][0],
                       buff_size,
                       MPI_UNSIGNED,
                       i,
                       tag,
                       Comm_pt->mpi_comm(),
                       &stat);

              // Receive actual data
              MPI_Recv(&buff_size,
                       1,
                       MPI_UNSIGNED,
                       i,
                       tag,
                       Comm_pt->mpi_comm(),
                       &stat);
              received_data[i - 1].resize(std::max(unsigned(1), buff_size));
              MPI_Recv(&received_data[i - 1][0],
                       buff_size,
                       MPI_DOUBLE,
                       i,
                       tag,
                       Comm_pt->mpi_comm(),
                       &stat);
            }

            // Analyse data for each point
            for (unsigned i = 0; i < Nplot_points; i++)
            {
              // Somebody has found it
              if (proc_point_found_plus_one[i] != 0)
              {
                // Root processor has found it
                if (proc_point_found_plus_one[i] == 1)
                {
                  // Copy directly from vec vector
                  coord_vec[i] = vec[i];
                }
                // Another (non-root) processor has found it
                else
                {
                  unsigned line_i = proc_point_found_plus_one[i] - 2;

                  // Resize data line
                  coord_vec[i].resize(received_size[line_i][counter_s[line_i]]);

                  // Copy values
                  for (unsigned j = 0;
                       j < received_size[line_i][counter_s[line_i]];
                       j++)
                  {
                    coord_vec[i][j] =
                      received_data[line_i][counter_d[line_i] + j];
                  }

                  // Increase counter
                  counter_d[line_i] += received_size[line_i][counter_s[line_i]];
                  counter_s[line_i]++;
                }
              } // end somebody has found it -- no output at all if nobody
              // has found the point (e.g. outside mesh)
            }
          }
          // Send data to root
          else
          {
            // Send the number of fields to the main process
            buff_size = size_values.size();
            MPI_Send(&buff_size, 1, MPI_UNSIGNED, 0, tag, Comm_pt->mpi_comm());

            // Send the sizes of fields to the main process
            if (buff_size == 0) size_values.resize(1);
            MPI_Send(&size_values[0],
                     buff_size,
                     MPI_UNSIGNED,
                     0,
                     tag,
                     Comm_pt->mpi_comm());

            // Send the number of data fields to the main process
            buff_size = local_values.size();
            MPI_Send(&buff_size, 1, MPI_UNSIGNED, 0, tag, Comm_pt->mpi_comm());

            // Send the data to the main process
            if (buff_size == 0) local_values.resize(1);
            MPI_Send(&local_values[0],
                     buff_size,
                     MPI_DOUBLE,
                     0,
                     tag,
                     Comm_pt->mpi_comm());
          }

#endif // Serial version
        }
      }
      else
      {
        get_local_plot_points_coordinates(coord_vec);
      }
    }

  private:
    ///  Max radius beyond which we stop searching the bin. Initialised
    /// to DBL_MAX so keep going until the point is found or until
    /// we've searched every single bin. Overwriting this means we won't search
    /// in bins whose closest vertex is at a distance greater than
    /// Max_search_radius from the point to be located.
    double Max_search_radius;

    ///  Pointer to communicator -- allows us to collect data on
    /// processor 0 if the mesh is distributed.
    OomphCommunicator* Comm_pt;

    /// Helper function to setup from file
    void setup_from_file(Mesh* mesh_pt,
                         const std::string file_name,
                         const double& scale)
    {
      // Open file use ifstream
      std::ifstream file_input(file_name.c_str(), std::ios_base::in);
      if (!file_input)
      {
        std::ostringstream error_message;
        error_message << "Cannot open file " << file_name << "\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      if (!file_input.is_open())
      {
        std::ostringstream error_message;
        error_message << "Cannot open file " << file_name << "\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Declaration of variables
      std::string line;
      Vector<Vector<double>> coord_vec_tmp; // Coord array

      // Loop over the lines of the input file
      while (getline(file_input, line) /* != 0*/)
      {
        // Test if the first char of the line is a number
        // using ascii enumeration of chars
        if (isdigit(line[0]))
        {
          Vector<double> tmp(3);

          // Read the 3 first doubles of the line
          // Return 3 if success and less if error
          int n =
            sscanf(line.c_str(), "%lf %lf %lf", &tmp[0], &tmp[1], &tmp[2]);

          if (n == 3) // success
          {
            // Rescaling
            for (unsigned i = 0; i < 3; i++)
            {
              tmp[i] *= scale;
            }

            // Add the new point to the list
            coord_vec_tmp.push_back(tmp);
          }
          else // error
          {
            oomph_info << "Line ignored \n";
          }
        }
      }

      // Call to the helper function
      setup(mesh_pt, coord_vec_tmp);
    }


    ///  Helper function to setup the output structures
    void setup(Mesh* mesh_pt, const Vector<Vector<double>>& coord_vec)
    {
      // Read out number of plot points
      Nplot_points = coord_vec.size();

      if (Nplot_points == 0) return;

      // Keep track of unlocated plot points
      unsigned count_not_found_local = 0;

      // Dimension
      unsigned dim = coord_vec[0].size();

      // Make space
      Plot_point.resize(Nplot_points);

      // Transform mesh into a geometric object
      MeshAsGeomObject mesh_geom_tmp(mesh_pt);

      // Limit the search radius
      mesh_geom_tmp.sample_point_container_pt()->max_search_radius() =
        Max_search_radius;

      // Loop over input points
      double tt_start = TimingHelpers::timer();
      for (unsigned i = 0; i < Nplot_points; i++)
      {
        // Local coordinate of the plot point with its element
        Vector<double> s(dim, 0.0);

        // Pointer to GeomObject that contains the plot point
        GeomObject* geom_pt = 0;

        // Locate zeta
        mesh_geom_tmp.locate_zeta(coord_vec[i], geom_pt, s);

        // Upcast GeomElement as a FiniteElement
        FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(geom_pt);

        // Another one not found locally...
        if (fe_pt == 0)
        {
          count_not_found_local++;
          /* oomph_info << "NOT Found the one at "  */
          /*            << coord_vec[i][0] << " " */
          /*            << coord_vec[i][1] << "\n"; */
        }
        else
        {
          /* oomph_info << "Found the one at "  */
          /*            << coord_vec[i][0] << " " */
          /*            << coord_vec[i][1] << "\n"; */
        }

        // Save result in a pair
        Plot_point[i] = std::pair<FiniteElement*, Vector<double>>(fe_pt, s);
      }


      oomph_info << "Number of points not found locally: "
                 << count_not_found_local << std::endl;

      // Global equivalent (is overwritten below if mpi)
      unsigned count_not_found = count_not_found_local;

      // Check communication pointer exists and problem is
      // distributed.
      if (Comm_pt != 0)
      {
        int nproc = Comm_pt->nproc();
        if (nproc > 1)
        {
#ifdef OOMPH_HAS_MPI

          // Declaration of MPI variables
          int my_rank = Comm_pt->my_rank();

          // Each processor indicates if it has found a given plot point.
          // Once this is gathered on the root processor we know
          // exactly which data we'll receive from where.
          Vector<unsigned> tmp_proc_point_found_plus_one(Nplot_points, 0);

          // Loop over the plot points
          for (unsigned i = 0; i < Nplot_points; i++)
          {
            // Found locally?
            if (Plot_point[i].first != 0)
            {
              tmp_proc_point_found_plus_one[i] = my_rank + 1;
            }
          }

          // Gather information on root

          // Find out who's found the points
          Vector<unsigned> proc_point_found_plus_one(Nplot_points, 0);
          MPI_Reduce(&tmp_proc_point_found_plus_one[0],
                     &proc_point_found_plus_one[0],
                     Nplot_points,
                     MPI_UNSIGNED,
                     MPI_MAX,
                     0,
                     Comm_pt->mpi_comm());


          // Main process analyses data
          if (my_rank == 0)
          {
            // Analyse data for each point
            count_not_found = 0;
            for (unsigned i = 0; i < Nplot_points; i++)
            {
              // Nobody has found it
              if (proc_point_found_plus_one[i] == 0)
              {
                count_not_found++;
              }
            }
          }

          // Now tell everybody about it
          MPI_Bcast(&count_not_found, 1, MPI_UNSIGNED, 0, Comm_pt->mpi_comm());

#endif
        }
      }
      oomph_info << "Number of plot points not found (with max search radius="
                 << Max_search_radius << ")]: " << count_not_found
                 << "\nTotal time for LineVisualiser setup [sec]: "
                 << TimingHelpers::timer() - tt_start << std::endl;
    }

    // Get coordinates of found points
    void get_local_plot_points_coordinates(Vector<Vector<double>>& data)
    {
      data.resize(Nplot_points);
      for (unsigned i = 0; i < Nplot_points; i++)
      {
        if (Plot_point[i].first != NULL)
        {
          unsigned dim = Plot_point[i].second.size();

          data[i].resize(dim);

          for (unsigned j = 0; j < dim; j++)
          {
            data[i][j] =
              Plot_point[i].first->interpolated_x(Plot_point[i].second, j);
          }
        }
      }
    }


    ///  Vector of pairs containing points to finite elements and
    /// local coordinates
    Vector<std::pair<FiniteElement*, Vector<double>>> Plot_point;

    /// Number of plot points
    unsigned Nplot_points;

  }; // end of class

} // namespace oomph

#endif
