/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_stream.h>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  class MoveWall : public OperatorNode
  {

    ADD_SLOT(Vec3d, init_normal, INPUT, Vec3d{1.0, 0.0, 0.0});
    ADD_SLOT(double, init_offset, INPUT, 0.0);
    ADD_SLOT(double, init_cutoff, INPUT, REQUIRED);
    ADD_SLOT(double, init_epsilon, INPUT, REQUIRED);
    ADD_SLOT(double, init_time, INPUT, REQUIRED);
    ADD_SLOT(double, init_velocity, INPUT, REQUIRED);
    ADD_SLOT(double, final_time, INPUT, OPTIONAL);
    ADD_SLOT(double, final_velocity, INPUT, OPTIONAL);
    ADD_SLOT(long, init_exponent, INPUT, 12);
    ADD_SLOT(double, physical_time, INPUT, REQUIRED);

    ADD_SLOT(MPI_Comm, mpi, INPUT, MPI_COMM_WORLD);
    ADD_SLOT(std::string, csv_filename, INPUT, OPTIONAL);
    ADD_SLOT(std::string, csv_separator, INPUT, ",");
    ADD_SLOT(bool, csv_append, INPUT, false);

    // outputs for wall
    ADD_SLOT(Vec3d, normal, OUTPUT);
    ADD_SLOT(double, offset, OUTPUT);
    ADD_SLOT(double, cutoff, OUTPUT);
    ADD_SLOT(double, epsilon, OUTPUT);
    ADD_SLOT(long, exponent, OUTPUT);

    std::ofstream m_csv_stream;

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Computes the (normal, offset, cutoff, epsilon, exponent) parameters of a moving wall
for a given instant, to be fed into the `wall` operator.

Before init_time, offset stays at init_offset. From init_time onward, offset advances
at constant init_velocity: offset = init_offset + (physical_time - init_time) * init_velocity.

final_time is optional. If set, motion instead follows constant acceleration from
init_velocity to final_velocity (reached at final_time), so non-constant velocity loading
(e.g. a ramped/accelerating piston) can be modeled; final_velocity defaults to init_velocity
if not given. Motion freezes once final_time is reached, and epsilon is then forced to 0,
so the `wall` operator it feeds has no more effect (wall "removed").
final_velocity requires final_time to be set. final_time must be greater than init_time.

If csv_filename is set, rank 0 appends one row per call (time, position, velocity,
acceleration) to that file, fields joined with csv_separator (default ","). The header
row is written once, only when the file is empty. csv_append (default false) controls
whether an existing file is truncated or appended to when (re)opened.

YAML example:

myoperator:
  - move_wall:
      init_normal: [1.0, 0.0, 0.0]
      init_offset: 0.0
      init_cutoff: 5.0 ang
      init_epsilon: 1.0e-19 J
      init_time: 10.0 ps
      init_velocity: 0.01 ang/ps
      final_time: 50.0 ps
      final_velocity: 0.05 ang/ps
      csv_filename: wall_trajectory.csv
      csv_separator: ";"
      csv_append: false
  - wall
)EOF";
    }

    inline void execute() override final
    {
      if (final_time.has_value() && *final_time <= *init_time)
      {
        fatal_error() << "move_wall: final_time (" << *final_time << ") must be greater than init_time (" << *init_time << ")" << std::endl;
      }
      if (final_velocity.has_value() && !final_time.has_value())
      {
        fatal_error() << "move_wall: final_velocity requires final_time to be set" << std::endl;
      }

      *normal = *init_normal;
      *cutoff = *init_cutoff;
      *epsilon = *init_epsilon;
      *offset = *init_offset;
      *exponent = *init_exponent;

      double velocity = 0.0;
      double acceleration = 0.0;

      if (*physical_time >= *init_time)
      {
        if (final_time.has_value())
        {
          const double v1 = final_velocity.has_value() ? *final_velocity : *init_velocity;
          const double T = *final_time - *init_time;
          double s = *physical_time - *init_time;
          const double accel = (v1 - *init_velocity) / T;
          if (s >= T) // motion frozen once final_time is reached (wall is disabled anyway)
          {
            s = T;
          }
          else
          {
            velocity = *init_velocity + accel * s;
            acceleration = accel;
          }
          *offset = *init_offset + (*init_velocity) * s + 0.5 * accel * s * s;
        }
        else
        {
          *offset = *init_offset + (*physical_time - *init_time) * (*init_velocity);
          velocity = *init_velocity;
        }
      }

      else
      {
        *offset = *init_offset;
      }

      if (final_time.has_value() && *physical_time >= *final_time)
      {
        *epsilon = 0.0;
      }

      ldbg << "offset=" << *offset << std::endl;

      if (csv_filename.has_value())
      {
        int rank = 0;
        MPI_Comm_rank(*mpi, &rank);
        if (rank == 0)
        {
          write_csv_row(*physical_time, *offset, velocity, acceleration);
        }
      }
    }

  private:
    inline void write_csv_row(double time, double position, double velocity, double acceleration)
    {
      if (!m_csv_stream.is_open())
      {
        m_csv_stream.open(*csv_filename, *csv_append ? std::ios::app : std::ios::trunc);
        if (m_csv_stream.tellp() == std::streampos(0))
        {
          m_csv_stream << "time" << *csv_separator << "position" << *csv_separator << "velocity" << *csv_separator << "acceleration" << "\n";
        }
      }
      m_csv_stream << time << *csv_separator << position << *csv_separator << velocity << *csv_separator << acceleration << "\n";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(move_wall)
  {
    OperatorNodeFactory::instance()->register_factory("move_wall", make_simple_operator<MoveWall>);
  }

}
