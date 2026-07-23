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

  class MoveCylinderWall : public OperatorNode
  {

    ADD_SLOT(Vec3d, init_origin, INPUT, Vec3d{0.0, 0.0, 0.0});
    ADD_SLOT(Vec3d, init_axis, INPUT, Vec3d{0.0, 0.0, 1.0});
    ADD_SLOT(double, init_radius, INPUT, REQUIRED);
    ADD_SLOT(double, init_cutoff, INPUT, REQUIRED);
    ADD_SLOT(double, init_epsilon, INPUT, REQUIRED);
    ADD_SLOT(double, init_time, INPUT, REQUIRED);
    ADD_SLOT(double, init_velocity, INPUT, REQUIRED);
    ADD_SLOT(double, final_time, INPUT, OPTIONAL);
    ADD_SLOT(double, final_velocity, INPUT, OPTIONAL);
    ADD_SLOT(long, init_exponent, INPUT, 12);
    ADD_SLOT(double, physical_time, INPUT, REQUIRED);

    // optional translation of the origin, sharing the same init_time/final_time window;
    // independent of axis (axis only orients the cylinder, it doesn't move)
    ADD_SLOT(Vec3d, direction, INPUT, Vec3d{0.0, 0.0, 0.0});
    ADD_SLOT(double, init_translation_velocity, INPUT, 0.0);
    ADD_SLOT(double, final_translation_velocity, INPUT, OPTIONAL);

    ADD_SLOT(MPI_Comm, mpi, INPUT, MPI_COMM_WORLD);
    ADD_SLOT(std::string, csv_filename, INPUT, OPTIONAL);
    ADD_SLOT(std::string, csv_separator, INPUT, ",");
    ADD_SLOT(bool, csv_append, INPUT, false);

    // outputs for cylinder_wall
    ADD_SLOT(Vec3d, origin, OUTPUT);
    ADD_SLOT(Vec3d, axis, OUTPUT);
    ADD_SLOT(double, radius, OUTPUT);
    ADD_SLOT(double, cutoff, OUTPUT);
    ADD_SLOT(double, epsilon, OUTPUT);
    ADD_SLOT(long, exponent, OUTPUT);

    std::ofstream m_csv_stream;

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Computes the (origin, axis, radius, cutoff, epsilon, exponent) parameters of a moving
cylindrical wall for a given instant, to be fed into the `cylinder_wall` operator.
`axis` is passed through unchanged (it only orients the cylinder, it never moves).

Before init_time, radius stays at init_radius and origin stays at init_origin. From
init_time onward, radius advances at constant init_velocity: radius = init_radius +
(physical_time - init_time) * init_velocity. A negative init_velocity shrinks the
cylinder (e.g. a converging cylindrical implosion).

Optionally, origin can translate along `direction` (a unit vector, independent of `axis`)
at init_translation_velocity, using the same init_time/final_time window as radius:
origin = init_origin + direction * displacement. direction defaults to {0,0,0} and
init_translation_velocity to 0, so no translation happens unless explicitly requested.

final_time is optional. If set, both radius and translation instead follow constant
acceleration (radius: init_velocity -> final_velocity, translation: init_translation_velocity
-> final_translation_velocity, both reaching their final value at final_time), so
non-constant velocity loading can be modeled; the "final_*" velocities default to their
"init_*" counterpart if not given. Motion freezes once final_time is reached, and epsilon
is then forced to 0, so the `cylinder_wall` operator it feeds has no more effect (wall
"removed"). final_velocity/final_translation_velocity require final_time to be set.
final_time must be greater than init_time.

If csv_filename is set, rank 0 appends one row per call (time, radius, radius_velocity,
radius_acceleration, origin_x, origin_y, origin_z, translation_velocity,
translation_acceleration) to that file, fields joined with csv_separator (default ",").
The header row is written once, only when the file is empty. csv_append (default false)
controls whether an existing file is truncated or appended to when (re)opened.

YAML example:

myoperator:
  - move_cylinder_wall:
      init_origin: [0.0, 0.0, 0.0]
      init_axis: [0.0, 0.0, 1.0]
      init_radius: 40.0 ang
      init_cutoff: 5.0 ang
      init_epsilon: 1.0e-19 J
      init_time: 10.0 ps
      init_velocity: -0.02 ang/ps
      final_time: 50.0 ps
      final_velocity: -0.1 ang/ps
      direction: [1.0, 0.0, 0.0]
      init_translation_velocity: 0.01 ang/ps
      csv_filename: cylinder_wall_trajectory.csv
  - cylinder_wall
)EOF";
    }

    inline void execute() override final
    {
      if (final_time.has_value() && *final_time <= *init_time)
      {
        fatal_error() << "move_cylinder_wall: final_time (" << *final_time << ") must be greater than init_time (" << *init_time << ")" << std::endl;
      }
      if (final_velocity.has_value() && !final_time.has_value())
      {
        fatal_error() << "move_cylinder_wall: final_velocity requires final_time to be set" << std::endl;
      }
      if (final_translation_velocity.has_value() && !final_time.has_value())
      {
        fatal_error() << "move_cylinder_wall: final_translation_velocity requires final_time to be set" << std::endl;
      }

      *origin = *init_origin;
      *axis = *init_axis;
      *cutoff = *init_cutoff;
      *epsilon = *init_epsilon;
      *radius = *init_radius;
      *exponent = *init_exponent;

      double radius_velocity = 0.0;
      double radius_acceleration = 0.0;
      double translation_velocity = 0.0;
      double translation_acceleration = 0.0;

      if (*physical_time >= *init_time)
      {
        if (final_time.has_value())
        {
          const double radius_v1 = final_velocity.has_value() ? *final_velocity : *init_velocity;
          const double translation_v1 = final_translation_velocity.has_value() ? *final_translation_velocity : *init_translation_velocity;
          const double T = *final_time - *init_time;
          double s = *physical_time - *init_time;
          const double radius_accel = (radius_v1 - *init_velocity) / T;
          const double translation_accel = (translation_v1 - *init_translation_velocity) / T;
          if (s >= T) // motion frozen once final_time is reached (wall is disabled anyway)
          {
            s = T;
          }
          else
          {
            radius_velocity = *init_velocity + radius_accel * s;
            radius_acceleration = radius_accel;
            translation_velocity = *init_translation_velocity + translation_accel * s;
            translation_acceleration = translation_accel;
          }
          *radius = *init_radius + (*init_velocity) * s + 0.5 * radius_accel * s * s;
          const double displacement = (*init_translation_velocity) * s + 0.5 * translation_accel * s * s;
          *origin = *init_origin + (*direction) * displacement;
        }
        else
        {
          const double s = *physical_time - *init_time;
          *radius = *init_radius + s * (*init_velocity);
          radius_velocity = *init_velocity;
          *origin = *init_origin + (*direction) * (s * (*init_translation_velocity));
          translation_velocity = *init_translation_velocity;
        }
      }

      else
      {
        *radius = *init_radius;
        *origin = *init_origin;
      }

      if (final_time.has_value() && *physical_time >= *final_time)
      {
        *epsilon = 0.0;
      }

      ldbg << "radius=" << *radius << " origin=" << *origin << std::endl;

      if (csv_filename.has_value())
      {
        int rank = 0;
        MPI_Comm_rank(*mpi, &rank);
        if (rank == 0)
        {
          write_csv_row(*physical_time, *radius, radius_velocity, radius_acceleration, *origin, translation_velocity, translation_acceleration);
        }
      }
    }

  private:
    inline void write_csv_row(double time, double radius_value, double radius_velocity, double radius_acceleration,
                               const Vec3d &origin_value, double translation_velocity, double translation_acceleration)
    {
      if (!m_csv_stream.is_open())
      {
        m_csv_stream.open(*csv_filename, *csv_append ? std::ios::app : std::ios::trunc);
        if (m_csv_stream.tellp() == std::streampos(0))
        {
          const std::string &sep = *csv_separator;
          m_csv_stream << "time" << sep << "radius" << sep << "radius_velocity" << sep << "radius_acceleration" << sep
                       << "origin_x" << sep << "origin_y" << sep << "origin_z" << sep
                       << "translation_velocity" << sep << "translation_acceleration" << "\n";
        }
      }
      const std::string &sep = *csv_separator;
      m_csv_stream << time << sep << radius_value << sep << radius_velocity << sep << radius_acceleration << sep
                   << origin_value.x << sep << origin_value.y << sep << origin_value.z << sep
                   << translation_velocity << sep << translation_acceleration << "\n";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(move_cylinder_wall)
  {
    OperatorNodeFactory::instance()->register_factory("move_cylinder_wall", make_simple_operator<MoveCylinderWall>);
  }

}
