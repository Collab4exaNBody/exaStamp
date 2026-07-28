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
#include <onika/physics/units.h>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  class MovePlaneWall : public OperatorNode
  {

    ADD_SLOT(Vec3d, init_normal, INPUT, Vec3d{1.0, 0.0, 0.0});
    ADD_SLOT(double, init_offset, INPUT, 0.0);
    ADD_SLOT(double, init_cutoff, INPUT, REQUIRED);
    ADD_SLOT(double, init_epsilon, INPUT, REQUIRED);
    ADD_SLOT(double, init_time, INPUT, OPTIONAL);
    ADD_SLOT(double, init_velocity, INPUT, OPTIONAL);
    ADD_SLOT(double, final_time, INPUT, OPTIONAL);
    ADD_SLOT(double, final_velocity, INPUT, OPTIONAL);
    ADD_SLOT(long, init_exponent, INPUT, 12);
    ADD_SLOT(double, physical_time, INPUT, REQUIRED);

    // alternative to init_time/init_velocity/final_time/final_velocity: an arbitrary
    // number of (time, velocity) waypoints, for multiple-shock (reshock) scenarios
    ADD_SLOT(std::vector<onika::physics::Quantity>, times, INPUT, OPTIONAL);
    ADD_SLOT(std::vector<onika::physics::Quantity>, velocities, INPUT, OPTIONAL);
    ADD_SLOT(std::string, velocity_profile, INPUT, "step");
    ADD_SLOT(bool, freeze_at_end, INPUT, true);

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

For multiple-shock (reshock) scenarios, `times` and `velocities` can be given instead of
init_time/init_velocity/final_time/final_velocity: two same-length lists of waypoints
(time[i], velocity[i]), at least 2 points, strictly increasing times. Before times[0], offset
stays at init_offset. times/velocities cannot be combined with
init_time/final_time/final_velocity.

Between waypoints, velocity_profile picks the shape:
  - "step" (default): velocity holds constant at velocity[i] over [time[i], time[i+1]) —
    i.e. it jumps to velocity[i] exactly at time[i] and stays there until time[i+1]. Each
    jump is a real velocity discontinuity, by design (that's the point of a reshock).
  - "linear": velocity ramps linearly from velocity[i] to velocity[i+1] over that interval,
    chaining the single init/final ramp above across all waypoints. Velocity is continuous
    across every waypoint (only the acceleration changes, from one segment's slope to the
    next's).

What happens once physical_time reaches times[last] is controlled by freeze_at_end:
  - true (default): motion freezes at the offset reached at times[last], and epsilon is
    forced to 0 (wall "removed"), same as final_time above. In "step" mode, velocity[last]
    is then never actually used to move the wall, since freezing happens before the
    interval starting at time[last] would begin; in "linear" mode velocity[last] is still
    used, as the ramp target the last segment reaches exactly at time[last].
  - false: the wall instead keeps moving at constant velocity[last] forever past
    times[last] (zero acceleration — this is a constant-velocity extrapolation, not
    another ramp segment), and the wall stays active. In "linear" mode this is a smooth
    continuation (velocity was already heading to velocity[last], so there's no jump,
    only the acceleration drops to 0). In "step" mode velocity[last] only takes effect at
    this point, and may jump relative to the previous step's velocity[last-1].

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

Multiple-shock example:

myoperator:
  - move_wall:
      init_offset: 0.0
      init_cutoff: 5.0 ang
      init_epsilon: 1.0e-19 J
      times: [10.0 ps, 30.0 ps, 50.0 ps, 80.0 ps]
      velocities: [0.01 ang/ps, 0.03 ang/ps, 0.06 ang/ps, 0.0 ang/ps]
      velocity_profile: step
  - wall
)EOF";
    }

    inline void execute() override final
    {
      const bool list_mode = times.has_value() || velocities.has_value();

      if (list_mode && (!times.has_value() || !velocities.has_value()))
      {
        fatal_error() << "move_wall: times and velocities must both be provided together" << std::endl;
      }
      if (list_mode && (init_time.has_value() || final_time.has_value() || final_velocity.has_value()))
      {
        fatal_error() << "move_wall: times/velocities cannot be combined with init_time/final_time/final_velocity" << std::endl;
      }
      if (list_mode && times->size() != velocities->size())
      {
        fatal_error() << "move_wall: times and velocities must have the same number of elements (" << times->size() << " vs " << velocities->size() << ")" << std::endl;
      }
      if (list_mode && times->size() < 2)
      {
        fatal_error() << "move_wall: times/velocities need at least 2 points" << std::endl;
      }
      if (list_mode && *velocity_profile != "step" && *velocity_profile != "linear")
      {
        fatal_error() << "move_wall: velocity_profile must be \"step\" or \"linear\", got \"" << *velocity_profile << "\"" << std::endl;
      }
      if (!list_mode && (!init_time.has_value() || !init_velocity.has_value()))
      {
        fatal_error() << "move_wall: init_time and init_velocity are required unless times/velocities are given" << std::endl;
      }
      if (!list_mode && final_time.has_value() && *final_time <= *init_time)
      {
        fatal_error() << "move_wall: final_time (" << *final_time << ") must be greater than init_time (" << *init_time << ")" << std::endl;
      }
      if (!list_mode && final_velocity.has_value() && !final_time.has_value())
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

      if (list_mode)
      {
        std::vector<double> t, v;
        t.reserve(times->size());
        v.reserve(velocities->size());
        for (auto &q : *times) t.push_back(q.convert());
        for (auto &q : *velocities) v.push_back(q.convert());
        for (size_t i = 1; i < t.size(); i++)
        {
          if (t[i] <= t[i - 1])
          {
            fatal_error() << "move_wall: times must be strictly increasing" << std::endl;
          }
        }

        const bool step = (*velocity_profile == "step");

        if (*physical_time < t.front())
        {
          *offset = *init_offset;
        }
        else if (*physical_time >= t.back())
        {
          *offset = *init_offset + waypoints_displacement(t, v, step, t.size() - 2, t.back());
          if (*freeze_at_end)
          {
            *epsilon = 0.0;
          }
          else
          {
            velocity = v.back();
            *offset += velocity * (*physical_time - t.back());
          }
        }
        else
        {
          size_t i = 0;
          while (i + 1 < t.size() && *physical_time >= t[i + 1]) i++;
          *offset = *init_offset + waypoints_displacement(t, v, step, i, *physical_time);
          if (step)
          {
            velocity = v[i];
          }
          else
          {
            acceleration = (v[i + 1] - v[i]) / (t[i + 1] - t[i]);
            velocity = v[i] + acceleration * (*physical_time - t[i]);
          }
        }
      }
      else if (*physical_time >= *init_time)
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

      if (!list_mode && final_time.has_value() && *physical_time >= *final_time)
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
    // displacement accumulated from t[0] up to tau, given tau is in [t[upto], t[upto+1]]
    static inline double waypoints_displacement(const std::vector<double> &t, const std::vector<double> &v, bool step, size_t upto, double tau)
    {
      double d = 0.0;
      for (size_t i = 0; i < upto; i++)
      {
        const double s = t[i + 1] - t[i];
        const double accel = step ? 0.0 : (v[i + 1] - v[i]) / s;
        d += v[i] * s + 0.5 * accel * s * s;
      }
      const double s = tau - t[upto];
      const double accel = step ? 0.0 : (v[upto + 1] - v[upto]) / (t[upto + 1] - t[upto]);
      d += v[upto] * s + 0.5 * accel * s * s;
      return d;
    }

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
  ONIKA_AUTORUN_INIT(move_plane_wall)
  {
    OperatorNodeFactory::instance()->register_factory("move_plane_wall", make_simple_operator<MovePlaneWall>);
  }

}
