// Mars lander simulator
// Version 1.11
// Mechanical simulation functions
// Gabor Csanyi and Andrew Gee, August 2019

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

#include "lander.h"

void autopilot (void)
  // Autopilot to adjust the engine throttle, parachute and attitude control
{
    //autopilot purpose set in lander.h
    //hovering is for info examples 1
    //NB: if not using landing, FUEL_RATE_AT_MAX_THRUST will be 0 not 0.5

    double lander_mass = UNLOADED_LANDER_MASS + FUEL_DENSITY * FUEL_CAPACITY * fuel;
    //only used in Examples 
    double F_eq = GRAVITY * MARS_MASS * lander_mass / pow(MARS_RADIUS + target_altitude, 2);

    //altitude
    double h = position.abs() - MARS_RADIUS;
    if (auto_purpose == Landing)
    {
        //Trial constant values; adjust depending on results for 100l
        //height constant, works at 0.018, optimised to 0.04
        const double K_h = 0.04;
        //controller gain, 0.34, 0.335 borderline
        const double K_p = 0.34;
        //offset, 0.08
        const double DELTA = 0.08;
        
        //error term, position.norm( ) is e_r
        double e = -(0.5 + K_h * h + velocity * position.norm());
        double P_out = K_p * e;

        //set throttle accordingly
        if (P_out <= -DELTA)
        {
            throttle = 0;
        }
        else if (-DELTA < P_out < 1 - DELTA)
        {
            throttle = DELTA + P_out;
        }
        else throttle = 1;

        //final part of task; values to file
        const bool write_to_file = true;
        if (write_to_file)
        {
            static ofstream fout;
            if (not fout.is_open())
                fout.open("descent_rates.txt");

            double target = 0.5 + K_h * h;
            double actual = -1 * velocity * position.norm();
            fout << h << ' ' << target << ' ' << actual << "\n";

            if (h == 0)
                fout.close();
        }
    }
    else if (auto_purpose == Examples1)
    {
        //hovering
        throttle = F_eq / MAX_THRUST;
        //hovering control system not implemented
    }
    else if (auto_purpose == Examples2 || auto_purpose == Examples3)
    {
        const float K_p = 0.091, K_d = 0.178;
        const bool part1 = false; //part 1 ex 2
        const bool no_feedback = true; //ex3 q5 ptd
        float error = target_altitude - h;
        float throttle_aim;
        if (part1 or no_feedback)
        {
            if (part1)
            {
                //K(s) = Kp therefore
                throttle_aim = K_p * error + F_eq / MAX_THRUST;

            }
            else if (no_feedback)
                throttle_aim = F_eq / MAX_THRUST;
            if (throttle_aim > 1)
                throttle = 1;
            else if (throttle_aim < 0)
                throttle = 0;
            else
                throttle = throttle_aim;
        }
        else
        {
            PDController(K_p, K_d, error, F_eq);
        }
    }
    else if (auto_purpose == Examples4)
        //delay or lag implemented
    {
        const float K_p = 0.091, K_d = 2*K_p;
        float error = target_altitude - h;
        PDController(K_p, K_d, error, F_eq);
    }
}
//sets throttle value based on PD controller with values given
void PDController(float K_p, float K_d, float error, float F_eq)
{
    static float last_error;
    if (simulation_time == 0)
            last_error = error;
    float edot = (error - last_error) / delta_t;
    float throttle_aim = K_d * edot + K_p * error + F_eq / MAX_THRUST;
    if (throttle_aim > 1)
        throttle = 1;
    else if (throttle_aim < 0)
        throttle = 0;
    else
        throttle = throttle_aim;

    last_error = error;

}

void numerical_dynamics (void)
  // This is the function that performs the numerical integration to update the
  // lander's pose. The time step is delta_t (global variable).
{
    // INSERT YOUR CODE HERE
    //variable to choose which integrator is being used
    const bool EULER_INTEGRATOR_CHOSEN = false; 
    static vector3d previous_position;
    vector3d thr = thrust_wrt_world();
    double density = atmospheric_density(position);
    vector3d drag = -0.5 * density * DRAG_COEF_LANDER *
        (M_PI * LANDER_SIZE * LANDER_SIZE) * velocity.abs2()
        * velocity.norm();
    if (parachute_status == DEPLOYED)
    {
        //each of the five parachute segments is a square with side length
        //2.0*LANDER_SIZE. Updates drag with drag from parachute
        drag += -0.5 * density * DRAG_COEF_CHUTE *
            (5 * 4 * LANDER_SIZE * LANDER_SIZE) * velocity.abs2() *
            velocity.norm();
    }
    double lander_mass = UNLOADED_LANDER_MASS + FUEL_DENSITY * FUEL_CAPACITY * fuel;
    vector3d gravity = -GRAVITY * MARS_MASS * lander_mass *
        position.norm() / position.abs2();

    if (auto_purpose == Examples3)
        thr += 100 * cos(0.1 * simulation_time) * position.norm(); //disturbance for examples 3
        //assumes thrust aligned with radial axis.

    vector3d net_force = drag + gravity + thr;

    if (EULER_INTEGRATOR_CHOSEN or simulation_time == 0)
    {
        //Euler integration
        previous_position = position;
        position += delta_t * velocity;
        velocity += delta_t * net_force / lander_mass;
    }
    else
    {
        //Verlet integration
        //When a new simulation is chosen, the time returns to 0
        //Verlet will not have a legit previous position, so 
        //Euler should be used, hence the condition
        vector3d new_position = 2 * position - previous_position +
            pow(delta_t, 2) * net_force / lander_mass;
        //velocity updated in sync with position for completeness
        velocity = 1 / delta_t * (new_position - position);
        previous_position = position;
        position = new_position;
    }

    

  // Here we can apply an autopilot to adjust the thrust, parachute and attitude
  if (autopilot_enabled) autopilot();

  // Here we can apply 3-axis stabilization to ensure the base is always pointing downwards
  if (stabilized_attitude) attitude_stabilization();
}

void initialize_simulation (void)
  // Lander pose initialization - selects one of 10 possible scenarios
{
  // The parameters to set are:
  // position - in Cartesian planetary coordinate system (m)
  // velocity - in Cartesian planetary coordinate system (m/s)
  // orientation - in lander coordinate system (xyz Euler angles, degrees)
  // delta_t - the simulation time step
  // boolean state variables - parachute_status, stabilized_attitude, autopilot_enabled
  // scenario_description - a descriptive string for the help screen

  scenario_description[0] = "circular orbit";
  scenario_description[1] = "descent from 10km";
  scenario_description[2] = "elliptical orbit, thrust changes orbital plane";
  scenario_description[3] = "polar launch at escape velocity (but drag prevents escape)";
  scenario_description[4] = "elliptical orbit that clips the atmosphere and decays";
  scenario_description[5] = "descent from 200km";
  scenario_description[6] = "altitude control test at 500m";
  scenario_description[7] = "altitude control test at 510m";
  scenario_description[8] = "altitude control test at 700m";
  scenario_description[9] = "";

  //default value, will be reassigned for scenarios that need it
  target_altitude = 0;

  switch (scenario) {

  case 0:
    // a circular equatorial orbit
    position = vector3d(1.2*MARS_RADIUS, 0.0, 0.0);
    velocity = vector3d(0.0, -3247.087385863725, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 1:
    // a descent from rest at 10km altitude
    position = vector3d(0.0, -(MARS_RADIUS + 10000.0), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = true; //changed from false
    break;

  case 2:
    // an elliptical polar orbit
    position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
    velocity = vector3d(3500.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 3:
    // polar surface launch at escape velocity (but drag prevents escape)
    position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE/2.0);
    velocity = vector3d(0.0, 0.0, 5027.0);
    orientation = vector3d(0.0, 0.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 4:
    // an elliptical orbit that clips the atmosphere each time round, losing energy
    position = vector3d(0.0, 0.0, MARS_RADIUS + 100000.0);
    velocity = vector3d(4000.0, 0.0, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 5:
    // a descent from rest at the edge of the exosphere
    position = vector3d(0.0, -(MARS_RADIUS + EXOSPHERE), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = true; //changed from false
    break;

    //6 to 8 for info examples 1 q3
  case 6:
    // 500m altitude control test at 500m altitude
    position = vector3d(0.0, -(MARS_RADIUS + 500), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.01;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = true; 
    target_altitude = 500;
    break;

  case 7:
    // 500m altitude control test at 510m altitude
    position = vector3d(0.0, -(MARS_RADIUS + 510), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.01;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = true;
    target_altitude = 500;
    break;

  case 8:
    // 500m altitude control test at 700m altitude
    position = vector3d(0.0, -(MARS_RADIUS + 700), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.01;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = true;
    target_altitude = 500;
    break;

  case 9:
    break;

  }
}
