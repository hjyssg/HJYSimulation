//
//  HJYPhysics.h
//  A simple physics engine
//  
//  Created by Junyang Huang on 1/15/14.
//  Copyright (c) 2014 www.hjyssg.com. All rights reserved.
//

#ifndef __HJYPhysics__
#define __HJYPhysics__

//The _OPENMP C preprocessor
#if _OPENMP
	//http://publib.boulder.ibm.com/infocenter/comphelp/v101v121/index.jsp?topic=/com.ibm.xlf121.aix.doc/proguide/openmpc.html
	#include <omp.h>

#endif 


#include <iostream>
#include <vector> 
#include <assert.h> 
#include <string>


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits>      
#include <stdlib.h>


//syntax candy for unit
#define M
#define MM *1e-3
#define NM *1e-9
#define PM *1e-12
#define S
#define MS *1e-3
#define NS *1e-9
#define PS *1e-12
#define FS *1e-15
#define KG 
#define N
#define J
#define C
#define MPS



#pragma mark some usefule constants
const int SECONDS_OF_HOUR = 3600;
const int SECONDS_OF_DAY = 24*SECONDS_OF_HOUR;
const int SECONDS_OF_YEAR = 365*SECONDS_OF_DAY;
#define GRAVIATIONAL_CONSTANT 6.67384e-11 //m3*kg-1*s-2
#define COULOMBS_CONSTANT 8.987552e+9
#define PI 3.14159265359
#define DEG_TO_RAD(a)  a/180*3.14159265359
#define RAD_TO_DEG(a)  a/3.14159265359*180




#pragma mark helper-functions and macro
//print function
#define  CLOG(a)  std::cout<<a<<"\n"
#define  DLOG(a)  std::cout<<__PRETTY_FUNCTION__<<": "<<a<<"\n"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


//tell if  min<num<max
bool is_in_range(double num, double min, double max);

//a small math function:  calculate error percentage from -100 to 100
double calculate_error_percentage(double estimate, double real);


//from http://stackoverflow.com/questions/2704521/generate-random-double-numbers-in-c
//Remember to call srand() with a proper seed each time your program starts.
double fRand(double fMin, double fMax);


#pragma mark some abbreviation used in the code:
/*
acl: acceleration
arr: array, list
num: number
pm: point mass
pos:position
spd: speed
*/

#pragma mark the engine




/*
represents the 3d math vertex. It has  x, y, z coordinate
*/
class Vertex
{

public:
	double x;
	double y;
	double z;

	//constructor
	Vertex();
	Vertex(double x_p, double y_p, double z_p);

	//convenient mehod, set x, y, z in one method call
	void set_xyz(double x_p, double y_p, double z_p);

	//calculate addition
	Vertex  add(Vertex v2);
	void  add_self(Vertex v2);

	//calculate substraction
	Vertex sub(Vertex v2);
	void sub_self(Vertex v2);

	//calculate crosss product
	Vertex crosss_product(Vertex v2);
	//calculate dot product
	double dot_product(Vertex v2);

	//same as add(), only difference is this one pass three doubles than vertex
	Vertex transfer( double x_p, double y_p, double z_p);
	inline void  transfer_self( double x_p, double y_p, double z_p);

	//mutliply number to x, y, z
	Vertex  scale( double x_s, double y_s, double z_s);
	//different numbers multiply with x, y,z
	inline void  scale_self( double x_s, double y_s, double z_s);  
	//one same number multiply with x, y,z
	inline void  scale_self( double p); 

	//return the distance/length to origin
	double length();

	//return the square of the distance/length to origin
	double length_square();

	//return normal vertex of this
	Vertex normalize();

	//normal self
	void normalize_self();

	//calculate distance to v
	double distance_to_v(Vertex v);

	//calculate distance^2 to v
	double distance_square_to_v(Vertex v);

	//tell if this is zero
	bool isZero();

	//tell if two vertexes equal
	bool isEqual(Vertex v);

	//tell if two vertexes are parallel
	bool isParrell(Vertex v);

	//rotate vertex along x axis
	//theta is in radians
	Vertex rotate_x(float theta);
	void rotate_x_self(float theta);

	//rotate vertex along y axis
	//theta is in radians
	Vertex rotate_y(float theta);
	void rotate_y_self(float theta);

	//rotate vertex along z axis
	//theta is in radians
	Vertex rotate_z(float theta);
	void rotate_z_self(float theta);

	//rotate vertex along arbitary axis specified by vtx
	//theta is in radians
	Vertex rotate_around_axis(Vertex vtx, float theta);

	//to string method
	std::string toString();

	//get identical vertex
	Vertex  get_copy();
};

/*
the most basic component in Newton Physics
represent a 3d point
*/
class PointMass
{
public:
	Vertex pos; //position     in meter
	Vertex spd; //speed        in meter/second
	Vertex acl; //acceleration in metet/second^2
	double collision_damping; //from 0 to 1, decide how much percent energy will lose when hitting others
	double mass;  //in kg
	double electric_charge; //in Coulomb, default is 0.0


	//default constructor that set all parameter to 0.0
	PointMass();
	PointMass(double p_x, double p_y, double p_z,
		double spd_x,double spd_y,double spd_z,
		double acl_x,double acl_y,double acl_z,
		double mass_p);

	PointMass(double p_x, double p_y, double p_z,
		double spd_x,double spd_y,double spd_z,
		double acl_x,double acl_y,double acl_z,
		double collision_damping, double mass_p);

	//calculete the gravitation acceleration that this pointMass poses on pm
	Vertex get_gravitation_acceleration_on_pm(PointMass *pm);

	//get electric field
	Vertex get_electric_field_at_position(Vertex pos);

	//aplly eletric field on acl
	void apply_electric_field_on_acl(Vertex EF);

	//get kinetic energy by KE = 0.5*M*V^2
	double get_kinetic_energy();

	//update the position, velocity and acceleration of point mass by simple euler integration
	//if the pointmass is added "theWorld", call "theWorld"'s update_state rather than this one
	//time step's unit is Second
	void update_state(double time_step);

	//same the update_state()
	//but this one use velocity_Verlet_integration
	void update_state_by_velocity_Verlet_integration(double time_step, Vertex old_acl);


	//to string method
	std::string toString();
};




/*
a spring connect two point mass objects
follow the Hookie's Law
FAQ: 
Q: if I want to connect a pointMass to fixed wall or something, how should I do?
A: set the mass of pm2 to a very huge value.a simpl e solution. 
*/
class SpringBond
{
public:
	//pointer to point masses
	PointMass *pm1;
	PointMass *pm2;

protected:
	//the zero-force length 
	double original_length;

	//according to Hookie's Law. F = K*x
	//K is the spring constant
	double spring_constant;

public:

	SpringBond();

	//constructor
	SpringBond(PointMass *pm1_p, PointMass *pm2_p, double original_length_p, double spring_constant_p);

	virtual Vertex get_force_on_pm1();
	virtual Vertex get_force_on_pm2();

	//get potential energy of the spring
	virtual double get_potential_energy();

	//tell if two spring equal
	virtual bool isEqual(SpringBond * other);
};



/*
a timer class used by "theWorld" class to keep track of time
*/
class TheWorldTimer
{
public:
	double elapsed_second;

	TheWorldTimer();

	//add time_step to elapsed_second
	void update_state(double time_step);

	int get_elapsed_femtosecond();

	int get_elapsed_nanosecond();

	int get_elapsed_picosecond();

	int get_elapsed_second();

	//tell how many minute passed
	int get_elapsed_minute();

	//tell how many days passed
	int get_elapsed_days();

	//tell how many year passed
	int get_elapsed_years();
};

enum IntergrationTypes
{
	simple_euler_intergration = 1, //also called "Rectangular Rule"
	velocity_Verlet_integration = 2, //Also called symplectic Euler, Leap-Frog
};


/*
a simulation to the real physics world where is not external force but the gravitation among point masses
FYI,"theWorld" is also the name of Dio's Stand from <JoJo's Bizarre Adventure> wwwwww
*/
class TheWorld
{
public:
	//a list containing point mass or rigid ball inside theWorld
	std::vector<PointMass *> point_mass_arr;

	//a list containing springBond inside theWorld
	std::vector<SpringBond *> spring_bond_arr;

	//internal timer to track how much time has passed
	TheWorldTimer timer;

	//defualt is false
	//on turn, object will affect each other according to the general gravitation
	//this function ,is especially for planet simulation
	bool gravitation_flag; 

	//defualt is false
	//on turn, object will affect each other by electriccharge
	bool electric_force_flag;

	//default is {0, 0, 0}
	//e.g the global acceleration of earth is {0,-9.8,0}
	Vertex global_acceleration;



	//a flag to enable/disable periodic box boundry
	//used to contrain objects in a box boundry
	//when an object collide or pass the boundries, it will bounce back
	bool periodic_box_boundry_flag;
	float x_negative_bnd, x_positive_bnd;
	float y_negative_bnd, y_positive_bnd;
	float z_negative_bnd, z_positive_bnd;


	//decide when update physical state which numeric intergration is used
	//different integrations have different performance cost, accuracy and stability
	IntergrationTypes intergration_type;


public:
	//constructor
	TheWorld();

	//default is {0, 0, 0}
	//e.g the global acceleration of earth is {0,-9.8,0}
	//set it to {0,-9.8,0}, all objects will have -9.8 m/s^2
	void set_global_acceleration(Vertex  gravitational_acceleration);

	//turn on or off gravitation calculation when "update_state(double time_step)"
	void set_gravitation_flag(bool on);



	//turn on or off the electric_force calculation   when "update_state(double time_step)"	
	void set_electric_force_flag(bool on);

	//turn on or off the periodic_box_boundry calculation   when "update_state(double time_step)
	void set_periodic_box_boundry_flag(bool on);

	//set the boundry of x, y, z
	//only work when periodic_box_boundry_flag is on
	//the defualt is from negative infinity to positive infinity
	void set_x_boundry(double x_negative_bnd_p, double x_positive_bnd_p);
	void set_y_boundry(double y_negative_bnd_p, double y_positive_bnd_p);
	void set_z_boundry(double z_negative_bnd_p, double z_positive_bnd_p);


	//default is rectangular rule
	//can set to others of enum IntergrationTypes
	//update_state will use the intergration type set
	void set_intergration_type(IntergrationTypes intergration_type_p);

	//add point mass to "point_mass_arr" of the world
	//so "update_state(double time_step)" will be effective on them
	//unexception will happen if call this after "update_state(double time_step)"
	void add_point_mass(PointMass * pm);


	//add spring bond
	//the two pms connected by sprintBond should be already added by "add_point_mass(PointMass * pm)"
	//unexception will happen if not added
	//unexception will happen if call this after "update_state(double time_step)"
	void add_spring_bond(SpringBond *spring_bond);

	//call this function  after adding point mass and etc 
	//will init correct acl of all objects
	//after init, must not change intergration_type
	void init();

	//update the position, velocity and acceleration of all object in theWorld
	//time step's unit is Second
	void update_state(double time_step);

	//to string method
	std::string toString();


	//get total kinetice energy of theworld
	double get_total_kinetic_energy();

	//get total spring_potential energy of theworld
	double get_total_spring_potential_energy();

	

private:

	//update state by using velocity verlet integration
	void update_state_by_velocity_Verlet_integration(double time_step);

	//update state by using rectangular method
	void update_state_by_rectangular(double time_step);

	//spring spring force on objects
	void apply_spring_force();

	//calculate gravitation and apply their acceleration
	void apply_gravitation_among_objects();


	//apply elertric force among object
	void apply_electric_among_objects();


	//apply all force method above according to the flags and setting
	//s.t change objects' acceleration
	void apply_force();
};




#endif /* defined(__HJYPhysics__) */
