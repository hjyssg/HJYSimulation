//
//  HJYPhysics.cpp
//  A simple physics engine
//
//  Created by Junyang Huang on 1/15/14.
//  Copyright (c) 2014 www.hjyssg.com. All rights reserved.
//

#include "HJYPhysics.h"

using namespace std;

bool is_in_range(double num, double min, double max)
{
	return num>min&&num<max;
}

double calculate_error_percentage(double estimate, double real)
{
	double difference = estimate - real;
	return difference/real*100.0;
}

double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

Vertex::Vertex()
{
	this->x = 0;
	this->y = 0;
	this->z = 0;
}

Vertex::Vertex(double x_p, double y_p, double z_p)
{
	this->x = x_p;
	this->y = y_p;
	this->z = z_p;
}

void Vertex::set_xyz(double x_p, double y_p, double z_p)
{
	this->x = x_p;
	this->y = y_p;
	this->z = z_p;
}

Vertex  Vertex::add(Vertex v2)
{
	Vertex result = Vertex(0,0,0);;
	result.x = x + v2.x;
	result.y = y + v2.y;
	result.z = z + v2.z;
	return result;
}

void   Vertex::add_self(Vertex v2)
{
	this->x +=  v2.x;
	this->y +=  v2.y;
	this->z +=  v2.z;
}


Vertex Vertex::sub(Vertex v2)
{
	Vertex result = Vertex(0,0,0);;
	result.x = x - v2.x;
	result.y = y - v2.y;
	result.z = z - v2.z;
	return result;
}


void  Vertex::sub_self(Vertex v2)
{
	this->x -=  v2.x;
	this->y -=  v2.y;
	this->z -=  v2.z;
}



Vertex Vertex::crosss_product(Vertex v2)
{
	Vertex result = Vertex(0,0,0);
	result.x =   y * v2.z - z*v2.y;
	result.y = -(x * v2.z - z*v2.x);
	result.z =   x * v2.y - y*v2.x;
	return result;
}

double Vertex::dot_product(Vertex v2)
{
	return  x * v2.x +  y * v2.y +  z * v2.z ;
}

Vertex   Vertex::transfer( double x_p, double y_p, double z_p)
{
	Vertex result = Vertex(0,0,0);;
	result.x = x + x_p;
	result.y = y + y_p;
	result.z = z + z_p;
	return result;
}

void  Vertex::transfer_self( double x_p, double y_p, double z_p)
{
	this->x += x_p;
	this->y += y_p;
	this->z += z_p;
}



//vtx is origined from the (0,0,0)
Vertex Vertex::rotate_around_axis(Vertex vtx, float theta)
{
	//if angle is 0 no need to calculate
	if(theta==0)
	{
		return get_copy();
	}

	static Vertex cache_v = Vertex();
	static float cache_theta = 0.0;

	static Vertex r1 = Vertex();
	static Vertex r2 = Vertex();
	static Vertex r3 = Vertex();

	//only calculate matrix once for each vertex and degress combination
	if (!(vtx.isEqual(cache_v)&&cache_theta==theta))
	{
		float cc = cos(theta);
		float ss = sin(theta);

		float _cc = 1-cc;

		float xy = vtx.x * vtx.y;
		float zy = vtx.z * vtx.y;
		float zx = vtx.z * vtx.x;

		float xx = vtx.x * vtx.x;
		float yy = vtx.y * vtx.y;
		float zz = vtx.z * vtx.z;

		r1 = Vertex( cc+ xx* _cc,
			xy*_cc - vtx.z*ss,
			zx*_cc + vtx.y*ss);

		r2 = Vertex( xy*_cc + vtx.z*ss,
			cc + yy *_cc,
			zy*_cc- vtx.x*ss
			);

		r3 = Vertex(zx *_cc - vtx.y*ss,
			zy *_cc + vtx.x*ss,
			cc +  zz*_cc
			);

		cache_theta = theta;
		cache_v = vtx;

	}

	Vertex result;
	result.x = this->dot_product(r1);
	result.y = this->dot_product(r2);
	result.z = this->dot_product(r3);

	return result;
}

Vertex  Vertex::scale( double x_s, double y_s, double z_s)
{
	Vertex result = Vertex(0,0,0);;
	result.x = x * x_s;
	result.y = y * y_s;
	result.z = z * z_s;
	return result;
}

void  Vertex::scale_self( double x_s, double y_s, double z_s)
{
	this->x *= x_s;
	this->y *= y_s;
	this->z *= z_s;
}

void  Vertex::scale_self( double p)
{
	this->x *= p;
	this->y *= p;
	this->z *= p;
}

Vertex Vertex::normalize()
{
	float len = this->length();
	float temp = 1/len;
	if (len==0)
	{
		temp = 1;
	}

	Vertex result;
	result.x = x * temp;
	result.y = y * temp;
	result.z = z * temp;

	return result;
}

void Vertex::normalize_self()
{
	float len = this->length();
	if (len==0)
	{
		return;
	}

	float temp = 1/len;

	x *= temp;
	y *= temp;
	z *= temp;
}

double Vertex::length()
{
	return  sqrt(x*x + y*y + z*z);
}

double  Vertex::length_square()
{
	return  x*x + y*y + z*z;
}

double Vertex::distance_to_v(Vertex v)
{
	double temp = (x-v.x)*(x-v.x) + (y-v.y)*(y-v.y)+(z-v.z)*(z-v.z);
	return sqrt(temp);
}



double Vertex::distance_square_to_v(Vertex v)
{
#if 0
	return sub(v).length_square();
#else
	double temp = (x-v.x)*(x-v.x) + (y-v.y)*(y-v.y)+(z-v.z)*(z-v.z);
	return temp;
#endif

}




bool Vertex::isZero()
{
	return x==0.0&&y==0.0&&z==0.0;
}

bool Vertex::isEqual(Vertex v)
{
	return (this->x==v.x && this->y==v.y && this->z==v.z);
}

bool Vertex::isParrell(Vertex v)
{
	//zero vector is parrell to any vector
	if (v.isZero()||this->isZero())
	{
		return true;
	}

	double r1 = this->x / v.x;
	double r2 = this->y / v.y;
	double r3 = this->z / v.z;

	return (r1==r2 && r2 ==r3&&r3 == r1);
}

Vertex Vertex::rotate_x(float theta)
{
	Vertex result;
	result.x = x;
	result.y = y*cos(theta) - sin(theta)*z;
	result.z = y*sin(theta) + cos(theta)*z;
	return result;
}

void Vertex::rotate_x_self(float theta)
{
	float tempy = y;
	float tempz = z;
	y = tempy*cos(theta) - sin(theta)*tempz;
	z = tempy*sin(theta) + cos(theta)*tempz;
}

Vertex Vertex::rotate_y(float theta)
{
	Vertex result;
	result.x = x*cos(theta) + sin(theta)*z;
	result.y = y;
	result.z = -x*sin(theta) + cos(theta)*z;
	return result;
}

void Vertex::rotate_y_self(float theta)
{
	float tempx = x;
	float tempz = z;

	x = tempx*cos(theta) + sin(theta)*tempz;
	z = -tempx*sin(theta) + cos(theta)*tempz;
}

Vertex Vertex::rotate_z(float theta)
{
	Vertex result;
	result.x = x*cos(theta) - sin(theta)*y;
	result.y = x*sin(theta) + cos(theta)*y;
	result.z = z;
	return result;
}

void Vertex::rotate_z_self(float theta)
{
	float tempx = x;
	float tempy = y;
	x = tempx*cos(theta) - sin(theta)*tempy;
	y = tempx*sin(theta) + cos(theta)*tempy;
}

Vertex  Vertex::get_copy()
{
	return Vertex(x,y,z);
}


ShapeType PointMass::get_shape()
{
	return none;
}

std::string Vertex::toString()
{
	static char tempStr[400];
	sprintf(tempStr, "[%.3e, %.3e,%.3e]", x,y,z);
	std::string result(tempStr);
	return result;
}

PointMass::PointMass()
{
	mass = 0.0;
	electric_charge = 0.0;
	pos = Vertex();
	spd = Vertex();
	acl = Vertex();
}


PointMass::PointMass(double p_x, double p_y, double p_z,
	double spd_x,double spd_y,double spd_z,
	double acl_x,double acl_y,double acl_z,
	double mass_p)
{
	this->mass = mass_p;
	this->collision_damping = 0.0;
	electric_charge = 0.0;
	pos = Vertex(p_x, p_y, p_z);
	spd = Vertex(spd_x, spd_y, spd_z);
	acl = Vertex(acl_x, acl_y, acl_z);
}


PointMass::PointMass(double p_x, double p_y, double p_z,
	double spd_x,double spd_y,double spd_z,
	double acl_x,double acl_y,double acl_z,
	double mass_p,
	double collision_damping)
{
	this->mass = mass_p;
	electric_charge = 0.0;
	this->collision_damping = collision_damping;
	pos = Vertex(p_x, p_y, p_z);
	spd = Vertex(spd_x, spd_y, spd_z);
	acl = Vertex(acl_x, acl_y, acl_z);
	this->collision_damping =  collision_damping;
}

//the gravitation that this poses on pm
//this  <-attract--- pm
Vertex PointMass::get_gravitation_acceleration_on_pm(PointMass *pm)
{
	Vertex direction = pos.sub(pm->pos);
	double distance_square = direction.length_square();
	double distance = sqrt(distance_square);

	double scalar_g = GRAVIATIONAL_CONSTANT*mass/distance_square;

	// normalized and then times scalar
	double temp = scalar_g/distance;
	direction.scale_self(temp, temp, temp);

	return direction;
}

Vertex PointMass::get_electric_field_at_position(Vertex pos)
{
	Vertex direction = this->pos.sub(pos);
	double distance_square = direction.length_square();
	double distance = sqrt(distance_square);

	double scalar_g = COULOMBS_CONSTANT*electric_charge/distance_square;

	// normalized and then times scalar
	double temp = -scalar_g/distance;
	direction.scale_self(temp, temp, temp);
	return direction;
}

void PointMass::apply_electric_field_on_acl(Vertex EF)
{
	double temp  = this->electric_charge/mass;
	EF.scale_self(temp,temp,temp);
	this->acl.add_self(EF);	
}

double PointMass::get_kinetic_energy()
{
	double temp = spd.x*spd.x + spd.y*spd.y + spd.z*spd.z;
	return 0.5*temp*mass;
}

void PointMass::update_state(double time_step)
{
	pos.x += spd.x*time_step;
	pos.y += spd.y*time_step;
	pos.z += spd.z*time_step;

	spd.x += acl.x*time_step;
	spd.y += acl.y*time_step;
	spd.z += acl.z*time_step;
}

void PointMass::update_state_by_velocity_Verlet_integration(double time_step, Vertex old_acl)
{
	pos.x += spd.x*time_step + 0.5*time_step*time_step*old_acl.x;
	pos.y += spd.y*time_step + 0.5*time_step*time_step*old_acl.y;
	pos.z += spd.z*time_step + 0.5*time_step*time_step*old_acl.z;

	spd.x += (old_acl.x+acl.x)*time_step/2;
	spd.y += (old_acl.y+acl.y)*time_step/2;
	spd.z += (old_acl.z+acl.z)*time_step/2;
}


std::string PointMass::toString()
{
	static char tt[100];
	sprintf(tt," mass: %.3e", mass);

	std::string temp;
	temp.append("pos: ");
	temp.append(pos.toString());
	temp.append(" spd: ");
	temp.append(spd.toString());
	temp.append(" acl: ");
	temp.append(acl.toString());
	temp.append(tt);
	return temp;
}

//use point mass parent constructor
RigidSphere::RigidSphere():PointMass()
{
	this->radius = 0.0;
}


ShapeType RigidSphere::get_shape()
{
	return sphere;
}

SpringBond::SpringBond()
{
}


SpringBond::SpringBond(PointMass *pm1_p, PointMass *pm2_p, double original_length_p, double spring_constant_p)
{
	this->pm1 = pm1_p;
	this->pm2 = pm2_p;
	this->original_length = original_length_p;
	this->spring_constant = spring_constant_p;
}

Vertex SpringBond::get_force_on_pm1()
{
	//pm1 -> pm2
	Vertex directionV = pm2->pos.sub(pm1->pos);
	double distance = directionV.length();

	//delta X
	double difference = distance - this->original_length;
	//double difference =   this->original_length - distance;


	//normaliza direction

	directionV.normalize_self();

	//F = K* delta X
	double scal = difference *  spring_constant;
	directionV.scale_self(scal, scal, scal);

	//CLOG("scal "<< scal << " difference"<<difference<<"spr const"<<spring_constant <<" force"<<directionV.toString());
	return directionV;
}

Vertex SpringBond::get_force_on_pm2()
{
	Vertex temp =  get_force_on_pm1();
	temp.scale_self(-1,-1,-1);
	return temp;
}

double SpringBond::get_potential_energy()
{
	double length = pm2->pos.distance_to_v(pm1->pos);
	double difference = length - original_length;

	double PE = 0.5*spring_constant*difference*difference;
	return PE;
}

bool SpringBond::isEqual( SpringBond * other)
{
	bool f1 = (other->original_length == this->original_length);
	bool f2 = (other->spring_constant == this->spring_constant);

	bool f3 = (other->pm1 == this->pm1)&&(other->pm2 == this->pm2);
	bool f4 = (other->pm1 == this->pm2)&&(other->pm2 == this->pm1);

	if ((f3&&f1&&f2)||(f4&&f1&&f2))
	{
		return true;
	}else
	{
		return false;
	}
}


PowerSpringBond::PowerSpringBond(PointMass *pm1_p, PointMass *pm2_p, double original_length_p, double spring_constant_p,int power_p)
{
	this->pm1 = pm1_p;
	this->pm2 = pm2_p;
	this->original_length = original_length_p;
	this->spring_constant = spring_constant_p;
	this->power = power_p;
}

Vertex PowerSpringBond::get_force_on_pm1()
{
	//pm1 -> pm2
	Vertex directionV = pm2->pos.sub(pm1->pos);
	double distance = directionV.length();

	//delta X
	double difference = distance - this->original_length;
	

	//normaliza direction

	directionV.normalize_self();

	//F = K* delta X
	long double temp = pow(difference, this->power);
	long double scal = temp* spring_constant;
	directionV.scale_self(scal, scal, scal);

	//CLOG(difference<<" "<<this->power<<" "<<spring_constant);
	//CLOG(temp<<" "<<scal<<directionV.toString());

	return directionV;
}


//get potential energy of the spring
double PowerSpringBond::get_potential_energy()
{
	double length = pm2->pos.distance_to_v(pm1->pos);
	double difference = length - original_length;

	double PE = 0.5*spring_constant* pow(difference, power)* pow(difference, power);
	return PE;
}

Molecule::Molecule()
{

}

void Molecule::add_atom(PointMass *atom)
{
	atoms.push_back(atom);
}

void Molecule::add_atom_bond(SpringBond * bond)
{
	atom_bonds.push_back(bond);
}

double  Molecule::get_mass()
{
	double mass = 0;
	for (int ii = 0; ii < atoms.size(); ++ii)
	{
		PointMass *atm = atoms[ii];
		mass  += atm->mass;
	}
	return mass;
}

Vertex Molecule::get_position()
{
	//http://hyperphysics.phy-astr.gsu.edu/hbase/cm.html
	Vertex pos = Vertex();
	for (int ii = 0; ii < atoms.size(); ++ii)
	{
		PointMass *atm = atoms[ii];
		double tm = atm->mass;
		Vertex * p = &atm->pos;
		pos.transfer_self(p->x*tm,p->y*tm,p->z*tm);
	}

	pos.scale_self(1.0/this->get_mass());
	return pos;
}

Vertex Molecule::get_acceleration()
{
	Vertex acl = Vertex();
	for (int ii = 0; ii < atoms.size(); ++ii)
	{
		PointMass *atm = atoms[ii];
		acl.add_self(atm->acl);
	}

	acl.scale_self(1.0/(double)atoms.size());
	return acl;
}

void Molecule::apply_force_on_molecule(Vertex force)
{
	for (int ii = 0; ii < atoms.size(); ++ii)
	{
		PointMass *atm = atoms[ii];
		force.scale_self(1/atm->mass);
		atm->acl.add_self(force);
	}
}
	

TheWorld::TheWorld()
{
	global_acceleration = Vertex(0,0,0);
	timer = TheWorldTimer();
	periodic_box_boundry_flag = false;
	gravitation_flag = false;
	electric_force_flag = false;
	rigid_body_collision_flag = false;
	lennard_jones_potential_flag = false;
	x_negative_bnd = -1.7e+308; 
	x_positive_bnd = 1.7e+308; 
	y_negative_bnd = -1.7e+308; 
	y_positive_bnd = 1.7e+308; 
	z_negative_bnd = -1.7e+308; 
	z_positive_bnd = 1.7e+308; 
	intergration_type = simple_euler_intergration;
}

void TheWorld::add_point_mass(PointMass * pm)
{
	point_mass_arr.push_back(pm);
}

void TheWorld::add_spring_bond(SpringBond *spring_bond)
{
	spring_bond_arr.push_back(spring_bond);
}

void TheWorld::add_molecule(Molecule * molecule)
{
	for (int ii = 0; ii < molecule->atoms.size(); ++ii)
	{
		PointMass *pm = molecule->atoms[ii];
		this->add_point_mass(pm);
	}

	for (int ii = 0; ii < molecule->atom_bonds.size(); ++ii)
	{
		SpringBond *bnd = molecule->atom_bonds[ii];
		this->add_spring_bond(bnd);
	}

	molecule_arr.push_back(molecule);
}




void TheWorld::set_intergration_type(IntergrationTypes intergration_type_p)
{
	intergration_type = intergration_type_p;
}

void TheWorld::set_electric_force_flag(bool on)
{
	electric_force_flag = on;
}

void  TheWorld::set_global_acceleration(Vertex  gravitational_acceleration)
{
	this->global_acceleration = gravitational_acceleration;
}

void TheWorld::set_gravitation_flag(bool on)
{
	this->gravitation_flag = on;
}
void TheWorld::set_rigid_body_collision_flag(bool on)
{
	this->rigid_body_collision_flag = on;
}


void TheWorld::set_lennard_jones_potential_flag(bool on)
{
	this->lennard_jones_potential_flag = on;
}

void TheWorld::set_periodic_box_boundry_flag(bool on)
{
	this->periodic_box_boundry_flag = on;
}

void TheWorld::set_x_boundry(double x_negative_bnd_p, double x_positive_bnd_p)
{
	this->x_negative_bnd = x_negative_bnd_p;
	this->x_positive_bnd = x_positive_bnd_p;
}

void TheWorld::set_y_boundry(double y_negative_bnd_p, double y_positive_bnd_p)
{
	this->y_negative_bnd = y_negative_bnd_p;
	this->y_positive_bnd = y_positive_bnd_p;
}

void TheWorld::set_z_boundry(double z_negative_bnd_p, double z_positive_bnd_p)
{
	this->z_negative_bnd = z_negative_bnd_p;
	this->z_positive_bnd = z_positive_bnd_p;
}

#if 0
void TheWorld::apply_gravitation_among_objects()
{
	//use (n^2)/2 algo now to firgure out graviation
	#pragma omp parallel for num_threads(4)
	for(int ii = 0; ii < point_mass_arr.size()-1;ii++)
	{
		PointMass *pm  = point_mass_arr.at(ii);

		for(int jj = ii+1; jj <  point_mass_arr.size();jj++)
		{
			PointMass *pm_2  =  point_mass_arr.at(jj);

			Vertex g =  pm_2->get_gravitation_acceleration_on_pm(pm);

			

			double temp = -pm->mass/pm_2->mass;
			Vertex g2 = g.scale(temp,temp, temp);
#pragma omp critical(apply_spring_force_lock)
		{
			pm->acl.transfer_self(g.x, g.y, g.z);
			pm_2->acl.transfer_self(g2.x, g2.y, g2.z); 
		}
		}
	}
}


#else
void TheWorld::apply_gravitation_among_objects()
{
	//CLOG("ssda");
	//use (n^2) algo now to firgure out graviation
#pragma omp parallel for num_threads(4)
	for(int ii = 0; ii < point_mass_arr.size();ii++)
	{
		PointMass *pm  = point_mass_arr.at(ii);

		for(int jj = 0; jj <  point_mass_arr.size();jj++)
		{
			if (ii==jj){continue;}

			PointMass *pm_2  =  point_mass_arr.at(jj);
			Vertex f =  pm_2->get_gravitation_acceleration_on_pm(pm);

			pm->acl.add_self(f);
		}
	}
}

#endif 

void TheWorld::apply_electric_among_objects()
{
	#pragma omp parallel for num_threads(4)
	for(int ii = 0; ii < point_mass_arr.size();ii++)
	{
		PointMass *pm  = point_mass_arr.at(ii);

		for(int jj = 0; jj <  point_mass_arr.size();jj++)
		{
			if (ii==jj){continue;}

			PointMass *pm_2  =  point_mass_arr.at(jj);

			Vertex fe =  pm_2->get_electric_field_at_position(pm->pos);

			fe.scale_self(1/120.0); //manual adjusting

			pm->apply_electric_field_on_acl(fe);
		}
	}
}

void TheWorld::apply_lennard_jones_potential()
{
	#pragma omp parallel for num_threads(4)
	for(int ii = 0; ii < molecule_arr.size();ii++)
	{
		Molecule *m1 = molecule_arr[ii];
		Vertex pos1 = m1->get_position();

		for(int jj = 0; jj < molecule_arr.size();jj++)
		{
			if (ii==jj){continue;}

			Molecule *m2 = molecule_arr[jj];
			Vertex pos2 = m2->get_position();

			Vertex force_by_pm2 = pos1.sub(pos2);
			double distance = force_by_pm2.length();

			//from http://www.numfys.net/examples/ex1_1_lennardjones.pdf
			const double epsilon  = 0.32e-9;
			const double sigma = 1.08e-21;
			

			//refer to http://www2.physics.umd.edu/~alaporta/Lennard-Jones.html
			

			//12σ12/r13 
			double temp1 = 12*(pow(sigma,12)/pow(distance,13));

			//6σ6/r7	
			double temp2 = 6*(pow(sigma,6)/pow(distance,7));
			//R(r) = 4ε [12σ12/r13 - 6σ6/r7]
			double R =4* epsilon * (temp1 - temp2); 


			R *=10; //manual adjusting



			force_by_pm2.scale_self(R);

			m1->apply_force_on_molecule(force_by_pm2);

		}
	}
}

void TheWorld::apply_spring_force()
{
	//apply spring bond to the acceleration
	//int size = ;
//#pragma omp parallel for num_threads(4)
	for (int ii = 0; ii < spring_bond_arr.size(); ++ii)
	{
		SpringBond * spd = spring_bond_arr[ii];
		PointMass *temp_pm1 = spd->pm1;
		PointMass *temp_pm2 = spd->pm2;

		Vertex f1 = spd->get_force_on_pm1();
		assert(temp_pm1->mass>0);
		double temp = 1/temp_pm1->mass;
		Vertex a1 = f1.scale(temp,temp,temp);

		//the force 2 is opposite direction
		assert(temp_pm2->mass>0);
		double temp2 = -1/temp_pm2->mass; 
		Vertex a2 = f1.scale(temp2,temp2,temp2);

//#pragma omp critical(apply_spring_force_lock)
		{
			temp_pm1->acl.add_self(a1);
			temp_pm2->acl.add_self(a2);
		}
	}
}

void TheWorld::apply_force()
{
	//init new acceleration
	for(int ii = 0; ii < point_mass_arr.size();ii++)
	{
		PointMass *pm  = point_mass_arr.at(ii);
		//pm->acl = global_acceleration;
		pm->acl.x = global_acceleration.x;
		pm->acl.y = global_acceleration.y;
		pm->acl.z = global_acceleration.z;
		//CLOG(pm->acl.toString());
	}

	//calculate gravitation acceleration
	if (gravitation_flag)
	{
		apply_gravitation_among_objects();
	}

	if (electric_force_flag)
	{
		apply_electric_among_objects();
	}

	if (lennard_jones_potential_flag)
	{
		apply_lennard_jones_potential();
	}

	apply_spring_force();
}

double TheWorld::get_total_kinetic_energy()
{
	double total_kinetic_energy = 0;
	for(int ii = 0; ii < point_mass_arr.size();ii++)
	{
		PointMass *pm  = point_mass_arr.at(ii);
		total_kinetic_energy += pm->get_kinetic_energy();
	}
	return total_kinetic_energy;
}


double TheWorld::get_total_spring_potential_energy()
{
	double total_spring_potential_energy = 0.0;

	for (int ii = 0; ii < this->spring_bond_arr.size(); ii++)
	{
		SpringBond * spd = this->spring_bond_arr[ii];
		double single_spe = spd->get_potential_energy();
		total_spring_potential_energy += single_spe;
	}
	return total_spring_potential_energy;
}

void TheWorld::update_state_by_velocity_Verlet_integration(double time_step)
{
	//refer to http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
	apply_force();

	double dt = time_step;
	double half_dt = 0.5*time_step;

#pragma omp parallel for num_threads(4)
	for(int ii = 0; ii < point_mass_arr.size();ii++)
	{
		PointMass *pm  = point_mass_arr.at(ii);

		//get v(t + 1/2*dt)
		pm->spd.transfer_self(pm->acl.x*half_dt,
			pm->acl.y*half_dt,
			pm->acl.z*half_dt);

		//x(t + dt)
		pm->pos.transfer_self(pm->spd.x*dt,
			pm->spd.y*dt,
			pm->spd.z*dt);
	}

	//get a(t+dt)
	apply_force();

#pragma omp parallel for num_threads(4)
	for(int ii = 0; ii < point_mass_arr.size();ii++)
	{
		PointMass *pm  = point_mass_arr.at(ii);
		//get v(t + dt)
		pm->spd.transfer_self(pm->acl.x*half_dt,
			pm->acl.y*half_dt,
			pm->acl.z*half_dt);
	}
}

void TheWorld::update_state_by_rectangular(double time_step)
{
	apply_force();

	//update velocity and velocity
	#pragma omp parallel for num_threads(4)
	for(int ii = 0; ii < point_mass_arr.size();ii++)
	{
		PointMass *pm  = point_mass_arr.at(ii);
		pm->update_state(time_step);		
	}
}

void TheWorld::update_state_by_simposon_rule(double time_step)
{
	//CLOG("HAVE NOT SUPPORT SIMPSON INTEGRATION YET!!");
	//refer to http://en.wikipedia.org/wiki/Simpson's_rule
	// = (b-a)/6[f(a)+4f((a+b)/2)+f(b)]

	apply_force();

	double one_sixth_dt = time_step/6.0;

	#pragma omp parallel for num_threads(4)
	for(int ii = 0; ii < point_mass_arr.size();ii++)
	{
		PointMass *pm  = point_mass_arr.at(ii);
		Vertex old_acl = pm->acl;
		Vertex old_spd = pm->spd;

		//calculate new speed
		//f(a)+f(b)
		Vertex spd_integral = pm_acl_one_time_step_ago[ii]->transfer(old_acl.x,
																	 old_acl.y,
																	 old_acl.z);
		//+4f((a+b)/2)
		spd_integral.transfer_self(pm_acl_half_time_step_ago[ii]->x*4,
								   pm_acl_half_time_step_ago[ii]->y*4,
								   pm_acl_half_time_step_ago[ii]->z*4);
		//*(b-a)/6
		spd_integral.scale_self(one_sixth_dt,one_sixth_dt,one_sixth_dt);

		//calculate new position
		//f(a)+f(b)
		Vertex pos_integral = pm_spd_one_time_step_ago[ii]->transfer(spd_integral.x,
																	 spd_integral.y,
																	 spd_integral.z);
		//+4f((a+b)/2)
		pos_integral.transfer_self(pm_spd_half_time_step_ago[ii]->x*4,
								   pm_spd_half_time_step_ago[ii]->y*4,
								   pm_spd_half_time_step_ago[ii]->z*4);
		//*(b-a)/6
		pos_integral.scale_self(one_sixth_dt,one_sixth_dt,one_sixth_dt);

		//save the old speed and old acl for next time using
		pm_spd_one_time_step_ago[ii] = pm_spd_half_time_step_ago[ii];
		pm_acl_one_time_step_ago[ii] = pm_acl_half_time_step_ago[ii];
		pm_spd_half_time_step_ago[ii] = &old_spd;
		pm_acl_half_time_step_ago[ii] = &old_acl;

		pm->spd.add_self(spd_integral);
		pm->pos.add_self(pos_integral);
	}
}


void TheWorld::init()
{
	this->apply_force();

	if (intergration_type == simpsons_rule_integration)
	{	
		//init fot later simplson update_state calculation
		for (int ii = 0; ii < point_mass_arr.size(); ++ii)
		{
			PointMass *pm = point_mass_arr[ii];
			pm_spd_one_time_step_ago.push_back(&pm->spd);
			pm_spd_half_time_step_ago.push_back(&pm->spd);

			pm_acl_one_time_step_ago.push_back(&pm->acl);
			pm_acl_half_time_step_ago.push_back(&pm->acl);

			//CLOG(pm_spd_one_time_step_ago[ii]->toString()<<" "<<pm_spd_half_time_step_ago[ii]->toString());
			//CLOG(pm_acl_one_time_step_ago[ii]->toString()<<" "<<pm_acl_half_time_step_ago[ii]->toString());
		}
	}
}

//assert at least two point mass
void TheWorld::update_state(double time_step)
{
	//update state according to the intergration specified
	if (intergration_type == simple_euler_intergration)
	{
		update_state_by_rectangular(time_step);
	}else if (intergration_type == velocity_Verlet_integration)
	{
		update_state_by_velocity_Verlet_integration(time_step);
	}else if (intergration_type ==  simpsons_rule_integration)
	{
		update_state_by_simposon_rule(time_step);
	}
	else
	{
		update_state_by_rectangular(time_step);
	}

	if (periodic_box_boundry_flag)
	{
		//CLOG("checking");
		for (int ii = 0; ii < point_mass_arr.size(); ++ii)
		{
			PointMass *pm = point_mass_arr[ii];

			//bounce back if out of boundry
			if (pm->pos.x < x_negative_bnd)
			{
				pm->spd.x = -pm->spd.x;
				pm->pos.x = x_negative_bnd;
			}
			else if (pm->pos.x > x_positive_bnd)
			{
				pm->spd.x = -pm->spd.x;
				pm->pos.x = x_positive_bnd;
			}


			if (pm->pos.y < y_negative_bnd)
			{
				pm->spd.y = -pm->spd.y;
				pm->pos.y = y_negative_bnd;
			}
			else if (pm->pos.y > y_positive_bnd)
			{
				pm->spd.y = -pm->spd.y;
				pm->pos.y = y_positive_bnd;
			}

			if (pm->pos.z < z_negative_bnd)
			{
				pm->spd.z = -pm->spd.z;
				pm->pos.z = z_negative_bnd;
			}
			else if (pm->pos.z > z_positive_bnd)
			{
				pm->spd.z = -pm->spd.z;
				pm->pos.z = z_positive_bnd;
			}
			
		}
	}

	//update internal timers
	this->timer.update_state(time_step);
}


void TheWorld::apply_collision_effect()
{
	for(int ii = 0; ii < point_mass_arr.size()-1;ii++)
	{
		PointMass *pm  = point_mass_arr.at(ii);
		ShapeType s1 = pm->get_shape();
		if (s1==none){continue;}

		for(int jj = ii+1; jj < point_mass_arr.size();jj++)
		{
			PointMass *pm_2  = point_mass_arr.at(jj);
			ShapeType s2 = pm_2->get_shape();
			if (s2==none){continue;}

			Vertex tempV = pm->pos.sub(pm_2->pos);
			double distanceSquare = tempV.length_square();
			double distance = sqrt(distanceSquare);
			if (s1==sphere&&s2==sphere)
			{
				RigidSphere * ball1 = (RigidSphere *)pm;
				RigidSphere * ball2 = (RigidSphere *)pm_2;

				if (distance < (ball1->radius+ball2->radius))
				{
					//wrong
				}
			}
		}
	}
}

std::string TheWorld::toString()
{
	std::string str;
	str.append("{ ");
	for (int ii = 0; ii < point_mass_arr.size(); ii++)
	{
		PointMass *pm  = point_mass_arr.at(ii);
		str.append(" [" );
		str.append(pm->toString());

		if (ii != point_mass_arr.size()-1)
		{
			str.append("], \n");
		}else
		{
			str.append("]");
		}
	}

	str.append(" }");
	return str;
}


TheWorldTimer::TheWorldTimer()
{
	elapsed_second = 0.0;
}

void TheWorldTimer::update_state(double time_step)
{
	elapsed_second += time_step;
}

int TheWorldTimer::get_elapsed_femtosecond()
{
	return (int)(elapsed_second/(1 FS));
}

int TheWorldTimer::get_elapsed_nanosecond()
{
	
	return (int)(elapsed_second/(1 NS));
}

int TheWorldTimer::get_elapsed_picosecond()
{
	CLOG(elapsed_second);
	return (int)(elapsed_second/(1 PS));
}


int TheWorldTimer::get_elapsed_second()
{
	int current_time = (int) elapsed_second;
	return current_time/(1 S);
}

int TheWorldTimer::get_elapsed_years()
{
	int current_time = (int) elapsed_second;
	return current_time/(31536000 S); //31536000 is 365*24*3600, second of one year;
}

int TheWorldTimer::get_elapsed_days()
{
	int current_time = (int) elapsed_second;
	return current_time/(86400 S); //86400 is 24*3600, second of one day
}

int TheWorldTimer::get_elapsed_minute()
{
	int current_time = (int) elapsed_second;
	return current_time/(60 S);
}