//
//  CrystalSimulation.cpp
//
//  Created by Junyang Huang on 02/04/14.
//  Copyright (c) 2014 www.hjyssg.com. All rights reserved.
//



#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>      
#include <time.h>  
#include <ctime>
#include <cstdlib>


#include "HJYPhysics.h"

#define VPYTHON_FLAG 1   //a flag to disable/enable vpython, so this program can run on machine without vpython
#define PRINT_FLAG 0     //a flag to disable/enable print 


#if VPYTHON_FLAG 
#include "\Python32\include\Python.h" 
char cmd[300];  //a resuable string for vpython  
std::string tempCMD;

bool vpy_animation = true;
bool vpy_graph = true;

#endif 


//time step of simulation
double time_step = 0.01;  

const int MOLECULE_SIDE_NUM = 4;
const double MOLECULE_SIDE_LENGTH =2.2 M;
const double MOLECULE_MASS = 0.2 KG;
const double MOCLECULE_RADIUS = 0.2 M;
const double MOLECULE_SPEED_RANGE = 2 MPS;
const double MOLECULE_BOND_SPRING_CONST = 50.0;
const double MOCLECULE_BOND_LENGTH = sqrt(MOLECULE_SIDE_LENGTH*MOLECULE_SIDE_LENGTH*2)+0.1;

TheWorld world;
double  initial_total_energy;
clock_t tt;


void init_physics()
{
	CLOG("init physics...");
	tt = clock();
	srand(NULL);

	world = TheWorld();
	world.set_intergration_type(velocity_Verlet_integration);

	//generate point mass
	for(int ii = -MOLECULE_SIDE_NUM/2; ii < MOLECULE_SIDE_NUM/2; ii++)
	{
		//inner to outer
		for(int jj = -MOLECULE_SIDE_NUM/2; jj < MOLECULE_SIDE_NUM/2;jj++)
		{
			//left to right
			for(int kk = -MOLECULE_SIDE_NUM/2; kk < MOLECULE_SIDE_NUM/2; kk++)
			{
				float px = kk * MOLECULE_SIDE_LENGTH;
				float pz = jj * MOLECULE_SIDE_LENGTH;
				float py = ii * MOLECULE_SIDE_LENGTH;

				float vx = fRand(-MOLECULE_SPEED_RANGE, MOLECULE_SPEED_RANGE);
				float vy = fRand(-MOLECULE_SPEED_RANGE, MOLECULE_SPEED_RANGE);
				float vz = fRand(-MOLECULE_SPEED_RANGE, MOLECULE_SPEED_RANGE);

				PointMass *pm = new PointMass(px, py, pz, vx,vy, vz, 0,0,0, MOLECULE_MASS);
				//CLOG(pm->spd.y);

				world.add_point_mass(pm);
			}
		}
	}

	//generate spring bond
	for (int ii = 0; ii <world.point_mass_arr.size(); ii++)
	{
		PointMass *pm1 = world.point_mass_arr[ii];

		for (int jj = ii+1; jj < world.point_mass_arr.size(); jj++)
		{
			PointMass *pm2 = world.point_mass_arr[jj];
			double distance2 = pm2->pos.distance_square_to_v(pm1->pos);
			//generate spring in allowed range
			if (distance2 <= MOCLECULE_BOND_LENGTH*MOCLECULE_BOND_LENGTH)
			{
				//use init distance as spring length so the cube will not collapse
				SpringBond *spd = new SpringBond(pm1, pm2, sqrtf(distance2), MOLECULE_BOND_SPRING_CONST);

				world.add_spring_bond(spd);				
			}
		}
	}

	world.init();
	initial_total_energy = world.get_total_kinetic_energy()+world.get_total_spring_potential_energy();

	tt = clock() - tt;
	CLOG("physics init "<<tt<<" millisecond");

}

void init_vpython()
{
#if VPYTHON_FLAG
	CLOG("init VPython...");
	Py_Initialize();

	//calling PyRun_SimpleString will slowdown the program
	//build a string and pass once will speedup the program
	tempCMD.clear();
	tempCMD.append("from visual import *\n");
	tempCMD.append("from visual.graph import *\n");

	sprintf(cmd,"scene1=display(autocenter = True, x =400, autoscale = True, width=900, height=650)\n" );  //PyRun_SimpleString(cmd);
	tempCMD.append(cmd);

	//sprintf(cmd,"label = label(yoffset=200,  line=4, text=\"%i ball\\n%i springs\\nspd range %.3f\\nspring const %.3f\\ntime step %f\")\n", world.point_mass_arr.size(), world.spring_bond_arr.size(), MOLECULE_SPEED_RANGE, MOLECULE_BOND_SPRING_CONST, time_step );	tempCMD.append(cmd);

	sprintf(cmd,"xaxis = curve(pos=[(0,0,0), (20,0,0)], color=color.red)\n");		tempCMD.append(cmd);
	sprintf(cmd,"Yaxis = curve(pos=[(0,0,0), (0,20,0)], color=color.green)\n");		tempCMD.append(cmd);
	sprintf(cmd,"Zaxis = curve(pos=[(0,0,0), (0,0,20)], color=color.blue)\n");		tempCMD.append(cmd);


	for (int ii = 0; ii < world.point_mass_arr.size(); ii++)
	{
		PointMass *pm = world.point_mass_arr[ii];
		sprintf(cmd,"pointMass%i=sphere(radius=%e,pos=(%e,%e,%e), color = color.orange, opacity = 0.9)\n", ii, MOCLECULE_RADIUS ,pm->pos.x,pm->pos.y,pm->pos.z); 
		tempCMD.append(cmd);
	}

	CLOG("there are "<<world.spring_bond_arr.size()<<" springs and "<<world.point_mass_arr.size()<<" molecules");
	for (int ii = 0; ii < world.spring_bond_arr.size(); ii++)
	{
		SpringBond *spd = world.spring_bond_arr[ii];
		Vertex pos1 = spd->pm1->pos;
		Vertex axis =  spd->pm2->pos.sub(pos1);

		sprintf(cmd,"spring_%i = cylinder(radius=%e,pos=(%f,%f,%f),axis=(%f,%f,%f), color = color.green, opacity = 0.3)\n", ii, 0.06 ,pos1.x,pos1.y,pos1.z, axis.x,axis.y,axis.z); 
		tempCMD.append(cmd);
	}


	sprintf(cmd, "energy_graph=gdisplay( title = 'energy: red is KE, green is PE, blue is total,white is init \
				 total, orange is error percentage', xtitle='second',x=0, y=660, width=1800, height=400, ytitle='energy')\n");
	tempCMD.append(cmd);
	sprintf(cmd, "PE_line = gcurve(color=color.green)\n");tempCMD.append(cmd);
	sprintf(cmd, "KE_line = gcurve(color=color.red)\n");	tempCMD.append(cmd);
	sprintf(cmd, "TE_line = gcurve(color=color.blue)\n");	tempCMD.append(cmd);
	sprintf(cmd, "init_TE_line = gcurve(color=color.white)\n");	tempCMD.append(cmd);
	sprintf(cmd, "TE_error_line = gcurve(color=color.orange)\n");	tempCMD.append(cmd);

	PyRun_SimpleString(tempCMD.c_str());

#endif
}

void update_vpython()
{
#if VPYTHON_FLAG
	tempCMD.clear();

	if (vpy_animation)
	{
		//updte point mass ui 
		for (int ii = 0; ii < world.point_mass_arr.size(); ii++)
		{
			PointMass *pm = world.point_mass_arr[ii];
			sprintf(cmd,"pointMass%i.pos=(%e,%e,%e)\n", ii, pm->pos.x,pm->pos.y,pm->pos.z); 
			tempCMD.append(cmd);
		}
		PyRun_SimpleString(tempCMD.c_str()); tempCMD.clear();

		//update spring ui
		for (int ii = 0; ii < world.spring_bond_arr.size(); ii++)
		{
			SpringBond * spd = world.spring_bond_arr[ii];
			Vertex pos1 = spd->pm1->pos;
			Vertex axis =  spd->pm2->pos.sub(pos1);

			sprintf(cmd,"spring_%i.pos = (%e,%e,%e)\n", ii, pos1.x,pos1.y,pos1.z);   tempCMD.append(cmd);
			sprintf(cmd,"spring_%i.axis = (%e,%e,%e)\n", ii, axis.x,axis.y,axis.z);   tempCMD.append(cmd);
		}
		PyRun_SimpleString(tempCMD.c_str());tempCMD.clear();
	}

	if (vpy_graph)
	{
		double KE = world.get_total_kinetic_energy();
		double SPE = world.get_total_spring_potential_energy();

		sprintf(cmd, "KE_line.plot(pos=(%f, %e))\n",world.timer.elapsed_second, KE); tempCMD.append(cmd);  
		sprintf(cmd, "PE_line.plot(pos=(%f, %e))\n",world.timer.elapsed_second, SPE); tempCMD.append(cmd);
		sprintf(cmd, "TE_line.plot(pos=(%f, %e))\n",world.timer.elapsed_second,KE+SPE ); tempCMD.append(cmd);
		sprintf(cmd, "init_TE_line.plot(pos=(%f, %e))\n",world.timer.elapsed_second,initial_total_energy ); tempCMD.append(cmd);
	}
	sprintf(cmd,"rate(10)\n");
	tempCMD.append(cmd);


	//reduce pyrun_simpleString call, improve performance 10 times.
	PyRun_SimpleString(tempCMD.c_str());tempCMD.clear();

#endif
}



int main(int argc, const char * argv[])
{
	init_physics();
	init_vpython();

	double second = world.timer.elapsed_second;

	while(true)
	{
		static long double time = 0;

		world.update_state(time_step);

		if (second < world.timer.elapsed_second-0.02)
		{
			update_vpython();
			second = world.timer.elapsed_second;
		}
	}

	system ("pause");
	return 0;
}