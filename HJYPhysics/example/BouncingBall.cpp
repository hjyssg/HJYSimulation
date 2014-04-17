//
//  main.cpp
//  Ball_Bouncing
//
//  Created by Junyang Huang on 1/15/14.
//  Copyright (c) 2014 www.hjyssg.com. All rights reserved.
//

#include <iostream>
#include "HJYPhysics.h"

const double HEIGHT = 80.0 M;  /* initial height of ball in meters*/
const double INITIAL_X_POS  = -100 M;  //initial x position of ball
const double INITIAL_X_SPEED =  15 MPS;


#define VPYTHON_FLAG 1   //a flag to disable/enable vpython, so this program can run on machine without vpython
#define PRINT_FLAG 0     //a flag to disable/enable print 

#if VPYTHON_FLAG 
    #include "\Python32\include\Python.h" 
    char cmd[100];  //a resuable string for vpython  
#endif 

float time_step = 0.05;  //time step of simulation


int main(int argc, const char * argv[])
{
	TheWorld world = TheWorld();
	world.set_global_acceleration(Vertex(0,-9.8,0));
	world.set_periodic_box_boundry_flag(true);
	world.set_y_boundry(0, 100);


	PointMass * pm1 = new PointMass(INITIAL_X_POS, HEIGHT, 10, //pos 
		INITIAL_X_SPEED, 0, 0,      //velocity
		0, 0, 0,    //acl
		1, 1);     //mass, hitting factor    



	world.add_point_mass(pm1);
	world.init();
	

#if VPYTHON_FLAG

	/* initialize graphics */
	Py_Initialize();
	PyRun_SimpleString("from visual import *\n");

	/* set the viewing position and range/scale */
	sprintf(cmd,"scene1=display(autocenter = 0, width=800, height=600,autoscale = 1, center=(0,%f,0),range=%f)\n",HEIGHT/2, (HEIGHT)+30.0);
	PyRun_SimpleString(cmd);

	//create visual ball for pointmass objects
	sprintf(cmd,"ball1=sphere(color=color.green,make_trail=True,radius=3.0,pos=(%f,%f,%f))\n", pm1->pos.x,pm1->pos.y,pm1->pos.z);
	PyRun_SimpleString(cmd); 

	/* put a landing pad at 0,0,0 */
	PyRun_SimpleString("box(length=10000.0, height=0.5, width=100.0, color=color.blue)\n");

#endif

	while(true)
	{
#if PRINT_FLAG
		CLOG(ii*time_step <<"," << world.toString() );
#endif

#if VPYTHON_FLAG
		sprintf(cmd,"ball1.pos=(%f,%f,%f)\n", pm1->pos.x,pm1->pos.y,pm1->pos.z);
		PyRun_SimpleString(cmd);

		PyRun_SimpleString("rate(50)");  
#endif

		world.update_state(time_step);
	}

	//system ("pause");
	return 0;
}

