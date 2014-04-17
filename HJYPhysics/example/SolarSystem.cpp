//
//  main.cpp
//  solarSystem
//
//  Created by Junyang Huang on 1/21/14.
//  Copyright (c) 2014 www.hjyssg.com. All rights reserved.
//

#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
 
#include "HJYPhysics.h"


#define VPYTHON_FLAG 1   //a flag to disable/enable vpython, so this program can run on machine without vpython
#define PRINT_FLAG 0     //a flag to disable/enable print 

#if VPYTHON_FLAG 
#include "\Python32\include\Python.h" 
char cmd[300];  //a resuable string for vpython  
#endif 


using namespace std;





TheWorld world;
PointMass * sun;
PointMass * mercury;
PointMass * venus;
PointMass * earth;
PointMass * mars;
PointMass * jupiter;
PointMass * saturn;
PointMass * uranus;
PointMass * neptune;
PointMass * moon;

vector<std::string> name_list;


double total_kinetic_energy = 0;

//time step of simulation
double time_step = 300;  

void init_physics()
{
	name_list.push_back("sun");
	name_list.push_back("mercury");
	name_list.push_back("venus");
	name_list.push_back("earth");
	name_list.push_back("mars");
	name_list.push_back("jupiter");
	name_list.push_back("saturn");
	name_list.push_back("uranus");
	name_list.push_back("neptune");

	world = TheWorld();
	world.set_gravitation_flag(true);
	world.set_intergration_type(velocity_Verlet_integration);

	//init based on http://ssd.jpl.nasa.gov/horizons.cgi 
	//from date 2013-01-01
	//solar system barcenter
	sun = new PointMass( -1.921258692976546E+05 KM, -3.672858221449571E+05 KM, -6.293939186492928E+03 KM,
		1.064112061805785E-02 KM, -2.282152601144673E-03 KM, -2.261768099104889E-04 KM,
		0,0,0,
		SUN_MASS);

	mercury = new PointMass(-2.542271688725740E+07 KM,-6.518047354640916E+07 KM, -2.986942437682569E+06 KM,
		3.561299513603979E+01 KM, -1.527095279264177E+01 KM, -4.514379444956798E+00 KM,
		0,0,0,
		MERCURY_MASS);
	venus = new PointMass(-6.905541781609924E+07 KM, -8.392401643421896E+07 KM,  2.823108036887005E+06 KM,
		2.679084974988556E+01 KM, -2.244015730475055E+01 KM, -1.853249972810995E+00 KM,
		0,0,0,
		VENUS_MASS );


	earth = new PointMass(-2.712311750453984E+07 KM, 1.442452249801158E+08 KM, -1.024866474659159E+04 KM,
		-2.975221878384466E+01 KM, -5.557630236369244E+00 KM,  1.166327378957226E-05 KM,
		0,0,0,
		EARTH_MASS);

	mars = new PointMass( 1.613831829365655E+08 KM, -1.299959988813151E+08 KM, -6.689343082176655E+06 KM,
		1.609480750581124E+01 KM,  2.096710871054539E+01  KM, 4.420934607818046E-02 KM,
		0,0,0,
		MARS_MASS);

	jupiter = new PointMass(2.130248609275347E+08 KM,  7.264621031147093E+08 KM, -7.796178050164084E+06 KM, 
		-1.269579412745502E+01 KM,  4.300631432663832E+00 KM,  2.662458765460191E-01 KM,
		0,0,0,
		JUPITER_MASS);

	saturn = new PointMass(-1.209082909091936E+09 KM, -8.253529571281528E+08 KM,  6.246931966183711E+07 KM, 
		4.923237564225631E+00 KM, -8.001281989714242E+00 KM, -5.660410531867004E-02 KM,
		0,0,0,
		SATURN_MASS);

	uranus = new PointMass(2.975359717954440E+09  KM,  3.846013966774312E+08  KM, -3.711880253105340E+07  KM, 
		-9.228235466197581E-01  KM,  6.436376882337913E+00  KM,  3.602154409843563E-02  KM,
		0,0,0,
		URANUS_MASS);

	neptune = new PointMass(3.973252374554337E+09 KM, -2.083244177306763E+09 KM, -4.866720990954759E+07 KM,  
		2.487210727347481E+00 KM,  4.845277256620428E+00 KM, -1.572810424785901E-01 KM,
		0,0,0,
		NEPTUNE_MASS);



	world.add_point_mass(sun);
	world.add_point_mass(mercury);
	world.add_point_mass(venus);
	world.add_point_mass(earth);
	world.add_point_mass(mars);
	world.add_point_mass(jupiter);
	world.add_point_mass(saturn);
	world.add_point_mass(uranus);
	world.add_point_mass(neptune);


	world.init();
	total_kinetic_energy = world.get_total_kinetic_energy();
}



void init_vpython()
{
#if VPYTHON_FLAG

	Py_Initialize();
	PyRun_SimpleString("from visual import *\n");
	PyRun_SimpleString("from visual.graph import *\n");

	sprintf(cmd,"scene1=display(autocenter = False, x =400,  width=900, height=650, scale = (7e-13, 7e-13, 7e-13)) \n" );PyRun_SimpleString(cmd);

	sprintf(cmd,"xaxis = curve(pos=[(0,0,0), (1e13,0,0)], color=color.red)\n");	PyRun_SimpleString(cmd);
	sprintf(cmd,"Yaxis = curve(pos=[(0,0,0), (0,1e13,0)], color=color.green)\n");	PyRun_SimpleString(cmd);
	sprintf(cmd,"Zaxis = curve(pos=[(0,0,0), (0,0,1e13)], color=color.blue)\n");	PyRun_SimpleString(cmd);

	//create planets
	//used the non-real visual size, otherwise some planet will be too small to see
	sprintf(cmd,"sun=sphere(color=color.red,make_trail=True,radius=%e,pos=(%e,%e,%e))\n", SUN_DIAMETER*10 ,sun->pos.x,sun->pos.y,sun->pos.z);PyRun_SimpleString(cmd);
	sprintf(cmd,"mercury=sphere(color=color.orange,make_trail=True,radius=%e,pos=(%e,%e,%e))\n", EARTH_DIAMETER*450,mercury->pos.x, mercury->pos.y,mercury->pos.z);PyRun_SimpleString(cmd);
	sprintf(cmd,"venus=sphere(color=color.blue,make_trail=True,radius=%e,pos=(%e,%e,%e))\n", EARTH_DIAMETER*450 ,venus->pos.x, venus->pos.y,venus->pos.z);	PyRun_SimpleString(cmd);
	sprintf(cmd,"earth=sphere(material = materials.earth, make_trail=True,radius=%e,pos=(%e,%e,%e))\n", EARTH_DIAMETER *500 ,earth->pos.x,earth->pos.y,earth->pos.z);PyRun_SimpleString(cmd);
	sprintf(cmd,"mars=sphere(color=color.green,make_trail=True,radius=%e,pos=(%e,%e,%e))\n", EARTH_DIAMETER*450,mars->pos.x, mars->pos.y,mars->pos.z);PyRun_SimpleString(cmd);
	sprintf(cmd,"jupiter=sphere(color=color.cyan,make_trail=True,radius=%e,pos=(%e,%e,%e))\n", EARTH_DIAMETER*450 ,jupiter->pos.x, jupiter->pos.y,jupiter->pos.z);PyRun_SimpleString(cmd);
	sprintf(cmd,"saturn=sphere(color=color.blue,make_trail=True,radius=%e,pos=(%e,%e,%e))\n", EARTH_DIAMETER*450,saturn->pos.x, saturn->pos.y,saturn->pos.z);PyRun_SimpleString(cmd);
	sprintf(cmd,"uranus=sphere(color=color.green,make_trail=True,radius=%e,pos=(%e,%e,%e))\n",EARTH_DIAMETER*450,uranus->pos.x, uranus->pos.y,uranus->pos.z);PyRun_SimpleString(cmd);
	sprintf(cmd,"neptune=sphere(color=color.white,make_trail=True,radius=%e,pos=(%e,%e,%e))\n", EARTH_DIAMETER*450 ,neptune->pos.x, neptune->pos.y,mars->pos.z);PyRun_SimpleString(cmd);

	//create label
	sprintf(cmd,"time_label = label(yoffset=280,  line=0, text=\"0 S\")" );	PyRun_SimpleString(cmd);

	//label for planets
	for (int ii = 0; ii< name_list.size();ii++)
	{
		const char *name = name_list[ii].c_str();
		sprintf(cmd,"%s_label = label(pos=%s.pos, text = \"%s\", yoffset=5, opacity=1, line=0, box = 0)\n",name,name,name );						PyRun_SimpleString(cmd);
		sprintf(cmd,"%s_varr = arrow(pos=%s.pos, axis=(0,0,0), color=color.yellow)\n",name,name);	PyRun_SimpleString(cmd);
	}

	//create graph
	sprintf(cmd, "energy_graph=gdisplay( title = 'energy', xtitle='day',x=0, y=660, width=900, height=300, ytitle='energy')\n");PyRun_SimpleString(cmd);
	sprintf(cmd, "simulated_kinetic_line = gcurve(color=color.green)\n");PyRun_SimpleString(cmd);
	sprintf(cmd, "real_kinetic_line = gcurve(color=color.blue)\n");PyRun_SimpleString(cmd);

#endif

}

void update_vpython()
{
#if VPYTHON_FLAG

	sprintf(cmd,"time_label.text='%i Years %i Days'\n",world.timer.get_elapsed_years(), world.timer.get_elapsed_days()%365 ); 
	PyRun_SimpleString(cmd);

	//synchronize ui with real data
	for (int ii = 0; ii< name_list.size();ii++)
	{
		const char *name = name_list[ii].c_str();
		PointMass *pm = world.point_mass_arr[ii];

		const double m = 1.5e5;
		Vertex temp =  pm->spd.normalize().scale(m,m,m);  //normilize and scale veloticy
		sprintf(cmd,"%s.pos=(%e,%e,%e)\n", name, pm->pos.x,pm->pos.y,pm->pos.z); PyRun_SimpleString(cmd);
		sprintf(cmd,"%s_label.pos=(%e,%e,%e)\n", name, pm->pos.x,pm->pos.y,pm->pos.z); PyRun_SimpleString(cmd);
		sprintf(cmd,"%s_varr.pos= (%e,%e,%e)\n",name,pm->pos.x,pm->pos.y,pm->pos.z);PyRun_SimpleString(cmd);
		sprintf(cmd,"%s_varr.axis= (%e,%e,%e)\n",name, temp.x*m,temp.y*m,temp.z*m);PyRun_SimpleString(cmd);
	}


	sprintf(cmd, "simulated_kinetic_line.plot(pos=(%i, %e))\n",
		world.timer.get_elapsed_days(), world.get_total_kinetic_energy()/1e+36 );   PyRun_SimpleString(cmd);
	sprintf(cmd, "real_kinetic_line.plot(pos=(%i, %e))\n",
		world.timer.get_elapsed_days(), total_kinetic_energy/1e+36);   PyRun_SimpleString(cmd);

	sprintf(cmd,"rate(100)");
	PyRun_SimpleString(cmd);

#endif
}



int main(int argc, const char * argv[])
{
	init_physics();
	init_vpython();

	int num_day = 0;

	while(true)
	{
#if PRINT_FLAG
		if (num_day < world.timer.get_elapsed_days())
		{
			CLOG("at day"<<elapse_day<<"s\nsun=" <<sun->toString()<<"\nearth="<<earth->toString());
		}
#endif

		//the ui update once a day
		if (num_day < world.timer.get_elapsed_days())
		{
			update_vpython();
			num_day = world.timer.get_elapsed_days();
		}

		world.update_state(time_step);
	}

	system ("pause");
	return 0;
}

