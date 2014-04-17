//
//  main.cpp
//  HJYAnalog
//
//  Created by Junyang Huang on 2/25/14.
//  Copyright (c) 2014 www.hjyssg.com. All rights reserved.
//

#include <iostream>
#include "HJYAnalog.h"
#include <math.h>
#include <algorithm>    // std::min

const float PI =  3.14159265;

bool  PRINT_FLAG = true;     //a flag to disable/enable print

//components will be used
AnalogCircuit * theCircuit = new AnalogCircuit();
Node *n1 = new Node(0);;
Node *n2 = new Node(1);
voltage_source *v1 = new voltage_source(5.0);
current_source *vs = new current_source(5.0);
resistor *r1 = new resistor();
Component *cc = NULL;


//time step of simulation
double time_step = 0.005 S;
int duration = 200;

int test = 3;




double linear_line(double x)
{
    double tt = 10*x;
    const double max = 2;
    
    return tt<max? tt:max;
}

double sin_wave(double x){
    return 5*sin(5*PI*x);
}

double step_wave(double x)
{
    if (x<0.01||x>0.02) {
        return 5;
    }else
    {
        return -5;
    }
}


void init_circuit()
{
    Node *grounded_node = new Node(1992);
    grounded_node->voltage = 0;
    theCircuit->ground_node = grounded_node;
    
    if(test==0)
    {
        //ohm's law RV circuit
        //one voltage with one resistor
        /*
         *     -------node1---------
         *     |              |     |
         *     V              R1    R2
         *     |              |     |
         *     |              |     |
         *     ----------------------
         *                    |
         *                  ----
         *
         */
        
        theCircuit->node_list.push_back(n1);
        
        //create voltage source
        theCircuit->add_component(v1);
        
        //connect v1 with n1
        v1->outputNode = n1;
        v1->inputNode = grounded_node;
        
        //use sine wave
        v1->set_voltage_function(sin_wave);
        
        
        //create resistor
        r1->resistance = 2;
        theCircuit->add_component(r1);
        
        //connect resistor with n1  and connect resistor with n2
        r1->outputNode = n1;
        r1->inputNode = grounded_node;
        
        cc = new resistor(3);
        theCircuit->add_component(cc);
        cc->outputNode = n1;
        cc->inputNode = grounded_node;
        
    }
    else if(test > 1)
    {
        /*
         *     node1-----R-----node2
         *     |                |
         *     V                other
         *     |                |
         *     |                |
         *     ------------------
         *                    |
         *                  -----
         *
         */
        
        theCircuit->node_list.push_back(n1);
        theCircuit->node_list.push_back(n2);
        
        
        
        
        if (false)
        {
            //create voltage source
            theCircuit->add_component(v1);
            //connect v1 with n1s
            v1->outputNode = n1;
            v1->inputNode = grounded_node;
            
            //use sine wave
            v1->set_voltage_function(sin_wave);
            
        }
        else
        {
            //create voltage source
            theCircuit->add_component(vs);
            //connect v1 with n1s
            vs->outputNode = n1;
            vs->inputNode = grounded_node;

            vs->set_current_function(sin_wave);
        }
        
       
        
        
        //create resistor
        r1->resistance = 500;
        theCircuit->add_component(r1);
        
        
        r1->inputNode = n2;
        
        r1->outputNode = n1;
        
        if (test == 2)
        {
            //create capacitor
            cc = new capacitor(1e-6);
            
        }else if (test == 3)
        {
            //create inductor
            cc = new inductor(1e-4);
        }
        else if (test == 4)
        {
            //create diode
            cc = new diode(0.025875, 1e-15);
            cc->voltage = 0.9;
        }
        else
        {
            cc = new resistor(2);
        }
        
        theCircuit->add_component(cc);
        
        
        cc->inputNode = grounded_node;
        
        cc->outputNode = n2;
    }

    //init will do some internal work to prepare the simulation
    theCircuit->init();
}



int main(int argc, const char * argv[])
{
    srand((unsigned int)time(NULL));
    
    
    //matrix_testing();
    
    CLOG("init component");
    init_circuit();
    theCircuit->log_component();
    theCircuit->log_nodes();
    
    
    //theCircuit->print_matrix_during_simulation = true;
    
    int count = 0;
    
    
    for (int ii = 0; ii < duration; ii++)
    {        theCircuit->update_state(time_step);
        
        
        if (PRINT_FLAG&&count == 0)
        {
            theCircuit->log_nodes();
            theCircuit->log_component();
        
            
            count = 1;
        }
        
        count--;
        
    }
    
    return 0;
}

