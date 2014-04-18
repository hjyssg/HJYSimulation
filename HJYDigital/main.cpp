/**
 * @file main.cpp
 *
 * @brief a example use HJYDigital
 *
 * @author Junyang Huang
 *
 * @date 2014/03/05
 *
 **/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>      
#include <time.h>  
#include <stdlib.h>  

//tech note: to use boost lib, add "path_to_boost_1_55_0" in VC++ directories -> include directory
//mine is "C:\boost_1_55_0"
// On Mac, add  "boost" folder path to   "Header Search Paths"
//e.g /Users/workingAccount/Dropbox/boost_1_55_0
//
#include<boost/tokenizer.hpp>
#include<boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "HJYDigital.h"

using namespace std;
using namespace boost;


DigitalCircuit *circuit;

//a flag to disable/enable print
bool PRINT_FLAG = true;     



 /**
 * @brief parse module string, generate module and add to circuit
 */
void parse_module_string(string line)
{
	char_separator<char> sep(", {}");
	tokenizer<char_separator<char>> tokens(line, sep);
	int ii = 0;
	Module *m = new Module();
	for ( tokenizer<char_separator<char>>::iterator it = tokens.begin(); it != tokens.end(); ++it)
	{   
        string token(*it);
        trim(token);
		if(token.size() == 0) {continue;}

		if (ii == 0) {
			//the first token is 'C' or 'c'
		}
		else 
		if (ii==1)
		{
			//specify type
			string temp(token);
			if (iequals(temp, "AND")){
				m->functionType = AND;
			}else if (iequals(temp, "OR")){
				m->functionType = OR;
			}else if (iequals(temp, "NOT")){
				m->functionType = NOT;
			}else if (iequals(temp, "NOR")){
				m->functionType = NOR;
			}
            		else if (iequals(temp, "NAND")){
				m->functionType = NAND;
			}
            		else{
				CLOG("SUPPORTED COMPONENT: ["<<(token)<<"]"<<line);
				return;
			}	
		}
		else if (ii == 2) {
			//specify id
			try{
				m->id = lexical_cast<int>(token);
			}
			catch(const bad_lexical_cast &){
				CLOG("WRONG MODULE ID: ["<<(token)<<"]"<<line);
				return;
			}
		}
		else if (ii == 3){
			//specify output
			try{
				int output_wire_id = lexical_cast<int>(token);
				Wire * w = circuit->getWireWithID(output_wire_id);

				//if wire exist, connect to
				//if not create new wire
				if (w == NULL){
					w = new Wire(output_wire_id, m, NULL);
					circuit->wireList.push_back(w);
				}
                
        			 //connect module and wire togher
        			 w->driver = m;
				m->output_wire = w;
			}
			catch(const bad_lexical_cast &){
				CLOG("WRONG MODULE OUTPUT WIRE ID: ["<<(token)<<"]"<<line);
				return;
			}
		}else 
		{
			//specify input
			try{
				int input_wire_id = lexical_cast<int>(token);
				Wire * w = circuit->getWireWithID(input_wire_id);

				//if wire exist, connect to
				//if not create new wire
				if (w == NULL){
					w = new Wire(input_wire_id, NULL,m);
					circuit->wireList.push_back(w);
				}

                 		//connect module and wire togher
        			 w->receiver = m;
				m->input_wires.push_back(w);
			}
			catch(const bad_lexical_cast &){
				CLOG("WRONG MODULE INPUT WIRE ID: ["<<(token)<<"]"<<line);
				return;
			}
		}

		ii++;
	}

	//add module to the circuit
	circuit->moduleList.push_back(m);
	return;
}


 /**
 * @brief parse event string, generate event and add to circuit
 */
void parse_event_string(string line)
{
	char_separator<char> sep(", ");   
	tokenizer<char_separator<char>> tokens(line, sep);
	int ii = 0;
	Event * event = new Event();
	for ( tokenizer<char_separator<char>>::iterator it = tokens.begin(); it != tokens.end(); ++it)
	{

		string token(*it);
		trim(token);
		if(token.size() == 0) {continue;}


		if (ii == 0) {
			//the first token is 'E' or 'e'
		}
		else if (ii==1)
		{
			try{
				int wire_id = lexical_cast<int>(token);
				Wire *w = circuit->getWireWithID(wire_id);

				if (w == NULL){
					w = new Wire(wire_id, NULL,NULL);
					circuit->wireList.push_back(w);
				}
				event->wire = w;
			}
			catch(const bad_lexical_cast &)
			{
				CLOG("WRONG EVENT ID: ["<<(token)<<"]"<<line);
				return;
			}	
		}
		else if (ii == 2) 
		{
			//specify time
			try{
				int tt = lexical_cast<int>(token);
				event->time = tt;
			}
			catch(const bad_lexical_cast &){
				CLOG("WRONG EVENT TIME: ["<<(token)<<"]"<<line);
				return;
			}

		}
		else if (ii == 3)
		{
			try{
				float vv = lexical_cast<float>(token);
				event->value = vv;
			}
			catch(const bad_lexical_cast &){
				CLOG("WRONG EVENT VALUE: ["<<(token)<<"]"<<line);
				return;
			}

		}else
		{
			CLOG("SUPPORTED EVENT: "<<line);
			return;
		}
		ii++;
	}


	//add event to the circuit
	circuit->add_event(event);
	return;
}

 /**
 * @brief init circuit with self-defined text file; add modules, wires, and events to it
 */
void init_circuit_with_file(std::string fileName)
{
	std::ifstream data_file(fileName);
	if (data_file)
	{	 
		std::string line;

		while (std::getline(data_file, line))
		{
			//ignore empty line
			if(line.size()==0){continue;}
			if(line[0] == '#') {continue;}

			if (line[0] == 'C'||line[0] == 'c'){
				//CLOG(line);
				parse_module_string(line);

			}else if (line[0] == 'E'||line[0] == 'e'){
				//CLOG(line);
				parse_event_string(line);
			}
		}  

		CLOG("init from "<< fileName);
	}
	else
	{
		CLOG("fail to load file name; run simple demo");
		Module *and1 = new Module(1, AND);
		circuit->moduleList.push_back(and1);

		Wire *input1 = new Wire(1, NULL, and1);
		and1->input_wires.push_back(input1);
		circuit->wireList.push_back(input1);

		Wire *input2 = new Wire(2, NULL, and1);
		and1->input_wires.push_back(input2);
		circuit->wireList.push_back(input2);

		Wire *outpu1 = new Wire(3, and1,NULL);
		circuit->wireList.push_back(outpu1);
		and1->output_wire = outpu1;

		Event * e2 = new Event(2, input2, 1.0);
		circuit->add_event(e2);

		Event * e1 = new Event(1, input1, 1.0);
		circuit->add_event(e1);
	}
}


void set_delay(int dd)
{
    for (int ii = 0; ii< circuit->moduleList.size(); ii++)
    {
        Module *m = circuit->moduleList[ii];
        
        m->delay = dd;
    }
}


 /**
 * @brief main function
 */
int main(int argc, char** argv)
{
	srand ((unsigned int) time(NULL) );

	circuit = new DigitalCircuit();
	init_circuit_with_file(std::string("test3.txt"));

    	set_delay(1);
    
	circuit->init();
	circuit->log_event();
    	circuit->log_wire();
    	circuit->log_module();
	
    	printf("=========Begin The Simulation=========\n");
	while (circuit->hasNextEvent())
	{
      
          int tt = circuit->get_current_time();
    
    	  if(PRINT_FLAG){
        	printf("at time %i ns:\n", tt);
		 circuit->log_wire();
    	  }

      		circuit->jump_to_next_event();
	}
    
        update_vpython();
        circuit->log_event();
	system("pause");
	return 0;
}

