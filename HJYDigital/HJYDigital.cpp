/*
 * File:   HJYDigital.cpp
 * Author: hjyssg
 *
 * Created on March 4, 2014, 11:12 AM
 */

#include "HJYDigital.h"
#include <assert.h>
#include <math.h>
#include <algorithm>    // std::sort


double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}


Wire::Wire()
{
    this->id = rand();
    this->driver = NULL;
    this->receiver = NULL;
    
    this->value = UNKNOWN;
}


Wire::Wire(int ID, Module * DRIVER, Module * RECEIVER)
{
    this->id = ID;
    this->driver = DRIVER;
    this->receiver = RECEIVER;
    
    this->value = UNKNOWN;
}


Module::Module()
{
	this->id = rand();
    this->functionType = OR;
    this->delay = 1;
}

Module::Module(int ID, ModuleFunctionType FUNCTION_TYPE)
{
    this->id = ID;
    this->functionType = FUNCTION_TYPE;
    
    this->delay = 1;
}

Event * Module::get_new_event_according_to_new_value(float newValue, int current_time)
{
    if (newValue != output_wire->value)
    {
        Event * event = new Event(current_time + delay, output_wire, newValue);
        return event;
    }else
    {
        return NULL;
    }
}


Event * Module::OR_function(int current_time)
{
    
    int new_value = 0;
    
    for (int ii = 0; ii < input_wires.size(); ii++)
    {
        Wire *w = input_wires[ii];
        
        //unknown input will cause unknown result
        if (w->value == UNKNOWN){
            new_value = UNKNOWN;
            break;
        }
        
        if (w->value == ONE){
            new_value = ONE;
        }
    }
    
    return get_new_event_according_to_new_value(new_value, current_time);
}

Event * Module::NAND_function(int current_time)
{
    int new_value = ZERO;
    int one_count = 0;
    int nan_count = 0;
    
    for (int ii = 0; ii < input_wires.size(); ii++)
    {
        Wire *w = input_wires[ii];
        
        if (w->value == ONE) {one_count++; }
        else if (w->value == UNKNOWN) {nan_count++;}
    }
    
    
    if(nan_count == input_wires.size()) {
        new_value = UNKNOWN;
    }else  if (one_count == input_wires.size()) {
        new_value = ZERO;
    }else  {
        new_value = ONE;
    }
    
    //printf("%i %f %f", current_time, new_value, output_wire->value);
    return get_new_event_according_to_new_value(new_value, current_time);

    

}

Event * Module::NOR_function(int current_time)
{
    //printf("NOR_function");
    
    int new_value = output_wire->value;
    int zero_count = 0;
    int nan_count = 0;
    
    for (int ii = 0; ii < input_wires.size(); ii++)
    {
        Wire *w = input_wires[ii];
        if (w->value == ZERO){ zero_count++; }
        if (w->value == UNKNOWN) {nan_count++;}
    }
    
    if (zero_count == input_wires.size())
    {
        new_value = ONE;
    }else if(nan_count == input_wires.size())
    {
        new_value = UNKNOWN;
    }else
    {
        new_value = ZERO;
    }
    
    return get_new_event_according_to_new_value(new_value, current_time);
}

Event * Module::AND_function(int current_time)
{
    double new_value = ONE;
    
    for (int ii = 0; ii < input_wires.size(); ii++)
    {
        Wire *w = input_wires[ii];
        
        //unknown input will cause unknown result
        if (w->value == UNKNOWN){
            new_value = UNKNOWN;
            break;
        }
        
        if (w->value != ONE){
            new_value = ZERO;
        }
    }
    
    //printf("%i %f %f", current_time, new_value, output_wire->value);
    return get_new_event_according_to_new_value(new_value, current_time);
}




Event * Module::NOT_function(int current_time)
{
    assert(input_wires.size()==1);
    int new_value;
    
    //unknown input will cause unknown result
    if (input_wires[0]->value == UNKNOWN){
        new_value = UNKNOWN;
    }
    else if (input_wires[0]->value == ZERO){
        new_value = ONE;
    }
    else{
        new_value = ZERO;
    }
    //printf("%i %f %f", current_time, new_value, output_wire->value);
    return get_new_event_according_to_new_value(new_value, current_time);
}

Event*  Module::function(int current_time)
{
    if (functionType == AND)
    {
        return AND_function(current_time);
    }
    else if (functionType == OR)
    {
        return OR_function(current_time);
        
    }else if (functionType == NOT)
    {
        return NOT_function(current_time);
    }else if(functionType == NAND)
    {
        return NAND_function(current_time);
    }
    else if(functionType == NOR)
    {
        
        return NOR_function(current_time);
    }else
    {
        return NULL;
    }
}

Event::Event()
{
	this->time = 0;
    this->wire = NULL;
    this->value = NAN;
}

Event::Event(int TIME, Wire * WIRE, int VALUE)
{
    this->time = TIME;
    this->wire = WIRE;
    this->value = VALUE;
}




DigitalCircuit::DigitalCircuit()
{
    this->current_time = 0;
    this->current_event_index = 0;
}


bool compare_event(Event * e1, Event *e2){return e1->time < e2->time;}
bool compare_wire(Wire *w1, Wire *w2){	return w1->id < w2->id;}

void DigitalCircuit::init()
{
	std::stable_sort (eventList.begin(), eventList.end(), compare_event);
	std::stable_sort (wireList.begin(), wireList.end(), compare_wire);
}

bool DigitalCircuit::hasNextEvent()
{
    if (this->current_event_index < (int)this->eventList.size())
    {
        return true;
    }else
    {
        return false;
    }
}

int DigitalCircuit::get_current_time()
{
#if 0
	if (current_event_index-1 <= 0)
	{
		return 0;
	}
    else if (current_event_index-1 < (int)eventList.size())
    {
        Event *current_event = eventList[current_event_index-1];
        return current_event->time;
    }else
    {
        return 100000000;
    }
    
#else
    
    return current_time;
#endif
    
}


Module * DigitalCircuit::getModuleWithID(int ID)
{
    for (int ii = 0; ii < moduleList.size(); ii++)
    {
        Module *m = moduleList[ii];
        if (m->id == ID) {
            return m;
        }
    }
    return NULL;
}
Wire * DigitalCircuit::getWireWithID(int ID)
{
    for (int ii = 0; ii < wireList.size(); ii++)
    {
        Wire *w = wireList[ii];
        if (w->id == ID) {
            return w;
        }
    }
    return NULL;
}

void DigitalCircuit::jump_to_next_event()
{
    bool jumped = true;
    
    //avoid out of index
    if (current_event_index < eventList.size())
    {
        Event *current_event = eventList[current_event_index];
        
        assert(current_event!= NULL);
        
        current_time = current_event->time;
        
        Wire *w = current_event->wire;
        
        //if the event does not make change, skip it
        if (w->value != current_event->value)
        {
            //printf("%f %f",w->value, current_event->value );
            
            w->value = current_event->value;
            
            //new trigget change on the receiver
            Module * m = w->receiver;
            
            if (m != NULL)
            {
                //module compute and generate new event
                Event *newEvent =  m->function(current_event->time);
                
                if (newEvent != NULL) {  add_event(newEvent);  }
            }
        }
        else
        {
            jumped = false;
        }
    }
    
    current_event_index++;
    
    if (!jumped) {
        //jump_to_next_event();
    }
}


void DigitalCircuit::add_event(Event * event)
{
#if  0
	//when the list is empty
    if (eventList.size()==0) {
        eventList.push_back(event);
        return;
    }
    
    //happen after the last event
    if ( eventList[eventList.size()-1]->time <= event->time )
    {
        eventList.push_back(event);
        return;
    }
    
    
    for (int ii = 0; ii < eventList.size(); ii++)
    {
        Event *tempE = eventList[ii];
        
        if (tempE->time <= event->time)
        {
            eventList.insert(eventList.begin()+ii+1, event);
            return;
        }
    }
    
#else
    
    for (int ii = 0; ii < eventList.size(); ii++)
    {
        Event *tempE = eventList[ii];
		if (tempE->time == event->time && tempE->wire->id == event->wire->id )
        {
            //identical event
            if(tempE->value == event->value)
            {
                return;
            }else
            {
                //try to change one to different values at the same time
                event->value = UNKNOWN;
            }
        }
    }
    
	eventList.push_back(event);
	std::stable_sort (eventList.begin(), eventList.end(), compare_event);
#endif
}


void DigitalCircuit::log_event()
{
	printf("events:\n");
    for (int ii = 0; ii < eventList.size(); ii++)
    {
        Event *E = eventList[ii];
		if (E->wire != NULL)
		{
			printf("	 time:%i wire:%i value:%i\n", E->time, E->wire->id,  E->value);
		}
    }
}


void DigitalCircuit::log_wire()
{
    printf("wires:\n");
    for (int ii = 0; ii < wireList.size(); ii++)
    {
        Wire *w = wireList[ii];
        
        if (w->value == 1.0 || w->value == 0.0){
            printf("    wire%i %i ",w->id, w->value);
        }else{
            printf("    wire%i X ",w->id);
        }        
        printf("\n");
    }
}


void DigitalCircuit::log_module()
{
    printf("modules:\n");
    for (int ii = 0; ii < moduleList.size(); ii++)
    {
        Module *m = moduleList[ii];

        if (m->functionType == AND){
            printf("    %i AND ",m->id);
        }
        else if (m->functionType == OR){
            printf("    %i OR ",m->id);
        }else if (m->functionType == NOT){
            printf("    %i NOT ",m->id);    
        }else if(m->functionType == NAND){
            printf("   %i NAND ",m->id);
        }
        else if(m->functionType == NOR){
            printf("   %i NOR ",m->id);  
        }else{
            printf("   %i UNKNOWN",m->id);
        }
        
        printf("\n");
    }
    
}