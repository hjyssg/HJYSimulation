/**
* @mainpage HJYDigital
* A simple digital circuit simulation engine
* @author Junyang Huang 黄竣杨
**/



/**
 * @file HJYDigital.h
 *
 * @brief A simple digital circuit simulation API
 *
 * @author Junyang Huang
 *
 * @date 2014/03/05
 *
 **/


#ifndef HJYDIGITAL_H
#define	HJYDIGITAL_H


#include <vector>


#define  CLOG(a)  std::cout<<a<<"\n"
#define  DLOG(a)  std::cout<<__PRETTY_FUNCTION__<<": "<<a<<"\n"

/**
* @brief genrate random float between fMin and fMax
*from http://stackoverflow.com/questions/2704521/generate-random-double-numbers-in-c
*Remember to call srand() with a proper seed each time your program starts.
*/
double fRand(double fMin, double fMax);



class Module;
class Wire;
class Event;


/**
 * @brief the wire used to connect two digital modules
 */
class Wire
{
public:
     /**
     * @brief the ID of wire must be unique
     */
    int id;
    
    /**
    * @brief the value of the value. default is nan
    **/
    int value;

    /**
    * @brief the module connected from 
    **/
    Module *driver;

    /**
    * @brief the module connected to
    **/
    Module *receiver;
    
public:
    Wire();
    Wire(int ID, Module * DRIVER, Module * RECEIVER);
    
};

/**
* @brief supported module functions
*/
enum ModuleFunctionType{
    NOT = -1,
    AND = 1,
    OR = 2,
    NAND = 5,
    NOR  = 6,
    OUTPUT = 3,
    INPUT = 4
};

/**
 * @brief Digital Module, e.g AND GATE, OR GATE
 */
class Module
{
public:
    /**
     * @brief the ID of module must be unique
     */
    int id;
    
    /**
     * @brief delay of signal propagation default is 1, change it if you need
     */
    int delay;
    
    /**
     * @brief function type of module defaut is OR
     */
    ModuleFunctionType functionType;
    
    /**
     * @brief the input wires
     */
    std::vector<Wire *> input_wires;
    
    /**
     * @brief the output wire
     */
    Wire * output_wire;
    
    
public:
	Module();
    Module(int ID, ModuleFunctionType FUNCTION_TYPE);

    /**
     * @brief given the current time, compute new event according to the inputs and function type
     * @param current_time
     * @return the new event If the output will not change, the event will be NULL
     */
    Event* function(int current_time);
    
private:
   
    /**
     * @brief OR function, used by  "Event* function(int current_time)"
     */
    Event * OR_function(int current_time);
    
    /**
     * @brief AND function, used by  "Event* function(int current_time)"
     */
    Event * AND_function(int current_time);
    
    /**
     * @brief NOT function, used by  "Event* function(int current_time)"
     */
    Event * NOT_function(int current_time);
    
    /**
     * @brief NAND function, used by  "Event* function(int current_time)"
     */
    Event * NAND_function(int current_time);
    
    /**
     * @brief NOR function, used by  "Event* function(int current_time)"
     */
    Event * NOR_function(int current_time);

    Event * get_new_event_according_to_new_value(float newValue, int current_time);
    
};


#define ONE 1
#define ZERO 0
#define UNKNOWN -1


/**
 * @brief The event is to tell a wire change its value to new value at certain time
 */
class Event
{
public:
    /**
     * @brief when the event happen. >= 0
     */
    int time;
    
    /**
     * @brief where the event happen
     */
    Wire *wire;
    
    /**
     * @brief how the value of the wire change to. 0.0<= value <= 1.0
     */
    int value;

    
public:
	Event();
    Event(int TIME, Wire * WIRE, int VALUE);


};

/**
 * @brief digital circuit contatins moduls, wires and events
 */
class DigitalCircuit
{

public:
	
	
    /**
     * @brief index points to the event in the eventList that is going to happen next
     */
    int current_event_index;
    
    
    int current_time;
    
    /**
     * @brief all wires in the circuit
     */
    std::vector<Wire *> wireList;
    
    /**
     * @brief all modules in the circuit
     */
    std::vector<Module *> moduleList;
    
    /**
     * @brief all events happened or going to happen in the circuit
     */
    std::vector<Event *> eventList;
    
public:
    DigitalCircuit();

    /**
     @brief after adding wires, modules and input events. init the circuit before simulation
     */
	void init();

    /**
     @brief if has next event
     */
    bool hasNextEvent();
    
    /**
     @brief get current time
     */
    int get_current_time();
    
    /**
     @brief get module from moduleList with id
     */
    Module * getModuleWithID(int ID);
    
    /**
     @brief get wire from wireList with id
     */
    Wire * getWireWithID(int ID);
    
    /**
     @brief jump to next event
     */
    void jump_to_next_event();
    
    /**
     @brief 
     */
    void add_event(Event * event);

    /**
     @brief print all events of eventlist for debugging
     */
	void log_event();
    
    /**
     @brief print all wires of wirelist for debugging
     */
    void log_wire();

    /**
     @brief print all modules of modulelist for debugging
     */
    void log_module();

    
    
};

#endif	/* HJYDIGITAL_H */

