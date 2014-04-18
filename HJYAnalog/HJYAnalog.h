
/**
 * @mainpage HJYAnalog.h
 * A simple digital circuit simulation engine
 * @author hjyssg
 **/


#ifndef __HJYAnalog__HJYAnalog__
#define __HJYAnalog__HJYAnalog__

#include <iostream>
#include <vector>

#include <string>
#include "Matrix.h"

#include <stdlib.h>



#pragma mark helper-functions and macro
//print function
#define  CLOG(a)  std::cout<<a<<"\n"
#define  DLOG(a)  std::cout<<__PRETTY_FUNCTION__<<": "<<a<<"\n"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


#pragma mark syntax candy for unit
#define S
#define MS *1e-3
#define NS *1e-9
#define PS *1e-12
#define FS *1e-15
#define J
#define OHM
#define H
#define F
#define V
#define A
#define UH *1e-6

#include <stdlib.h>
#include <stdbool.h>

typedef double **  matrix_t;

///init a empty matrix
matrix_t matrix_init(int row_num, int col_num);

///apply gauss elimination on a matrix
bool matrix_gauss(matrix_t M, int row_num, int col_num);

///return the transpose of matrix
matrix_t matrix_transpose(matrix_t m, int row_num, int col_num);

///print a matrix for debugging
void matrix_print(matrix_t M, int row_num, int col_num);

///free matrix memory
void matrix_free(matrix_t M, int row_num);

///simple unit testing for methods above
void matrix_testing();



class Node;
class Timer;

///supported components
enum ComponentType{
    None,
    Voltage_Source,
    Current_Source,
    Resistor
};

enum TerminalType {
    PositiveTerminal = 1,
    NegativeTerminal = -1
};


///root component for specific component
class Component
{
public:
    //a arbitary id
    //must be unique among all component
    //will be generated randomly by constructor
    //or set by other classes
    int id;
    
    ///in F
    double capacitance;
    
     // in ohm
    double resistance;
    
    ///in  henries
    double inductance;
    
    /// in A
    double current;
    
    /// in V
    double voltage;
    
    //flag used by graph traversal algorithm
    bool PT_visited;
    
    //flag used by graph traversal algorithm
    bool NT_visited;
    
    ///the negative terminal
    Node * inputNode;
    
    ///the positive terminal
    Node * outputNode;
    
public:
    ///constructor
    Component();
    
    ///get current
    virtual double get_current();
    
    ///get component type
    virtual ComponentType get_type();
    
    ///to string method
    virtual std::string toString();

    virtual double get_conductance(double time_step);
};

///voltage source
class voltage_source:public Component
{
    ///if on, used lambda function passed by the user to generate voltage
    ///default is false
    bool use_voltage_function;
    
private:
    std::function<double(double)> voltage_function;
    
public:
    voltage_source();
    
    ///construct with the voltage that will not change during the simulation
    voltage_source(double const_voltage);
    
    void set_voltage_function(std::function<double(double)> voltage_function_p);
    
    double get_voltage(double time);
    
    ComponentType get_type();
    

    ///to string method
    std::string toString();
};

///current source
class current_source:public Component
{
public:
    
    ///if on, used lambda function passed by the user to generate voltage
    ///default is false
    bool use_current_function;
    
private:
    std::function<double(double)> current_function;

public:
    current_source();
    
    ///construct with the current that will not change during the simulation
    current_source(double const_current);
    
    void set_current_function(std::function<double(double)> current_function_p);
    
    double get_current(double time);
    
    ComponentType get_type();

    ///to string method
    std::string toString();
};


///resisotr
class resistor:public Component
{
public:

public:
    resistor();
    
    ///construct with the resistance that will not change during the simulation
    resistor(double const_resis);
    
    ComponentType get_type();

    double get_conductance(double time_step);
    
    ///to string method
    std::string toString();
};


/// node that connects component toghter
class Node
{
public:
    
    ///in V
    double voltage;
    
    ///should numbered from 0, 1, 2, 3...
    int id;

    Node(int ID);
    
    ///to string method
    std::string toString();
};


/// a timer class used by "TheAnalogCircuit" class to keep track of time
class Timer
{
public:
	double elapsed_second;
    
	Timer();
    
	//add time_step to elapsed_second
	void update_state(double time_step);
    
	int get_elapsed_femtosecond();
    
	int get_elapsed_nanosecond();
    
	int get_elapsed_picosecond();
    
	int get_elapsed_second();
    
	int get_elapsed_minute();
};


///Analog Circuit
class AnalogCircuit
{
public:
    ///track the time
    Timer time;
    
    ///the list of all components
    std::vector<Component *> component_list;
    
    ///the list of all nodes
    std::vector<Node *> node_list;
    
    
    ///you must set one node to grounded. but do not add it
    Node *ground_node;
    
    ///a flag to decide if to print the matrix during matrix
    bool print_matrix_during_simulation;
    
    ///for internal algorithm used
    ///only independent voltage source
    std::vector<voltage_source *> voltage_sources;
    
    ///for internal algorithm used
    ///only independent current source
    std::vector<current_source *> current_sources;
    
private:
     ///for internal algorithm used
    //a matrix created by appending A+z
    //size :
    matrix_t old_matrix_Q;
    
    int row_num_q;
    int col_num_q;
    
    
    int row_num_a;
    int col_num_a;
    
    int row_num_z;
    int col_num_z;
    
    ///number of nodes
    int NN;
    
    ///node of indepent voltage sources
    int MM;
    
public:
    
    AnalogCircuit();
    
    void add_component(Component * component);
    
    ///check if the circuit is correct
    bool check_if_circuit_correct();
    
    ///after adding all component, call this function
    void init();
    
 
    void update_state(double time_step);
    
    ///print component condition
    void log_component();
    
    ///print node for debugging
    void log_nodes();
    
    ~AnalogCircuit();
    
private:
    
    ///get matrix_A for MNA
    ///return values by pointer
    /// referto http://qucs.sourceforge.net/tech/node14.html
   matrix_t MNA_get_A_Matric( double time_step);
    
    ///get matrix_Z for MNA
    ///return values by pointer
    matrix_t MNA_get_Z_Matric(double time_step);
    
    
    matrix_t MNA_get_Q_matrix(matrix_t MATRIX_A, matrix_t MATRIX_Z);
    
    void apply_matrix_B_and_C_to_Matrix_A(matrix_t MATRIX_A);
    
    void add_current_to_matrix_Z(matrix_t MATRIX_Z, Component * comp,double time_step);
    
    
    ///after solving the linear equation
    ///result to components
    void apply_matrix_Q_to_component(matrix_t new_matrix_q , double time_step);
    
    
    ///find voltage_source index in voltage_sources
    ///if return find, return index
    ///if not find, return a negative number
    int get_voltage_num(voltage_source *vv);
    
    
    ///find voltage_source index in voltage_sources
    ///if return find, return index
    ///if not find, return a negative number
    int get_current_num(current_source *cc);
};



#endif /* defined(__HJYAnalog__HJYAnalog__) */
