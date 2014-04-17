//
//  HJYAnalog.cpp
//  HJYAnalog
//
//  Created by Junyang Huang on 2/25/14.
//  Copyright (c) 2014 www.hjyssg.com. All rights reserved.
//

#include "HJYAnalog.h"
#include <math.h>
#include <limits>
#include <assert.h>
#include <algorithm>





Component::Component()
{
    static int count = 0;
    PT_visited = false;
    NT_visited = false;
    id = count;
    voltage = 0.0;
    current = 0.0;
    count ++;
}

double Component::get_current()
{
    return current;
}

ComponentType Component::get_type()
{
    return None;
}

double Component::get_conductance(double time_step)
{
    return 0.0;
}

std::string Component::toString()
{
    char s[100];
    sprintf(s, "id:%i type:unknown %.2eV %.2eA", id,voltage,current);
    std::string temp(s);
    return s;
}


#pragma mark voltage_source
voltage_source::voltage_source():Component()
{
    use_voltage_function = false;
}


voltage_source::voltage_source(double const_voltage)
{
    this->voltage = const_voltage;
}

double voltage_source::get_voltage(double time)
{
    if (!use_voltage_function) {
        return voltage;
    }else
    {
        this->voltage = voltage_function(time);
        return voltage;
    }
}

void voltage_source::set_voltage_function(std::function<double(double)> voltage_function_p)
{
    use_voltage_function = true;
    voltage_function = voltage_function_p;
}

ComponentType voltage_source::get_type()
{
    return Voltage_Source;
}

std::string voltage_source::toString()
{
    char s[100];
    sprintf(s, "id:%i type:VS %.2eV %.2eA", id,voltage,current);
    std::string temp(s);
    return s;
}


#pragma mark current_source
current_source::current_source():Component()
{
    
}


current_source::current_source(double const_current)
{
    this->current = const_current;
}

ComponentType current_source::get_type()
{
    return Current_Source;
}

void current_source::set_current_function(std::function<double(double)> current_function_p)
{
    use_current_function = true;
    current_function = current_function_p;
}


double current_source::get_current(double time)
{
    if (!use_current_function) {
        return current;
    }else
    {
        return current_function(time);
    }
}

std::string current_source::toString()
{
    char s[100];
    sprintf(s, "id:%i type:CS %.2eV %.2eA", id,voltage,current);
    std::string temp(s);
    return s;
}


#pragma mark resistor_source
resistor::resistor():Component()
{
    
}

resistor::resistor(double const_resis)
{
    this->resistance = const_resis;
}

ComponentType resistor::get_type()
{
    return Resistor;
}

double resistor::get_conductance(double time_step)
{
    return 1/resistance;
}


std::string resistor::toString()
{
    char s[100];
    sprintf(s, "id:%i type:RST %.3eV %.3eA %.3eΩ", id,voltage,current, resistance);
    std::string temp(s);
    return s;
}


#pragma mark Timer
Timer::Timer()
{
	elapsed_second = 0.0;
}

void Timer::update_state(double time_step)
{
	elapsed_second += time_step;
}

int Timer::get_elapsed_femtosecond()
{
	return (int)(elapsed_second/(1 FS));
}

int Timer::get_elapsed_nanosecond()
{
    
	return (int)(elapsed_second/(1 NS));
}

int Timer::get_elapsed_picosecond()
{
	CLOG(elapsed_second);
	return (int)(elapsed_second/(1 PS));
}


int Timer::get_elapsed_second()
{
	int current_time = (int) elapsed_second;
	return current_time/(1 S);
}



int Timer::get_elapsed_minute()
{
	int current_time = (int) elapsed_second;
	return current_time/(60 S);
}

#pragma mark node
Node::Node(int ID)
{
    id = ID;
}

std::string Node::toString()
{
    char s[100];
    sprintf(s, "Node_%i %.2eV", id,voltage);
    std::string temp(s);
    return s;
}




#pragma mark analog circuit
AnalogCircuit::AnalogCircuit()
{
    print_matrix_during_simulation = false;
    ground_node = NULL;
}

void AnalogCircuit::add_component(Component * component)
{
    this->component_list.push_back(component);
}


bool AnalogCircuit::check_if_circuit_correct()
{
    if (ground_node == NULL) {
        printf("NO GROUND NODE!!\n");return false;
    }
    
    if (component_list.size() == 0) {
        printf("NO ELECTRIC COMPONENTS!!\n");return false;
    }
    
    if(node_list.size() == 0)  {
        printf("ZERO WIRE CIRCUIT!!\n");  return false;
    }
    
    for (int ii = 0; ii < node_list.size(); ii++)
    {
        if(node_list[ii]->id != ii){
            printf("NODE ID SHOULD START FROM 0, 1, 2, 3....\n");  return false;
        }else if (node_list[ii] == ground_node){
            printf("DO NOT ADD GROUND NODE TO NODE LIST....\n");  return false;
        }
    }
    
    if (ground_node == NULL){
        printf("NO GROUND NODE ....\n");  return false;
    }
    
    if (ground_node->voltage != 0.0){
        printf("GROUND NODE VOLTAGE MUST BE ZERO....\n");  return false;
    }
    
    
    //generate nodelist
    for (int ii = 0; ii < component_list.size();ii++)
    {
        if (component_list[ii] == NULL) {
            printf("NULL ELECTRIC COMPONENTS!!\n");return false;
        }else if (component_list[ii]->inputNode == NULL)
        {
            CLOG(component_list[ii]->toString()<<" HAS NO INPUT NODE");return false;
        }
        else if (component_list[ii]->outputNode == NULL)
        {
            CLOG(component_list[ii]->toString()<<" HAS NO OUT NODE");return false;
        }
    }
    
    return true;
}

bool compare_node(Node *n1, Node *n2){ return n1->id < n2->id;}

void AnalogCircuit::init()
{
    //sort node list
    std::sort(node_list.begin(), node_list.end(), compare_node);
    
    if(check_if_circuit_correct())
    {
        for (int ii = 0; ii< this->component_list.size(); ii++)
        {
            Component *c = (Component *)component_list[ii];
            
            if (c->get_type() == Voltage_Source) {
                voltage_sources.push_back( (voltage_source *) c);
            }
            else if (c->get_type() == Current_Source) {
                current_sources.push_back( (current_source *) c);
            }
        }
        
        ///init old matrix z
        NN = (int)node_list.size();
        MM = (int)voltage_sources.size();
        row_num_q = NN+MM;
        col_num_q = NN+MM+1;
        
        row_num_a = NN+MM;
        col_num_a = NN+MM;
        
        row_num_z = NN+MM;
        col_num_z = 1;
        
        this->old_matrix_Q = matrix_init(row_num_q, col_num_q);
        
    }else
    {
        printf("FAIL TO INIT CIRCUIT");
        assert(1==2);
    }
}


#pragma mark matrix-A related algorithm

void AnalogCircuit::apply_matrix_B_and_C_to_Matrix_A(matrix_t MATRIX_A)
{
    matrix_t matrixA = MATRIX_A;
    
    //then generate matrix B and put it into A
    //B has a column for each voltage source, V
    for (int ii = 0; ii < voltage_sources.size(); ii++)
    {
        voltage_source *vv = voltage_sources[ii];
        
        /// If V has a positive terminal on node j, put 1 in row j.
        int jj = vv->outputNode->id;
        
        assert(jj>=0 && NN+ii < NN+MM);
        matrixA[jj][NN+ii] = 1;
        
        //If other end of V is on ground, stop
        if (vv->inputNode != this->ground_node)
        {
            ///Else, if V has a negative terminal on node k put -1 on row k
            int kk = vv->inputNode->id;
            matrixA[kk][NN+ii] = -1;
        }
    }
    
    
    ///transpose B into C region
    for (int ii = 0; ii < NN+MM; ii++)
    {
        for (int jj = NN; jj < NN+MM; jj++)
        {
            matrixA[jj][ii] = matrixA[ii][jj];
        }
    }
}


matrix_t  AnalogCircuit::MNA_get_A_Matric( double time_step)
{
    
    ///init empty matrix A
    matrix_t matrixA = matrix_init((MM+NN), (MM+NN));
    
    
    //put in matrix G into A
    for (int ii = 0; ii< this->component_list.size(); ii++)
    {
        Component * c =  component_list[ii];         assert(c != NULL);
        
        double g1 = c->get_conductance(time_step);
        
        if (g1 == 0.0) {
            continue;
        }
        
        Node * node_k = (Node *) c->outputNode;      assert(node_k != NULL);
        Node * node_j = (Node *) c->inputNode;       assert(node_j != NULL);
        
        //if node_k is ground, swap them
        if (node_k == this->ground_node)
        {
            std::swap(node_j, node_k);
        }
        
        int kk = node_k->id;   assert(kk>=0&&kk<NN);
        
        //add value of 1/R on diagonal (k,k)
        matrixA[kk][kk] += + g1;
        
        // If other end of R is on ground, stop
        if (node_j != this->ground_node)
        {
            int jj = node_j->id;   assert(jj>=0&&jj<NN);
            /// Else, if R is on node j add 1/R to other diagonal (j,j) and add (-1/R) on elements (k,j) and (j,k)
            matrixA[jj][jj] += g1;
            matrixA[jj][kk] -= g1;
            matrixA[kk][jj] -= g1;
        }
        
    }
    
    apply_matrix_B_and_C_to_Matrix_A(matrixA);
    
    
    return matrixA;
}

#pragma mark matrix-Z related algorithm
void AnalogCircuit::add_current_to_matrix_Z(matrix_t MATRIX_Z, Component * comp, double time_step)
{
    matrix_t matrix_z  = MATRIX_Z;
    
    assert(comp->outputNode!=NULL);  assert(comp->inputNode!=NULL);
    
    if (comp->outputNode == ground_node) {
        std::swap(comp->outputNode,comp->inputNode);
    }
    
    int kk = comp->outputNode->id; assert(kk>=0 && kk<NN);
    
    double old_voltage = comp->voltage;
    double old_current = comp->get_current();

    double g = comp->get_conductance(time_step);
        
    //calculate new current
    double new_current = 0.0;
    if (comp->get_type() == Current_Source)
    {
        current_source * cs = (current_source*)comp;
        new_current = cs->get_current(time.elapsed_second);
    }
    else
    {
        return;
    }
    
    
    // there is a current source at node k. if so, the value is +current positive node
    matrix_z[kk][0] = new_current;
    
    //If other end is not on ground, then value is –current at other node.
    if (comp->inputNode != this->ground_node)
    {
        int kk = comp->inputNode->id;
        
        assert(kk>=0 && kk<NN);
        matrix_z[kk][0] = -new_current;
    }
}

matrix_t  AnalogCircuit::MNA_get_Z_Matric( double time_step)
{
    
    matrix_t matrix_z = matrix_init(row_num_z,col_num_z);
    
    for (int ii = 0; ii < current_sources.size(); ii++){
        add_current_to_matrix_Z(matrix_z, current_sources[ii], time_step);
    }
    

    
    /// If V has a positive terminal on node j, put 1 in row j.
    for (int ii = 0; ii < voltage_sources.size(); ii++)
    {
        voltage_source *vv = voltage_sources[ii];
        assert(NN+ii<row_num_z);
        matrix_z[NN+ii][0] = vv->get_voltage(time.elapsed_second);
    }
    
    return  matrix_z;
}



#pragma mark matrix-Q related algorithm

matrix_t AnalogCircuit::MNA_get_Q_matrix(matrix_t MATRIX_A,matrix_t MATRIX_Z)
{
    matrix_t matrix_Q = matrix_init(row_num_q, col_num_q);
    
    //first, copy A
    for (int ii = 0; ii< row_num_a; ii++){
        for (int jj = 0; jj < col_num_a; jj++){
            matrix_Q[ii][jj] = MATRIX_A[ii][jj];
        }
    }
    
    //then, copy Z
    for (int ii = 0; ii< row_num_z; ii++){
        for (int jj = 0; jj < col_num_z; jj++){
            matrix_Q[ii][jj+col_num_a] = MATRIX_Z[ii][jj];
        }
    }
    
    return matrix_Q;
}



void AnalogCircuit::apply_matrix_Q_to_component(matrix_t new_matrix_q, double time_step)
{
    //the top n elements are either zero or the sum and difference of independent current sources in the circuit.
    for (int ii = 0; ii < NN; ii++)
    {
        Node * node = node_list[ii];
        node->voltage = new_matrix_q[ii][col_num_q-1];
       
        
        //printf("node_%i %.3e\n", ii, node->voltage);
    }
    

    //the bottom m elements represent the m independent voltage sources in the circuit.
    for (int ii = 0; ii < voltage_sources.size(); ii++)
    {
        voltage_source *vs = voltage_sources[ii];
        vs->current = new_matrix_q[NN+ii][col_num_q-1];
        
        //printf("source_%i %.3e A\n", ii, vs->current);
    }
    
    
    for (int ii = 0; ii < component_list.size(); ii++)
    {
        Component *comp = component_list[ii];
        
        comp->voltage = comp->outputNode->voltage - comp->inputNode->voltage;
        comp->current = comp->get_conductance(time_step) * comp->voltage;
    }
}


void AnalogCircuit::update_state(double time_step)
{
    time.update_state(time_step);
    
    matrix_t matrixA = MNA_get_A_Matric(time_step);;

    matrix_t matrixZ = MNA_get_Z_Matric(time_step);

    matrix_t matrix_Q = MNA_get_Q_matrix(matrixA,matrixZ);
    
    
    if (print_matrix_during_simulation){
        printf("======unsolved Q matrix========\n");
        matrix_print(matrix_Q, row_num_q, col_num_q);
    }
    
    matrix_gauss(matrix_Q, row_num_q, col_num_q);
    
    if (print_matrix_during_simulation){
        printf("======solved Q matrix========\n");
        matrix_print(matrix_Q, row_num_q, col_num_q);
    }
    
    apply_matrix_Q_to_component(matrix_Q, time_step);
    
    matrix_free(matrixA, row_num_a);
    matrix_free(matrixZ, row_num_z);
    matrix_free(old_matrix_Q, row_num_q);
    old_matrix_Q = matrix_Q;
    
}



AnalogCircuit::~AnalogCircuit()
{
    matrix_free(old_matrix_Q, row_num_q);
}

void AnalogCircuit::log_component()
{
    //generate nodelist
    for (int ii = 0; ii < component_list.size();ii++){
        CLOG(component_list[ii]->toString());
    }
}


void AnalogCircuit::log_nodes()
{
    for (int ii = 0; ii < this->node_list.size(); ii++)
    {
        Node *nd = node_list[ii];
        CLOG(nd->toString());
    }
}