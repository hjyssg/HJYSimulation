//
//  main.cpp
//  HJYAnalog
//
//  Created by Junyang Huang on 2/25/14.
//  Copyright (c) 2014 www.hjyssg.com. All rights reserved.
//

#include <iostream>
#include <math.h>
#include <algorithm>   
#include "HJYAnalog.h"

const float PI =  3.14159265;

//a flag to disable/enable print
bool  PRINT_FLAG = true;   

matrix_t matrix_init(int row_num, int col_num)
{
    double** matrix = (double**) malloc(row_num*sizeof(double*));
    for (int ii = 0; ii < row_num; ii++)
    {
        matrix[ii] = (double*) malloc(col_num*sizeof(double));
        
        //init all element to zero
        //C would not do for you
        for (int jj = 0; jj < col_num; jj++)
        {
            matrix[ii][jj] = 0.0;
        }
        
    }
    return matrix;
}

bool matrix_gauss(matrix_t M, int row_num, int col_num)
{
    for (int diag = 0; diag < row_num; diag++)
    {
        /* for each row = variable = diagonal element */
        double pivot = M[diag][diag];
        if (pivot == 0) {
            //printf("Gauss: Singular Matrix!\n");
            //return false;
            continue;
        }
        
        /* from diagonal to right, divide out by the pivot diagonal element */
        for (int col = diag; col < col_num; col++)
        {
            M[diag][col] /= pivot;
        }
        
        /* for all other rows subtract scalled diag row, zeroing pivot column */
        for(int row = 0; row < row_num; row++)
        {
            if (row != diag)
            {
                /* scale all elements of diag row by value which will zero diag element of this row,
                 and subtract from this row */
                double pivot2 = M[row][diag];
                for (int col = diag; col < col_num; col++)
                {
                    M[row][col] -= M[diag][col]*pivot2;
                }
            }
        }
    }
    
    return true;
}


matrix_t matrix_transpose(matrix_t m, int row_num, int col_num)
{
    double ** m2 =  matrix_init(col_num, row_num);
    
    for (int ii = 0; ii< row_num; ii++)
    {
        for (int jj = 0; jj < col_num; jj++)
        {
            m2[jj][ii] = m[ii][jj];
        }
    }
    
    return m2;
}



void matrix_print(matrix_t M, int row_num, int col_num)
{
#if 0
    char * str =(char *) malloc(col_num * 15 * sizeof(char));
    char temp[15];
    
    for (int ii = 0; ii < row_num; ii++)
    {
        for (int jj = 0; jj< col_num; jj++)
        {
            if (col_num == 1)
            {
                sprintf(temp,"|%.3e|",  M[ii][jj]);
                strcat(str, temp);
            }
            else if (jj==0)
            {
                sprintf(temp,"|%.3e, ",  M[ii][jj]);
                strcat(str, temp);
            }
            else if (jj!=col_num-1)
            {
                sprintf(temp,"%.3e, ",  M[ii][jj]);
                strcat(str, temp);
            }
            else
            {
                sprintf(temp,"%.3e|\n",  M[ii][jj]);
                strcat(str, temp);
            }
            
            printf("%s",str);
            str[0] ='\0'; // erase the string
        }
    }
    
    free(str);
    
#else
    for (int ii = 0; ii < row_num; ii++) {
        printf("|");
        for (int jj = 0; jj <col_num; jj++)
        {
            printf("%.2e ",M[ii][jj]);
        }
        printf("|\n");
    }
    
#endif
}


void matrix_free(matrix_t M, int row_num)
{
    for(int ii = 0; ii < row_num; ii++)
    {
        free(M[ii]);
    }
    
    free(M);
}

void matrix_testing()
{
    int row_num = 3;
    int col_num = 4;
    
    double ** m1 = matrix_init(row_num, col_num);
    
    for (int ii = 0; ii < row_num; ii++)
    {
        for (int jj = 0; jj< col_num; jj++)
        {
            m1[ii][jj] = (ii+1) * (jj+1) * (rand()%50);
        }
    }
    
    printf("init matrix\n");
    matrix_print(m1, row_num, col_num);
    
    
    printf("its transpose\n");
    double ** m2 = matrix_transpose(m1, row_num, col_num);
    matrix_print(m2, col_num, row_num);
    
    matrix_gauss(m1, row_num, col_num);
    
    printf("after gauss elimination\n");
    matrix_print(m1, row_num, col_num);
}

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
        

            cc = new resistor(2);
        
        
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

