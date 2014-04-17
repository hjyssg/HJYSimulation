/* 
 * File:   VertexAndPointmassUnitTest.cpp
 * Author: junyang huang
 *
 *
 * Created on Mar 3, 2014, 3:49:03 PM
 */

#include <stdlib.h>
#include <iostream>
#include "HJYPhysics.h"

/*
 * so far, only test "Vertex" and "PointMass" class
 */



void testPointMass() {
   PointMass pointmass = PointMass(0,0,0, 
                             1,1,1,
                             1,1,1,
                             1);

    bool t1 = pointmass.pos.isEqual(Vertex(0,0,0));
    bool t2 = pointmass.spd.isEqual(Vertex(1,1,1));
    bool t3 = pointmass.acl.isEqual(Vertex(1,1,1));

   
      if (!t1||!t2||!t3) {
        std::cerr << "--TEST_FAILED--  testname=testPointMass()" << std::endl;
    }
}



void testGet_kinetic_energy() {
      PointMass p1 = PointMass(0,0,0,
                             1,1,1,
                             1,1,1,
                             1);
      bool t1 = (float)p1.get_kinetic_energy() == 0.5*p1.mass*p1.spd.length_square();
    
    if (!t1) {
        std::cerr << "--TEST_FAILED--  testname=testGet_kinetic_energy()" << std::endl;
    }
}



void testUpdate_state() {
    double time_step = 1.0;
    PointMass p1= PointMass(0,0,0,
                             1,1,1,
                             1,1,1,
                             1);
    p1.update_state(time_step);
    bool t1 = p1.pos.isEqual(Vertex(1,1,1));
    bool t2 = p1.spd.isEqual(Vertex(2,2,2));
    if (!t1||!t2 ) {
        std::cerr << "--TEST_FAILED--  testname=testUpdate_state()" << std::endl;
    }
}



void testAdd() {
    
    Vertex v1 = Vertex(1,1,1);
    Vertex v2 = Vertex(1,1,1);
    Vertex v3 = v1.add(v2);
    bool t1 = v3.x == 2&& v3.y ==2 &&v3.z == 2;

    if (!t1) {
        std::cerr << "--TEST_FAILED--  testname=testAdd()" << std::endl;
    }
}


void testAdd_self() {
    
    Vertex v1 = Vertex(1,1,1);
    Vertex v2 = Vertex(1,1,1);
    v1.add_self(v2);
    bool t1 = v1.x == 2&& v1.y ==2 &&v1.z == 2;

    if (!t1) {
        std::cerr << "--TEST_FAILED--  testname=testAdd_self()" << std::endl;
    }
}



void testCrosss_product() {
    Vertex v2 = Vertex(312, 231, 132);
    Vertex vertex =    Vertex(12,33,4142);
    Vertex result = vertex.crosss_product(v2);
   
    if (!result.isEqual( -952446,1290720,-7524 )) {
        std::cerr << "--TEST_FAILED--  testname=testCrosss_product()" << std::endl;
    }
}


void testDistance_square_to_v() {
    Vertex v = Vertex(1,2,3);
    Vertex vertex = Vertex(10,3,0);
    double result = vertex.distance_square_to_v(v);
    if ((float)result != (9*9+1+9)) {
        std::cerr << "--TEST_FAILED--  testname=testDistance_square_to_v()" << std::endl;
    }
}



void testDistance_to_v() {
    Vertex v = Vertex(1,2,3);
    Vertex vertex = Vertex(10,3,0);
    double result = vertex.distance_to_v(v);
    if ((float)result != sqrtf(9*9+1+9)) {
        std::cerr << "--TEST_FAILED--  testname=testDistance_to_v()" << std::endl;
    }
}


void testDot_product() {
    Vertex v2 = Vertex(3, 4, 5);
    Vertex vertex = Vertex(2,4,6);
    double result = vertex.dot_product(v2);
    if (result != (6+16+30)) {
        std::cerr << "--TEST_FAILED--  testname=testDot_product()" << std::endl;
    }
}


void testGet_copy() {
    Vertex vertex = Vertex(1.2, 1.3, 1.2);
    Vertex result = vertex.get_copy();
    if (!result.isEqual(vertex)) {
        std::cerr << "--TEST_FAILED--  testname=testGet_copy()" << std::endl;
    }
}


void testIsEqual() {
    Vertex v = Vertex(1,2,3);
    Vertex vertex = Vertex(1,2,3);
    bool result = vertex.isEqual(v);
    if (!result) {
        std::cerr << "--TEST_FAILED--  testname=testIsEqual()" << std::endl;
    }
}



void testIsParrell() {
    Vertex v = Vertex(10, 10, 10);
    Vertex vertex = Vertex(2, 2, 2);
    bool result = vertex.isParrell(v);
    if (!result) {
        std::cerr << "--TEST_FAILED--  testname=testIsParrell()" << std::endl;
    }
}



void testIsZero() {
    Vertex vertex;
    bool result = vertex.isZero();
    if (!result) {
        std::cerr << "--TEST_FAILED--  testname=testIsZero()" << std::endl;
    }
}



void testLength() {
    Vertex vertex = Vertex(3,4,0);
    double result = vertex.length();
    if ((float)result != 5.0) {
        std::cerr << "--TEST_FAILED--  testname=testLength()" << std::endl;
    }
}

void testLength_square() {
    Vertex vertex = Vertex(3,4,0);
    double result = vertex.length_square();
    if ((float)result != 3*3+4*4) {
        std::cerr << "--TEST_FAILED--  testname=testLength_square()" << std::endl;
    }
}



void testNormalize() {
    Vertex vertex = Vertex(100,-200,200);
    Vertex result = vertex.normalize();
    if ((float)result.length()!= 1.0) {
        std::cerr << "--TEST_FAILED--  testname=testNormalize()" << std::endl;
    }
}


void testNormalize_self() {
    Vertex vertex = Vertex(100,-200,200);
    vertex.normalize_self();
    if ((float)vertex.length()!= 1.0) {
        std::cerr << "--TEST_FAILED--  testname=testNormalize_self()" << std::endl;
    }
}



void testScale() {
    Vertex vertex = Vertex(-1, 0.1, 4);
    Vertex result = vertex.scale(2, 5, 8);
    if (!result.isEqual(-2, 0.5,4*8)) {
        std::cerr << "--TEST_FAILED--  testname=testScale()" << std::endl;
    }
}



void testScale_self() {
    double x_s = 10;
    double y_s = 20;
    double z_s = 20;
    Vertex vertex = Vertex(1,1,0);
    vertex.scale_self(x_s, y_s, z_s);
    if (!vertex.isEqual(10,20,0)) {
        std::cerr << "--TEST_FAILED--  testname=testScale_self()" << std::endl;
    }
}



void testScale_self2() {
    double t = 10;
    Vertex vertex = Vertex(1,1,0);
    vertex.scale_self(t);
    if (!vertex.isEqual(10,10,0)) {
        std::cerr << "--TEST_FAILED--  testname=testScale_self2()" << std::endl;
    }
}



void testSet_xyz() {
    double x_p = -10;
    double y_p = 10000.0;
    double z_p = 1e123;
    Vertex vertex;
    vertex.set_xyz(x_p, y_p, z_p);
    if (!vertex.isEqual(x_p,y_p,z_p)) {
        std::cerr << "--TEST_FAILED--  testname=testSet_xyz()" << std::endl;
    }
}


void testSub() {
    Vertex v2 = Vertex(1,2,3);
    Vertex vertex = Vertex();
    Vertex result = vertex.sub(v2);
    if (!result.isEqual(-1, -2, -3)) {
        std::cerr << "--TEST_FAILED--  testname=testSub()" << std::endl;
    }
}



void testSub_self() {
    Vertex v2 = Vertex(10,16,-20);
    Vertex vertex;
    vertex.sub_self(v2);
    if (!vertex.isEqual(-10,-16, 20)) {
        std::cerr << "--TEST_FAILED--  testname=testSub_self()" << std::endl;
    }
}


void testTransfer() {
    double x_p = 1.1;
    double y_p = -200;
    double z_p = 300;
    Vertex vertex;
    Vertex result = vertex.transfer(x_p, y_p, z_p);
    if (!result.isEqual(x_p, y_p, z_p)) {
        std::cerr << "--TEST_FAILED--  testname=testTransfer()" << std::endl;
    }
}



void testTransfer_self() {
    double x_p = 1.43214;
    double y_p = -200.12312;
    double z_p = 300.213124234;
    Vertex vertex;
    vertex.transfer_self(x_p, y_p, z_p);
    if (!vertex.isEqual(x_p, y_p, z_p)) {
        std::cerr << "--TEST_FAILED--  testname=testTransfer_self()" << std::endl;
    }
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% HJYPhysicsUnitTest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% testPointMass()" << std::endl;
    testPointMass();
    std::cout << "%TEST_FINISHED%  testPointMass()" << std::endl;

    std::cout << "%TEST_STARTED% testGet_kinetic_energy()" << std::endl;
    testGet_kinetic_energy();
    std::cout << "%TEST_FINISHED%  testGet_kinetic_energy()" << std::endl;

    std::cout << "%TEST_STARTED% testUpdate_state()" << std::endl;
    testUpdate_state();
    std::cout << "%TEST_FINISHED%  testUpdate_state()" << std::endl;

    std::cout << "%TEST_STARTED% testAdd()" << std::endl;
    testAdd();
    std::cout << "%TEST_FINISHED%  testAdd()" << std::endl;

    std::cout << "%TEST_STARTED% testAdd_self()" << std::endl;
    testAdd_self();
    std::cout << "%TEST_FINISHED%  testAdd_self()" << std::endl;

    std::cout << "%TEST_STARTED% testCrosss_product()" << std::endl;
    testCrosss_product();
    std::cout << "%TEST_FINISHED%  testCrosss_product()" << std::endl;

    std::cout << "%TEST_STARTED% testDistance_square_to_v()" << std::endl;
    testDistance_square_to_v();
    std::cout << "%TEST_FINISHED%  testDistance_square_to_v()" << std::endl;

    std::cout << "%TEST_STARTED% testDistance_to_v()" << std::endl;
    testDistance_to_v();
    std::cout << "%TEST_FINISHED%  testDistance_to_v()" << std::endl;

    std::cout << "%TEST_STARTED% testDot_product()" << std::endl;
    testDot_product();
    std::cout << "%TEST_FINISHED%  testDot_product()" << std::endl;

    std::cout << "%TEST_STARTED% testGet_copy()" << std::endl;
    testGet_copy();
    std::cout << "%TEST_FINISHED%  testGet_copy()" << std::endl;

    std::cout << "%TEST_STARTED% testIsEqual()" << std::endl;
    testIsEqual();
    std::cout << "%TEST_FINISHED%  testIsEqual()" << std::endl;

    std::cout << "%TEST_STARTED% testIsParrell()" << std::endl;
    testIsParrell();
    std::cout << "%TEST_FINISHED%  testIsParrell()" << std::endl;

    std::cout << "%TEST_STARTED% testIsZero()" << std::endl;
    testIsZero();
    std::cout << "%TEST_FINISHED%  testIsZero()" << std::endl;

    std::cout << "%TEST_STARTED% testLength()" << std::endl;
    testLength();
    std::cout << "%TEST_FINISHED%  testLength()" << std::endl;

    std::cout << "%TEST_STARTED% testLength_square()" << std::endl;
    testLength_square();
    std::cout << "%TEST_FINISHED%  testLength_square()" << std::endl;

    std::cout << "%TEST_STARTED% testNormalize()" << std::endl;
    testNormalize();
    std::cout << "%TEST_FINISHED%  testNormalize()" << std::endl;

    std::cout << "%TEST_STARTED% testNormalize_self()" << std::endl;
    testNormalize_self();
    std::cout << "%TEST_FINISHED%  testNormalize_self()" << std::endl;

    std::cout << "%TEST_STARTED% testScale()" << std::endl;
    testScale();
    std::cout << "%TEST_FINISHED%  testScale()" << std::endl;

    std::cout << "%TEST_STARTED% testScale_self()" << std::endl;
    testScale_self();
    std::cout << "%TEST_FINISHED%  testScale_self()" << std::endl;

    std::cout << "%TEST_STARTED% testScale_self2()" << std::endl;
    testScale_self2();
    std::cout << "%TEST_FINISHED%  testScale_self2()" << std::endl;

    std::cout << "%TEST_STARTED% testSet_xyz()" << std::endl;
    testSet_xyz();
    std::cout << "%TEST_FINISHED%  testSet_xyz()" << std::endl;

    std::cout << "%TEST_STARTED% testSub()" << std::endl;
    testSub();
    std::cout << "%TEST_FINISHED%  testSub()" << std::endl;

    std::cout << "%TEST_STARTED% testSub_self()" << std::endl;
    testSub_self();
    std::cout << "%TEST_FINISHED%  testSub_self()" << std::endl;

    std::cout << "%TEST_STARTED% testTransfer()" << std::endl;
    testTransfer();
    std::cout << "%TEST_FINISHED%  testTransfer()" << std::endl;

    std::cout << "%TEST_STARTED% testTransfer_self()" << std::endl;
    testTransfer_self();
    std::cout << "%TEST_FINISHED%  testTransfer_self()" << std::endl;

    std::cout << "%SUITE_FINISHED% " << std::endl;

    return (EXIT_SUCCESS);
}

