HJYPhysics
==========

##a simple C++ OOP 3D physics engine


This engine is for simple scientific simulation. It only consists one header file and one source file. So it is very  easy to import and use in your project.  
HJYPhysics uses **double** for all computation to offer high precision.    
All computation use the **metric system**: Meter, Second.  
If your compiler supports OpenMp, please enable it. This engine supports openMp multithreading and will increase performance. If not, don't worry, the code will run without problem.



####It has following class:  

* **Vertex**:  It represents a mathematical point in the 3D space. It has x, y, z coordinates and related methods: add,substract, dot-product, cross-product and etc. It is used to represent position, velocity and acceleration bt the following classes.
*  **PointMass**: the most basic component in Newton Physics. It has position, velocity, acceleration, mass. Point mass does not have volume or radius. It is a point with mass. 
* **RigidSphere**: subclassing from PointMass. Different than PointMass, It has radius and can collide with other rigidSphere.
* **SpringBond**: The representation of a spring connecting two point mass objects and follows the Hooke's Law. 
*  **Molecule**: A simplified representation of molecule. It used “Pointmass” and “springbond” to represent atom and interatomic bonds. 
*  **TheWorldTimer**: A small timer used by “theWorld” to track the time.
*  **TheWorld**: The representation of Newton physics world. Point mass and other objects will be added to it and interact each the time steps. By specifying flags, it can integrate by using different integration methods and enable/disable different force computation. Please don't turn on all force computaion flag. Only keep on nessary flags, otherwise the computaion will be slow.

####pseudo-code: how a simulation program looks like:
	initilize a theWorld instance and specify computation flags
	create pointMass and other object
	add them to theWorld engine
	theWorld init
	
	while(true)
	{
	   theWorld.update_state()
	   
	   visualize the state 
	   or log the state
	}



####Visualization:  
The examples use VPython to visualize the simulation result.  
In each example, there is a VPython_flag. By setting it to 1, the visualization will be turned on. If set to 0, the program will run without visualization.  
To set up on windows machine,First, download and install the 32-bit Python-3.2.2 from python.org. Second, download and install VPython-Win-Py3.2-5.74. In Visual studio, add VPython directory to the c++ linker directory. 
For more detial, you can also refer to this [page](http://kona.ee.pitt.edu/1180wiki/doku.php?id=how_to_mix_c_and_python). 


####unit test
In the "unittest" folder, I wrote the a simple unit test for "Vertex" and "Pointmass". If I have more time, I will cover more in my unit test.

####screenshots
![screenshot1](./screenshot/Bouncing_ball_SS.png?raw=true) 
[bouncing ball](http://youtu.be/mDeFVCZtp5Y "bouncing ball")
![screenshot1](./screenshot/Matrix_system_SS.png?raw=true)
[spring and ball matrix](http://youtu.be/Qg-D4nlz-s8 "spring and ball matrix")
![screenshot1](./screenshot/Solar_system_SS.png?raw=true)
[solar system](http://youtu.be/HIRNhrII4ho "solar system")








 
