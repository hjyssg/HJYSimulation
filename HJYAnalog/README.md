###HJYAnalog

A simple analog circuit simulation API.

It supports three analog components: 

* Voltage Source
* Current Source
* Resistor


It uses [MNA](http://www.swarthmore.edu/NatSci/echeeve1/Ref/mna/MNA5.html) internally to solve the circuit.   
You need to create components and connect components with nodes to form a circuit. The engine will solve the circuit by time-step.  
For more detail, check the main file.  


![screenshot1](./screenshot/1.png?raw=true)