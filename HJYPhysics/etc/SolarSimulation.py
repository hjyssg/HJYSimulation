#!/usr/bin/python
# -*- coding: utf-8 -*-

#########################
#
#  created by hjyssg
#
#########################

from visual import *
from visual.graph import *

SUN_MASS = 1.99e+30  # kg
SUN_DIAMETER = 1391e+6  # m
MERCURY_MASS = 0.33e24  # kg
MERCURY_DIAMETER = 4897e3  # m
VENUS_MASS = 4.87e24  # kg
VENUS_DIAMETER = 12104e3  # m
EARTH_MASS = 5.972e+24  # kg
EARTH_DIAMETER = 12756e3  # m
MARS_MASS = 0.642e24  # kg
MARS_DIAMETER = 6792e3  # m
JUPITER_MASS = 1898e24  # kg
JUPITER_DIAMETER = 142984e3
SATURN_MASS = 568e24  # kg
SATURN_DIAMETER = 120536e3  # m
URANUS_MASS = 86.8e24  # kg
URANUS_DIAMETER = 51118e3  # m
NEPTUNE_MASS = 102e24  # kg
NEPTUNE_DIAMETER = 49528e3  # m

scene1 = display(autocenter=False, x=400, width=900, height=650,
                 scale=(7e-13 * 2, 7e-13 * 2, 7e-13 * 2))

xaxis = curve(pos=[(0, 0, 0), (1e13, 0, 0)], color=color.red)
Yaxis = curve(pos=[(0, 0, 0), (0, 1e13, 0)], color=color.green)
Zaxis = curve(pos=[(0, 0, 0), (0, 0, 1e13)], color=color.blue)

KM = 1000.0

##we do not display the real size of planets, otherwise you will see nothing
r1 = SUN_DIAMETER * 20
r2 = SUN_DIAMETER * 2

#init based on the data from http://ssd.jpl.nasa.gov/horizons.cgi
# date 2013-01-01
# solar system barcenter
sun = sphere(color=color.red, make_trail=True, radius=r1)
sun.pos = vector(-1.921258692976546E+05 * KM, -3.672858221449571E+05
                 * KM, -6.293939186492928E+03 * KM)
sun.spd = vector(1.064112061805785E-02 * KM, -2.282152601144673E-03
                 * KM, -2.261768099104889E-04 * KM)
sun.mass = SUN_MASS
mercury = sphere(color=color.orange, make_trail=True, radius=r2)
mercury.pos = vector(-2.542271688725740E+07 * KM,
                     -6.518047354640916E+07 * KM,
                     -2.986942437682569E+06 * KM)
mercury.spd = vector(3.561299513603979E+01 * KM, -1.527095279264177E+01
                     * KM, -4.514379444956798E+00 * KM)
mercury.mass = MERCURY_MASS
venus = sphere(color=color.blue, make_trail=True, radius=r2)
venus.pos = vector(-6.905541781609924E+07 * KM, -8.392401643421896E+07
                   * KM, 2.823108036887005E+06 * KM)
venus.spd = vector(2.679084974988556E+01 * KM, -2.244015730475055E+01
                   * KM, -1.853249972810995E+00 * KM)
venus.mass = VENUS_MASS
earth = sphere(material=materials.earth, make_trail=True, radius=r2)
earth.pos = vector(-2.712311750453984E+07 * KM, 1.442452249801158E+08
                   * KM, -1.024866474659159E+04 * KM)
earth.spd = vector(-2.975221878384466E+01 * KM, -5.557630236369244E+00
                   * KM, 1.166327378957226E-05 * KM)
earth.mass = EARTH_MASS
mars = sphere(color=color.green, make_trail=True, radius=r2)
mars.pos = vector(1.613831829365655E+08 * KM, -1.299959988813151E+08
                  * KM, -6.689343082176655E+06 * KM)
mars.spd = vector(1.609480750581124E+01 * KM, 2.096710871054539E+01
                  * KM, 4.420934607818046E-02 * KM)
mars.mass = MARS_MASS
jupiter = sphere(color=color.cyan, make_trail=True, radius=r2)
jupiter.pos = vector(2.130248609275347E+08 * KM, 7.264621031147093E+08
                     * KM, -7.796178050164084E+06 * KM)
jupiter.spd = vector(-1.269579412745502E+01 * KM, 4.300631432663832E+00
                     * KM, 2.662458765460191E-01 * KM)
jupiter.mass = JUPITER_MASS
saturn = sphere(color=color.blue, make_trail=True, radius=r2)
saturn.pos = vector(-1.209082909091936E+09 * KM, -8.253529571281528E+08
                    * KM, 6.246931966183711E+07 * KM)
saturn.spd = vector(4.923237564225631E+00 * KM, -8.001281989714242E+00
                    * KM, -5.660410531867004E-02 * KM)
saturn.mass = SATURN_MASS
uranus = sphere(color=color.green, make_trail=True, radius=r2)
uranus.pos = vector(2.975359717954440E+09 * KM, 3.846013966774312E+08
                    * KM, -3.711880253105340E+07 * KM)
uranus.spd = vector(-9.228235466197581E-01 * KM, 6.436376882337913E+00
                    * KM, 3.602154409843563E-02 * KM)
uranus.mass = URANUS_MASS
neptune = sphere(color=color.white, make_trail=True, radius=r2)
neptune.pos = vector(3.973252374554337E+09 * KM, -2.083244177306763E+09
                     * KM, -4.866720990954759E+07 * KM)
neptune.spd = vector(2.487210727347481E+00 * KM, 4.845277256620428E+00
                     * KM, -1.572810424785901E-01 * KM)
neptune.mass = NEPTUNE_MASS

sun.name = 'sun'
mercury.name = 'mercury'
venus.name = 'venus'
earth.name = 'earth'
mars.name = 'mars'
jupiter.name = 'jupiter'
saturn.name = 'saturn'
uranus.name = 'uranus'
neptune.name = 'neptune'

solar_system_objects = list()
solar_system_objects.append(sun)
solar_system_objects.append(mercury)
solar_system_objects.append(venus)
solar_system_objects.append(earth)
solar_system_objects.append(mars)
solar_system_objects.append(jupiter)
solar_system_objects.append(saturn)
solar_system_objects.append(uranus)
solar_system_objects.append(neptune)

# //create label

time = 0
time_label = label(yoffset=280, line=0)
solar_object_labels = list()

for ii in xrange(0, len(solar_system_objects)):
    obj = solar_system_objects[ii]
    obj.acl = vector(0.0, 0.0, 0.0)
    temp_l = label(yoffset=5, opacity=1, line=0, box=0)
    temp_l.text = obj.name
    solar_object_labels.append(temp_l)

init_total_ke = 0.0
for obj in solar_system_objects:
    init_total_ke += mag2(obj.spd) * obj.mass / 2.0

# create graph
energy_graph = gdisplay(
    title='kinectic energy',
    xtitle='day',
    x=0,
    y=660,
    width=900,
    height=300,
    ytitle='energy',
    )
simulated_kinetic_line = gcurve(color=color.green)
real_kinetic_line = gcurve(color=color.blue)

GRAVITATIONAL_CONSTANT = 6.67384e-11  # m3*kg-1*s-2
time_interval = 1000.0

while true:
    for obj in solar_system_objects:
        obj.spd += obj.acl * time_interval / 2
        obj.pos += obj.spd * time_interval

    # calculal acl for each object
    for ii in xrange(0, len(solar_system_objects)):
        obj1 = solar_system_objects[ii]
        obj1.acl = vector(0.0, 0.0, 0.0)
        for jj in xrange(0, len(solar_system_objects)):
            if jj == ii:
                continue
            obj2 = solar_system_objects[jj]
            direction = obj1.pos - obj2.pos
            dist = mag2(direction)
            norm_v = norm(direction) * -1
            gravitation_acl = obj2.mass * norm_v \
                * GRAVITATIONAL_CONSTANT / dist
            obj1.acl += gravitation_acl

        # #update position and velocity
        obj1.spd += obj1.acl * time_interval
        obj1.pos += obj1.spd * time_interval

        # also update label
        solar_object_labels[ii].pos = obj1.pos

    for obj in solar_system_objects:
        obj.spd += obj.acl * time_interval / 2.0

    total_ke = 0.0
    for obj in solar_system_objects:
        total_ke += mag2(obj.spd) * obj.mass / 2.0

    # update plot

    simulated_kinetic_line.plot(pos=(time, total_ke * 1e-35  ))
    real_kinetic_line.plot(pos=(time, init_total_ke * 1e-35 ))

    time += time_interval
    time_label.text = '%d day ' % (time / 3600 / 24)
    rate(200)


