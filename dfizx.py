import copy
from vectors import *
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import numpy as np
import map_handler

#####Coordinates
  #x,y,z coordinates in meters      
#####Velocity 
  #x,y,z velocities and rotaitonal velocity in m/s  
#####Angles
  #Everything in radians
  #AoA = angle of attack
  #aoar: angle of attack in radians

#####Unit Vectors
  #vel_unit:  unit vector of total disc velocity
  #lift_unit: unit vector 90 degrees from vel_unit in line with discs 'up' vector
"""
         ^  lift_unit
         |
         |
         |/
         /--vel_unit--->  direction of travel
        /     
       
  disc is facing up in this diagrams, 
  lift unit would to/away from the screen if the disc were sideways 
  or down if the disc were upside down
"""    
  #disc_normal_unit:  unit vector of disc which points 'up' relative to the disc
  #disc_unit_y:  unit vector of velocity vector projected onto discs planes (angle between this and vel_unit is angle of attack)
  #disc_unit_x:  unit vector of disc which points perpendicular to direction of travel but lays on the discs plane
    
"""   Disc units        z is disc_normal_unit which is normal to the disc's plane, x and y are both contained by the plane of the disc
      _____
     /     \
    (   z y )    ---->  direction of travel
     \__x__/
"""  
  #rotation direction goes from disc_unit_x to disc_unit_y
  #spin is in rads/s
  
  
# a Disc contains all of the aerodynamic properties of a specific type of disc
temporaryadjustment = 57.3 *10
class Disc (object):
  def __init__(self,Cl0,Cla,Cd,Cm0,Cma,Acp,diameter = 0.21,):
    #Zero pitch lifting coefficient
    self.Cl0=Cl0
    #Alpha (AoA) lift coefficient
    self.Cla=np.deg2rad(Cla) * temporaryadjustment
    #Non stalled Drag Coefficient of disc
    self.Cd=Cd
    #CMo zero lift pitching moment coefficient 
    self.Cm0=Cm0
    #Alpha pitching moment coefficient
    self.Cma=np.deg2rad(Cma) * temporaryadjustment
    #Chord point (currently unused)
    self.Acp = Acp
    #Disc diameter
    self.diameter = diameter
    self.radius = diameter/2
    #Disc area
    self.area = math.pi * (diameter/2) ** 2


    
####Disc List          first set is from a thesis experiment, Cd appeared to be too high for them
#disc_ =        Disc( Cl0 , Cla , Cd  ,  Cm0 , Cma)         disc_ = Disc()
disc_Aviar    = Disc(0.152,0.044,0.083,-0.018,0.002,0.05)# 2|3| 0|1
disc_Buzz     = Disc(0.099,0.041,0.061,-0.033,0.004,0.10)# 5|4|-1|1
disc_Roc      = Disc(0.053,0.043,0.067,-0.015,0.003,0.07)# 4|4| 0|3
disc_Flick    = Disc(0.100,0.038,0.076,-0.007,0.008,0.21)# 9|4| 1|4
disc_Storm    = Disc(0.107,0.045,0.057,-0.026,0.004,0.09)# 8|4| 0|3     %%%grain of salt
disc_Wraith   = Disc(0.143,0.040,0.055,-0.020,0.006,0.15)#11|5|-1|3
disc_QuarterK = Disc(0.138,0.039,0.065,-0.038,0.005,0.13)

#Custom made for simulation discs
disc_Driver_OS   = Disc(0.05,0.04,0.020,-0.002,0.001,0.21)
disc_Driver_S    = Disc(0.11,0.04,0.020,-0.003,0.002,0.21)
disc_Driver_US   = Disc(0.15,0.04,0.020,-0.004,0.003,0.21)

disc_Mid_OS      = Disc(0.10,0.042,0.040,-0.005,0.0025,0.21)
disc_Mid_S      = Disc(0.14,0.042,0.040,-0.010,0.0040,0.21)
disc_Mid_US      = Disc(0.16,0.042,0.040,-0.015,0.0050,0.21)

disc_Putter      = Disc(0.15,0.044,0.055,-0.010,0.0010,0.05)




#The location and orientation of a disc at a moment in time
class DiscFlight (object):
  def __init__(self,disc,throw_hyzer_angle,velocity,time_coefficiant):
    self.pos_x = np.array([])
    self.pos_y = np.array([])
    self.pos_z = np.array([])
    self.orient_x = np.array([])
    self.orient_y = np.array([])
    self.orient_z = np.array([])
    self.disc = copy.copy(disc)
    self.throw_hyzer_angle = copy.copy(throw_hyzer_angle)
    self.velocity = copy.copy(velocity)
    self.time_coefficiant = copy.copy(time_coefficiant)
    
  def record_state (self,coordinates,orientation):
    self.pos_x = np.append(self.pos_x, coordinates.x)
    self.pos_y = np.append(self.pos_y, coordinates.y)
    self.pos_z = np.append(self.pos_z, coordinates.z)
    self.orient_x = np.append(self.orient_x, orientation.x)
    self.orient_y = np.append(self.orient_y, orientation.y)
    self.orient_z = np.append(self.orient_z, orientation.z)

  def finish_flight (self):
    self.seconds_of_travel = self.pos_x.size * self.time_coefficiant
    self.max_height = np.amax(self.pos_z) 
    self.touch_times = (np.where(self.pos_z<= 0,self.pos_z*0+1,self.pos_z*0))
    self.touch_times = np.nonzero(self.touch_times)
    index_on_first_touch = self.touch_times[0][0]
    self.throw_hang_time = index_on_first_touch * self.time_coefficiant
    self.air_distance_thrown = math.sqrt((self.pos_x[index_on_first_touch]-self.pos_x[0])**2+(self.pos_y[index_on_first_touch]-self.pos_y[0])**2+(self.pos_z[index_on_first_touch]-self.pos_z[0])**2)
    self.total_distance_thrown = math.sqrt((self.pos_x[-1]-self.pos_x[0])**2+(self.pos_y[-1]-self.pos_y[0])**2+(self.pos_z[-1]-self.pos_z[0])**2)
    print ("Calculating flight stats")
    print ("Seconds of travel: "+ str(self.seconds_of_travel) )
    print ("Maximum height   : "+ str(self.max_height))
    print ("Airtime on throw : "+ str(self.throw_hang_time))
    print ("Distance of throw: "+ str(self.air_distance_thrown))
    print ("Total distance   : "+ str(self.total_distance_thrown))
#####Main function, simulates a disc throw 
#sim time in ms    should be adjustable without affecting results....should be < 10m when collisions are introduced
delta_t =  50
#returns a Trajectory  (list of DiscStates)         %%%%Todo, add hand and throw type variable LHBH, RHFH...%%%%
def do_throw (coordinates,velocity,throw_hyzer_angle= 0,disc = disc_Putter):
  
  ####Flags
  checking_for_aerial_collisions = True
  
  
  
  
  
  
  #Where the disc is thrown from as a vector in meters
  coords = copy.copy(coordinates)
  #The direction of travel for the disc in m/s
  vel = copy.copy(velocity)
  #Disc being thrown
  disc = disc
  #Throw angle is the amount of hyzer/anhyzer on the disc. Positive number are degrees of hyzer, negative numbers are degress of anhyzer
  
  
  
  
  
  #some required initializing and constants
  global_up_vector = Vector (0,0,1)
  global_down_vector = Vector (0,0,1)
  Wr = 0
  #rotations_per_meter when one edge is stationary
  rotations_per_meter = 1/ (disc.diameter*math.pi)
  #radians per second of rotational speed         starts at approx perfect spin (zero rim velocity on release side) 
  rps = vel.magnitude() * rotations_per_meter
  #air density
  p = 1.225
  #disc masss
  mass = 0.170
  a = disc.area 
  #moment of inertia      %%%%0.9 is approx k value, update and add to discs%%%%%%
  Im = mass * (0.9*disc.radius)**2
  
  #amount of time between each simulation in fractional seconds
  time_coefficiant = delta_t / 1000
  #useful for accelerations
  seconds_over_mass = time_coefficiant / mass
  # G = 9.8 m/s^2
  GRAV = 9.80665

  gravitational_force_magnitude = -GRAV * mass
  gravitational_force_vector = global_down_vector.multiply(gravitational_force_magnitude)
  
  
  
  
  
  
  #####Make disc_normal_unit from the velocity vector of the disc
  
  ##Pre Hyzer Angle
  #this is only disc_x_unit until hyzer is added later
  disc_unit_x = vel.cross(global_up_vector)
  disc_unit_x = make_unit_vector (disc_unit_x)
  
  disc_normal_unit= disc_unit_x.cross(vel)
  disc_normal_unit= make_unit_vector (disc_normal_unit)
  
  #adds hyzer/anhyzer angle to the orientation of the disc
  #angle of hyzer
  throw_hyzer_angle = np.deg2rad(throw_hyzer_angle)
  #disc_normal_unit and disc_unit_x are perpendicular unit vectors which describe a plane. 
  #We can substitue these vectors for 2D x and y trig calculations 
  orientation_component_y = disc_normal_unit.multiply(math.cos (throw_hyzer_angle))
  orientation_component_x =      disc_unit_x.multiply(math.sin (throw_hyzer_angle))
  #adding the x and y components to create the resultant vector
  disc_normal_unit= orientation_component_y.sum (orientation_component_x)
  
  
  
  
  
  # Trajectory 
  #is a list of DiscStates that describes the flight of a disc in constant time increments
  disc_flight = DiscFlight (disc,throw_hyzer_angle,velocity,time_coefficiant)
  disc_flight.record_state(coords,disc_unit_x)
  
  #####Main Simulation Loop
  #loop_count is for emergency exit conditions and tests
  loop_count = 0
  #                      %%%Todo: find better exit conditon%%%
  print ("Starting Speed: " + str(vel.magnitude()) + "      starting sping in rps: "+ str(rps))
  while loop_count<400 and vel.magnitude()>2:
  
    
    
    
    
    
    
    ###Disc State Detection
    ####Skipping | Sliding | Rolling
    if coords.z <0:
      
      
      #contact point between the disc and ground (disc orientation normal projected onto the ground plane)
      collision_vector = Vector (disc_normal_unit.x,disc_normal_unit.y,0)
      
      if disc_normal_unit.angle(global_up_vector) < 0.15 or disc_normal_unit.angle(Vector(0,0,-1)) < 0.15 or rps < 10:
        #     %%%%needs work maybe?%%%
        print ("slide   on frame: " + str(loop_count))
        
        vel.z *= -0.1
        disc_normal_unit= global_up_vector
        vel =vel.divide(1+ 2 * time_coefficiant)
      elif disc_unit_x.angle(global_up_vector) < 1.57:    #roll
        print ("roll   on frame: " + str(loop_count))
        
        #flag variable to be checked during moments secti
        rolling_moment = 1
        rps = vel.magnitude() * rotations_per_meter
        # Stop lateral sliding is an intersection of 2 planes problem...the line they create is the crossproduct of the normal of both planes
        vel_unit = make_unit_vector (disc_normal_unit.cross(global_up_vector))
        
        old_speed =  vel.magnitude()
        new_speed =old_speed * math.cos(vel_unit.angle(vel))
        
        vel = vel_unit.multiply(new_speed)  
      else:                                               #skip
        print ("skip   on frame: " + str(loop_count)+ "    speed in: "+str(vel.magnitude()))
        
        
        collision_normal = collision_vector.cross (disc_normal_unit)
      
        collision_normal = collision_normal.divide(collision_normal.magnitude())
        rps += vel.z/5 
      
        #striking force on edge     force = mass * acceleration
        force = mass * vel.z * 0.8  #assuming 1 second "bounce" time    %%%%probably needs work%%%
        ####friction stuff         force will be negative
        new_speed =vel.magnitude() + 0.4 * force / mass
        vel = vel_unit.multiply (new_speed)
        #torque on edge
        pitching_moment = math.cos(disc_normal_unit.angle(global_up_vector)) * disc.radius * force
        #Wp is rads/sec of off axis torque
        Wp = pitching_moment / (Im * rps)
        print ("angle change of :" + str(np.rad2deg(Wp))+"   while spinning at rps: " +str (rps))
        orientation_change = collision_normal.multiply (-math.tan(Wp ))
    
        disc_normal_unit= disc_normal_unit.sum(orientation_change)  
        disc_normal_unit= disc_normal_unit.divide(disc_normal_unit.magnitude())
        ## amount of bounce in the skip
        vel.z *= -0.8

      coords.z = 0.0
    else :
      #flag to not do rolling effects
      rolling_moment = 0
    
    
    
    
    """
    #####Unit Vectors
    #vel_unit:  unit vector of total disc velocity
    #lift_unit: unit vector 90 degrees from vel_unit in line with discs 'up' vector
    
    #disc_normal_unit:  unit vector of disc which points 'up' relative to the disc
    #disc_unit_y:  unit vector of velocity vector projected onto discs planes (angle between this and vel_unit is angle of attack)
    #disc_unit_x:  unit vector of disc which points perpendicular to direction of travel but lays on the discs plane
    
    #rotation direction goes from disc_unit_x to disc_unit_y
    """
    
    #####Unit Vectors for orientation of forces
    disc_normal_unit= make_unit_vector (disc_normal_unit)
    #unit vector of velocity
    vel_unit = make_unit_vector (vel)
    if not disc_normal_unit.parallel (vel_unit):
      #unit vector perp to vel and orient
      disc_unit_x = vel_unit.cross(disc_normal_unit)
      disc_unit_x = make_unit_vector (disc_unit_x)
      disc_unit_y = disc_unit_x.cross(disc_normal_unit)
      disc_unit_y = make_unit_vector (disc_unit_y)
      #direction that lift force will be applied in
      lift_unit = disc_unit_x.cross(vel_unit)  
    else:
      disc_unit_x = Vector (0,0,0)
      disc_unit_y = Vector (0,0,0)
      lift_unit = Vector (0,0,0)
    
    
    
    
    
    #####Multiuse Variables
    #AoAr angle of attack (radians)   
    aoar = vel_unit.angle(disc_normal_unit)-np.deg2rad(90)
    #velocity squared
    V2 = (vel.magnitude()) ** 2
    #0.5 * pressure * area * velocity^2
    pav2by2 = p * a * V2 / 2
    
    
    
    
    
    #####Coefficients
    #stall detection and lift coefficiant    
    if aoar > -0.52 and aoar < 0.52:
      #normal flight conditions
      coefficient_curve = (0.5*math.sin(6*aoar)+math.sin(2*aoar))
      Cl = disc.Cl0 + disc.Cla * coefficient_curve
      Cm = disc.Cm0 + disc.Cma * coefficient_curve
      Cds = 0
    else:
      #stall conditions
      stall_curve = math.sin(2*aoar)
      Cl = disc.Cl0 + disc.Cla * stall_curve
      Cm = disc.Cm0 + disc.Cma * stall_curve
      #stall induced drag   %%%might need work%%%%
      Cds = -math.cos(2*aoar)+0.55
      print ("Stalling!")
    ##print (str(aoar) + "  " +str(coefficient_curve)+ "   " + str(Cm))  
    ##parasitic drag  |  lift induced drag   |   stall drag    
    Cd =   disc.Cd    +  Cl**2/(3.1415*0.8)  +     Cds
    
    
    
    
    
    #####Forces
    #magnitudes
    #lift_force_magnitude = 0.5 * p * a * V2 * cl
    lift_force_magnitude = pav2by2 * Cl
    drag_force_magnitude = pav2by2 * Cd
    #gravitational_force_magnitude = -9.81 * mass      is a throw constant
    
    #with direction
    lift_force_vector = lift_unit.multiply(lift_force_magnitude)
    drag_force_vector = vel_unit.multiply(-drag_force_magnitude)
    #gravitational_force_vector = global_down_vector.multiply(gravitational_force_magnitude)
    
    total_force_vector = lift_force_vector.sum(drag_force_vector)
    total_force_vector = total_force_vector.sum (gravitational_force_vector)
    total_force_magnitude = total_force_vector.magnitude()
    
    
    
    
    
    ######Moments
    #pitching moment (lift induced torque caused by lifting force coming from off center)
    pitching_moment = pav2by2 * Cm * disc.diameter 
    #Disc State Detection section sets rolling moment to 1 if disc is rolling
    if rolling_moment:
      disc_normalizing_angle = total_force_vector.angle(disc_normal_unit)
      disc_normalized_force_magnitude =-total_force_magnitude * math.cos (disc_normalizing_angle)
      rolling_moment = disc_normalized_force_magnitude * disc.diameter / 2
      #total_force_vector = total_force_vector.sum (drag_force_vector.divide(10))
    else:
      rolling_moment = 0
      
      
      
      
      
    ######Gyroscopic Precession
    #W is rads/sec of off axis rotation
    #Wp is pitch     Wr is roll
    Wp = pitching_moment / (Im * rps)
    Wr = rolling_moment / (Im * rps)
  
    #Pitch and Roll rotations
    #####       %%%%might need work regarding summing of moment pre gyroscopic precession
    orientation_change = disc_unit_x.multiply (math.tan(Wp * time_coefficiant ))
    disc_normal_unit= disc_normal_unit.sum(orientation_change)
    
    orientation_change = disc_unit_y.multiply (-math.tan(Wr * time_coefficiant ))
    disc_normal_unit= disc_normal_unit.sum(orientation_change)

    disc_normal_unit= make_unit_vector (disc_normal_unit)
    
  
  
  
    #####Integration of acceleration and displacement
    total_acceleration_vector = total_force_vector.multiply (seconds_over_mass)
    
    
    #####Collision checking for trees, bushes, and  %%%%future todo%%% baskets
    """
    Collision detection will check for collision with bushes and tree in movement steps of 1/10 the velocity
    -given a max 30m/s release this would mean checking every 0.15m
    This 1/10 application of velocity will be called a substep
    
    1 - check to see what collision_map coordinate the substep will move into
    2 - perform a collision on the predicted movement
    3a- if no collision is detected then move the disc this substep
    3b- if a collision is detected then move the disc to the point of the collision and change the velocity based on the collision
    4 - exit the loop when 10 substep have been completed
    
    some displacement will be lost on every 3b, but this should be uncommon enough to be negligable
    """
    if checking_for_aerial_collisions:
      #vel_substep is one tenth the displacement of a normal simulation step
      vel_substep = vel.multiply(time_coefficiant/10)
      #used as the disc coordinates through the collision substeps  
      coords_substep = coords
      for substep in range (10):
        #move the disc
        coords_substep = coords_substep.sum (vel_substep)
        #check to see if there is a collision in the new location
        #detect_collision will return a list ["collision type as string",var1,var2(sometimes)]
        collision_type = map_handler.detect_collision (coords_substep)
        if collision_type[0] == "none":
          pass
        #collision with a tree trunk
        elif collision_type [0] == "trunk":
          #%%%todo%%% add pitching / rolling moment to tree collision
          trunk_center = collision_type [1]
          trunk_radius = collision_type [2]
          #x and y distance of disc to trunk centers
          x_distance = coords_substep.x - trunk_center.x
          y_distance = coords_substep.y - trunk_center.y
          center_center_distance = math.sqrt((x_distance)**2+(y_distance)**2)
          #center to center distance required for a collision
          collision_distance = trunk_radius + disc.radius
          
          #if there is a collision
          if center_center_distance < collision_distance: 
            print ("hit on frame " + str(loop_count))
            
            #flattening the vectors in order to do the collision in 2 space    %%%todo: look into 3 space solutions%%%
            flat_vel = Vector (vel.x,vel.y,0)
            trunk_normal = Vector (x_distance,y_distance,0)
            trunk_normal = make_unit_vector (trunk_normal)
            
            #angle between velocity vector and the trunks normal at the collision point
            collision_angle = flat_vel.angle(trunk_normal)
            
            #speed relative to the trunks normal
            collision_vector = vel.multiply( math.cos (collision_angle))
            collision_speed = collision_vector.magnitude ()
            bounce_acceleration = trunk_normal.multiply(collision_speed*2)
            
            vel = vel.sum (bounce_acceleration)
            #speed absorbed by bouncing   %%%todo: maximum force of bounce, based on rps maybe%%%^
            vel = vel.multiply (0.8)
            #recalculate the velocity substep after bounce
            vel_substep = vel.multiply(time_coefficiant/10)
            
            #Place disc outside of trunk
            correction_factor = collision_distance / center_center_distance
            coords_substep = Vector (trunk_center.x+correction_factor*x_distance,trunk_center.y+correction_factor*y_distance,coords_substep.z)
          else:
            print ("no hit")
        elif collision_type [0] == "foliage":
          density = collision_type [1]
        elif collision_type [0] == "branchy foliage":
          density = collision_type [1]
        elif collision_type [0] == "basket":
          basket_center = collision_type [1]
      coords = coords_substep
    else:
      coords = coords.sum(vel.multiply(time_coefficiant))
      #coords = coords.sum(total_acceleration_vector.divide(2))
      #^above is more accurate, but for the sake of simplicity in collision detection it is abstracted out
    #accelerate the disc
    vel = vel.sum (total_acceleration_vector)
    
    
    
      
    #rotational drag %%%%%%%needs work%%%%%%
    rps -= rps*0.02 * time_coefficiant  
    #trajectory creation     
    disc_flight.record_state(coords,disc_unit_x)
    
    loop_count += 1  
  else:
    if loop_count == 400:
      print ("Simulation did not finish properly")
    disc_flight.finish_flight ()  
  return (disc_flight)
  
  

  






####Throw Constants
coordinates = Vector(65,0,1)
ground_roll = Vector(65,0,0)

flat = 0
hyzer = 15
hyzer_30 = 30
hyzer_45 = 45
hyzer_spike = 75
anhyzer = -15
anhyzer_45 = -45
tomahawk = -90


########################################################################Throw Done Here#######################################################

fig = plt.figure(figsize=(20, 10))
ax = Axes3D(fig)#fig.add_subplot(111, projection='3d')
map_handler.load_map_from_filepath ("test_map.bmp")

#nominal speeds and angles for disc flights!!!
#      direction, loft  ,  speed
#Drivers   0    ,  9    ,   25
#Mids      0    ,  11   ,   20
#Putter    0    ,  13   ,   15


#big air shot
test_throw = do_throw (coordinates,vector_from_direction_loft_speed (1.5,9,28),0,disc_Driver_S)
#skip shot
#test_throw = do_throw (coordinates,vector_from_direction_loft_speed (0,3,25),15,disc_Driver_OS)
#test_throw = do_throw (coordinates,vector_from_direction_loft_speed (0,5,25),15,disc_Driver_OS)
#big roller
#test_throw = do_throw (coordinates,vector_from_direction_loft_speed (15,5,28),6,disc_Driver_US)
#tomahawk
#test_throw = do_throw (coordinates,vector_from_direction_loft_speed (0,45,30),-90,disc_Driver_US)

#animate plot in realtime using dt and FuncAnimation
animate_plot = 1



#####Display Stuff

# build arrays of coords
display_x = test_throw.pos_x + test_throw.orient_x
display_y = test_throw.pos_y + test_throw.orient_y
display_z = test_throw.pos_z + test_throw.orient_z
display_x2 = test_throw.pos_x - test_throw.orient_x
display_y2 = test_throw.pos_y - test_throw.orient_y
display_z2 = test_throw.pos_z - test_throw.orient_z
pos_x = test_throw.pos_x
pos_y = test_throw.pos_y
pos_z = test_throw.pos_z
array_length = pos_x.size



_k = array_length-1

ax.set_xlabel('Distance in Meters')
ax.set_ylabel('Distance in Meters')
ax.set_zlabel('Height in Meters')
ax.set_aspect('equal', 'box')
ax.set_xlim(0,175)
ax.set_ylim(0,175)
ax.set_zlim(0,175)


###map display stuff
greenx = []
greeny = []
greenz = []
brownx = []
browny = []
brownz = []

for x in range(map_handler.collision_map.shape[0]):
    for y in range(map_handler.collision_map.shape[1]):
      for z in range(map_handler.collision_map.shape[2]):
        if map_handler.collision_map [x][y][z] in range (71,80):
          greenx = np.append (greenx,x)
          greeny = np.append (greeny,y)
          greenz = np.append (greenz,z)
        elif map_handler.collision_map [x][y][z] in range (81,90):
          greenx = np.append (greenx,x)
          greeny = np.append (greeny,y)
          greenz = np.append (greenz,z)
        elif map_handler.collision_map [x][y][z] in range (91,100):
          brownx = np.append (brownx,x)
          browny = np.append (browny,y)
          brownz = np.append (brownz,z)

ax.scatter(greenx, greeny, greenz, c='g', marker='o')
ax.scatter(brownx, browny, brownz, c='brown', marker='o')

#normal ol' boring plot
if not animate_plot:
  #plot coords and shadow

  ax.plot(np.flipud(pos_x), np.flipud(pos_y), np.flipud(pos_z), linewidth=3, linestyle='-', color=[0.0, 0.8, 0.0], alpha=1, \
        marker='o', markersize=15, markerfacecolor=[0.8, 0.8, 0.0], markevery=array_length+1)
  ax.plot(np.flipud(pos_x), np.flipud(pos_y), np.flipud(pos_z)*0, linewidth=4, linestyle='-', color=[0.4, 0.4, 0.4], alpha=0.4, \
        marker='o', markersize=10, markerfacecolor=[0.8, 0.8, 0.8], markeredgewidth=0.0, markevery=array_length+1)
        
  ax.plot(np.flipud(display_x), np.flipud(display_y), np.flipud(display_z), linewidth=2, linestyle='-', color=[0.0, 0.6, 0.0], alpha=1, \
        markerfacecolor=[0.8, 0.8, 0.0], markevery=array_length+1)  
  ax.plot(np.flipud(display_x2), np.flipud(display_y2), np.flipud(display_z2), linewidth=2, linestyle='-', color=[0.0, 0.6, 0.0], alpha=1, \
        markerfacecolor=[0.8, 0.8, 0.0], markevery=array_length+1)  
        
  plt.show()

else:
  #use FuncAnimation and a handler function to update the datasets and animate in real time
  #animation handler function
  def update_disc_pose(frameNum, signals, px, py, pz, px2, py2, pz2, px3, py3, pz3, al):

    global _k

    #plot the flipped arrays backwards so we can have the disc 'marker' on the right side

    #disc
    signals[0][0].set_data(px[_k:al-1], py[_k:al-1])
    signals[0][0].set_3d_properties(pz[_k:al-1])
    #signals[0][0].set_markerevery(_k-1) #always show 'disc' marker at first and last point
    
    #disc_orientation
    signals[2][0].set_data(px2[_k:al-1], py2[_k:al-1])
    signals[2][0].set_3d_properties(pz2[_k:al-1])
    signals[3][0].set_data(px3[_k:al-1], py3[_k:al-1])
    signals[3][0].set_3d_properties(pz3[_k:al-1])
    #signals[0][0].set_markerevery(_k-1) #always show 'disc' marker at first and last point

    #shadow
    signals[1][0].set_data(px[_k:al-1], py[_k:al-1])
    signals[1][0].set_3d_properties(pz[_k:al-1]*0)

    #increment count
    _k = _k-1

    if(_k <= 0 ):
      _k = 0


  signals = []

  signals.append(ax.plot([], [], [], linewidth=3, linestyle='-', color=[0.0, 1, 0.0], alpha=1, markevery=array_length+1))
  signals.append(ax.plot([], [], [], linewidth=5, linestyle='-', color=[0.4, 0.4, 0.4], alpha=0.5, markevery=array_length+1))

  signals.append(ax.plot([], [], [], linewidth=2, linestyle='-', color=[0.0, 0.6, 0.0], alpha=1, markevery=array_length+1))
  signals.append(ax.plot([], [], [], linewidth=2, linestyle='-', color=[0.0, 0.6, 0.0], alpha=1, markevery=array_length+1))
  #play animation real-time (use delta_t)
  anim = animation.FuncAnimation(fig, update_disc_pose, \
    fargs=([signals, np.flipud(pos_x), np.flipud(pos_y), np.flipud(pos_z), np.flipud(display_x), np.flipud(display_y), \
        np.flipud(display_z), np.flipud(display_x2), np.flipud(display_y2), np.flipud(display_z2), array_length]), interval=(delta_t))

  plt.show()


















