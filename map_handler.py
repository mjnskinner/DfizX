from PIL import Image 
import numpy as np


#####Overview
"""
Map handler  
imports pictures and turn them into map data
map data has a few stages

>>image_array
2d array - has rgb values which are mapped into objects in
>>object_map
2d arraay - has int values which are mapped to different object types (grass type, bushes, trees)
this is how maps are stored, but not used. To be used they must first be converted into
>>collision_map
3d array - has int values which are mapped to different collision types
each value represents a 1x1x1 meter cube of airspace
z at 0 represents ground skippiness
"""    
####Main array types              contains                                     notes
#image_array   [x][y][rgb]      |rgb int values
#object_map    [x][y]           |object_id int values
#collision_map [x][y][z][value] |collision_id int values              z is max 25 for memory considerations


####Main Routines
#
#       import_image (image_filepath)
# take a filepath of an image and returns an image_array
#
#       make_object_map (image_array)
# takes image_array and returns an object_map filled with object identities
#
#       object_identities
# a list of object identities and their values
#
#       make_collision_map (object_map)
# takes object_map and creates a global collision_map from it, filled with collision_id value 


#routines to make    %%%todo%%%%
#load_map_from_filepath ()
#filepath > loaded collision map
#
#create and maintain an object_map database?
#
# load_object_map ()    | load a saved object map from file
# save_object_map ()    | save an object map to file
# 



#####Image filepath to image_array
def import_image (image_filepath):   
  map_image = Image.open (image_filepath)
  image_array = np.array(map_image)
  image_array = np.rot90 (image_array,-1)
  #map_image.show()
  #image_array = np.flip(image_array,0)
  return image_array
  



#####object_identities
#all trees and 
fairway     = 1
rough       = 2
deep_rough  = 3
packed_dirt = 4
pavement    = 5
water       = 6
sand        = 7
small_bush    = 51
medium_bush   = 52
large_bush    = 53
tall_bush     = 54   #can make hedges
baby_tree         = 100
skinny_tree       = 101
tall_skinny_tree  = 102
medium_tree       = 103
tall_tree         = 104
thick_tree        = 105
tall_thick_tree   = 106
giant_tree        = 107






######Image to object map
#object map is a flat 2d map with integers representing different objects


#returns an object identity based on a pixels rgb values
def rgb_to_object_map (rgb):
  if rgb [0] == 0 and rgb [1] == 100 and rgb [2] == 0:
    return small_bush
  elif rgb [0] == 0 and rgb [1] == 110 and rgb [2] == 0:
    return tall_bush
  elif rgb [0] == 0 and rgb [1] == 255 and rgb [2] == 0:
    return grass
  elif rgb [0] == 100 and rgb [1] == 50 and rgb [2] == 0:
    return medium_tree
  else:
    return 0

#takes an image array and converts it to an object map    
def make_object_map (image_array):
  global object_map
  object_map = []
  object_map = np.zeros ([image_array.shape[0],image_array.shape[1]])
  for x in range(image_array.shape[0]):
    for y in range(image_array.shape[1]):
      object_map [x][y] = rgb_to_object_map (image_array [x][y])
      #print (rgb_to_object_map (image_array [x][y]))
  return object_map



 
######object map to collision map

#collision map is a 3d array of ints that represent certain shapes and sizes of collidable objects and ground
#z = 0  represents ground data 
#z above zero represents collision data from z-1 meters to z meters

#when z > 0
# int value of collision_map [x][y][z]  | defines
#                                     0 | nothing
                              #  91-100 | diameter of tree trunk in 0.1m increments
                              #  81-90  | density (hit chance) of dense branchy foliage
                              #  71-80  | desntiy (hit chance) of loose leafy foliage
                              
def add_collision_type (x,y,z,add_collision_type):
  if add_collision_type > collision_map[x][y][z]:
    collision_map[x][y][z] = add_collision_type
    
def object_to_collision_mapping (x,y,object_type):
  if object_type == small_bush:
    add_collision_type (x,y,1,89)
  elif object_type == tall_bush:
    add_collision_type (x,y,1,89)
    add_collision_type (x,y,2,88)
    add_collision_type (x,y,3,87)
    add_collision_type (x,y,4,86)
  elif object_type == medium_tree:
    add_collision_type (x,y,1,95)
    add_collision_type (x,y,2,95)
    add_collision_type (x,y,3,95)
    add_collision_type (x,y,4,95)
    add_collision_type (x,y,5,95)
    add_collision_type (x,y,6,90)
    add_collision_type (x-1,y,6,85)
    add_collision_type (x+1,y,6,85)
    add_collision_type (x,y-1,6,85)
    add_collision_type (x,y+1,6,85)
    add_collision_type (x,y,7,85)
    add_collision_type (x-1,y-1,6,75)
    add_collision_type (x-1,y+1,6,75)
    add_collision_type (x+1,y+1,6,75)
    add_collision_type (x+1,y-1,6,75)
    add_collision_type (x-1,y,7,73)
    add_collision_type (x+1,y,7,73)
    add_collision_type (x,y-1,7,73)
    add_collision_type (x,y+1,7,73)
    
def make_collision_map (object_map):
  #collision map is 5m wider to make room for tree growth without error handling
  global collision_map
  collision_map = []
  collision_map = np.zeros ([object_map.shape[0]+10,object_map.shape[1]+10,25])
  
  for x in range(object_map.shape[0]):
    for y in range(object_map.shape[1]):
      object_to_collision_mapping (x+5,y+5,object_map [x][y])
      #print (rgb_to_object_map (object_map [x][y]))
      
def load_map_from_filepath (filepath):
  make_object_map (import_image (filepath))
  make_collision_map (object_map)
  
#load_map_from_filepath ("test_map.bmp")