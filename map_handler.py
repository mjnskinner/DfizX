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

collision_map = []


#####Image filepath to image_array
def import_image (image_filepath):   
  map_image = Image.open (image_filepath)
  image_array = np.array(map_image)
  map_image.show()
  image_array = np.flip(image_array,0)
  return image_array
  
image_array = import_image ("test.bmp")  


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
    return bush
  elif rgb [0] == 0 and rgb [1] == 255 and rgb [2] == 0:
    return grass
  else:
    return 0

#takes an image array and converts it to an object map    
def make_object_map (image_array):
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
def object_to_collision_mapping (x,y):
  pass

def make_collision_map (object_map):
  #collision map is 5m wider to make room for tree growth without error handling
  global collision_map
  collision_map = np.zeros ([object_map.shape[0]+10,object_map.shape[1]+10,25])
  
  for x in range(object_map.shape[0]):
    for y in range(object_map.shape[1]):
      object_to_collision_mapping (x+5,y+5,object_map [x][y])
      #print (rgb_to_object_map (object_map [x][y]))