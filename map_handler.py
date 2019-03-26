from PIL import Image 
import numpy as np





    
map_image = Image.open ("test.bmp")
image_array = np.array(map_image)

map_image.show()

image_array = np.flip(image_array,0)


#####object identities
bush = 1
grass = 2






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
def object_to_collision_mapping (x,y):
  pass

def make_collision_map (image_map):
  collision_map = np.zeros ([image_array.shape[0]+6,image_array.shape[1]+6,25])
  for x in range(image_array.shape[0]):
    for y in range(image_array.shape[1]):
      create (image_array [x][y])
      #print (rgb_to_object_map (image_array [x][y]))