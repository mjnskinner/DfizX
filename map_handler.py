from PIL import Image 
import numpy as np
map_image = Image.open ("test.bmp")
image_array = np.array(map_image)

map_image.show()
print (image_array[1][0])
image_array = np.flip(image_array,0)
print (image_array[2][0])