from PIL import Image, ImageChops 
import os
parent_path=os.getcwd()
for path, subdirs, files in os.walk(os.getcwd()):
    for filename in files:
        if '.jpeg' not in filename: continue
        # Opening the test image, and saving it's object 
        img = Image.open(filename) 
          
        # Passing the image object to invert()   
        #inv_img = ImageChops.invert(img) 
        #inv_img.save(filename) 
        img.save(new_filename)
