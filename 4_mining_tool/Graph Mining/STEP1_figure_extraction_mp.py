from tqdm import tqdm
import fitz # PyMuPDF
import io
from PIL import Image
import os
import cv2
import numpy as np
from multiprocessing import Process, Manager
from multiprocessing import Pool
from functools import partial
import subprocess
from os import listdir
from os.path import isfile, join
from tqdm import tqdm


def check_and_inv(file_path):

	# Reading an image in default mode:
	inputImage = cv2.imread(file_path)
	
	# Convert RGB to grayscale:
	originalGrayscale = cv2.cvtColor(inputImage, cv2.COLOR_BGR2GRAY)
	
	# Equalize histogram
	grayscaleImage = cv2.equalizeHist(originalGrayscale)
	
	# It might be interesting to you to check out the image equalization:
	cv2.waitKey(0)
	
	# Binarize the image with a fixed threshold:
	minThresh = 128
	_, binaryImage = cv2.threshold(grayscaleImage, minThresh, 255, cv2.THRESH_BINARY)
	
	# Compute the percent of white pixels:
	(imageHeight, imageWidth) = binaryImage .shape[:2]
	whitePercent = cv2.countNonZero(binaryImage)/(imageHeight * imageWidth)

	if whitePercent < 0.2:
		#print("Correcting images...")
	
		# Correct the equalized image:
		grayscaleImage = 255 - grayscaleImage
		cv2.imwrite(file_path, grayscaleImage )


def image_extraction(database_path,result_folder,output_list,file_list,index):
	# get the file path
	file_name=file_list[index]
	file_path=database_path+'/'+file_name
	# make a folder to store result
	article_DOI=file_name.replace('.pdf','')
	with fitz.open(file_path) as my_pdf_file:
		#loop through every page
		count=0
		for page_number in range(0, len(my_pdf_file)):
	
			# acess individual page
			page = my_pdf_file[page_number]
	
	
			# check if images are there
			#$if images:
			#	print(f"There are {len(images)} image/s on page number {page_number}[+]")
			#else:
			#	print(f"There are No image/s on page number {page_number}[!]")
	
			# loop through all images present in the page 
			for image_number, image in enumerate(page.get_images(), start=0):
	
				#access image xerf
				xref_value = image[0]
				
				#extract image information
				base_image = my_pdf_file.extract_image(xref_value)
	
				# access the image itself
				image_bytes = base_image["image"]
	
				#get image extension
				ext = base_image["ext"]
	
				#load image
				image = Image.open(io.BytesIO(image_bytes))
	
				#save image locally
				count+=1
				out_im_path=f"{result_folder}/{article_DOI}_{count}.{ext}"
				image.save(open(out_im_path, "wb"))
				check_and_inv(out_im_path)
	output_list.append('.')
	print(len(output_list))
# file path you want to extract images from
#file_name = "0-0040403996013020.pdf"
if __name__=='__main__':
	
	

	#cwd = os.getcwd()
	
	#database_path='D:/DATACHEM/co_activation_database'
	#database_path='D:/DATACHEM/picture_extraction/broken pdf'
	
	database_path='PDF_folders'
	file_list = [f for f in listdir(database_path) if isfile(join(database_path, f))]
	

	#out_path='pictures from pdf'
	#out_path='broken pdf figure extraction'
	
	out_path='extracted_figs'

	try:
		os.mkdir(out_path)
	except:
		print('result folder already exist, carefull of overwriting')
	#file_name = "acs.organomet.5b00874.pdf"
	
	# this part for single-processing
	#for i in tqdm(range(len(file_list))):
	#for index,file_name in enumerate(file_list):
		#print(index,len(file_list))
		#file_name=file_list[i]
		#image_extraction(file_name,database_path,out_path)
		#print(file_name)

	# this part for multi-processing
	output_list = Manager().list()
	input_list=range(0,len(file_list))
	print(input_list)
	pool = Pool(processes=4)
	func = partial(image_extraction,database_path,out_path,output_list,file_list)
	pool.map(func, input_list)
	pool.close()
	pool.join()
	#

