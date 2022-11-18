from tqdm import tqdm
import fitz # PyMuPDF
import io
from PIL import Image
from IPython.display import Image as ImageShow
import os
import cv2
import numpy as np
from os import listdir
from os.path import isfile, join
from matplotlib import pyplot as plt
from scipy.spatial import distance
import pandas as pd
plt.rcParams['figure.dpi'] = 300 # to increase the solution in jupyter notebook


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

def line_extraction(picture_path,output_folder):
	#try:
		#output_list.append('.')
		#print(len(output_list),'|', end =" ")
		#picture_path=file_list[index]
	#print('picture_path',picture_path)
	file_name=picture_path.split('/')[-1]
	#print(file_name)
	# reduce maxLineGap_set untill you get one arrow line
	def get_arrow(img):
		height, width = img.shape
		arrow_count=0 # intital the value
		maxLineGap_set=0 # inital the value
		#intialize value
		x1a=0
		x2a=0
		y1a=0
		y2a=0

		while arrow_count!=1 and maxLineGap_set<height:
			lines = cv2.HoughLinesP(img, 2, np.pi/180,
									5, #
									minLineLength=width*0.2, # minimum length
									maxLineGap=maxLineGap_set)
			arrow_count=0
			min_pl_arrow=0.2 # minimum percentage length of arrow
			max_pl_arrow=0.6 # maximum percentage length of arrow
			pl_value_list=[] #pl : percentage of length compare to width
			if lines is None:
				return 0,0,0,0,0
			for index,line in enumerate(lines):
				max_length=0
				for x1, y1, x2, y2 in line:
					point1=np.array((x1,y1))
					point2=np.array((x2,y2))
					pl=np.linalg.norm(point1 - point2)/width
					if y1!=y2: continue
					pl=round(pl,2)
					if pl in pl_value_list:
						continue
					if min_pl_arrow<pl<max_pl_arrow:
						arrow_count+=1
						x1a=x1
						x2a=x2
						y1a=y1
						y2a=y2
						pl_value_list.append(pl)
			#print('mlg',maxLineGap_set,'#a',arrow_count,'|', end =" ")
			maxLineGap_set+=height*0.1 # increase maxLineGap_set to make openCV more sensitive
		return x1a,x2a,y1a,y2a,arrow_count

	# pre_processing_the_figure

	file_path=picture_path
	img = cv2.imread(file_path, cv2.IMREAD_GRAYSCALE)
	#print('file_path',	file_path)
	# Initialize output
	out = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
	# Median blurring to get rid of the noise; invert image
	img = 255 - cv2.medianBlur(img, 3)
	cv2.imwrite('out_temp.png', out)
	height, width = img.shape



	# Get the list of line coordinate'

	x1a,x2a,y1a,y2a,arrow_count=get_arrow(img)
	if arrow_count==0 or arrow_count>1: return False # no arrow -> end function here
	cv2.line(out, (x1a, y1a), (x2a, y2a), (0, 0, 255), 2)
	# start to crop image here

	img_total=out
	print('y1_arrow',y1a,'height',height)
	img_rxn_scheme=img_total[0:y1a*2, 0:] # the mainscheme
	img_rxn_substrate=img_total[y1a*2:height, 0:] # the one with optimized substrate
	img_rxn_cond_top=img_total[0:y1a, x1a:x2a]
	img_rxn_cond_bot=img_total[y1a:y1a*2, x1a:x2a]
	img_rxn_cond_all=img_total[0:y1a*2, x1a:x2a]
	img_rxn_reactant=img_total[0:y1a*2, 0:x1a]
	img_rxn_product=img_total[0:y1a*2, x2a:width]


	# output the file

	file_name_wo_png=file_name.replace('.png','')

	# make output folder
	output_path=output_folder+'/'+file_name_wo_png
	try:
		os.mkdir(output_path)
	except:
		pass

	# generating outputpath
	op7=output_path+'/cond.png'
	op6=output_path+'/cond_bot.png'
	op5=output_path+'/cond_top.png'
	op4=output_path+'/scheme.png'
	op3=output_path+'/substrate.png'
	op2=output_path+'/reactant.png'
	op1=output_path+'/product.png'
	op0=output_path+'/0riginal.png'

	# write output file
	try:

		cv2.imwrite(op0,img)
		check_and_inv(op0)

		cv2.imwrite(op4, img_rxn_scheme)
		check_and_inv(op4)


		cv2.imwrite(op5, img_rxn_cond_top)
		check_and_inv(op5)

		cv2.imwrite(op6, img_rxn_cond_bot)
		check_and_inv(op6)

		cv2.imwrite(op7, img_rxn_cond_all)
		check_and_inv(op7)


		cv2.imwrite(op2, img_rxn_reactant)
		check_and_inv(op2)

		cv2.imwrite(op1, img_rxn_product)
		check_and_inv(op1)

		cv2.imwrite(op3, img_rxn_substrate) # it may not have substrate
		check_and_inv(op3)

	except:
		pass
		#except:
		#	pass

def do_task(output_folder,count_list,path_list,index):

	count_list.append('.')
	print(len(count_list),'_____________________________________')
	picture_path=path_list[index]
	line_extraction(picture_path,output_folder)
	# this part for multi-processing

	#

if __name__ == '__main__':
	output_folder='segmented_figs'
	try:
		os.mkdir(output_folder)
	except:
		pass
	database_path='extracted_figs'
	file_list = [f for f in listdir(database_path) if isfile(join(database_path, f))]
	path_list = [database_path+'/'+file_name for file_name in file_list]

	count_list = Manager().list()
	input_list=range(0,len(path_list))
	print(input_list)
	pool = Pool(processes=4)
	func = partial(do_task,output_folder,count_list,path_list)
	pool.map(func, input_list)
	pool.close()
	pool.join()
