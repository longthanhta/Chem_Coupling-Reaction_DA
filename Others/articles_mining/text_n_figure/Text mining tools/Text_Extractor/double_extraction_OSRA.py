import subprocess
import cv2

file = "0329test.png"
# Define command as string and then split() into list format
cmd = 'osra -c '+file
# Use shell to execute the command
sp = subprocess.Popen(cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True)

# Separate the output and error.
# This is similar to Tuple where we store two values to two different variables
out,err=sp.communicate()

# Store the return code in rc variable
rc=sp.wait()

print('Return Code:',rc,'\n')
print('output is: \n', out)
print('error is: \n', err)  

print("--------------Finished First SMILES and coordinate Extraction------------------")
text_file = open("SMILES&corr_record.txt", "w")
text_file.write(out)
text_file.close()
text_file = open("SMILES&corr_record.txt", "r")
Lines = text_file.readlines()
print("Number of SMILES detected:"+ str(len(Lines)))
#print(Lines)
print("--------------Restoring coordinate record----------------------------------------")
box=[]
for x in Lines:
    empty_position=x.find(" ")
    #print(x[empty_position+1:])
    box.append(x[empty_position+1:])
print(box)
XY1_corr = []
XY2_corr = []
for i in box:
    splitposition=i.find("-")
    XY1_corr.append(i[:splitposition])
    XY2_corr.append(i[splitposition+1:-1])   
print(XY1_corr)
print(XY2_corr)
X1_corr = []
X2_corr = []
Y1_corr = []
Y2_corr = []
for i in XY1_corr:
    xposition=i.find("x")
    X1_corr.append(i[:xposition])
    Y1_corr.append(i[xposition+1:])   
for i in XY2_corr:
    xposition=i.find("x")
    X2_corr.append(i[:xposition])
    Y2_corr.append(i[xposition+1:])   
print(X1_corr)
print(X2_corr)
print(Y1_corr)
print(Y2_corr)
print("--------------Rescalling coordinate by 5%----------------------------------------")
Scalling_Factor=0.05
X1_ncorr = []
X2_ncorr = []
Y1_ncorr = []
Y2_ncorr = []
for i in X1_corr:
    i = int(int(i)-int(i)*(Scalling_Factor))
    X1_ncorr.append(i)
for i in Y1_corr:
    i = int(int(i)-int(i)*(Scalling_Factor))
    Y1_ncorr.append(i)
for i in X2_corr:
    i = int(int(i)+int(i)*(Scalling_Factor))
    X2_ncorr.append(i)
for i in Y2_corr:
    i = int(int(i)+int(i)*(Scalling_Factor))
    Y2_ncorr.append(i)
print(X1_ncorr)
print(X2_ncorr)
print(Y1_ncorr)
print(Y2_ncorr)
print("--------------Finished Rescalling coordinate by 5%----------------------------------------")
print("--------------Croping individual Structure----------------------------------------")
img = cv2.imread(file)
COUNT=0
for i in X1_ncorr:
    x1=X1_ncorr[COUNT]
    x2=X2_ncorr[COUNT]
    y1=Y1_ncorr[COUNT]
    y2=Y2_ncorr[COUNT]
    Y =y1
    h = y2-y1
    X=x1
    w=x2-x1
    crop_img = img[Y:Y+h, X:X+w]
    cv2.imwrite('crop_img'+str(COUNT)+".jpg", crop_img)
    COUNT+=1
print("--------------Finished Croping individual Structure----------------------------------------")
print("--------------Individual SMILES Extraction----------------------------------------")
COUNT=0

for i in X1_ncorr:
    individual_picture = "crop_img"+str(COUNT)+".jpg"
    cmd = 'osra '+individual_picture
    # Use shell to execute the command
    sp = subprocess.Popen(cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True)

# Separate the output and error.
# This is similar to Tuple where we store two values to two different variables
    out,err=sp.communicate()

# Store the return code in rc variable
    rc=sp.wait()


#print('Return Code:',rc,'\n')
    print(out)
    text_file2 = open("SMILES_second_extraction.txt", "a")
    text_file2.write(out)
    text_file2.close()
#print('error is: \n', err) 
    COUNT+=1

