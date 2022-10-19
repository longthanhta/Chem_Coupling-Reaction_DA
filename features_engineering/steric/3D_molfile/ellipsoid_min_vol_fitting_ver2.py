import pandas as pd
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdchem
from random import randrange
import numpy as np
from numpy import random, cos, sin, sqrt, pi
import subprocess
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA 
from numpy import linalg

#use conda install -c conda-forge openbabel for installing openbabel on Mac

def main():
    result_list_list=[[],[],[]]
    # conditions can be val_e, elemt, label, atommap, dist_to_rc
    # all in AND relationship e.g. {'elmt':'C','atommap':888} = both a C and has amap 888
    rc_cond = {'elmt':'C','atommap':888}
    SMILES_list = [
                    '[C:888][C:4]#[N:5]',
                    ]
    for i in range(len(SMILES_list)):
        SMILES = SMILES_list[i]
        print(i,SMILES)
        if i > 100:
            continue
        '''r1,r2,r3 = ellip_info(SMILES,rc_cond,max_scan_dist=100,plot_graph=True,mode="radii")
        result_list_list[0].append(r1)
        result_list_list[1].append(r2)
        result_list_list[2].append(r3)'''
        # round to closest 0.2
        (axes_zylong,axes_zyshort), deviation_displ = ellip_info(SMILES,rc_cond,max_scan_dist=100,plot_graph=False,mode="2D_ellipse")
        result_list_list[0].append(axes_zylong)
        result_list_list[1].append(axes_zyshort)
        result_list_list[2].append(deviation_displ)

    for z in range(len(result_list_list)):
        print(z)
        for value_ in result_list_list[z]:
            print(value_)
    return


def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def set_Xaxis(coord_df):

    # set rc coordinate as 0, 0, 0
    rc_coord_copy = coord_df.loc[(coord_df['dist_to_rc'] == 0)]['coord'].values[0]
    for i, row in coord_df.iterrows():
        new_value = row['coord'] - rc_coord_copy
        coord_df.at[i,'coord'] = new_value

    
    neigh_coord_vector = coord_df.loc[(coord_df['dist_to_rc'] == 1)].head(1)['coord'].values[0]
    neigh_coord_scalar = sum([i**2 for i in neigh_coord_vector])**(1/2)

    

    # find the proper rotational matrix
    # if the axis is already x-axis than dont do rotation
    if list(neigh_coord_vector[1:]) != [0. , 0.]:
        rotational_matrix = rotation_matrix_from_vectors(neigh_coord_vector,np.array([neigh_coord_scalar,0,0]))
    else:
        rotational_matrix = [   [1.0, 0.0, 0.0], 
                                [0.0, 1.0, 0.0], 
                                [0.0, 0.0, 1.0],  ]

    # set rc-neigh to x-axis and change other vectors
    for i, row in coord_df.iterrows():
        new_value = np.matmul(rotational_matrix,row['coord'])
        coord_df.at[i,'coord'] = new_value

    return coord_df


def ellip_info(SMILES,
                rc_cond,
                file_name='test',
                max_scan_dist = 10000,
                print_detail=False,
                plot_graph=False,
                random_count=200,
                rc_atommap = 888,
                extra_atommap = 999,
                radical_rc = False,
                radical_atom = 'C',
                mode = "radii"
                ):

    if radical_rc:
        SMILES = SMILES.replace(str(rc_atommap)+']',str(rc_atommap)+']([{0}:{1}])'.format(radical_atom,str(extra_atommap)))
    
    # generate mol file
    generate_molfile(SMILES,file_name)

    
    atom_str_2d_list, _ = parse_molfile('test_pre3d.mol')
    atommap_list = [int(atom_str[-10:-7]) for atom_str in atom_str_2d_list]


    with open(file_name+'.mol', "r") as myfile:
        molfile_str=myfile.readlines()

    # First three line doesn't have information
    line_3_str=molfile_str[3]
    no_of_atom=int(line_3_str[0:3])
    no_of_bond=int(line_3_str[3:6])
    
    up_limit=4+no_of_atom

    # Create data table from mol file
    coord_dict={'label':[],'elmt':[],'coord':[],'neigh':[],'val_e':[],'atommap':[]}
    for indx,line_str in enumerate(molfile_str):

        # skip first three line
        if indx <4: continue

        # skip last line
        if indx >= up_limit: continue
        label=indx-3
        line_inf=line_str.split(' ')
        for entry in line_inf:
            if len(entry)==0:
                line_inf.remove(entry)

        # Get coordinate
        x_c=line_inf[0]
        y_c=line_inf[1]
        z_c=line_inf[2]
        ele=line_inf[3]
        if indx-4 < len(atommap_list):
            atommap = atommap_list[indx-4]
        else:
            atommap = 0
        # and other information
        coord_dict['label'].append(label)
        coord_dict['coord'].append(np.array([x_c,y_c,z_c],dtype=float))
        coord_dict['elmt'].append(ele)
        coord_dict['atommap'].append(atommap)
        val_e=0
        neigh_lst=[]

        # now loop in bond data (the later half in mol file)
        for indx2,line_str in enumerate(molfile_str):
            if indx2 <up_limit: continue
            if indx2 >=up_limit+no_of_bond: continue
            
            
            begin_a=line_str[0:3]
            end_a=line_str[3:6]
            bond_type=line_str[6:9]

            #Get neighbour information and numbe of valance electron
            if int(begin_a)==label:
                neigh_lst.append(int(end_a))
                val_e+=int(bond_type)
            if int(end_a)==label:
                neigh_lst.append(int(begin_a))
                val_e+=int(bond_type)
        coord_dict['val_e'].append(val_e)
        coord_dict['neigh'].append(neigh_lst)

    # Show the table, these information can be double-checked with Avogadro
    mol_df=pd.DataFrame(coord_dict)
    #print('mol_df',mol_df)

    # Get reaction center (C with 3 val electron)

    RC_df=mol_df
    for _label, _value in rc_cond.items():
        RC_df=RC_df[(RC_df[_label]==_value)]
    #print('RC_df',RC_df)
    rc_label=int(RC_df['label'].values[0])
    #print('rc_label',rc_label)
    
    # Now calculte bond distance to reaction center
    bond_dist_lst=[]
    for indx3,row in mol_df.iterrows():

        if row.label==rc_label:
            bond_dist_lst.append(0)
            continue
        start_point=row['label']

        scanned_atoms=list(row['neigh'])

        bond_distance=1
        while rc_label not in scanned_atoms:
            bond_distance+=1

            # Warning, use the variable directly, otherwise
            # those two variables just like two link to same list
            # try to delete "list(", result will be different
            scanned_atoms_temp=list(scanned_atoms)
            for label in scanned_atoms_temp:
                neigh_lst=mol_df[(mol_df['label']==label)]['neigh'].values
                scanned_atoms.extend(neigh_lst[0])
                scanned_atoms=list(set(scanned_atoms))
        bond_dist_lst.append(bond_distance)

    # Add bond distance to reaction center
    #print('bond_dist_lst',bond_dist_lst)
    mol_df['dist_to_rc']=bond_dist_lst
    #print('add distance')
    #print(mol_df)


    ########################
    # REMOVE STRUCTURE HERE#
    ########################

    # simple modification: remove atom far from reacion center
    sub_mol_df=mol_df[(mol_df['dist_to_rc']<=max_scan_dist) & (mol_df['atommap']!=extra_atommap)]
    print('after remove\n',sub_mol_df)

    # set one of the rc-neigh bond as x-axis
    sub_mol_df = set_Xaxis(sub_mol_df)
    #print('after remove\n',sub_mol_df)

    # add radius of sphere
    ar_df=pd.read_excel('ar.xlsx')
    radii_list,(axes_zylong,axes_zyshort), deviation_displ  = get_vol(sub_mol_df,ar_df,file_name,print_detail=print_detail,plot_graph=plot_graph,random_count=random_count)
    if print_detail:
        print('3 radius are',radii_list)

    if mode == "radii":
        return radii_list
    elif mode == "2D_ellipse":
        # round to closest 0.2 to make result more consistent
        axes_zylong = round(0.2*round(axes_zylong/0.2),1)
        axes_zyshort = round(0.2*round(axes_zyshort/0.2),1)
        deviation_displ = round(0.2*round(deviation_displ/0.2),1)
        return (axes_zylong,axes_zyshort), deviation_displ


def parse_molfile(filename):
    with open(filename, "r") as myfile:
        molfile_str_list=myfile.readlines()[3:]
    info_str = molfile_str_list[0]

    no_of_atom=int(info_str[0:3])
    no_of_bond=int(info_str[3:6])

    atom_str_list = molfile_str_list[1:1+no_of_atom]
    bondstr_list = [str[:9] for str in molfile_str_list[1+no_of_atom:1+no_of_atom+no_of_bond]]
    
    return (atom_str_list,bondstr_list)

def generate_molfile(sub_smi,file_name):
    sub_mol=Chem.MolFromSmiles(sub_smi)

    # Write 2D mol file
    molfile_str=Chem.MolToMolBlock(sub_mol)
    with open(file_name+'_pre3d.mol', 'w') as mf:
        mf.write(molfile_str)
    command_line=["obabel",\


    # Use obabel to convert it to 3D
    file_name+'_pre3d.mol',\
    '--gen3D','-O',file_name+'.mol','-h']
    print('run command line',' '.join(command_line))
    subprocess.call(command_line)

    _, bondstr_2d_list = parse_molfile(file_name+'_pre3d.mol')
    _, bondstr_3d_list = parse_molfile(file_name+'.mol')

    #check if there is changes in index
    if not all(i in bondstr_3d_list for i in bondstr_2d_list):
        print('There is some error with the conversion, please check the 2d and 3d mol files')

def get_vol(coord_df,ar_df,file_name,print_detail=True,plot_graph=True,random_count=200,skip_rc = False):
    coord_lst=[]
    for index,row in coord_df.iterrows():
        if skip_rc and row['dist_to_rc'] == 0:
            continue
        x_org,y_org,z_org=row['coord']
        x_org=float(x_org)
        y_org=float(y_org)
        z_org=float(z_org)

        # Add the sphere with the size equal to atomic radius
        coord_lst.append([x_org,y_org,z_org])
        elemt=row['elmt']
        ar=get_ar(elemt,ar_df)/100 # unit is angstrong
        x_sph,y_sph,z_sph=rand_sphere(x_org,y_org,z_org,ar,random_count=random_count)
        for i_sph in range(0,len(x_sph)):
            new_coord=[x_org+x_sph[i_sph],y_org+y_sph[i_sph],z_org+z_sph[i_sph]]
            coord_lst.append(new_coord)



    # Find point with maximum distance to the center (the radius of sphere)
    x_p = [p[0] for p in coord_lst]
    y_p = [p[1] for p in coord_lst]
    z_p = [p[2] for p in coord_lst]

    if print_detail:
        print('total number of point',len(x_p))

        print('max according to x',max(x_p))
        print('min according to x',min(x_p))

        print('max according to y',max(y_p))
        print('min according to y',min(y_p))

        print('max according to z',max(z_p))
        print('min according to z',min(z_p))

    # plot it


    P=np.array(coord_lst)


    
    # old version that uses a PCA component to align with y-axis
    '''yzcoord = [i[1:] for i in P]
    pca = PCA(n_components=2)
    pca.fit_transform(yzcoord)

    rotation_3dmatrix = rotation_matrix_from_vectors(np.array([0]+list(pca.components_[0])),np.array([0,1,0]))

    rotation_2dmatrix = np.array([ rotation_3dmatrix[i][1:3] for i in range(1,3)])


    yzcoord_rotated = []
    for row in yzcoord:
        new_value = np.matmul(rotation_2dmatrix,row)
        yzcoord_rotated.append(new_value)

    new_P = np.array([np.matmul(np.array(i),np.array(rotation_3dmatrix)) for i in P])
    P = new_P

    yzcoord_rotated_y = [i[0] for i in yzcoord_rotated]
    yzcoord_rotated_z = [i[1] for i in yzcoord_rotated]
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    ax.plot(yzcoord_rotated_y,yzcoord_rotated_z,'o',color='k')
    ax.set_title ('Principal Component Analysis',fontsize=16,fontweight='bold',family='sans-serif')
    ax.set_xlabel ('PC1',fontsize=14,fontweight='bold')
    ax.set_ylabel ('PC2',fontsize=14,fontweight='bold')

    plt.tick_params ('both',width=2,labelsize=12)
    plt.tight_layout()
    plt.show()'''

    

    # find the ellipsoid
    ET = EllipsoidTool()
    (center, radii, rotation) = ET.getMinVolEllipse(P, .01)
    print('radii',radii)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_box_aspect([1,1,1])
    ax.set_xlim3d(-10, 10)
    ax.set_ylim3d(-10, 10)
    ax.set_zlim3d(-10, 10)
    # plot points
    ax.scatter(P[:,0], P[:,1], P[:,2], color='g', s=100)

    # 2d center (ignore x coord) of deviation from the origin (rc atomic center)
    deviation_coord = center[1:]
    deviation_displ = (deviation_coord[0]**2 + deviation_coord[1]**2 )**(1/2)

    # plot ellipsoid
    axes = ET.plotEllipsoid(center, radii, rotation, ax=ax, plotAxes=True)
    # get axes projection on x-axis
    axes_list = np.abs(axes).tolist()
    axes_zy_list = [ (i[1]**2 + i[2]**2)**(1/2) for i in axes_list ]
    axes_zy_list.sort()
    print(axes_zy_list)
    axes_zylong = axes_zy_list[-1]
    axes_zyshort = axes_zy_list[-2]



    if plot_graph:
        plt.show()


    return(radii),(axes_zylong,axes_zyshort), deviation_displ

class EllipsoidTool:

    # original link: https://github.com/minillinim/ellipsoid/blob/master/ellipsoid.py

    """Some stuff for playing with ellipsoids"""
    def __init__(self): pass
    
    def getMinVolEllipse(self, P=None, tolerance=0.01):
        """ Find the minimum volume ellipsoid which holds all the points
        
        Based on work by Nima Moshtagh
        http://www.mathworks.com/matlabcentral/fileexchange/9542
        and also by looking at:
        http://cctbx.sourceforge.net/current/python/scitbx.math.minimum_covering_ellipsoid.html
        Which is based on the first reference anyway!
        
        Here, P is a numpy array of N dimensional points like this:
        P = [[x,y,z,...], <-- one point per line
             [x,y,z,...],
             [x,y,z,...]]
        
        Returns:
        (center, radii, rotation)
        
        """
        (N, d) = np.shape(P)
        d = float(d)
    
        # Q will be our working array
        Q = np.vstack([np.copy(P.T), np.ones(N)]) 
        QT = Q.T
        
        # initializations
        err = 1.0 + tolerance
        u = (1.0 / N) * np.ones(N)

        # Khachiyan Algorithm
        while err > tolerance:
            V = np.dot(Q, np.dot(np.diag(u), QT))
            M = np.diag(np.dot(QT , np.dot(linalg.inv(V), Q)))    # M the diagonal vector of an NxN matrix
            j = np.argmax(M)
            maximum = M[j]
            step_size = (maximum - d - 1.0) / ((d + 1.0) * (maximum - 1.0))
            new_u = (1.0 - step_size) * u
            new_u[j] += step_size
            err = np.linalg.norm(new_u - u)
            u = new_u

        # center of the ellipse 
        center = np.dot(P.T, u)
    
        # the A matrix for the ellipse
        A = linalg.inv(
                       np.dot(P.T, np.dot(np.diag(u), P)) - 
                       np.array([[a * b for b in center] for a in center])
                       ) / d
                       
        # Get the values we'd like to return
        U, s, rotation = linalg.svd(A)
        radii = 1.0/np.sqrt(s)
        
        return (center, radii, rotation)

    def getEllipsoidVolume(self, radii):
        """Calculate the volume of the blob"""
        return 4./3.*np.pi*radii[0]*radii[1]*radii[2]

    def plotEllipsoid(self, center, radii, rotation, ax=None, plotAxes=False, cageColor='b', cageAlpha=0.2):
        """Plot an ellipsoid"""
        make_ax = ax == None
        if make_ax:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            
        u = np.linspace(0.0, 2.0 * np.pi, 100)
        v = np.linspace(0.0, np.pi, 100)
        
        # cartesian coordinates that correspond to the spherical angles:
        x = radii[0] * np.outer(np.cos(u), np.sin(v))
        y = radii[1] * np.outer(np.sin(u), np.sin(v))
        z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
        # rotate accordingly
        for i in range(len(x)):
            for j in range(len(x)):
                [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotation) + center
    
        # make some purdy axes
        axes = np.array([[radii[0],0.0,0.0],
                            [0.0,radii[1],0.0],
                            [0.0,0.0,radii[2]]])
        # rotate accordingly
        for i in range(len(axes)):
            axes[i] = np.dot(axes[i], rotation)
    
    
        if plotAxes:
            # plot axes
            for p in axes:
                X3 = np.linspace(-p[0], p[0], 100) + center[0]
                Y3 = np.linspace(-p[1], p[1], 100) + center[1]
                Z3 = np.linspace(-p[2], p[2], 100) + center[2]
                ax.plot(X3, Y3, Z3, color=cageColor)
    
        # plot ellipsoid
        ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color=cageColor, alpha=cageAlpha)
        
        if make_ax:
            plt.show()
        
        return axes

def rand_sphere(x,y,z,ar,random_count = 200):
  # n points distributed evenly on the surface of a unit sphere
  z = 2 * random.rand(random_count) - 1   # uniform in -1, 1
  t = 2 * pi * random.rand(random_count)   # uniform in 0, 2*pi
  x = sqrt(1 - z**2) * cos(t)
  y = sqrt(1 - z**2) * sin(t)
  return ar*x, ar*y, ar*z

def get_ar(sym,ar_df):
    match_row=ar_df[(ar_df['sym']==sym)]
    return match_row['ar'].values[0]

if __name__ == '__main__':
    main()
