import numpy as np
import MDAnalysis as mda




#Calculating Area per Lipid (APL) for a membrane
def APL(pdb,dcd,out):
    # Define the input parameters
    trajectory_file = dcd
    pdb_path = pdb
    f =  open(out,'w')
    f.write("frame  APL_upper APL_lower\n")  
    # Define the distance cutoff (10 angstroms in this case)
    
    protein_cutoff = 10
    
    max_points = 100
    cutoff_distance = 10
    
    # Load the PDB file using MDAnalysis
    universe = mda.Universe(pdb_path,trajectory_file)
    
    # Select the lipid atoms in the system
    lipid_atoms = universe.select_atoms("resname POPC and name P")
    protein_atoms = universe.select_atoms("protein")
    t = 0
    for i in universe.trajectory[::100]:
      print(t)
      box = universe.dimensions[:3]
      
      
      COM_z = lipid_atoms.center_of_mass()[2]
      
      
      # Define a function to check if a point is within the cutoff distances
      def is_point_valid(point):
          # Check if the point is within the protein cutoff distance of any protein atom
          for atom in protein_atoms:
              distance = np.linalg.norm(atom.position - point)
              if distance <= protein_cutoff:
                  return False
          # If the point is not within either cutoff distance, it is valid
          return True
      
      nr_points = 0
      # Initialize the counters for the upper and lower leaflets
      upper_leaflet_count = 0
      lower_leaflet_count = 0
      
      
      while nr_points < max_points:
          # Define the coordinates of the chosen point
          chosen_point = np.array([np.random.random_sample()*(box[0]-40)+20, np.random.random_sample()*(box[1]-40)+20, COM_z])
          if is_point_valid(chosen_point):
              #print(chosen_point)
              nr_points += 1
              
              # Iterate over the lipid atoms
              for atom in lipid_atoms:
                  # Calculate the distance between the atom and the chosen point
                  distance = np.linalg.norm(atom.position[:2] - chosen_point[:2])
              
                  # Check if the distance is within the cutoff distance
                  if distance <= cutoff_distance:
                      # Detemine if the atom is in the upper or lower leaflet 
                      if atom.position[2] > COM_z:
                          upper_leaflet_count += 1
                      else:
                          lower_leaflet_count += 1
              
      # Print the results
      f.write(str(t)+ ' '+str(3.1415/(upper_leaflet_count/max_points))+' '+str(3.1415/(lower_leaflet_count/max_points))+'\n')
      #print("area per lipid (upper leaflet):", 3.1415/(upper_leaflet_count/max_points),  'nm2')
      #print("area per lipid (lower leaflet):", 3.1415/(lower_leaflet_count/max_points), 'nm2')
      t += 1
    f.close()

def min_distance(point, points_array):
    dists = np.linalg.norm(points_array - point, axis=1)
    return np.min(dists)

def Membrane_thickness(pdb,dcd,out):
    # Define the input parameters
    trajectory_file = dcd
    pdb_path = pdb
    f =  open(out,'w')
    f.write("frame  average_thickness std\n")

    # Load the PDB file using MDAnalysis
    universe = mda.Universe(pdb_path,trajectory_file)

    # Select the lipid atoms in the system
    lipid_atoms = universe.select_atoms("resname POPC and name P")
    t = 0
    for i in universe.trajectory[::10]:
      COM_z = lipid_atoms.center_of_mass()[2]
      print(t)
      mem_width = []

      for atomu in lipid_atoms: 
        if atomu.position[2] > COM_z:
          atomsl = [atoml.position for atoml in lipid_atoms if atoml.position[2] < COM_z]
          mem_width += [min_distance(atomu.position, np.array(atomsl))]
          #print(mem_width)
      
      # Print the results
      f.write(str(t)+ ' '+str(np.mean(mem_width))+' '+str(np.std(mem_width))+'\n')
      t += 1
      print(t)#,str(t)+ ' '+str(np.mean(mem_width))+' '+str(np.std(mem_width)))
      #print("area per lipid (upper leaflet):", 3.1415/(upper_leaflet_count/max_points),  'nm2')
    f.close()
#Membrane_thickness('../Coel_Kir61_closed/Runs/Run1/md.gro','../Coel_Kir61_closed/Runs/Run1/md_trjconv.xtc','tmp.dat')

from multiprocessing import Process
import os



p=[]

aliases = ['../Coel_Kir61_closed/Runs/Run1/','../Coel_Kir61_closed/Runs/Run2/','../Coel_Kir61_closed/Runs/Run3/','../Coel_Kir61_open/Runs/Run1/', '../Coel_Kir61_open/Runs/Run2/','../Coel_Kir61_open/Runs/Run3/','../Coel_Kir62_closed/Runs/Run1/', '../Coel_Kir62_closed/Runs/Run2/','../Coel_Kir62_closed/Runs/Run3/','../Coel_Kir62_open/Runs/Run1/','../Coel_Kir62_open/Runs/Run2/','../Coel_Kir62_open/Runs/Run3/','../Coel_Kir63_closed/Runs/Run1/','../Coel_Kir63_closed/Runs/Run2/','../Coel_Kir63_closed/Runs/Run3/','../Coel_Kir63_open/Runs/Run1/','../Coel_Kir63_open/Runs/Run2/','../Coel_Kir63_open/Runs/Run3/']

for a in aliases:
	
	pdb = a + 'md.gro'
	#new = a + 'apl.dat'
	new = a + 'mem_thickness.dat'
	dcd = a + 'md_trjconv.xtc'
	#p.append(Process(target=APL, args=(pdb,dcd,new)))
	p.append(Process(target=Membrane_thickness, args=(pdb,dcd,new)))
for p_ in p:
	p_.start()
for p_ in p:
	p_.join()

