#!/usr/bin/env /Users/mon/opt/anaconda3/envs/instadeep_env/bin/python3.10

#Author: M. Penaloza-Amion
#Date : 11.12.2023


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import warnings
# suppress some MDAnalysis warnings when writing PDB files
warnings.filterwarnings('ignore')


class MyPDBTool:

    """
    Description:

    Super class to manage data of the PDB file using MDAnalysis package.

    Input: 
    - pdb_code: pdb code or name of pdb file it will be used for analysis.


    """

    def __init__(self, pdb_code: int):
        self.pdb_code = pdb_code

    def load_pdb(self):
        """ 
        Load PDB file after intializing the instance of class.
        
        """
        self.pdb_filename = f'{self.pdb_code}.pdb'.lower()
        self.pdb_file = mda.Universe(self.pdb_filename)

    def describe(self):

        """ 
        Provides informacion of chains, total atoms and total residues in PDB file.
        
        """

        self.total_atoms = len(self.pdb_file.atoms)
        self.chains = list(self.pdb_file.atoms.segments)
        self.total_aa = len(list(self.pdb_file.residues))

        print(f'*** PDB file {self.pdb_filename} ***')
        print(f'Total of chains: {len(self.chains)}')
        print(f'Chains available: {self.chains}')
        print(f'Total atoms : {self.total_atoms}')
        print(f'total residues : {self.total_aa}')
        print('*******************************')

        for chain_id in range(0, len(self.pdb_file.atoms.segments)) :

            chain = self.pdb_file.atoms.segments[chain_id]
            print(chain)
            print(f'Chain id : {chain_id}')
            print('Total atoms : ' + str(len(chain.atoms)))
            print('Total residues : ' + str(len(chain.residues)))
            print('*******************************')

    def select_chain(self, chain_id: int):
        """ 
        From instance object creates a object according to chain (Segment) selection.
        Input:

        chain_id : id that makes reference to chain ids displayed from discribe() method.
        
        """

        self.my_chain = self.pdb_file.atoms.segments[chain_id]
        return self.my_chain
    
    def dimensions(self):
        """ 
        Calculates the dimension as vectors and numpy array from PDB file system. Needed for distance calculations.
        
        """
        self.system_dimension = self.pdb_file.dimensions


    def get_list_id_from_pdb_id(self, pdb_id:int , protein_obj: object):

        """ 
        PDB ids doens't match with list ids from python. Here is a method to obtain list ids from PDB ids.
        
        Input:

        - pdb_id: residue id from PDB file
        - protein_obj: Object from mda instance or protein selection.

        Output:

        Int value of list id
        """
        self.pdb_id = pdb_id
        self.protein_obj = protein_obj

        for list_id, resid in enumerate(list(self.protein_obj.residues.resids)):

            if resid == self.pdb_id:
                return list_id
                #print(list(self.protein_obj.residues)[list_id])
                #print(f'list_id : {list_id}, pdb_residue : {resid}')

    def get_pdb_id_from_list_id(self, list_id: int, protein_obj: object):

        """ 
        PDB ids doens't match with list ids from python. Here is a method to obtain PDB ids from python list ids.
        
        Input:

        - list_id: residue id from list id

        - protein_obj: Object from mda instance or protein selection

        Output:

        Int value of PDB id


        """

        self.list_id = list_id
        self.protein_obj = protein_obj

        for idx, resid in enumerate(list(self.protein_obj.residues.resids)):

            if idx == self.list_id:
                return resid
                #print(list(self.protein_obj.residues)[idx])
                #print(f'list_id : {idx}, pdb_residue : {resid}')


class ProteinAProteinB(MyPDBTool):
    """

    Description:

    Usage:

    - prot_A: mda object with prot_A selection
    - prot_B: mda_object with prot_B selection
    - box_dimension : np.array of box dimension from mdAnalysis (MyPDBTools.dimension())
    - cutoff: max distance between residues from prot_A and residues prot_B
    """
     
    def __init__(self, prot_A, prot_B, box_dimension):
         #super().__init__(self.pdb_code)
         self.prot_A = prot_A
         self.prot_B = prot_B
         self.box_dimension = box_dimension
         #self.cutoff = cutoff

    def get_com_residue_residue_distance(self, resid1: int, resid2: int):

        """
        Obtain dintance between centrer of mass (com) of two residues
        res1 : residue id number
        res2 : residue id number        
        """
        self.resid1 = list(self.prot_A.residues)[resid1]
        self.resid1_com = self.resid1.atoms.center_of_mass(compound='residues')

        self.resid2 = list(self.prot_B.residues)[resid2]
        self.resid2_com = self.resid2.atoms.center_of_mass(compound='residues')

        distance = distances.distance_array(self.resid1_com, self.resid2_com, box=self.box_dimension)[0][0]

        return distance        

    def get_com_protein_residue_distance(self, resid: int, prot_A_com=True):

        """
        Obtain dintance between centrer of mass (com) of protein and residue

        prot_A_com: determines which protein will be taken as COM reference. First protein argument (True), second protein argument (False).

        resid : residue id from the protein is not taken as COM reference.
        
        """
        distance = 0
        if prot_A_com == True:

            total_residues_B = len(self.prot_B.residues)

            assert resid <= total_residues_B, f'resid : {resid} is biger than the total number of residues ({total_residues_B})'
            
            self.protein_com = self.prot_A.atoms.center_of_mass()
            self.resid = list(self.prot_B.residues)[resid]
            self.resid_com = self.resid.atoms.center_of_mass(compound='residues')
            distance = distances.distance_array(self.protein_com, self.resid_com, box=self.box_dimension)[0][0]

        else:


            total_residues_A = len(self.prot_A.residues)
            assert resid <= total_residues_A, f'resid : {resid} is biger than the total number of residues ({total_residues_A})'

            self.protein_com = self.prot_B.atoms.center_of_mass()
            self.resid = list(self.prot_A.residues)[resid]
            self.resid_com = self.resid.atoms.center_of_mass(compound='residues')
            distance = distances.distance_array(self.protein_com, self.resid_com, box=self.box_dimension)[0][0]

        return distance
    
    def SORT(self, sublist:list, position: int):

        """
        Method to sort nested lists. 

        sublist: Nested list to sort
        position: list id reference to sort
        
        """
        sublist.sort(key = lambda x: x[position])
        return sublist
    
    def search_protein_A_protein_B_interface(self, cutoff: float, second_prot_name: str, surface_stimate=0.2):

        """ 
        Method 1 to search for prot_A and prot_B interface residues

        Custom parameters:

        - cutoff: 8 (distances are between COM of diff residues)
        - surface_stimate: 20% as default.

        Search Algorithm:
        smaller protein is chosen, eg. Cov-2 receptor is smaller than ACE2
        1. COM of smaller protein is calculated and coordinates is taken as reference. Distances from all residues of bigger protein are meassured and sorted. Sublist1 is created to stored [resid, distance] and then sorted (small to big) and a % (surface stimate) of the closest residues to COM of reference is considered.
        2. Sublist1 residues ids (sublist1[0]) is used to calculated distances with all the residues of big protein. Finally, Sublist2 = [resid_A, resid_B , distance] is stored.
        
        3. Generate PDB files with custom occupancy values to visualize (eg. VMD) contact residues.

        4. Return nested list with [resid_A, resid_B, distance] of interface with PDB ids.
        
        """
        self.surface_estimate_res = int(len(self.prot_A.residues) * surface_stimate) #20 %
        self.cutoff = cutoff

        total_res_A = len(self.prot_A.residues)
        total_res_B = len(self.prot_B.residues)
        sublist1 = []
        self.contact_list = []
        
        if total_res_A < total_res_B :
            #1. scanning prot_A, prot_B_COM
            list_resid_A = list(self.prot_A.residues)
            for idx in range(len(list_resid_A)):
                dist_value = self.get_com_protein_residue_distance(idx, prot_A_com=False)
                sublist1.append([idx, dist_value]) # prot_A_idx, dist_to_COM_B
            subl_sorted = self.SORT(sublist1, 1)[:self.surface_estimate_res]

            #2. scanning prot_B with sublist from prot_A
            list_resid_B = list(self.prot_B.residues)
            for resid_B in range(len(list_resid_B)): #source of prot_B_ids
                for id_A in subl_sorted: #source of prot_A_ids
                    resid_A = id_A[0]
                    dist_res_res = self.get_com_residue_residue_distance(resid_A, resid_B)
                    if dist_res_res < self.cutoff:
                        #print([resid_A, resid_B, dist_res_res])
                        self.contact_list.append([resid_A, resid_B, round(dist_res_res,2)])
        else:
            #1. scanning prot_B, prot_A_COM
            list_resid_B = list(self.prot_B.residues)
            for idx in range(len(list_resid_B)):
                dist_value = self.get_com_protein_residue_distance(idx, prot_A_com=True)
                sublist1.append([idx, dist_value]) # prot_B_idx, dist_to_COM_A
            subl_sorted = self.SORT(sublist1, 1)[:self.surface_estimate_res]
            #print('2 ')

            #2. scanning prot_A with sublist from prot_B
            list_resid_A = list(self.prot_A.residues)
            for resid_A in range(len(list_resid_A)): #source of prot_A_ids

                for id_B in subl_sorted: #source of prot_B_ids
                    resid_B = id_B[1]

                    dist_res_res = self.get_com_residue_residue_distance(resid_A, resid_B)
                    if dist_res_res < self.cutoff:
                        self.contact_list.append([resid_A, resid_B, round(dist_res_res,2)])

        #create PDB files with occupancy value for visualization
        #for prot_A
        label_0_A = []
        label_1_A =[]
        tmp_1_A = []
        #one list first, then second, then merge.T_T
        for ai in range(total_res_A):
            label_0_A.append([ai,0])        

        for ai in self.SORT(self.contact_list,0):
            label_1_A.append([ai[0],ai[2]]) 
            tmp_1_A.append(ai[0])

        check_list = list(set(tmp_1_A))
        for ax in label_0_A:
            if ax[0] in check_list:
                label_0_A.remove(ax)
                label_0_A.insert(ax[0], [ax[0], 1])
        total_res_A_ids_labels = label_0_A

        #now from residues labels to atom labels, create list for atom_labels
        self.A_atom_labels = []
        for res in total_res_A_ids_labels:
            for atom in self.prot_A.residues[res[0]].atoms:
                self.A_atom_labels.append("%.2f" % res[1])

        #add occupancy lables to PDB
        occ_A = np.array(self.A_atom_labels)
        self.prot_A.atoms.write('tmp_prot_A_chain.pdb')
        p_A = mda.Universe('tmp_prot_A_chain.pdb')
        p_A.add_TopologyAttr('occupancies', values=occ_A)
        p_A.atoms.write('prot_A_chain.pdb')

        #for prot_B
        label_0_B = []
        label_1_B =[]
        tmp_1_B = []
        #one list first, then second, then merge.T_T
        for ai in range(total_res_B):
            label_0_B.append([ai,0])        

        for ai in self.SORT(self.contact_list,0):
            label_1_B.append([ai[1],ai[2]]) 
            tmp_1_B.append(ai[1])

        check_list2 = list(set(tmp_1_B))
        for ax in label_0_B:
            if ax[0] in check_list2:
                label_0_B.remove(ax)
                label_0_B.insert(ax[0], [ax[0], 1])
        total_res_B_ids_labels = label_0_B

        #now from residues labels to atom labels, create list for atom_labels
        self.B_atom_labels = []
        for res in total_res_B_ids_labels:
            for atom in self.prot_B.residues[res[0]].atoms:
                self.B_atom_labels.append("%.2f" % res[1])

        #add occupancy lables to PDB
        occ_B = np.array(self.B_atom_labels)
        self.prot_B.atoms.write('tmp_prot_B_chain.pdb')
        p_B = mda.Universe('tmp_prot_B_chain.pdb')
        p_B.add_TopologyAttr('occupancies', values=occ_B)
        p_B.atoms.write(f'prot_{second_prot_name}_chain.pdb')

        print('*** PDB files with modified occupancies were generated ***')

        os.system('rm tmp*')

        #convert list_dis to PDB ids
        self.contact_list_PDB_ids = []
        self.contact_list_PDB_ids_to_display = []

        for ax in self.contact_list:
            new_x = self.get_pdb_id_from_list_id(ax[0], self.prot_A)
            name_x = self.prot_A.residues[ax[0]].resname
            new_y = self.get_pdb_id_from_list_id(ax[1], self.prot_B)
            name_y = self.prot_B.residues[ax[1]].resname
            dist = ax[2]
            self.contact_list_PDB_ids.append([new_x, new_y, dist])
            self.contact_list_PDB_ids_to_display.append([f'{name_x}-{new_x}', f'{name_y}-{new_y}', dist])

        #sort by distances!
        self.contact_list_PDB_ids_to_display = self.SORT(self.contact_list_PDB_ids_to_display,2)

        #Save it as CSV file
        contact_list_array = np.array(self.contact_list_PDB_ids_to_display, dtype="O")
        contact_list_DF = pd.DataFrame(contact_list_array, columns=['Prot_A res_id', f'Prot_{second_prot_name} res_id', 'Distance (Angstrom)'])
        contact_list_DF.to_csv(f'contact_matrix_prot_A_prot_{second_prot_name}.csv')
        print('*** Contact matrix csv file was generated ***')
        
        return contact_list_DF

    def interface_contacts_plot(self, second_prot_name: str, color_cmap='gist_heat', PDB_ids=True):

        # 1. all the points
        #res values for id_list
        x_l= []
        y_l= []
        x_lowest= []
        y_lowest= []
        #res value for pdb_ids
        x_pdb_id = []
        y_pdb_id = []
        x_pdb_id_lowest = []
        y_pdb_id_lowest = []

        #this value is the same for both cases
        cmap_l= []
        cmap_lowest= []
        for idx in self.contact_list:
            xi = idx[0]
            yi = idx[1]
            zi = idx[2]
            x_l.append(xi)
            y_l.append(yi)
            #for pdb_ids
            xi_pdb_id = self.get_pdb_id_from_list_id(xi, self.prot_A)
            yi_pdb_id = self.get_pdb_id_from_list_id(yi, self.prot_B)
            x_pdb_id.append(xi_pdb_id)
            y_pdb_id.append(yi_pdb_id)
            #color map values
            cmap_l.append(zi)

        #2. lowest distance points cutoff = 5 to highlight the lowest points
            if idx[2] <=5.0:
                x_lowest.append(xi)
                y_lowest.append(yi)
                x_pdb_id_lowest.append(self.get_pdb_id_from_list_id(xi, self.prot_A))
                y_pdb_id_lowest.append(self.get_pdb_id_from_list_id(yi, self.prot_B))
                cmap_lowest.append(zi)
        #arrays        
        x_array= np.array(x_l)
        y_array= np.array(y_l)
        x_pdb_id_array= np.array(x_pdb_id)
        y_pdb_id_array= np.array(y_pdb_id)
    
        x_array_lowest= np.array(x_lowest)
        y_array_lowest= np.array(y_lowest)
        x_pdb_id_lowest_array= np.array(x_pdb_id_lowest)
        y_pdb_id_lowest_array= np.array(y_pdb_id_lowest)

        cmap= np.array(cmap_l)

        if PDB_ids == True: 

            ###PLOT settings
            plt.figure()
            plt.title(f'Contact map interface prot_A // prot_{second_prot_name} \n', fontsize=14)
            plt.xlabel(r'Prot_A res_id')
            plt.ylabel(rf'Prot_{second_prot_name} res_id')
            plt.scatter(x_pdb_id_array, y_pdb_id_array, c=cmap, cmap=color_cmap, s=50, marker=(','), alpha=0.85, linewidth=0)
            clb = plt.colorbar()
            clb.ax.set_title('Distance \n (Angstrom) \n')
            plt.scatter(x_pdb_id_lowest_array, y_pdb_id_lowest_array, color='black', s=50, marker=(','), linewidth=0)
            plt.rcParams["figure.figsize"] = (10,3.33)
            plt.tight_layout()
            #for i, ax in enumerate(x_lowest):
            #    plt.annotate([x_lowest[i], y_lowest[i]], (x_lowest[i]+5*i, y_lowest[i]-5*i))
            #plot_name = "prot_A_prot_B_interface.png"
            #plt.savefig(plot, dpi=300, bbox_inches='tight')
            plt.show()
        
        else:
            ###PLOT settings
            plt.figure()
            plt.title(f'Contact map interface prot_A // prot_{second_prot_name} \n', fontsize=14)
            plt.xlabel(r'Prot_A res_id')
            plt.ylabel(rf'Prot_{second_prot_name} res_id')
            plt.scatter(x_array, y_array, c=cmap, cmap=color_cmap, s=50, marker=(','), alpha=0.85, linewidth=0)
            clb = plt.colorbar()
            clb.ax.set_title('Distance \n (Angstrom) \n')
            plt.scatter(x_array_lowest, y_array_lowest, color='black', s=50, marker=(','), linewidth=0)
            plt.rcParams["figure.figsize"] = (10,3.33)
            plt.tight_layout()
            #for i, ax in enumerate(x_lowest):
            #    plt.annotate([x_lowest[i], y_lowest[i]], (x_lowest[i]+5*i, y_lowest[i]-5*i))
            #plot_name = f"prot_A_prot_{second_prot_name}_interface.png"
            #plt.savefig(plot, dpi=300, bbox_inches='tight')
            plt.show()

    def aminoacid_content(self, second_prot_name: str):
        """ 
        Calculate the frequency of each aminoacid for prot_A and prot_B.

        Input:
        second_prot_name : B or C to differenciate from other protein systems.
        
        """

        #for prot_A
        tmp_list1 = []
        for ai in self.SORT(self.contact_list,0):
            tmp_list1.append(ai[0])

        check_list1 = list(set(tmp_list1))
        aa_list_A = []
        for ax in check_list1:
            resname = self.prot_A.residues[ax].resname
            aa_list_A.append(resname)
        aa_list_compact_A = list(set(aa_list_A))

        #couting aa in A
        self.aa_count_A = []
        self.aa_count_A_percent = []
        total_aa_count_A = len(aa_list_A)
        for aa in aa_list_compact_A:
            aa_count = []
            for aa_item in aa_list_A:
                if aa_item == aa:
                    aa_count.append(aa_item)

            self.aa_count_A.append([aa, len(aa_count)])
            aa_percent = (len(aa_count)/total_aa_count_A)*100
            self.aa_count_A_percent.append([aa, round(aa_percent,1)])

        #for prot_B
        tmp_list2 = []
        for ai in self.SORT(self.contact_list,0):
            tmp_list2.append(ai[1])

        check_list2 = list(set(tmp_list2))
        aa_list_B = []
        for ax in check_list2:
            resname = self.prot_B.residues[ax].resname
            aa_list_B.append(resname)
        aa_list_compact_B = list(set(aa_list_B))

        #couting aa in B
        self.aa_count_B = []
        self.aa_count_B_percent = []
        total_aa_count_B = len(aa_list_B)
        for aa in aa_list_compact_B:
            aa_count = []
            for aa_item in aa_list_B:
                if aa_item == aa:
                    aa_count.append(aa_item)

            self.aa_count_B.append([aa, len(aa_count)])
            aa_percent = (len(aa_count)/total_aa_count_B)*100
            self.aa_count_B_percent.append([aa, round(aa_percent,1)])

        #polar, non-polar aa content

        non_polar = ['ALA','GLY','VAL', 'ILE', 'LEU', 'MET','CYS', 'PHE', 'TYR', 'TRP','PRO']
        polar_charged = ['SER', 'THR', 'ASN', 'GLN']
        positive = ['ARG', 'HIS','LYS']
        negative = ['ASP', 'GLU']

        non_polar_count_A =[]
        polar_charged_count_A = []
        positive_count_A = []
        negative_count_A = []

        for resname in self.aa_count_A_percent:
            if resname[0] in non_polar:
                non_polar_count_A.append(resname[1])
            elif resname[0] in polar_charged:
                polar_charged_count_A.append(resname[1])
            elif resname[0] in positive:
                positive_count_A.append(resname[1])
            elif resname[0] in negative:
                negative_count_A.append(resname[1])

        print('***** Aminoacid restype % prot_A ****')
        print(f'non-polar count : {round(sum(non_polar_count_A),2)} %')
        print(f'polar charged count : {round(sum(polar_charged_count_A),2)} %')
        print(f'positive count : {round(sum(positive_count_A),2)} %')
        print(f'negative count : {round(sum(negative_count_A),2)} %')
        print('********* \n')


        non_polar_count_B =[]
        polar_charged_count_B = []
        positive_count_B = []
        negative_count_B = []

        for resname in self.aa_count_B_percent:
            if resname[0] in non_polar:
                non_polar_count_B.append(resname[1])
            elif resname[0] in polar_charged:
                polar_charged_count_B.append(resname[1])
            elif resname[0] in positive:
                positive_count_B.append(resname[1])
            elif resname[0] in negative:
                negative_count_B.append(resname[1])

        print(f'***** Aminoacid restype % prot_{second_prot_name} ****')
        print(f'non-polar count : {round(sum(non_polar_count_B),2)} %')
        print(f'polar charged count : {round(sum(polar_charged_count_B),2)} %')
        print(f'positive count : {round(sum(positive_count_B),2)} %')
        print(f'negative count : {round(sum(negative_count_B),2)} %')
        print('********* \n')


        #Bar Plot
        #merge labels in the same order as all the input lists
        labels = sorted(list(set(aa_list_A + aa_list_B)))
        A_B_difference = list(set(aa_list_A).difference(aa_list_B)) #to add in B
        B_A_difference = list(set(aa_list_B).difference(aa_list_A)) # to add in A

        #1. add diff aa to lists and asign value 0

        for aa_diff in B_A_difference:
            self.aa_count_A_percent.append([aa_diff, 0.0])

        for aa_diff in A_B_difference:
            self.aa_count_B_percent.append([aa_diff, 0.0])
 
        self.sorted_aa_count_A_percent = sorted(self.aa_count_A_percent, key=lambda x:x[0])
        self.sorted_aa_count_B_percent = sorted(self.aa_count_B_percent, key=lambda x:x[0])

        #get list to plot
        labels = []
        prot_A_bar_values = []
        prot_B_bar_values = []

        for ax in self.sorted_aa_count_A_percent:
            labels.append(ax[0])
            prot_A_bar_values.append(ax[1])

        for ax in self.sorted_aa_count_B_percent:
            prot_B_bar_values.append(ax[1])

        x=np.arange(len(labels))
        my_font_size = 10
        width = 0.35
        fig, ax = plt.subplots()
        bar1 = ax.bar(x - width/2, prot_A_bar_values, width, label='prot_A',color='c', ecolor= 'black', capsize=5, tick_label=my_font_size)

        bar2 = ax.bar(x + width/2, prot_B_bar_values, width, label=f'prot_{second_prot_name}',color='m', ecolor= 'black', capsize=5, tick_label=my_font_size)
        #plot settings
        ax.set_title(f'Aminoacid content in interface prot_A // prot_{second_prot_name}', fontsize=my_font_size)
        ax.set_ylabel('Frequency (%)', fontsize=my_font_size)
        plt.yticks(size=12)
        plt.xticks(x, labels, rotation=45, fontsize=my_font_size)
        ax.set_ylim([0, 100])
        ax.set_xticks(x)
        ax.tick_params(direction='out', length=4, width=2, colors='black', grid_color='black', grid_alpha=0.5)
        ax.legend(loc='upper right',prop={'size': my_font_size})
        fig.tight_layout()
        fig.set_size_inches(6.66, 3.33)
        #plt.savefig(f'aminoacid_content_prot_A_prot_{second_prot_name}.png', dpi=300, bbox_inches='tight')
        plt.show()

        #save data in csv file

        aa_content_dict = {'Aminoacid':labels, 'prot_A %': prot_A_bar_values, 'prot_B %': prot_B_bar_values}
        aa_contet_DF = pd.DataFrame(aa_content_dict)
        aa_contet_DF.to_csv(f'aa_content_interface_A{second_prot_name}.csv')