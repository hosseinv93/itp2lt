#important: in the PDB file the x,y, and z are 6, 7, and 8th columns. If they are ither columns in your file, modify the code in the PDB section.

tt=0.238846 #kj2kcal
def convert_itp_to_lt(itp_file):
    lt_lines = []
    atom_section = False
    bond_section = False
    angle_section = False
    dihedral_section = False
    improper_section = False
    counter=False
    already_added_atomtype = set()
    
    # Parse GROMACS ITP file
    bonds = []
    angles = []
    dihedrals = []
    impropers = []
    masses=[]
    data=[]
    bond_coeffs=False
    angle_coeffs=False
    dihedral_coeffs=False
    improper_coeffs=False
    pair_coeffs=False

    with open(itp_file, 'r') as f, open('ffnonbonded.itp', 'r') as g:

        lines = f.readlines()

    for line in lines:
        line = line.strip()
        if line.startswith('#define dfTPP_bon'):
            if not bond_coeffs:
                lt_lines.append(f'\n\t# Bond coefficients')
                lt_lines.append(f"\twrite_once(\"In Settings\") {{\n\n\t# Bond coefficients")
            bond_coeffs=True
            tokens = line.split()
            definition = tokens[1]
            value1 = float(tokens[3])
            value2 = float(tokens[4])
            lt_lines.append(f'\t bond_coeff @bond:{definition} harmonic {value2*0.00238846:.6f} {value1*10:.6f}')

        elif line.startswith('#define dfTPP_ang'):
            if not angle_coeffs:
                lt_lines.append(f'\t}}\n\n \t# Angle coefficients')
                lt_lines.append(f"\twrite_once(\"In Settings\") {{\n\n\t# Angle coefficients")
            angle_coeffs=True
            tokens = line.split()
            definition = tokens[1]
            value1 = float(tokens[3])
            value2 = float(tokens[4])
            lt_lines.append(f'\t angle_coeff @angle:{definition} harmonic {value2*0.238846:.6f} {value1:.6f}')

        elif line.startswith('#define dfTPP_dih'):
            if not dihedral_coeffs:
                lt_lines.append(f'\t}}\n\n \t# Dihedral coefficients')
                lt_lines.append(f"\twrite_once(\"In Settings\") {{\n\n\t# Dihedral coefficients")
            dihedral_coeffs=True
            tokens = line.split()
            definition = tokens[1]
            value1 = float(tokens[3])
            value2 = float(tokens[4])
            value3 = float(tokens[5])
            value4 = float(tokens[6])
            lt_lines.append(f'\t dihedral_coeff @dihedral:{definition} opls {value1*tt:.6f} {value2*tt:.6f} {value3*tt:.6f} {value4*tt:.6f}')

        elif line.startswith('#define improper'):
            if not improper_coeffs:
                lt_lines.append(f'\t}}\n\n \t# Improper coefficients')
                lt_lines.append(f"\twrite_once(\"In Settings\") {{\n\n\t# Improper coefficients")
            improper_coeffs=True
            tokens = line.split()
            definition = tokens[1]
            value1 = float(tokens[2])
            value2 = float(tokens[3])
            lt_lines.append(f'\t improper_coeff @improper:{definition} harmonic {value2*tt:.6f} {value1:.6f}')
            
        # Parse atom information
        elif line.startswith('[ atoms ]'):
                atom_section = True
                # Skip headers

        elif atom_section:
                tokens = line.split()
                if len(tokens) >= 8:
                    atom_index = tokens[0]
                    atom_type = tokens[1]
                    residue_index = tokens[2]
                    residue_name = tokens[3]
                    atom_name = tokens[4]
                    charge = tokens[6]
                    mass = tokens[7]

                    # data atoms
                    with open('mol.pdb', 'r') as pdb_file:
                        lines_pdb = pdb_file.readlines()
                        # Find the matching lines in the PDB
                        matching_lines = []
                        for line_pdb in lines_pdb:
                            if line_pdb.startswith('ATOM') and line_pdb.split()[1] == atom_index:
                                data.append(f'\t $atom:{atom_index} $mol:m1 @atom:{atom_type} {charge} {line_pdb.split()[6]} {line_pdb.split()[7]} {line_pdb.split()[8]} ')  


                    #Finding pair coeffs
                    # Read ffnonbonded.itp
                    with open('ffnonbonded.itp', 'r') as file2:
                        lines_file2 = file2.readlines()
                        # Find the matching lines in the second file
                        matching_lines = []
                        if not pair_coeffs:
                            lt_lines.append(f'\t}}\n\n \t# Pair coefficients')
                            lt_lines.append(f"\twrite_once(\"In Settings\") {{\n\n\t# Pair coefficients")
                        pair_coeffs=True
                        for line_file2 in lines_file2:
                            if line_file2.split()[0] in atom_type and line_file2.split()[0] not in already_added_atomtype:
                                lt_lines.append(f'\t pair_coeff @atom:{line_file2.split()[0]} @atom:{line_file2.split()[0]} lj/cut/coul/long {float(line_file2.split()[7])*tt:.5f} {float(line_file2.split()[6])*10:.6f}')

                                # masses
                                masses.append(f'\t @atom:{atom_type} {mass}')

                                already_added_atomtype.add(line_file2.split()[0])
                
                elif line.strip() == '':
                    atom_section = False
                    lt_lines.append(f'\t}}')
                    #break  # End of atom section
        #bonds
        elif line.startswith('[ bonds ]'):
                bond_section = True
                i=1
        elif bond_section:
            tokens = line.split()
            if len(tokens) >= 3:
                bond_index = f"b{i}"
                atom1_index = tokens[0]
                atom2_index = tokens[1]
                bond_type = tokens[2]
                i+=1
                bonds.append(f'\t $bond:{bond_index} @bond:{bond_type} $atom:{atom1_index} $atom:{atom2_index}')

            elif line.strip() == '':
                    bond_section = False
                    #break  # End of bond section
        #angles
        elif line.startswith('[ angles ]'):
                angle_section = True
                i=1
        elif angle_section:
            tokens = line.split()
            if len(tokens) >= 4:
                angle_index = f"ang{i}"
                atom1_index = tokens[0]
                atom2_index = tokens[1]
                atom3_index = tokens[2]
                angle_type = tokens[3]
                i+=1
                angles.append(f'\t $angle:{angle_index} @angle:{angle_type} $atom:{atom1_index} $atom:{atom2_index} $atom:{atom3_index}')

            elif line.strip() == '':
                    angle_section = False
                    #break  # End of bond section

        #dihedrals
        elif line.startswith('[ dihedrals ]') and not counter:
                counter = True
                dihedral_section = True
                i=1
        elif dihedral_section:
            tokens = line.split()
            if len(tokens) >= 5:
                dihedral_index = f"d{i}"
                atom1_index = tokens[0]
                atom2_index = tokens[1]
                atom3_index = tokens[2] 
                atom4_index = tokens[3]
                dihedral_type = tokens[4]
                i+=1
                dihedrals.append(f'\t $dihedral:{dihedral_index} @dihedral:{dihedral_type} $atom:{atom1_index} $atom:{atom2_index} $atom:{atom3_index} $atom:{atom4_index}')

            elif line.strip() == '':
                    dihedral_section = False
                    #break  # End of bond section

        #dihedrals
        elif line.startswith('[ dihedrals ]') and counter:
                improper_section = True
                i=1
        elif improper_section:
            tokens = line.split()
            if len(tokens) >= 6:
                improper_index = f"I{i}"
                atom1_index = tokens[0]
                atom2_index = tokens[1]
                atom3_index = tokens[2] 
                atom4_index = tokens[3]
                improper_type = tokens[5]
                i+=1
                impropers.append(f'\t $improper:{improper_index} @improper:{improper_type} $atom:{atom1_index} $atom:{atom2_index} $atom:{atom3_index} $atom:{atom4_index}')

            elif line.strip() == '':
                    improper_section = False
                    #break  # End of bond section

        
    return lt_lines, masses, bonds, angles, dihedrals, impropers, data


# Usage example
itp_file = './mol.itp'
lt_lines = convert_itp_to_lt(itp_file)[0]
masses = convert_itp_to_lt(itp_file)[1]
bonds = convert_itp_to_lt(itp_file)[2]
angles = convert_itp_to_lt(itp_file)[3]
dihedrals = convert_itp_to_lt(itp_file)[4]
impropers = convert_itp_to_lt(itp_file)[5]
data = convert_itp_to_lt(itp_file)[6]

g = open("mol.lt", "w")
# Print the Moltemplate input
g.write('molecule {\n\n \t# Interaction coefficients\n')

for line in lt_lines:
    g.write('{}\n'.format(line))

g.write('\n\n\t# Masses')
g.write(f"\n\twrite_once(\"Data Masses\") {{\n")
for line in masses:
    g.write('{}\n'.format(line))
g.write('\t}')

g.write('\n\n\t# Data Files : atom-id mol-id atom-type charge  X     Y      Z ')
g.write(f"\n\twrite(\"Data Atoms\") {{\n")
for line in data:
    g.write('{}\n'.format(line))
g.write('\t}')

g.write('\n\n\t# Bonds')
g.write(f"\n\twrite(\"Data Bonds\") {{\n")
for line in bonds:
    g.write('{}\n'.format(line))
g.write('\t}')

g.write('\n\n\t# Angles')
g.write(f"\n\twrite(\"Data Angles\") {{\n")
for line in angles:
    g.write('{}\n'.format(line))
g.write('\t}')

g.write('\n\n\t# Dihedrals')
g.write(f"\n\twrite(\"Data Dihedrals\") {{\n")
for line in dihedrals:
    g.write('{}\n'.format(line))
g.write('\t}')

g.write('\n\n\t# Impropers')
g.write(f"\n\twrite(\"Data Impropers\") {{\n")
for line in impropers:
    g.write('{}\n'.format(line))
g.write('\t}')


g.write('\n}')
g.close()
        

