import gemmi
import os

bound_occ = 0.3

ground = gemmi.read_structure('init.pdb')
bound = gemmi.read_structure('pandda.pdb')

ground_model = ground[0]
bound_model = bound[0]

for cra in ground_model.all():
	current_occ = cra.atom.occ
	cra.atom.occ = current_occ - current_occ * 0.3

for cra in bound_model.all():
	current_occ = cra.atom.occ
	cra.atom.occ = current_occ * 0.3

# the bound state should be the first since it is the scientifically most relevant
bound.add_model(ground_model, pos=2)
bound.write_pdb('ensemble.pdb')
os.system('pdb_extract -iPDB ensemble.pdb -o ensemble.mmcif -NMR > pdb_extract.log')

cmd = ( 'refmac5 hklin free.mtz hklout ensemble.mtz xyzin ensemble.pdb xyzout refmac.pdb libin VT00025.cif << eof > refmac.log\n'
	'ncyc 0\n'
	'end\n'
	'eof\n'
	)
os.system(cmd)

# remove loop/ table
s = gemmi.cif.read_file('refmac.mmcif')
#b = s.sole_block()
for b in s:
	if b.name == "XXXX":
		coordinate_block = b
		header = b.find_loop('_atom_site.group_PDB').get_loop().tags
		table = b.find(header)
		table.erase()
		break

e = gemmi.cif.read_file('ensemble.mmcif')
eb = e.sole_block()

atom_site_header = eb.find_loop('_atom_site.group_PDB').get_loop().tags
# need to do this, otherwise end up with _atom_site._atom_site.group_PDB in xyz_table below
atom_site_labels = []
for i in atom_site_header:
	atom_site_labels.append(i.replace('_atom_site.', ''))

xyz = eb.find(atom_site_header)
xyz_table = coordinate_block.init_loop('_atom_site.', atom_site_labels)
for row in xyz:
	if not row['_atom_site.type_symbol'] == 'H':
		xyz_table.add_row(row)

s.write_file('test.mmcif')


