#! /usr/bin/awk -f
#
#   Calculate RMSD of atoms with the same name in two PDB files		James Holton 9-5-11
#
#	The PDB feild:
#           |<--- here    -->|
#ATOM      1  N   ALA A 327      40.574  34.523  43.012  1.00 34.04
#
#	is used to determine if two atoms are a "pair"
#
BEGIN {
if(! atom) atom = "CA"
maxXYZ = maxdB = maxdO = 0
max_atom_XYZ = max_atom_dB = max_atom_dO = 0
maxXYZ_ID = maxdB_ID = maxdO_ID = "-"
max_atom_XYZ_ID = max_atom_dB_ID = max_atom_dO_ID = "-"
}

/^ATOM|^HETATM/{
    # read in values (watching for duplications)
    ID = substr($0,12,18)
    ++count[ID]

    if(count[ID] == 1)
    {
	# not initialized yet
	X[ID] = substr($0, 31, 8)+0
	Y[ID] = substr($0, 39, 8)+0
	Z[ID] = substr($0, 47, 8)+0
	O[ID] = substr($0, 55, 6)+0
	B[ID] = substr($0, 61, 6)+0
    }
    
    if(count[ID] == 2)
    {
	++pairs
    
	# seen this before, subtract values
	dX     = X[ID] - substr($0, 31, 8)
	dY     = Y[ID] - substr($0, 39, 8)
	dZ     = Z[ID] - substr($0, 47, 8)
	dO[ID] = O[ID] - substr($0, 55, 6)
	dB[ID] = B[ID] - substr($0, 61, 6)
	
	# get drift (and add up squares of drifts)
	sqrD   = dX*dX + dY*dY + dZ*dZ
	dXYZ[ID] = sqrt(sqrD)

	# remember maximum shifts
	if(dXYZ[ID] > maxXYZ) {maxXYZ  = dXYZ[ID]; maxXYZ_ID = ID }
	if(dO[ID]*dO[ID] > maxdO*maxdO) {maxdO = dO[ID]; maxdO_ID = ID }
	if(dB[ID]*dB[ID] > maxdB*maxdB) {maxdB = dB[ID]; maxdB_ID = ID }

	# maintain mean-square sums	
	sumXYZ += sqrD
	sumO   += dO[ID]*dO[ID]
	sumB   += dB[ID]*dB[ID]

	# separate stats for special atom type
	split(ID,word)
	if(word[1] == atom)
	{
	    ++atom_pairs

	    # maintain separate mean-square sums
	    sum_atom_XYZ += sqrD
	    sum_atom_O   += dO[ID]*dO[ID]
	    sum_atom_B   += dB[ID]*dB[ID]
	    
	    # remember maximum drifts too
	    if(dXYZ[ID] > max_atom_XYZ) {max_atom_XYZ  = dXYZ[ID]; max_atom_XYZ_ID = ID }
	    if(dO[ID]*dO[ID] > max_atom_dO*max_atom_dO) {max_atom_dO = dO[ID]; max_atom_dO_ID = ID }	    
	    if(dB[ID]*dB[ID] > max_atom_dB*max_atom_dB) {max_atom_dB = dB[ID]; max_atom_dB_ID = ID }	    
	}
	# debug output
	if(debug)
	{
	    printf("%s moved %8.4f (XYZ) %6.2f (occ) %6.2f (B) at %s\n", ID, dXYZ[ID],dO[ID], dB[ID], "cen_x " substr($0,31,24))
	}	
    }

    if(count[ID] > 2)
    {
	print "WARNING: " ID " appeared more than twice! "
    }
}


END{
    
    if(pairs+0 == 0) 
    {
	print "no atom pairs found"
	exit
    }
    rmsXYZ = sqrt(sumXYZ/pairs)
    rmsO = sqrt(sumO/pairs)
    rmsB = sqrt(sumB/pairs)
    if(atom_pairs+0 != 0)
    {
	rms_atom_XYZ = sqrt(sum_atom_XYZ/atom_pairs)
	rms_atom_O = sqrt(sum_atom_O/atom_pairs)
	rms_atom_B = sqrt(sum_atom_B/atom_pairs)
    }
    

    if(! xlog) 
    {
	print pairs " atom pairs found"
	print "RMSD("atom" )= " rms_atom_XYZ " ("atom_pairs, atom " pairs)"
	print "RMSD(all)= " rmsXYZ " ("pairs" atom pairs)"
	print "RMSD(Bfac)= " rmsB
	if(maxdO_ID != "-") print "RMSD(occ)= " rmsO

	print "MAXD(all)= " maxXYZ "\tfor " maxXYZ_ID
	print "MAXD(Bfac)= " maxdB "\tfor " maxdB_ID
	if(maxdO_ID != "-") print "MAXD(occ)= " maxdO "\tfor " maxdO_ID
	
	# final check for orphan atoms 
	for(id in count)
	{
	    if(count[id]<2) print "WARNING: " id " only found once"
	}
    }
    else
    {
#	printf "%10.8f %10.8f %10.5f %8.4f %10.4f %8.3f %8.3f    %s %s %s\n", rms_atom_XYZ, rmsXYZ, rmsO, rmsB, maxXYZ, maxO,maxdB, maxXYZ_ID,maxdB_ID,maxdO_ID
	printf "%10.8f %10.8f %10.5f %8.4f %10.4f %8.3f %8.3f\n", rms_atom_XYZ, rmsXYZ, rmsO, rmsB, maxXYZ, maxdO,maxdB
    }
}
