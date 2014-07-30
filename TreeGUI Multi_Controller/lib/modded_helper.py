import re
import filecmp
import math
import os
# Reads in the raw output of GTP (or the two file version)
# and creates a dissimilarity matrix file of the distances.
# Arguments: raw distance file output from GTP
#            output file
#            "symmetric" or "non" : whether matrix should be symmetric or not
# Returns:  nothing (writes to output file)

for geo_file_name in os.listdir(os.getcwd()):
    print(geo_file_name)
    sym_flag = "symmetric"
    if(geo_file_name.endswith("_GTPrawOut.txt")):
        diss_file_name = geo_file_name.replace("_GTPrawOut.txt","") + "_symMatrix.txt"
        print(diss_file_name)
        max_row = 0
        max_col = 0
        # Loop through the geo file the first time to get the dimensions
        # of diss_matrix
        geo_dist_file = open(geo_file_name,'r')
        for line in geo_dist_file:
            match_line = re.match(r'([0-9]+)\s([0-9]+)\s([0-9\.E-]+)',line)
            if (match_line):
                if ( int(match_line.group(1)) > max_row):
                    max_row = int(match_line.group(1))
                if ( int(match_line.group(2)) > max_col):
                    max_col = int(match_line.group(2))
        geo_dist_file.close()
     
        # Create and populate the diss_matrix
        geo_dist_file = open(geo_file_name,'r')
     
        # add one to max_row if symmetric matrix
        # (because last row is not output by GTP because determined by symmetry)
        if (sym_flag == "symmetric"):
            max_row= max_row + 1
            if (max_row != max_col):
                print "Error:  symmetric matrix but read in %d rows and %d columns.  Exiting" % (max_row + 1, max_col +1)
                exit(1)
     
        # add one because rows/cols were labelled 0, 1, 2, ... by GTP
        diss_matrix = [ [0]*(max_col+1) for x in xrange(max_row+1)]
     
        # Loop through the geo file the second time to populate diss_matrix
        for line in geo_dist_file:
            match_line = re.match(r'([0-9]+)\s([0-9]+)\s([0-9\.E-]+)',line)
     
            if (match_line):
                diss_matrix[int(match_line.group(1))][int(match_line.group(2))] = match_line.group(3)
                if (sym_flag == "symmetric"):
                    diss_matrix[int(match_line.group(2))][int(match_line.group(1))] = match_line.group(3)
     
        geo_dist_file.close()
     
        # print the geodesic distances as a dissimilarity matrix to the output file.
        o = open(diss_file_name,'w')
     
        for row in range(max_row + 1):
            for col in range(max_col + 1):
                o.write(str(diss_matrix[row][col]) + ' ')
            o.write('\n')
        o.close()
 

    
