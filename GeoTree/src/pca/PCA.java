package pca;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import polyAlg.PolyMain;
import distanceAlg1.*;
import distances.XORShiftRandom;
import polyAlg.Tools;
import static polyAlg.PolyMain.getGeodesic;
import static polyAlg.PolyMain.calcGeoDist;

public class PCA {

	
	/** Help message (ie. which arguments can be used, etc.)
	 * 
	 */
	public static void displayHelp() {
		System.out.println("Command line syntax:");
		System.out.println("java -jar pca.jar [options] treefile");
		System.out.println("Optional arguments:");	
		System.out.println("\t -a <algorithm> \t specifies what to compute");
		System.out.println("\t -e <epsilon> \t specifies precision for projections");
		System.out.println("\t -i <iterations> \t specifies numnber of iterations");
		System.out.println("\t -o <outfile> \t store the output in the file <outfile>.  Default is output.txt");
		System.out.println("\t -r \t outputs number of orthants the geodesic passes through");
		System.out.println("\t -u \t trees are unrooted. Default is trees are rooted.");
		System.out.println("\n");
		System.out.println("Algorithms are:");
		System.out.println("\t random \t randomly selects two trees, computes the geodesic between them and the projection score onto that geodesic");
	}
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Boolean rooted = true;   // default
		Boolean orthants = false;   // default
		String treeFile = null;
		String otherTreeFile = null;
		String algorithm = null;
		int numIter = 100;   // default
		String outfile = "output.txt";
		Double epsilon = 0.5; //default
		PhyloTree[] otherTrees = null;

		/* Parse command line arguments */

		if (args.length < 1) {
			System.out.println("Error: Missing input file name"); displayHelp(); System.exit(1);
		}
		treeFile = args[args.length-1];
		for (int i = 0; i < args.length - 1; i++) {
			
			if (!args[i].startsWith("-")) { System.out.println("Invalid command line option"); displayHelp(); System.exit(1); }
			// specify algorithm
			else if (args[i].equals("-a")) {
				if (i < args.length -2) { algorithm = args[i+1]; i++; }
				else { System.err.println("Error: algorithm not specified"); displayHelp(); System.exit(1); }
			}
			// epsilon, if desired
			else if (args[i].equals("-e")) {
				if (i < args.length -2) { epsilon = Double.valueOf(args[i+1]); i++; }
				else { System.err.println("Error: value for epsilon not specified"); displayHelp(); System.exit(1); }
			}
			// other tree file, if desired
			else if (args[i].equals("-f")) {
				if (i < args.length -2) { otherTreeFile = args[i+1]; i++; }
				else { System.err.println("Error: name of second tree file not specified"); displayHelp(); System.exit(1); }
			}
			// number of iterations, if desired
			else if (args[i].equals("-i")) {
				if (i < args.length -2) { numIter = Integer.valueOf(args[i+1]); i++; }
				else { System.err.println("Error: number of iterations not specified"); displayHelp(); System.exit(1); }
			}
			// output file
			else if (args[i].equals("-o")) {
				if (i < args.length -2) { outfile = args[i+1]; i++; }
				else { System.err.println("Error: output file not specified");  displayHelp(); System.exit(1); }
			}	
				
			// all other arguments.  Note we can have -vn
			else { 
				for (int j = 1; j<args[i].length(); j++) {
					switch(args[i].charAt(j)) {						
					// display help
					case 'h': displayHelp(); System.exit(0); break;
					// output number of orthants the geodesic goes through
					case 'r': orthants = true; break;
					case 'u': rooted = false; break;
					
					default: System.out.println("Illegal command line option.\n"); displayHelp(); System.exit(1); break;
					} // end switch
				} // end j loop (arguments without parameter)
			} // end parsing an individual argument
		}  // end for i (looping through arguments)
		
		
		// Read in trees from file.
		PhyloTree[] trees = PolyMain.readInTreesFromFile(treeFile,rooted);

		if (otherTreeFile != null) {
			otherTrees = PolyMain.readInTreesFromFile(otherTreeFile,rooted);
		}
	
		// Algorithms that don't write to outFile.
		if (algorithm.equals("random")) {
    		// XXX: Should check that numIter, etc make sense
			// outfile is used a prefix for a group of output files.
    		projToRandomGeodesics(trees,orthants,numIter,epsilon,outfile);
    		System.exit(0);
    	}
		
		// Algorithms that write to outfile.
		// Open output file
		PrintWriter outfileStream = null;

		try {
        	outfileStream = new PrintWriter(new FileWriter(outfile));
		
        	if (algorithm.equals("projection_indices")) {
        		Geodesic geo = getGeodesic(otherTrees[0], otherTrees[1],null);
        		double[] indices = getProjIndices(geo,trees,epsilon);
        		for (int i = 0; i < indices.length; i++ ) {
        			outfileStream.println(indices[i]);
        		}
        	}
        	else if (algorithm.equals("distances_to_geo")) {
        		Geodesic geo = getGeodesic(otherTrees[0], otherTrees[1],null);
        		for (PhyloTree t: trees) {
        			outfileStream.println(calcGeoDist(t,t.projectToGeo(geo, epsilon)));
        		}
        	}
        	else {
				System.out.println("Error:  no command specified.\n");
				System.exit(1);
        	}
		
        	if (outfileStream != null) {
    			outfileStream.close();
    		}
		} catch (FileNotFoundException e) {
			System.out.println("Error opening or writing to " +outfile + ": "+ e.getMessage());
			System.exit(1);
		} catch (IOException e) {
			System.out.println("Error opening or writing to " + outfile + ": " + e.getMessage());
			System.exit(1);
		}
		
		System.exit(0);		
	}
		
		
	/**  Choose 2 trees at random from trees, compute the geodesic between them,
	 *  compute the sum of square distances of each tree in trees, projected onto
	 *  this geodesic. 
	 *  Repeat this for numIter iterations, and for each iteration write the following
	 *  to a new file called outFile_<iteration number>:  
	 *  first tree index in trees, second tree index in trees, 
	 *  num of orthants geodesic passes through (if orthants is true), sum of square distances
	 * 
	 * @param trees
	 * @param orthants
	 * @param numIter
	 * @param outFile
	 */
	public static void projToRandomGeodesics(PhyloTree[] trees, Boolean orthants, int numIter, double epsilon,String outFile) {
		XORShiftRandom r = new XORShiftRandom();
		
		for (int i = 0; i < numIter; i++) {
		
			// print projected tree to file
			PrintWriter outputStream = null;
	  
			try {
				outputStream = new PrintWriter(new FileWriter(outFile + "_" + i));
		
				if (orthants) {
					outputStream.println("Index 1\tIndex 2\tNum Orthants\tGeo Distance\tScore");
				}
				else {
					outputStream.println("Index 1\tIndex 2\tGeo Distance\tScore");
				}
		
	       		// pick two trees at random
	       		int i1 = r.nextInt(trees.length);
	       		int i2 = r.nextInt(trees.length);
	       		while (i2 == i1) {
	       			i2 = r.nextInt(trees.length);
	       		}
	       		
	       		// print tree indices to file
	       		outputStream.print(i1 + "\t" + i2 + "\t");
			
	       		Geodesic geo  = getGeodesic(trees[i1],trees[i2],null);
	       		// print number of orthants geo passes through to file, if desired
	       		if (orthants) {
		       		// get number of orthants that the geodesic passes through
	       			outputStream.print(geo.numTopologies() + "\t");
	       		}
	       		//print geodesic distance to file
	       		outputStream.print(geo.getDist() + "\t");
	       		
	       		// print sum of squares distance to file
	       		outputStream.println(Tools.roundSigDigits(scoreProjsToGeodesic(geo,trees, epsilon),6));
			
		
	       		if (outputStream != null) {
	       			outputStream.close();
	       		}
			} 
			catch (FileNotFoundException e) {
				System.out.println("Error opening or writing to " + outFile + ": "+ e.getMessage());
	    		System.exit(1);
			} 
			catch (IOException e) {
				System.out.println("Error opening or writing to " + outFile + ": " + e.getMessage());
				System.exit(1);
			} // end try/catch
		} // end for
	}	// end method
		
	
	/** Computes the sum of square distances of the trees in trees to
	 * their projections onto the geodesic between trees t1 and t2.
	 * epsilon is the accuracy to which we compute the projected tree.
	 * 
	 * @param t1
	 * @param t2
	 * @param trees
	 * @param epsilon
	 * @return
	 */
	public static double scoreProjsToGeodesic (Geodesic geo, PhyloTree[] trees, double epsilon) {
		double score = 0; 

		for (PhyloTree t : trees ) {
			score = score + Math.pow(calcGeoDist(t,t.projectToGeo(geo, epsilon)),2);
		}
		return score;
	}
	
	
	/** Computes the indices of the projections of tree onto the geodesic geo.
	 * i.e. for tree i returns the number between 0 and 1 representing the location of the projection of tree i
	 * onto the geodesic geo (geodesic goes from 0 to 1)
	 * 
	 * @param geo
	 * @param trees
	 * @param epsilon
	 * @return
	 */
	public static double[] getProjIndices (Geodesic geo, PhyloTree[] trees, double epsilon) {
		double[] indices = new double[trees.length]; 
		
		// for each tree, projection onto geo and store index
		for (int i = 0; i < trees.length; i++) {
			indices[i] = trees[i].projectToGeoIndex(geo,epsilon);
 		}
		
		return indices;
	}
}
