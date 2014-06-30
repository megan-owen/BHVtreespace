package compGeo;

import polyAlg.PolyMain;
import distanceAlg1.PhyloTree;

import java.util.*;

import distanceAlg1.Geodesic;

public class CompGeoMain {

	/** Returns the tree representing the weighted combination of the elements
	 * of trees, as given by s.
	 * s has the following form: "t_1 w_1 t_2 w_2 t_3 ... w_{i-1} t_i"
	 * and corresponds to 
	 * (1-{w_{i-1})...((1-w_2)*((1 - w_1)*T_{t_1} + w_1*T_{t_2}) + w_2*T_{t_3}) ... + w_{i-1}*T_{t_i} 
	 * where (1-w)*T_{t_1} + w*T_{t_2}, 0<= w <= 1 means the tree w along the geodesic 
	 * from T_{t_1} to T_{t_2} 
	 *
	 * Note:  this doesn't give all weighted combinations.  In particular, it misses
	 * ones like (T_1 + T_2) + (T_3 + T_4)
	 * @param trees
	 * @param s
	 * @return
	 */
	public static PhyloTree restrictedWeightedComb(PhyloTree[] trees, String s) {
		PhyloTree comb;
		
		// split up string s
		String[] tokens = s.split(" ",-1);	
		
		System.out.println("s split up is ");
		for(int i = 0;i < tokens.length; i++) {
			System.out.println("i = " + i + ": " + tokens[i]);
		}
		return trees[0];
	}
	
	
	/** Return the trees on the pairwise geodesics between the inputTrees
	 *  that are on the boundaries of orthants.
	 *  Each tree in the boundary should only appear once.
	 * @param inputTrees
	 * @return
	 * @throws Exception 
	 */
	public static Vector<PhyloTree> getBoundaryTrees (PhyloTree[] inputTrees) throws Exception {
		Vector<PhyloTree> boundaryTrees = new Vector<PhyloTree>();
		
		for (int i = 0; i < inputTrees.length; i++) {
			for (int j = i + 1; j < inputTrees.length; j++) {
				Vector<PhyloTree> treesToAdd = Geodesic.getBoundaryTrees(inputTrees[i],inputTrees[j]);
				// don't add any duplicate trees
				treesToAdd.removeAll(boundaryTrees);
				boundaryTrees.addAll(treesToAdd);
			}
		}
		return boundaryTrees;
	}
	
	
	/** Same function as above, but takes in a vector instead of an array.
	 * 
	 * @param inputTrees
	 * @return
	 * @throws Exception 
	 */
	public static Vector<PhyloTree> getBoundaryTrees (Vector<PhyloTree> inputTrees) throws Exception {
		Vector<PhyloTree> boundaryTrees = new Vector<PhyloTree>();
		
		for (int i = 0; i < inputTrees.size(); i++) {
			for (int j = i + 1; j < inputTrees.size(); j++) {
				Vector<PhyloTree> treesToAdd = Geodesic.getBoundaryTrees(inputTrees.get(i),inputTrees.get(j));
				// don't add any duplicate trees
				treesToAdd.removeAll(boundaryTrees);
				boundaryTrees.addAll(treesToAdd);
			}
		}
		return boundaryTrees;
	}
	
	
	/** Returns the original set of trees, plus all trees on the intersection of geodesics between them
	 *  and orthant boundaries.  All of these trees make up the set of input trees for the next repetition.
	 * @param inputTrees
	 * @param numReps
	 * @return
	 * @throws Exception 
	 */
	public static Vector<PhyloTree> getExtremalTrees(Vector<PhyloTree> inputTrees, int numReps ) throws Exception {
		Vector<PhyloTree> boundaryTrees = inputTrees;
		Vector<PhyloTree> newBoundaryTrees;
		
		Boolean newTrees = true;
		
		int i = 0;
		while ((i < numReps) && newTrees)  {
			System.out.println("boundaryTrees are: " + boundaryTrees);
			newBoundaryTrees = getBoundaryTrees(boundaryTrees);
			System.out.println("newBoundaryTrees are: " + newBoundaryTrees);
			newBoundaryTrees.removeAll(boundaryTrees);
			if (newBoundaryTrees.size() ==0) {
				newTrees = false;
			}
			else {
				boundaryTrees.addAll(newBoundaryTrees);
			}
			i++;
		}
		return boundaryTrees; 
	}
	

	
	
	
	/** Help message (ie. which arguments can be used, etc.)
	 * 
	 */
	public static void displayHelp() {
		System.out.println("Command line syntax:");
		System.out.println("java -jar analysis.jar [options] treefile");
		System.out.println("Optional arguments:");	
		System.out.println("\t -a <algorithm> \t specifies what to compute");
		System.out.println("\t -f <otherTreeFile> \t reads in an additional tree file");
		System.out.println("\t -o <outfile> \t store the output in the file <outfile>.  Default is output.txt");
		System.out.println("\t -s <string> \t specifies the string needed for algorithm weighted_comb");
		System.out.println("\t -u \t trees are unrooted. Default is trees are rooted.");
		System.out.println("\n");
		System.out.println("Algorithms are:");
		System.out.println("\t weighted_comb \t computes the weighted combination of the input trees.");
		System.out.println("\t\t Specify the combination as a string \"T_1 w_1 T_2 w_2 T_3 ... w_{i-1} T_i\"");
		System.out.println("\t\t which corresponds to");
		System.out.println("\t\t(1-{w_{i-1})...((1-w_2)*((1 - w_1)*T_1 + w_1*T_2) + w_2*T_3) ... + w_{i-1}*T_i");
		System.out.println("\t\twhere (1-w)*T_1 + w*T_2, 0<= w <= 1 means the tree w along the geodesic from T_1 to T_2");
	}
	
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String algorithm = "";
		String otherTreeFile = null;
		String treeFile;
		String outfile = "output.txt";
		String s = "";
		Boolean rooted = true;
		PhyloTree [] otherTrees;
		int numIter = 0;	// specify number of iterations
		
		/* Parse command line arguments */

		if (args.length < 1) {
			System.out.println("Error: Missing input file name"); displayHelp(); throw new Exception();
		}
		treeFile = args[args.length-1];
		for (int i = 0; i < args.length - 1; i++) {
			
			if (!args[i].startsWith("-")) { System.out.println("Invalid command line option"); displayHelp(); throw new Exception(); }
			// other tree file, if desired
			else if (args[i].equals("-a")) {
				if (i < args.length -2) { algorithm = args[i+1]; i++; }
				else { System.err.println("Error: algorithm not specified"); displayHelp(); throw new Exception(); }
			}
			// other tree file, if desired
			else if (args[i].equals("-f")) {
				if (i < args.length -2) { otherTreeFile = args[i+1]; i++; }
				else { System.err.println("Error: name of other tree file not specified"); displayHelp(); throw new Exception(); }
			}
			// number of iterations
			else if (args[i].equals("-i")) {
				if (i < args.length -2) { numIter = Integer.valueOf(args[i+1]); i++; }
				else { displayHelp(); throw new Exception(); }
			}
			// output file
			else if (args[i].equals("-o")) {
				if (i < args.length -2) { outfile = args[i+1]; i++; }
				else { System.err.println("Error: output file not specified");  displayHelp(); throw new Exception(); }
			}
			// string, if needed
			else if (args[i].equals("-s")) {
				if (i < args.length -2) { s = args[i+1]; i++; }
				else { System.err.println("Error: string not specified");  displayHelp(); throw new Exception(); }
			}
				
			// all other arguments.  Note we can have -vn
			else { 
				for (int j = 1; j<args[i].length(); j++) {
					switch(args[i].charAt(j)) {						
					// display help
					case 'h': displayHelp(); System.out.println("Finished"); break;
					case 'u': rooted = false; break;
					
					default: System.out.println("Illegal command line option.\n"); displayHelp(); throw new Exception();
					} // end switch
				} // end j loop (arguments without parameter)
			} // end parsing an individual argument
		}  // end for i (looping through arguments)
		
		/* Read in the trees  */
		
		PhyloTree[] trees = PolyMain.readInTreesFromFile(treeFile,rooted);
		
		if (otherTreeFile != null) {
			otherTrees = PolyMain.readInTreesFromFile(otherTreeFile,rooted);
		}
		
		if (algorithm.equals("weighted_comb")) {
			// otherTreeFile should contain a single tree
			// treeFile contains the list of trees to compute sum of square distance to
			PhyloTree weighted_comb = restrictedWeightedComb(trees,s);
			
			System.out.println("Weighted combination tree is " + weighted_comb);
			
		}
		else if (algorithm.equals("boundary_trees")) {
			Vector<PhyloTree> inputTrees = new Vector<PhyloTree>(Arrays.asList(trees));
			Vector<PhyloTree> boundaryTrees = getExtremalTrees(inputTrees,numIter);
			
			for (PhyloTree tree: boundaryTrees) {
				System.out.println("" + tree.getNewick(true));
			}
		}

	}
	
	

}
