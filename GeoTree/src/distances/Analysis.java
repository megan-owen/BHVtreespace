package distances;

import polyAlg.PolyMain;
import static polyAlg.PolyMain.getGeodesic;
import static polyAlg.PolyMain.calcGeoDist;
import static centroid.CentroidMain.getCentroidViaRandPermCauchy;
import polyAlg.Tools;
import distanceAlg1.PhyloTree;
import distanceAlg1.Geodesic;
import distanceAlg1.PhyloTreeEdge;
import distanceAlg1.EdgeAttribute;
import distanceAlg1.Bipartition;
import distances.Presentation;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.*;

import centroid.*;

public class Analysis {
	
	public static int verbose = 0;
	
	/**  Returns the average geodesic distance for all the geodesics in geos.
	 * 
	 * @param geos
	 * @return
	 */
	public static double averageGeoDist(Geodesic[][] geos) {
		int numTrees = geos.length;
		double average = 0;
		int numDists = 0;
		
		for (int i = 0; i < numTrees - 1; i++) {
			for (int j = i + 1; j < numTrees; j++) { 
				average = average + (geos[i][j]).getDist();
				numDists++;
			}
		}
		return average / numDists;
	}
	
	/** Returns the average distance of the trees from the origin.
	 *  i.e. the average lengths/norms of the trees
	 * @param trees
	 * @return
	 */
	public static double avgDistanceFromOrigin(PhyloTree [] trees) {
		int numTrees = trees.length;
		double avg = 0;
		
		for(int i = 0; i < trees.length; i++) {
			avg = avg + trees[i].getDistanceFromOrigin();
		}
		return avg/numTrees;
	}
	
	
	/** Given a set of trees, removes all edges common to these trees, and returns the new trees.
	 *  TODO:  Finish method
	 *  
	 * @param inTrees
	 * @return
	 */
	public static PhyloTree[] removeCommonSplits(PhyloTree[] inTrees) {
		PhyloTree[] outTrees = new PhyloTree[inTrees.length];
		
		/* Fill in  */  // TODO
		Vector<Bipartition> commonSplits = getCommonSplits(inTrees);
		
		if (commonSplits == null) {
			return inTrees;
		}
		
		if (commonSplits.size() == 0) {
			return inTrees;
		}
		
		for (int i = 0; i < inTrees.length;i++) {
			outTrees[i] = new PhyloTree(inTrees[i]);
			outTrees[i].removeSplits(commonSplits);
		}
		
		return outTrees;
	}	
		
	
	/** Returns a vector of all edges common to the trees in trees.  Returns null if no common edges.
	 * 
	 * @param trees
	 * @return
	 */
		
	// XXX: This method assumes that there are no 0 length edges.  However, this is currently (June 2011)
	// not compatible with the rest of the code (but shouldn't cause a problem in most cases.)
	// TODO:  assumes PhyloTree contains the same number of trees as its length
	public static Vector<Bipartition> getCommonSplits(PhyloTree[] trees) {
		Vector<Bipartition> commonSplits = new Vector<Bipartition>();
				
		if (trees == null) {
			System.err.println("No trees for finding common edges: returning null.");
			return null;
		}
		
		if (trees[0] == null) {
			System.err.println("No trees for finding common edges:  returning null");
			return null;
		}
		
		commonSplits = trees[0].getSplits();
		
		if (trees.length ==1) {
			System.out.println("One tree:  Returning all its edges as common");
			return commonSplits; 
		}
		
		for (int i = 1; i < trees.length; i++) {
			commonSplits.retainAll(trees[i].getSplits());
		}
		return commonSplits;
	}
		
		
		
	/**  Do numRep bootstraps of finding the centroid of a set of trees sampled randomly with replacement from trees.
	 * 	 Do at most numIter iterations of Sturm's algorithm to determine the centroid at each iteration.
	 *   If the algorithm hasn't converged by then, centroid will be null;
	 * @param trees
	 * @param numRep
	 * @param numIter
	 * @return
	 */
	public static PhyloTree[] bootstrap(PhyloTree[] trees, int numRep, long numIter) {
		PhyloTree[] centroids = new PhyloTree[numRep];
		int numTrees = trees.length;
		
		PhyloTree[] sampledTrees = new PhyloTree[numTrees];
		
		int counter = 0;  //  remove all statements with // next to them.
		
		for (int i = 0; i< numRep; i++) {
			sampledTrees = sampleWithReplacement(trees);
			
			double epsilon = CentroidMain.getEpsilon(sampledTrees,5000);
			
			centroids[i] = CentroidMain.getCentroidViaRandomCauchy(sampledTrees, numIter,5, epsilon, null,0,0);	
			
			double var = CentroidMain.variance(centroids[i],sampledTrees);  //
			
			PhyloTree starTree = getStarTree(trees); //  remove all statements with // next to them.
			
			double starVar = CentroidMain.variance(starTree, sampledTrees);  //
			
			if (var <= starVar) {  //
				System.out.println("------------------------------------  STAR NOT THE MEAN ---------------- ");  //
				System.out.println("Star is " + starTree.getNewick(true));  //
				System.out.println("Centroid is " + centroids[i].getNewick(true));   //
				System.out.println("Star tree variance is " + starVar);    //
				System.out.println("Centroid variance is " + var);    //
				
				double dist = calcGeoDist(centroids[i], starTree); //
				System.out.println("Distance to star tree is " + dist);  //
				counter++;
			}    //
		}
		
		System.out.println(counter + " bootstraps where star didn't have a lower variance");
		
		return centroids;
	}
	
	
	
	/**  Returns the proportion of geodesics that are the cone path.
	 * 
	 * @param geos
	 * @return
	 */
	public static double proportionConePath(Geodesic[][] geos) {
		int numTrees = geos.length;
		int numConePath = 0;
		int numDists = 0;
		
		for (int i = 0; i < numTrees - 1; i++) {
			for (int j = i + 1; j < numTrees; j++) { 
				if (geos[i][j].getRS().size() == 1) {
					numConePath++;
				}
				numDists++;
			}
		}
		return (double) numConePath / numDists;
	}
	

	public static void geodesicStats(PhyloTree[] trees) {
		int numTrees = trees.length;
		double averageDist = 0;
		int numDists = 0;
		int numLeaves = trees[0].getLeaf2NumMap().size();
		
		double[][] dist = new double[numTrees][numTrees];
		
		/*  Initialize an array to hold the number of geodesics with the given number of ratios in its ratio sequence. */
		int[] numRSOfEachSize = new int[numLeaves];
		Arrays.fill(numRSOfEachSize,0);
		
		/* Initialize an array to hold the number of geodesics with the given number of common edges. */
		int[] numCommonEdgesOverGeos = new int[numLeaves];
		Arrays.fill(numCommonEdgesOverGeos,0);
		
		for (int i = 0; i < numTrees - 1; i++) {
			for (int j = i + 1; j < numTrees; j++) { 
				Geodesic geo = PolyMain.getGeodesic(trees[i], trees[j], null);
				dist[i][j] = geo.getDist();
				
				numRSOfEachSize[geo.getRS().size()]++;
				numCommonEdgesOverGeos[geo.numCommonEdges()]++;
				
				averageDist = averageDist + geo.getDist();
				numDists++;
			}
		}
		
		double[] starDist = distancesToStarTree(trees);
		
				
		printDistancesToFile(dist,starDist, "montpellier_analysis0to19_2.txt");
		
		System.out.println("Average length of geodesic is " + averageDist/numDists);
		
		for (int i = 1; i < numTrees; i++) {
			System.out.println("# of geodesics with " + i + " support pairs: " + numRSOfEachSize[i]);
		}
		
		System.out.println();
		
		for (int i = 0; i < numTrees; i++) {
			System.out.println("# of geodesics with " + i + " common edges: " + numCommonEdgesOverGeos[i]);
		}
	}
	
	
	
	/**  Returns the star tree with averaged leaf edge lengths.
	 * 
	 * @param trees
	 * @return
	 */
	public static PhyloTree getStarTree(PhyloTree[] trees) {
		if (trees == null) {
			System.err.println("Warning:  tree set null when getting star tree; returning null");
			return null;
		}
		
		if (trees[0] == null) {
			System.err.println("Warning:  getting star tree for empty set of trees; returning null");
			return null;
		}
		
		int numLeaves = trees[0].getLeaf2NumMap().size();
		int numTrees = trees.length;
		
		Vector<PhyloTreeEdge> starEdges = new Vector<PhyloTreeEdge>();
		
		PhyloTree star = new PhyloTree(starEdges,trees[0].getLeaf2NumMap(), trees[0].isRooted()); 
		
		EdgeAttribute[] starLeafAttribs = new EdgeAttribute[numLeaves];
		
		
		for (int leaf= 0; leaf < numLeaves; leaf++) {
			starLeafAttribs[leaf] = trees[0].getLeafEdgeAttribs()[leaf].clone();
			for (int tree = 1; tree < numTrees; tree++) {
				starLeafAttribs[leaf] = EdgeAttribute.add(starLeafAttribs[leaf], trees[tree].getLeafEdgeAttribs()[leaf].clone());
			}
			starLeafAttribs[leaf].scaleBy(1.0/numTrees);
		}
		star.setLeafEdgeAttribs(starLeafAttribs);
		
		return star;
	}
	
	
	public static double[] distancesToStarTree(PhyloTree[] trees) {
		int numTrees = trees.length;
		
		PhyloTree star = getStarTree(trees);
		
		double[] dist = new double[numTrees];
		
		for (int i = 0; i < numTrees; i++) {
			dist[i] = calcGeoDist(star,trees[i]);
		}
		return dist;
	}
	
	
	
	public static void printDistancesToFile(double[][] dist, double[] specialDist, String fileName) {
		int numTrees = dist.length;
		
		PrintWriter outputStream = null;
  
        /* Outputs the distances in a matrix, with 0 entries on the diagonals, and the rows labelled from 1 to numTrees. */
        try {
        	outputStream = new PrintWriter(new FileWriter(fileName));
 
    		for (int i = 0; i < numTrees; i++) {
    			outputStream.print(i);
    			
    			if (i != 0) {
    				outputStream.print("\t" + specialDist[i-1]);
    			}
    			
    			for (int j = 1; j< i; j++) {
    				outputStream.print("\t" + dist[j-1][i-1]);
    			}
    			outputStream.println("\t0");
    		}
    		if (outputStream != null) {
    			outputStream.close();
    		}
        } catch (FileNotFoundException e) {
        	System.out.println("Error opening or writing to " +fileName + ": "+ e.getMessage());
        	System.exit(1);
        } catch (IOException e) {
        	System.out.println("Error opening or writing to " + fileName + ": " + e.getMessage());
        	System.exit(1);
        }
	}
	

	
	public static void writeLengths(PhyloTree[] trees, PrintWriter out) {		
    		for (int i = 0; i < trees.length; i++) {
    			out.println(trees[i].getDistanceFromOrigin());
    		}
    		if (out != null) {
    			out.close();
    		}
	}
	
	
	/** Returns an array of the trees of the same size as trees, with the entries drawn uniformly from trees with replacement.
	 * 
	 * 
	 * @param trees
	 * @return
	 */
	public static PhyloTree[] sampleWithReplacement(PhyloTree[] trees) {
		int numTrees = trees.length;
		PhyloTree[] sampledTrees = new PhyloTree[numTrees];
		
		XORShiftRandom rng = new XORShiftRandom();
		
		
		for (int i = 0; i < numTrees; i++) {
			sampledTrees[i] = trees[rng.nextInt(numTrees)].clone();
		}
		
		return sampledTrees;		
	}
	
	
	public static void main(String[] args) {
		String outfile = "output.txt";
		String treeFile = "";
		String otherTreeFile = null;
		Boolean rooted = true;
		String algorithm = "";
		PhyloTree[] otherTrees = null;
		double epsilon = -1;
		double samplePt = -1;
		
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
				else { System.err.println("Error: name of other tree file not specified"); displayHelp(); System.exit(1); }
			}
			// output file
			else if (args[i].equals("-o")) {
				if (i < args.length -2) { outfile = args[i+1]; i++; }
				else { System.err.println("Error: output file not specified");  displayHelp(); System.exit(1); }
			}	
			// sample point on geodesic
			else if (args[i].equals("-s")) {
				if (i < args.length -2) { samplePt = Double.valueOf(args[i+1]); i++; }
			}
			// all other arguments.  Note we can have -vn
			else { 
				for (int j = 1; j<args[i].length(); j++) {
					switch(args[i].charAt(j)) {						
					// display help
					case 'h': displayHelp(); System.exit(0); break;
					case 'u': rooted = false; break;
					
					default: System.out.println("Illegal command line option.\n"); displayHelp(); System.exit(1); break;
					} // end switch
				} // end j loop (arguments without parameter)
			} // end parsing an individual argument
		}  // end for i (looping through arguments)
		
		/* Read in the trees  */
		if (algorithm.equals("verify_treefile")) {
			PolyMain.readInTreesFromFile(treeFile, rooted);
			System.out.println("Done reading in trees from file " + treeFile);
			System.exit(0);
		}
		
		
		PhyloTree[] trees = PolyMain.readInTreesFromFile(treeFile,rooted);
		
		if (otherTreeFile != null) {
			otherTrees = PolyMain.readInTreesFromFile(otherTreeFile,rooted);
		}
		
		// Algorithms that don't write to an output file
		if (algorithm.equals("compute_SSD")) {
			// otherTreeFile should contain a single tree
			// treeFile contains the list of trees to compute sum of square distance to
			if (otherTreeFile == null) {
				System.out.println("Error:  second tree file not specified");
				System.exit(1);
			}
				
			double ssd = CentroidMain.sumOfSquareDist(otherTrees[0], trees);
			System.out.println("Sum of squares distance is " + ssd);
			
			System.exit(0);
		}

		// Algorithms that write to an outfile	
		// Open output file
		PrintWriter outfileStream = null;

		try {
        	outfileStream = new PrintWriter(new FileWriter(outfile));
		
    		if (algorithm.equals("gtp_twofiles")) {
    			// otherTreeFile contains one list of trees to compare against the trees in treeFile
    			if (otherTreeFile == null) {
    				System.out.println("Error:  second tree file not specified");
    				System.exit(1);
    			}
    			
    			computeAllGeodesicsBtwLists(trees,otherTrees,outfile);
    			
    			System.exit(0);
    		}
    		
    		else if (algorithm.equals("normalize")) {
    			for(int i = 0; i < trees.length; i++) {
    				trees[i].normalize();
    				outfileStream.println(trees[i].getNewick(true));
    			}
    		}
    		
    		//  Compute projection of first tree in otherTreeFile onto the geodesic between
    		// the first two trees in treeFile
    		else if (algorithm.equals("project")) {
    			projectTreeToGeo(trees,otherTrees,epsilon,outfile);
    			System.exit(0);
    		}
    		// Compute the distance of the projection of the first tree in otherTreeFile
    		// on the geodesic between the first two trees in treeFile
    		else if ( algorithm.equals("project_distance")) {
    			Geodesic geo = getGeodesic(trees[0],trees[1],null);
    			PhyloTree projection = otherTrees[0].projectToGeo(geo,epsilon);
    			double dist = calcGeoDist(otherTrees[0],projection);
    			Presentation.printStringToFile(""+ dist, outfile);
    			System.exit(0);
    		}
    		// Writes the lengths of the trees in treefile to outfile, 
    		// one per line
    		else if (algorithm.equals("lengths")) {
    			writeLengths(trees, outfileStream);
    			System.exit(0);
    		}
    		
    		//  Nicely prints out the splits that differ in the trees in trees from the trees in otherTrees
    		else if (algorithm.equals("diff_splits")) {
    			printIncompatibleSplits(trees,otherTrees,outfile, false);
    			System.exit(0);
    		}
    		
    	//  Nicely prints out the splits that differ between tree i in trees and tree i in otherTrees
    		else if (algorithm.equals("diff_splits_paired")) {
    			printIncompatibleSplits(trees,otherTrees,outfile, true);
    			System.exit(0);
    		}
    		
    		// returns the tree at point samplePt along the geodesic starting at tree 0 in treefile
    		// and ending at tree 1 in treefile.
    		// samplePt must be between 0 and 1
    		else if (algorithm.equals("sample_point")) {
    			if ((samplePt < 0) || (samplePt > 1)) {
    				System.err.println("Error: sample point is either missing or not between 0 and 1");
    				System.exit(1);
    			}
    			if (trees.length < 2) {
    				System.err.println("Error: need two trees in treefile");
    				System.exit(1);
    			}
    			PhyloTree sample = (getGeodesic(trees[0],trees[1],null)).getTreeAt(samplePt,trees[0].getLeaf2NumMap(),rooted); 
    			Presentation.printStringToFile(sample.getNewick(true), outfile);
    			System.exit(0);
    		}
    		// returns a count of the different topologies in the tree file
    		else if (algorithm.equals("topology_count")) {
    			
    			countTopologies(trees, outfileStream);
    			System.exit(0);
    		}
    		// returns a count of the different splits in the tree file
    		else if (algorithm.equals("split_count")) {
    			
    			countSplits(trees, outfileStream);
    			System.exit(0);
    		}
    		// computes the log map coordinates for each tree in trees,
    		// relative to the first tree in otherTrees
    		else if (algorithm.equals("log_map")) {
    			computeAllLogMaps(trees,otherTrees[0],outfileStream);
    		}
    		// returns 4 lines to the output file corresponding to the 4 angles
    		else if (algorithm.equals("endray_angles")) {
    			endRayAngles(trees, otherTrees, outfileStream);
    			System.exit(0);
    		}
    		else {
    			System.out.println("Error:  no algorithm specified.\n");
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
	
	/** Help message (ie. which arguments can be used, etc.)
	 *  XXX:  rename sample_point to getTreeAt
	 */
	public static void displayHelp() {
		System.out.println("Command line syntax:");
		System.out.println("java -jar analysis.jar [options] treefile");
		System.out.println("Optional arguments:");	
		System.out.println("\t -a <algorithm> \t specifies what to compute");
		System.out.println("\t -e <epsilon> \t specifies epsilon, if needed.");
		System.out.println("\t -f <otherTreeFile> \t reads in an additional tree file");
		System.out.println("\t -o <outfile> \t store the output in the file <outfile>.  Default is output.txt");
		System.out.println("\t -s <sample point> \t a number between 0 and 1, inclusive.");
		System.out.println("\t -u \t trees are unrooted. Default is trees are rooted.");
		System.out.println("\n");
		System.out.println("Algorithms are:");
		System.out.println("\t compute_SSD \t computes the sum of square distance from the tree in otherTreeFile to the trees in treefile");
		System.out.println("\t gtp_twofiles \t computes the geodesic distance between all trees in treefile and all trees in otherTreeFile");
		System.out.println("\t verify_treefile \t verifies that all trees in the treefile can be read in without errors");
		System.out.println("\t project \t projects the first tree in otherTreeFile onto the geodesic between the first two trees in treefile");
		System.out.println("\t project_distance \t projects the first tree in otherTreeFile onto the geodesic between the first two trees in treefile and outputs the distance");
		System.out.println("\t lengths \t computes the lengths (norm of edge vectors) of all trees in treefile");
		System.out.println("\t diff_splits \t nicely prints out (to the output file) the splits that differ between trees in treefile and trees in otherTreeFile");
		System.out.println("\t diff_splits_paired \t nicely prints out (to the output file) the splits the differ between tree i in treefile and tree i in otherTreeFile");
		System.out.println("\t sample_point \t returns the tree on the geodesic between the trees in treefile at <sample point>");
		System.out.println("\t topology_count \t returns a file containing information about the topologies in treefile");
		System.out.println("\t split_count \t returns a file containing information about the splits appearing in the trees in treefile");
		System.out.println("\t log_map \t returns a file containing the log map coordinates of the trees in treefile relative to the first tree in otherTreeFile (centre tree)");
		System.out.println("\t endray_angles \t returns a file with four lines corresponding to the angles made by the endrays of the two geodesics given (one in treefile another in otherTreeFile)");
	}
	
	
	/** Counts the different topologies and writes the info to a file.
	 * 
	 * @param trees
	 * @param outfile
	 */

	public static void countTopologies(PhyloTree[] trees, PrintWriter out) throws IOException {
		// for each distinct topology, we want to store the number of trees 
		// with this topology and their indices
		//
		// In choosing the data structure, we have two concerns:
		// 1) being able to search and find each topology as we encounter it
		// 2) when we are done, sort the list by number of trees of each topology
		ArrayList<ArrayList> topologies = new ArrayList<ArrayList>(trees.length);
		Boolean foundTopology = false;
		int[] topNumOfTree = new int[trees.length];		// store topology number of each tree, to output at end
														// of human readable part for use by computer
		
		// for each tree:
		for (int i = 0; i < trees.length; i++) {
			foundTopology = false;
			// search through already encountered topologies to see if it matches any of them
			for (int j = 0; j < topologies.size(); j++ ) {
				ArrayList entry = (ArrayList) topologies.get(j);
				if (trees[i].hasSameTopology( (PhyloTree) entry.get(0) )) {
					// add in this index ( + 1 to get file row)
					entry.add(i + 1);
					foundTopology = true;
					break;	// we are done
				}
			}
			// if we haven't found the topology, create an entry
			if (!foundTopology) {
				ArrayList newEntry = new ArrayList();
				newEntry.add(trees[i]);
				newEntry.add(i + 1);   // ( + 1 to get file row)
				topologies.add(newEntry);
			}
		}
		
		// now we want to write the data structure to a file

 
		// put a vector of the tree counts in the first line of the file
		String vCounts = "Raw topology counts: ";
		String topCounts = ""; 	// at the end of the file, put the topology of the tree on that line
		
        int counter = 1;
    	while(topologies.size() > 0) {
    		int maxLength = ((ArrayList) topologies.get(0)).size();
    		int maxIndex = 0;
    		for (int i = 1; i < topologies.size(); i++) {
    			int thisLen = ((ArrayList) topologies.get(i)).size();
    			if (thisLen > maxLength) {
    				maxLength = thisLen;
    				maxIndex = i;
    			}	
    		}
    		// now we have found the max length and index, print this entry
    		// to the file and delete from topologies array
    		ArrayList entry = (ArrayList) topologies.remove(maxIndex);
    		out.println(counter + ". " + ( (PhyloTree) entry.get(0)).getNewick(false));
    		String s = "";
    		if (maxLength == 2) {
    			s = "\t" + (maxLength-1) + " tree at file row: ";
    		}
    		else {
    			s = "\t" + (maxLength-1) + " trees at file rows: ";
    		}
    		vCounts = vCounts + " " + (maxLength-1);
    		for (int j = 1; j < entry.size(); j++) {
    			s = s + entry.get(j) + ",";
    			// store entry in appropriate box in topNumOfTree array
    			topNumOfTree[(Integer) entry.get(j) - 1] = counter;  			
    		}
    		s = s.substring(0,s.length()-1);
    		out.println(s);
    		counter++;
    	}
    	out.println(vCounts);
    	
    	for (int i = 0; i < topNumOfTree.length; i++) {
    		out.println(topNumOfTree[i]);
    	}
    		
    	if (out != null) {
    		out.close();
    	}
	}
	
	public static void countSplits(PhyloTree[] trees, PrintWriter out) throws IOException {
		// for each distinct split, we want to store the number of trees 
		// with this split and their indices
		//
		// In choosing the data structure, we have two concerns:
		// 1) being able to search and find each split as we encounter it
		// 2) when we are done, sort the list by number of trees of each split
		ArrayList<ArrayList> splits = new ArrayList<ArrayList>(trees.length);
		Boolean foundSplit = false;
		
		// for each tree:
		for (int i = 0; i < trees.length; i++) {
			Vector<Bipartition> treeSplitSet = trees[i].getSplits();
			// for each split in the tree
			for (int k = 0; k < treeSplitSet.size(); k++ ) {
				foundSplit = false;
				// search through already encountered splits to see if it matches any of them
				for (int j = 0; j < splits.size(); j++ ) {
					ArrayList entry = (ArrayList) splits.get(j);
					if (((Bipartition) treeSplitSet.get(k)).equals((Bipartition) entry.get(0) )) {
						// add in this index ( + 1 to get file row)
						entry.add(i + 1);
						foundSplit = true;
						break;	// we are done
					}
				}
				// if we haven't found the split, create an entry
				if (!foundSplit) {
					ArrayList newEntry = new ArrayList();
					newEntry.add(treeSplitSet.get(k));
					newEntry.add(i + 1);   // ( + 1 to get file row)
					splits.add(newEntry);
				} 
			}
		}
 
		// put a vector of the tree counts in the first line of the file
		String v_counts = "Raw split counts: ";
		
        int counter = 1;
    	while(splits.size() > 0) {
    		int maxLength = ((ArrayList) splits.get(0)).size();
    		int maxIndex = 0;
    		for (int i = 1; i < splits.size(); i++) {
    			int thisLen = ((ArrayList) splits.get(i)).size();
    			if (thisLen > maxLength) {
    				maxLength = thisLen;
    				maxIndex = i;
    			}	
    		}
    		// now we have found the max length and index, print this entry
    		// to the file and delete from splits array
    		ArrayList entry = (ArrayList) splits.remove(maxIndex);
    		out.println(counter + ". " + Bipartition.toStringVerbose( ((Bipartition) entry.get(0)).getPartition(), trees[0].getLeaf2NumMap()));
    		String s = "";
    		if (maxLength == 2) {
    			s = "\t" + (maxLength-1) + " split in tree at file row: ";
    		}
    		else {
    			s = "\t" + (maxLength-1) + " splits in trees at file rows: ";
    		}
    		v_counts = v_counts + " " + (maxLength-1);
    		for (int j = 1; j < entry.size(); j++) {
    			s = s + entry.get(j) + ",";
    		}
    		s = s.substring(0,s.length()-1);
    		out.println(s);
    		counter++;
    	}
    	out.println(v_counts);

    	if (out != null) {
    		out.close();
    	}	
	}
	
	/**  Computes all geodesics between the two lists.
	 *   Prints the results to file.
	 * 
	 * @param trees
	 * @param otherTrees
	 * @param rooted
	 * @param outfile
	 */
	public static void computeAllGeodesicsBtwLists(PhyloTree[] trees,PhyloTree[] otherTrees,String outFileName) {
		int numTrees = trees.length;
		int numOtherTrees = otherTrees.length;
	    
	    // print distances to file
	    PrintWriter outputStream = null;
	  
	    // Outputs the distances in a column, with the first two columns being the trees numbers and the third
	    // number the geodesic distance between those trees
	    try {
	       	outputStream = new PrintWriter(new FileWriter(outFileName));
	 
	    	for (int i = 0; i < numTrees ; i++) {
	    		for (int j = 0; j< numOtherTrees; j++) {
	    			double dist = PolyMain.getGeodesic(trees[i], otherTrees[j], "geo_" + i + "_" + j).getDist();
	    			outputStream.println(i + "\t" + j + "\t" + Tools.roundSigDigits(dist, 6));
				}
				outputStream.println();
			}
			if (outputStream != null) {
	            outputStream.close();
	        }
	    } catch (FileNotFoundException e) {
	        System.out.println("Error opening or writing to " + outFileName + ": "+ e.getMessage());
	        System.exit(1);
	    } catch (IOException e) {
	    	System.out.println("Error opening or writing to " + outFileName + ": " + e.getMessage());
	    	System.exit(1);
	    }
	}
	
	/** Wrapper for code that projects one tree onto the geodesic between two others.
	 * Computes the projection of otherTrees[0] onto the geodesic between trees[0] and trees[1].
	 * Writes the projected tree to outFileName.
	 * 
	 * @param trees
	 * @param otherTrees
	 * @param outFileName
	 */
	public static void projectTreeToGeo(PhyloTree[] trees, PhyloTree[] otherTrees, double epsilon, String outFileName) {
		if ((trees.length < 2) || (otherTrees.length < 1)) {
			System.out.println("Error: not enough trees in input files");
			System.exit(1);
		}
		
		if (epsilon <= 0) {
			epsilon = 0.05;
			System.out.println("Using epsilon = 0.05 when projecting tree onto geodesic");
		}
		
		PhyloTree treeToProject = otherTrees[0];
		PhyloTree t1 = trees[0];
		PhyloTree t2 = trees[1];
		
		PhyloTree projectedTree = treeToProject.projectToGeo(getGeodesic(t1,t2,null),epsilon); 
		
		 // print projected tree to file
	    PrintWriter outputStream = null;
	  
	    try {
	       	outputStream = new PrintWriter(new FileWriter(outFileName));
	 
			outputStream.println(projectedTree.getNewick(true));
			
			if (outputStream != null) {
	            outputStream.close();
	        }
	    } catch (FileNotFoundException e) {
	        System.out.println("Error opening or writing to " + outFileName + ": "+ e.getMessage());
	        System.exit(1);
	    } catch (IOException e) {
	    	System.out.println("Error opening or writing to " + outFileName + ": " + e.getMessage());
	    	System.exit(1);
	    }
	}

	
	/** Prints nicely to a file info about which splits in the trees in otherTree in are incompatible
	 *  with the trees in goodTrees (considered individually).
	 * If goodTrees contains only one tree (i.e. an original tree), then this prints 
	 * information about which splits in the trees in otherTrees are incompatible with it.
	 * 
	 */
	public static void printIncompatibleSplits(PhyloTree[] goodTrees, PhyloTree[] otherTrees, String outFilename, Boolean paired) {
		String output = "";
		Vector<String> leaf2NumMap = goodTrees[0].getLeaf2NumMap();
		int startIndex, endIndex;
		
		for (int i = 0; i < goodTrees.length; i++) {
			output = output + "------------------------------------\n";
			output = output + "------------------------------------\n";
			output = output + "Tree " + i + " in input file:\n";
			
			
			if (paired) {
				startIndex = i;
				endIndex = i + 1;
			}
			else {
				startIndex = 0;
				endIndex = otherTrees.length;
			}
			
			for (int j = startIndex; j < endIndex; j++) {
				output = output + "\n";
				
				Vector<PhyloTreeEdge> incompEdges = goodTrees[i].getEdgesIncompatibleWith(otherTrees[j]);
				
//				System.out.println("otherTree edges incomp with good tree: ");
				
				if (incompEdges.size() >0) {
					output = output + "\tOther tree " + j + " is missing splits:\n";
					// display edges in verbose output, in a numbered list
					for (int k = 0; k < incompEdges.size(); k++) {
						output = output + "\t" + k + ". " + incompEdges.get(k).toStringReroot(leaf2NumMap, "Opossum") + "\n";
//						System.out.println(incompEdges.get(k).toStringVerbose(leaf2NumMap));
					}
					
				}
				
				incompEdges = otherTrees[j].getEdgesIncompatibleWith(goodTrees[i]);
				
				if (incompEdges.size() > 0) {
					output = output + "\n\thas extra splits:\n";
					
					// display edges in verbose output, in a numbered list
					for (int k = 0; k < incompEdges.size(); k++) {
						output = output + "\t" + k + ". " + incompEdges.get(k).toStringVerbose(leaf2NumMap) + "\n";
					}
				}
				
			}
			output = output + "\n";
		}
		
		
		// write the output to the file
		// write verbose output to geofile, if in verbose mode
	    PrintWriter outputStream = null;
	        
	    try {
	        outputStream = new PrintWriter(new FileWriter(outFilename));
	        	
	        outputStream.println(output);
	        
	    	if (outputStream != null) {
	               outputStream.close();
	        }
	    } catch (FileNotFoundException e) {
	        System.out.println("Error opening or writing to " + outFilename + ": "+ e.getMessage());
	        System.exit(1);
	    } catch (IOException e) {
	        System.out.println("Error opening or writing to " + outFilename + ": " + e.getMessage());
	        System.exit(1);
	    }
	}
	
	/**  Computes all log map coordinates for the trees in Trees,
	 *   relative to centreTree.
	 * @param nums
	 * @return
	 */
	public static void computeAllLogMaps(PhyloTree[] trees,PhyloTree centreTree, PrintWriter outfileStream) {
		double[] coords;
		
		for (PhyloTree tree: trees) {
			coords = centreTree.getLogMap(tree);
			outfileStream.println(Tools.doubleArray2String(coords));
		}
		
	}
	
	
	public static double mean(double [] nums) {
		double sum = 0;
		
		for (double n: nums) {
			sum = sum + n;
		}
		return sum/nums.length;
	}
	
	public static double var(double [] nums) {
		double squaresSum = 0;
		double sum = 0;
		
		for (double n: nums) {
			squaresSum = squaresSum + Math.pow(n,2);
			sum = sum + n;
		}
		
		return squaresSum/nums.length - Math.pow(sum/nums.length,2);
	}

	
	public static double scatter(double [] nums) {
		double sum = 0;
		double avg = 0;
		double intScatter = 0;
		
		for (double n: nums) {
			sum = sum + n;
		}
		
		avg = sum/nums.length;
		
		for (double n: nums) {
			intScatter = intScatter + Math.pow(n - avg, 2);			
		}		
		
		return intScatter;
	}

	//Gets 4 angles to print out
	public static void endRayAngles(PhyloTree[] trees, PhyloTree[] otherTrees, PrintWriter outFileStream) {
		Geodesic g1 = getGeodesic(trees[0], trees[1],null);
		Geodesic gA = getGeodesic(otherTrees[0], otherTrees[1], null);
		
		for (double e: Geodesic.getEndpointAngles(g1, gA, trees[0].getLeaf2NumMap(), trees[0].isRooted())) {
			double ne = e*180/Math.PI;
			
			if (outFileStream!=null) outFileStream.println(ne);
			System.out.println(ne);	
		}
		if (outFileStream!=null) outFileStream.close();
	}
}
