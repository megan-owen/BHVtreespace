package distances;

import distanceAlg1.EdgeAttribute;
import distanceAlg1.PhyloTree;
import distanceAlg1.PhyloTreeEdge;
import polyAlg.Tools;
import polyAlg.PolyMain;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class Distances {
	/** Help message (ie. which arguments can be used, etc.)
	 * 
	 */
	public static void displayHelp() {
		System.out.println("Command line syntax:");
		System.out.println("distance [options] treefile");
		System.out.println("Optional arguments:");	
//		System.out.println("\t -c \t double check results, by computing each distance with the target tree as the starting tree and vice versa; default is false");
		System.out.println("\t -d <distance> \t uses the given distance for the computations, where the possible distances are given below");
		System.out.println("\t -f <otherTreeFile> \t reads in an additional tree file and computes distances between all trees in treefile and all trees in otherTreeFile");
		System.out.println("\t -h || --help \t displays this message");
		System.out.println("\t -n \t normalize (vector of the lengths of all edges has length 1)");
		System.out.println("\t -o <outfile> \t store the output in the file <outfile>");
		System.out.println("\t -u \t unrooted trees (default is rooted trees)");
		System.out.println("\n");
		System.out.println("Distances are:");
		System.out.println("\t RF \t Robinson-Foulds distance");
		System.out.println("\t weightedRF \t weighted Robinson-Foulds distance");
	}
	
	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String treeFile = "";
		String otherTreeFile = null;
		String outFile = "dist_"; // default
		String distance = "RF"; // default
//		boolean doubleCheck = false;
		boolean rooted = true;
		PhyloTree[] otherTrees = null;
		

		if (args.length < 1) {
			displayHelp();
			System.out.println("Finished");
		}
		treeFile = args[args.length-1];
		for (int i = 0; i < args.length - 1; i++) {
			
			if (!args[i].startsWith("-")) { System.out.println("Invalid command line option"); displayHelp(); System.out.println("Finished"); }

			// specify distance
			else if (args[i].equals("-d")) {
				if (i < args.length -2) { distance = args[i+1]; i++; }
				else { System.err.println("Error: distance not specified"); displayHelp(); throw new Exception(); }
			}
			
			// other tree file, if desired
			else if (args[i].equals("-f")) {
				if (i < args.length -2) { otherTreeFile = args[i+1]; i++; }
				else { System.err.println("Error: name of other tree file not specified"); displayHelp(); throw new Exception(); }
			}
			
			// output file
			else if (args[i].equals("-o")) {
				if (i < args.length -2) { outFile = args[i+1]; i++; }
				else { displayHelp(); System.out.println("Finished"); }
			}
			
			// all other arguments.  Note we can have -vn
			else {
				for (int j = 1; j<args[i].length(); j++) {
					
					switch(args[i].charAt(j)) {						
					// doublecheck distances
//					case 'c':  doubleCheck = true;  break;
						
					// display help
					case 'h':  displayHelp();  System.out.println("Finished");  break;				
												
					// unrooted trees?
					case 'u': rooted = false;  break;
							
					default:
						System.out.println("Illegal command line option.\n");  displayHelp();  System.out.println("Finished");  break;
					} // end switch
				} // end for j
			} // end parsing an individual argument
		}  // end for i (looping through arguments)

		/* Read in the trees  */
		PhyloTree[] trees = PolyMain.readInTreesFromFile(treeFile,rooted);
		
		if (otherTreeFile != null) {
			otherTrees = PolyMain.readInTreesFromFile(otherTreeFile,rooted);
		}	
		
		// TODO:  routine to normalize tree lengths
		
		
		// print distances to file
	    PrintWriter outputStream = null;
	  
	    // Outputs the distances in a column, with the first two columns being the trees numbers and the third
	    // number the geodesic distance between those trees
	    try {
	       	outputStream = new PrintWriter(new FileWriter(outFile));
	 
	       	// computing distances between trees in two files
	       	if (otherTreeFile != null) {
	       		for (int i = 0; i < trees.length ; i++) {
	       			for (int j = 0; j< otherTrees.length; j++) {
	       				
	       				if (distance.equals("RF")) {
	       					int dist = rf(trees[i], otherTrees[j]);
	       					outputStream.println(i + "\t" + j + "\t" + dist);
	       				}
	       				else if (distance.equals("weightedRF")) {
	       					double dist = weightedRF(trees[i], otherTrees[j]);
	    	    			outputStream.println(i + "\t" + j + "\t" + Tools.roundSigDigits(dist, 6));
	       				}
	       				else {
	       					System.out.println("Error:  invalid distance.\n");
	       					throw new Exception();
	       				}
	       			}
				}
	       	}
	       	// computing distances between all trees in file
	       	else {	
	        	for (int i = 0; i < trees.length -1 ; i++) {
	        		for (int j = i + 1; j< trees.length; j++) {
		       				
		       			if (distance.equals("RF")) {
		       				int dist = rf(trees[i], trees[j]);
		       				outputStream.println(i + "\t" + j + "\t" + dist);
		       			}
		       			else if (distance.equals("weightedRF")) {
		       				double dist = weightedRF(trees[i], trees[j]);
		    	    		outputStream.println(i + "\t" + j + "\t" + Tools.roundSigDigits(dist, 6));
		       			}
		       			else {
		       				System.out.println("Error:  invalid distance.\n");
		       				throw new Exception();
		       			}
		       		}
	       		}
			}
			if (outputStream != null) {
	            outputStream.close();
	        }
	    } catch (FileNotFoundException e) {
	        System.out.println("Error opening or writing to " + outFile + ": "+ e.getMessage());
	        throw new Exception();
	    } catch (IOException e) {
	    	System.out.println("Error opening or writing to " + outFile + ": " + e.getMessage());
	    	throw new Exception();
	    }
		
	}

	/**  Returns the Robinson-Foulds distance between trees t1 and t2.  
	 *   That is, returns the number of splits in t1 but not t2, 
	 *   plus the number of splits in t2 but not t1.
	 *   This means that a split that appears in one trees 
	 *   and is compatible with the other tree, but does not appear in it,
	 *   contributes 1 to the RF distance.
	 *   Since 0 length edges are not added to trees when they are constructed,
	 *   we don't have to check for this.
	 *   	 
	 * @param t1
	 * @param t2
	 * @return
	 * @throws Exception 
	 */
	public static int rf(PhyloTree t1, PhyloTree t2) throws Exception {
		
		if (!(t1.getLeaf2NumMap().equals(t2.getLeaf2NumMap()))){
			System.out.println("Error: the two trees do not have the same leaves");
			System.out.println("Tree 1 leaves: " + t1.getLeaf2NumMap() );
			System.out.println("Tree 2 leaves: " + t2.getLeaf2NumMap() );
			throw new Exception();
		}
		
		int rf = 0;
		// Count the number of splits which are in t1, but not in t2.
		for (PhyloTreeEdge e : t1.getEdges()) {
			if (!t2.getSplits().contains(e.asSplit() ) ){
				rf = rf + 1;
			}
		}
		
		// Count the number of splits which are in t2, but not in t1.
		for (PhyloTreeEdge e : t2.getEdges()) {
			if (!t1.getSplits().contains(e.asSplit() ) ){
				rf = rf + 1;
			}
		}
		return rf;
	}
	
	/** Returns the contribution of the leaves to the geodesic distance.  Namely sum_over_all_leaves (t1_leaf_edge_length - t2_leaf_edge_length)^2.
	 * 
	 * @param t1
	 * @param t2
	 * @return
	 * @throws Exception 
	 */
	public static double getLeafSpaceL2Dist(PhyloTree t1, PhyloTree t2) throws Exception {
		double sumOfSquares = 0;
		
		if (!(t1.getLeaf2NumMap().equals(t2.getLeaf2NumMap()))){
			System.out.println("Error: the two trees do not have the same leaves");
			System.out.println("Tree 1 leaves: " + t1.getLeaf2NumMap() );
			System.out.println("Tree 2 leaves: " + t2.getLeaf2NumMap() );
			throw new Exception();
		}
		
		for(int i = 0; i < t1.getLeafEdgeAttribs().length; i++) {
			sumOfSquares = sumOfSquares + Math.pow(EdgeAttribute.difference(t1.getLeafEdgeAttribs()[i],t2.getLeafEdgeAttribs()[i]).norm(),2);
		}
		return Math.sqrt(sumOfSquares);
	}
	
	/** Computes the weighted Robinson-Foulds distance between trees t1 and t2.
	 *  
	 * @param t1
	 * @param t2
	 * @return
	 * @throws Exception 
	 */
	public static double weightedRF(PhyloTree t1, PhyloTree t2) throws Exception {
		double rf = 0.0;
				
		if (!(t1.getLeaf2NumMap().equals(t2.getLeaf2NumMap()))){
			System.out.println("Error: the two trees do not have the same leaves");
			System.out.println("Tree 1 leaves: " + t1.getLeaf2NumMap() );
			System.out.println("Tree 2 leaves: " + t2.getLeaf2NumMap() );
			throw new Exception();
		}
		
		// leaf contributions
		for(int i = 0; i < t1.getLeafEdgeAttribs().length; i++) {
			rf = rf + EdgeAttribute.difference(t1.getLeafEdgeAttribs()[i],t2.getLeafEdgeAttribs()[i]).norm();
		}
		
		// common edge contributions
		for (PhyloTreeEdge e : PhyloTree.getCommonEdges(t1, t2) ) {
			rf = rf + e.getNorm();
		}
			
		// edges only in t1 contribution
		for (PhyloTreeEdge e : t1.getEdgesIncompatibleWith(t2) ) {
			rf = rf + e.getNorm();
		}
		
		// edges only in t2 contribution
		for (PhyloTreeEdge e : t2.getEdgesIncompatibleWith(t1) ) {			
			rf = rf + e.getNorm();
		}
		return rf;
	
	}
	
	/** Returns the length of the cone path from t1 to t2.
	 *  This length is the square root of:  (norm of non-common edges of t1 + norm of non-common edges of t2)^2 + sum_{e is common split}(length of e in t1 - length of e in t2)^2 + sum_{e is a leaf split}(length of e in t1 - length of e in t2)^2
	 * @param t1
	 * @param t2
	 * @return
	 * @throws Exception 
	 */
	public static double conePath(PhyloTree t1, PhyloTree t2) throws Exception {
		if (!(t1.getLeaf2NumMap().equals(t2.getLeaf2NumMap()))){
			System.out.println("Error: the two trees do not have the same leaves");
			System.out.println("Tree 1 leaves: " + t1.getLeaf2NumMap() );
			System.out.println("Tree 2 leaves: " + t2.getLeaf2NumMap() );
			throw new Exception();
		}
		
		double t1norm = l2norm(t1.getEdgesIncompatibleWith(t2));		// norm of the vectors of the edges only in t1
		double t2norm = l2norm(t2.getEdgesIncompatibleWith(t1));		// norm of the vectors of the edges only in t2
		double commonNorm = l2norm(PhyloTree.getCommonEdges(t1,t2));
		
		return Math.sqrt(Math.pow(getLeafSpaceL2Dist(t1,t2),2) + Math.pow(t1norm + t2norm,2) + Math.pow(commonNorm,2));
		
	}
	
	/** Returns the length of the cone path from t1 to t2, assuming disjoint and not including leaf edges.
	 *  This length is the square root of:  (norm of non-common edges of t1 + norm of non-common edges of t2)^2 
	 * @param t1
	 * @param t2
	 * @return
	 * @throws Exception 
	 */
	public static double conePathDisjointNoLeaves(PhyloTree t1, PhyloTree t2) throws Exception {
		if (!(t1.getLeaf2NumMap().equals(t2.getLeaf2NumMap()))){
			System.out.println("Error computing the cone path between disjoint trees: the two trees do not have the same leaves");
			System.out.println("Tree 1 leaves: " + t1.getLeaf2NumMap() );
			System.out.println("Tree 2 leaves: " + t2.getLeaf2NumMap() );
			throw new Exception();
		}
		
		double t1norm = l2norm(t1.getEdges());		// norm of the vectors of the edges only in t1
		double t2norm = l2norm(t2.getEdges());		// norm of the vectors of the edges only in t2
		
		return t1norm + t2norm;
		
	}
	
	
	/**  Returns the l2 norm of the gives split lengths.
	 * 
	 * @param edges
	 * @return
	 */
	public static double l2norm(Vector<PhyloTreeEdge> edges) {
		double norm = 0;
		
		for (int i = 0; i < edges.size(); i++ ) {
			norm = norm + Math.pow(edges.get(i).getNorm(),2);
		}
		
		return Math.sqrt(norm);
	}
	
	/**  Returns the similarly index (or whatever) between two disjoint trees:  geodesic distance/cone path length.
	 * Ignores the leaf split lengths.
	 * 
	 * @param t1
	 * @param t2
	 * @return
	 * @throws Exception 
	 */
	public static double siDisjointNoLeaves(PhyloTree t1, PhyloTree t2) throws Exception {
		// Check that trees are disjoint and exit with error if they have common edges.
		Vector<PhyloTreeEdge> commonEdges = PhyloTree.getCommonEdges(t1, t2);
		if (commonEdges.size() > 0) { 
			System.out.println("Error:  trees need to be disjoint to compute siDisjointNoLeaves:  common edges are " + commonEdges);
			throw new Exception();
		}
		
		double conePathDist = conePathDisjointNoLeaves(t1,t2);
		if (conePathDist == 0) {
			return 0;
		}
		double si = PolyMain.getGeodesicNoCommonEdges(t1, t2).getDist()/conePathDist;
		return (2 + Math.sqrt(2))*(si - (Math.sqrt(2)/2));		// scale based on range of possible values
	}
	
	
	
	/**  Returns the similarly index (or whatever):  geodesic distance/cone path length.
	 * 
	 * @param t1
	 * @param t2
	 * @return
	 * @throws Exception 
	 */
	public static double si(PhyloTree t1, PhyloTree t2) throws Exception {
		return PolyMain.getGeodesic(t1, t2, null).getDist()/conePath(t1,t2);
	}
	
	
	public static double si2(PhyloTree t1, PhyloTree t2) throws Exception {
		return PolyMain.getGeodesic(t1,t2,null).getDist()/(t1.getDistanceFromOrigin() + t2.getDistanceFromOrigin());
	}
}


// Old main method.  Removed from being main method on Sept. 28, 2012.

//PhyloTree[] trees = PolyMain.readInTreesFromFile(treeFile,rooted);
//
//int[][] rfDists = new int[trees.length][trees.length];
//double[][] wrfDists = new double[trees.length][trees.length];
//double[][] geoDists = new double[trees.length][trees.length];
//
///* i = row, j = col */
//for( int i = 0; i < trees.length; i++) {
//	rfDists[i][i] = 0;
//	wrfDists[i][i] = 0;
//	geoDists[i][i] = 0;
//	for (int j = 0; j < i; j++) {
//		System.out.println("Tree " + i + " vs tree " + j);
//		rfDists[i][j] = rf(trees[i],trees[j]);
//		if (rfDists[i][j] != rf(trees[j],trees[i])) {
//			System.out.println("#### Error?  " + i + " -> " + j + " = " + rfDists[i][j] + " but " + j + " -> " + i + " = " + rf(trees[j],trees[i]));
//		}
//		System.out.println("RF distance: " + rfDists[i][j]); 
//		
//		wrfDists[i][j] = weightedRF(trees[i],trees[j]);
//		System.out.println("Weighted RF distance: " + wrfDists[i][j]);
//		if (Tools.round(wrfDists[i][j],3) != Tools.round(weightedRF(trees[j],trees[i]),3)) {
//			System.out.println("#### Error?  " + i + " -> " + j + " = " + wrfDists[i][j] + " but " + j + " -> " + i + " = " + weightedRF(trees[j],trees[i]));
//		}
//		
//		geoDists[i][j] =PolyMain.getGeodesic(trees[i],trees[j],null).getDist();
//		System.out.println("Geodesic distance: " + geoDists[i][j] + "\n");
//		if (Tools.round(geoDists[i][j],3) != Tools.round(PolyMain.getGeodesic(trees[j],trees[i],null).getDist(),3)) {
//			System.out.println("#### Error?  " + i + " -> " + j + " = " + geoDists[i][j] + " but " + j + " -> " + i + " = " + PolyMain.getGeodesic(trees[j],trees[i],null).getDist());
//		}
//		
////		System.out.println("RF distance reversed: " + rf(trees[j],trees[i])); 
////		System.out.println("Weighted RF distance reversed: " + weightedRF(trees[j],trees[i]));
////		System.out.println("Geodesic distance: " + PolyMain.getGeodesic(trees[j],trees[i],null).getDist() + "\n");
//	}
//}
//
////print centroid to outfile
////   PrintWriter outputStream = null;
//
//try {
//	PrintWriter outputStreamRF = new PrintWriter(new FileWriter(outFile + "_rf.txt"));
//	PrintWriter outputStreamWRF = new PrintWriter(new FileWriter(outFile + "_weightRF.txt"));
//	PrintWriter outputStreamGeo = new PrintWriter(new FileWriter(outFile + "_geo.txt"));
//
//	String[] treeLabels = {"Orig", "ML", "Maj_0", "Maj_1", "Maj_2", "Maj_3", "Maj_4", "Cen_r1_0", "Cen_r1h_0","Cen_r2_0", "Cen_r2h_0", "Cen_r1_1", "Cen_r1h_1", "Cen_r2_1", "Cen_r2h_1", "Cen_r1_2", "Cen_r1h_2", "Cen_r2_2", "Cen_r2h_2", "Cen_r1_3", "Cen_r1h_3", "Cen_r2_3", "Cen_r2h_3", "Cen_r1_4", "Cen_r1h_4", "Cen_r2_4", "Cen_r2h_4", "MAP_0", "MAP_1", "MAP_2", "MAP_3", "MAP_4"};
//
//	/*  Print the labels on the first line. */
//	for (int i = 0; i < treeLabels.length; i++) {
//		outputStreamRF.print(treeLabels[i] + "\t");
//		outputStreamWRF.print(treeLabels[i] + "\t");
//		outputStreamGeo.print(treeLabels[i] + "\t");
//	}
//	outputStreamRF.println();
//	outputStreamWRF.println();
//	outputStreamGeo.println();
//	
//	for( int i = 0; i < trees.length; i++) {
//		outputStreamRF.print(treeLabels[i] + "\t");
//		outputStreamWRF.print(treeLabels[i] + "\t");
//		outputStreamGeo.print(treeLabels[i] + "\t");
//		
//		for (int j = 0; j < i; j++) {
//				
//			outputStreamRF.print(rfDists[i][j] + "\t");
//			outputStreamWRF.print(wrfDists[i][j] + "\t");
//			outputStreamGeo.print(geoDists[i][j] + "\t");
//		}
//		outputStreamRF.println(rfDists[i][i]);
//		outputStreamWRF.println(wrfDists[i][i]);
//		outputStreamGeo.println(geoDists[i][i]);
//	}
//
//	if (outputStreamRF != null) {
//		outputStreamRF.close();
//	}
//	if (outputStreamWRF != null) {
//		outputStreamWRF.close();
//	}
//	if (outputStreamGeo != null) {
//		outputStreamGeo.close();
//	}
//} catch (FileNotFoundException e) {
//	System.out.println("Error opening or writing to " + outFile + ": "+ e.getMessage());
//	throw new Exception();
//} catch (IOException e) {
//	System.out.println("Error opening or writing to " + outFile + ": " + e.getMessage());
//	throw new Exception();
//}
//
//System.out.println("Finished"); 
