/** This file is part of sturmMean, a program form computing the
 * Frechet mean between phylogenetic trees.
    Copyright (C) 2008-2012  Megan Owen, Scott Provan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

package centroid;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;

import distanceAlg1.*;
import distances.Analysis;
import static polyAlg.PolyMain.getGeodesic;
import static polyAlg.PolyMain.calcGeoDist;
import polyAlg.PolyMain;

public class CentroidMain {
	
	/** Help message (ie. which arguments can be used, etc.)
	 * 
	 */
	public static void displayHelp() {
		System.out.println("Command line syntax:");
		System.out.println("java -jar sturmMean.jar [options] treefile");
		System.out.println("Optional arguments:");	
		System.out.println("\t -a <algorithm> \t algorithm to choose the next tree.  Choices are \"random\" and \"rand_perm\". ");
		System.out.println("\t -c <length> \t length of Cauchy sequence for determining convergence.  Default is 5.");
		System.out.println("\t -e <epsilon> \t set the value of epsilon used to determine convergence.");
		System.out.println("\t -f <epsilonFactor> \t calculate the value of epsilon automatically, using the specified factor. Default for the factor is 5000");
		System.out.println("\t -h \t display this message.");
		System.out.println("\t -i <interval> \t display the current mean tree and its standard deviation every <interval> iterations.  Default does not display intermediate trees.");
		System.out.println("\t -n <numIterations> \t max number of iterations for algorithm.  Default is 100 000.");
		System.out.println("\t -o <outfile> \t store the output in the specified file.  Default is output.txt");
		System.out.println("\t -p \t randomly and independently permute leaves of each tree");
		System.out.println("\t -s start> \t display the current mean tree and its standard deviation starting at iteration <start>.  Default is 1.");
		System.out.println("\t -u \t set input trees as unrooted. Default is trees are rooted.");
	}
	
	/**
	 * @param args
	 */
	public static void main (String[] args) {		
		/* Variable initialization */
		String treeFile = "";
		String outfile = "output.txt"; // default
		boolean rooted = true;
		String algorithm = "rand_perm";
		int numIter = 100000;	// default; cannot have more than 2,147,483,647 iterations.
		int cauchyLength = 5;
		int displayIter = 0;
		int displayStart = 1;
		Boolean permute = false;
		Boolean twoRuns = false;
		double epsilonFactor = 5000;
		double epsilon = 0;
		int bootstrapRep = 0;
		int p = 2;
		
		/* Parse command line arguments */

		if (args.length < 1) {
			System.out.println("Error: Missing input file name"); displayHelp(); System.exit(1);
		}
		treeFile = args[args.length-1];
		
		// check that the user hasn't inputted just -h
		if (treeFile.equals("-h")) {
			displayHelp(); System.exit(0);
		}
		
		for (int i = 0; i < args.length - 1; i++) {
			
			if (!args[i].startsWith("-")) { System.out.println("Invalid command line option"); displayHelp(); System.exit(1); }

			// number of iterations
			else if (args[i].equals("-n")) {
				if (i < args.length -2) {
					numIter = Integer.valueOf(args[i+1]); i++;
					if (numIter < 2) { System.err.println("Error:  number of iterations must be greater than 1"); System.exit(1); }
				}
				else { System.err.println("Error: number of iterations not specified"); displayHelp(); System.exit(1); }
			}
			// output file
			else if (args[i].equals("-o")) {
				if (i < args.length -2) { outfile = args[i+1]; i++; }
				else { System.err.println("Error: output file not specified");  displayHelp(); System.exit(1); }
			}	
			// algorithm for choosing next tree
			else if (args[i].equals("-a")) {
				if (i < args.length -2) { algorithm = args[i+1]; i++; }
				else { System.err.println("Error: algorithm not specified"); displayHelp(); System.exit(1); }
			}
			// use a set of trees bootstrapped from the input trees; 
			// parameter is the number of bootstrap replicates
			else if (args[i].equals("-b")) {
				if (i < args.length - 2) { bootstrapRep = Integer.valueOf(args[i+1]); i++; }
				else { System.err.println("Error: number of bootstrap replicates not specifie"); displayHelp(); System.exit(1); }
			}
			// display states every i iterations
			else if (args[i].equals("-i")) {
				if (i < args.length -2) { displayIter = Integer.valueOf(args[i+1]); i++; }
				else { displayHelp(); System.exit(1); }
			}
			// start displaying the state after s iterations
			else if (args[i].equals("-s")) {
				if (i < args.length -2) { displayStart = Integer.valueOf(args[i+1]); i++; }
				else { displayHelp(); System.exit(1); }
			}
			// length of Cauchy sequence determining convergence
			else if (args[i].equals("-c")) {
				if (i < args.length - 2) { cauchyLength = Integer.valueOf(args[i+1]); i++; }
				else { System.err.println("Error: length of Cauchy sequence for determining convergence not specified"); displayHelp(); System.exit(1); }
			}
			// epsilon 
			else if (args[i].equals("-e")) {
				if (i < args.length - 2) { epsilon = Double.valueOf(args[i+1]); i++; }
				else { System.err.println("Error: epsilon for determining convergence not specified"); displayHelp(); System.exit(1); }
			}
			// epsilon factor (used in computing epsilon)
			else if (args[i].equals("-f")) {
				if (i < args.length - 2) { epsilonFactor = Double.valueOf(args[i+1]); i++; }
				else { System.err.println("Error: epsilon factor for determining convergence not specified"); displayHelp(); System.exit(1); }
			}	
			
			// compute the tree minimizing the sum of the geodesic distances to the p-th power
			else if (args[i].equals("-q")) {
				if (i < args.length - 2) { p = Integer.valueOf(args[i+1]); i++; }
				else { System.err.println("Error: p value for finding the tree minimizing the sum of the distances to the p-th power not specified"); displayHelp(); System.exit(1); }
			}
				
			// all other arguments.  Note we can have -vn
			else { 
				for (int j = 1; j<args[i].length(); j++) {
					switch(args[i].charAt(j)) {						
					// display help
					case 'h': displayHelp(); System.exit(0); break;
					// permute leaves	
					case 'p': permute = true; break;
					// set trees to be unrooted
					case 'u': rooted = false; break;
					
					default: System.out.println("Illegal command line option.\n"); displayHelp(); System.exit(1); break;
					} // end switch
				} // end j loop (arguments without parameter)
			} // end parsing an individual argument
		}  // end for i (looping through arguments)
		
		/* Read in the trees  */
		
		PhyloTree[] trees = PolyMain.readInTreesFromFile(treeFile,rooted);
		int numTrees = trees.length;
		System.out.println("Num. trees = " + numTrees);
		
		/*  Permute leaves in all trees, if desired. */
		if (permute) {
			for (int i = 0; i < trees.length; i++) {
				trees[i].permuteLeaves();
			}
			System.out.println("Permuting leaves of input trees.");
		}
		
		/* Use bootstrap replicates instead, if desired. */
		if (bootstrapRep > 0) {
			trees = PhyloTree.resample(trees,bootstrapRep);
			numTrees = trees.length;
			System.out.println("Using bootstrap replicates as tree set.");
		}
		
		// if there was only one tree in the file, that tree is the mean tree
		// TODO:  write to file
		if (numTrees == 1) {
			System.out.println("Mean tree is " + trees[0]);
			System.exit(0);
		}

		
		PhyloTree star = Analysis.getStarTree(trees);	
		System.out.println("Star tree is " + star.getNewick(true));
		
		double stdDevStar = stdDev(star,trees);	
		System.out.println("Standard deviation if star tree is mean: " + stdDevStar);
		
		/* Pre-processing to determine epsilon for convergence test. */
		if (epsilon <= 0) {
			epsilon = getEpsilon(trees,epsilonFactor);  // info about this printed in getEpsilon method
		}
			
		if (p >= 2) {
			getPMeanViaRandPermCauchy(trees,p,numIter,cauchyLength,epsilon,outfile,displayIter,displayStart);
			System.exit(0);
		}
		
		
		// algorithm chooses the next tree at random; convergence tested by Cauchy sequence
		if (algorithm.equals("random") && !twoRuns) {
			getCentroidViaRandomCauchy(trees,numIter,cauchyLength,epsilon,outfile,displayIter,displayStart);
		}
		else if (algorithm.equals("random") && twoRuns) {
			getCentroidViaRandomTwoRuns(trees,numIter,cauchyLength,epsilon,outfile,displayIter,displayStart);
		}
		// algorithm orders all trees randomly, chooses them in this order, 
		// then chooses a new random order of the trees for the next round
		else if (algorithm.equals("rand_perm") && twoRuns) {
			getCentroidViaRandPermTwoRuns(trees,numIter,cauchyLength,epsilon,outfile,displayIter,displayStart);
		}
		else if( algorithm.equals("rand_perm") && !twoRuns ) {
			getCentroidViaRandPermCauchy(trees,numIter,cauchyLength,epsilon,outfile,displayIter,displayStart);
		}
		else {
			System.out.println("Error:  unknown algorithm");
			System.exit(1);
		}

		System.exit(0);
	}
	
	
	/**  Uses Sturm's algorithm to get the centroid of trees.
	 *   For each "round" of Sturm's algorihtm, permutes the array trees and uses the trees in that order.  Repeat for each round.
	 *   Convergence iff there are cauchyLength consecutive trees within geodesic distance epsilon of each other.
	 *   Returns null if no convergence in numIter iterations.  Otherwise returns the centroid.
	 *   
	 * 
	 * @param trees
	 * @param epsilon
	 * @param numIter
	 * @return
	 */
	public static PhyloTree getCentroidViaRandPermCauchy(PhyloTree[] trees, int numIter, int cauchyLength, double epsilon,String outfile, int displayIter, int displayStart) {
		DecimalFormat d6o = new DecimalFormat("#0.######");
		Boolean converge = false;
		int numTrees = trees.length;
		Vector<PhyloTree> oldCentroids = new Vector<PhyloTree>();	  // stores the last 5 centroids form for checking for cauchy sequence		
			
		/*  Permute the trees for this round. */
		PhyloTree [] shuffledTrees = permuteTrees(trees);
			
		// find the first centroid candidate
		PhyloTree centroid = shuffledTrees[0];
		oldCentroids.add(centroid.clone());
		int cauchyIndex = 0;  // last position in the vector that there is a old centroid at
			
		int i = 1;
		while ((i < numIter ) && (!converge)) {
			// Shuffle trees if necessary
			if (i % numTrees == 0) {
				shuffledTrees = permuteTrees(trees);
			}
			
			centroid =  getGeodesic(centroid, shuffledTrees[i % numTrees],null).getTreeAt((double)1/(i+1),centroid.getLeaf2NumMap(),centroid.isRooted());	


			if (displayIter > 0) {
				if ((i%displayIter ==0) &&  (i >= displayStart)) {
					System.out.println("Iteration " + i + ". Variance, tree: " + d6o.format(variance(centroid,trees)) + ", " + centroid.getNewick(true));
				}
			}
			
			/* Check for Cauchy sequence to determine convergence. */
			
			for(int j = cauchyIndex; j >= 0; j--) {
				double dist = calcGeoDist(oldCentroids.get(j), centroid);
				if (dist > epsilon) {
					// don't have a Cauchy sequence.  Move all already checked centroids to the front of the queue.
					for (int k = j; k >= 0; k--) {
						oldCentroids.remove(k);
					}
					break;
				}
			}
			oldCentroids.add(centroid.clone());
			cauchyIndex = oldCentroids.size() - 1;
			
			// Check if we have a Cauchy sequence of the desired length.
			if (oldCentroids.size() == cauchyLength ) {
				converge = true;
			}
			i++;
		}
		
		displayCauchy("random permutations of the tree set.",i,cauchyLength,epsilon,converge,centroid,trees,outfile);

		if (converge) {
			return centroid;
		}
		else {
			return null;
		}
	}
	
	
	
	public static PhyloTree getCentroidViaRandomCauchy(PhyloTree[] trees, long numIter, int cauchyLength, double epsilon, String outfile, int displayIter, int displayStart) {
		DecimalFormat d6o = new DecimalFormat("#0.######");
		Boolean converge = false;
		int numTrees = trees.length;

		int numNextTree;
		Random r = new Random();
		PhyloTree centroid = null;
		
		/* Variables for testing for a Cauchy sequence to determine convergence */
		int cauchyIndex = -1;  // last position in the vector that there is a old centroid at
		Vector<PhyloTree> oldCentroids = new Vector<PhyloTree>();	  // stores the last 5 centroids form for checking for cauchy sequence		

		
		// choose a tree at random to start
		centroid = trees[r.nextInt(numTrees)];			
	
		int i = 1;
		while ((i < numIter) && (!converge)) {
			numNextTree = r.nextInt(numTrees);
			centroid =  getGeodesic(centroid, trees[numNextTree],null).getTreeAt((double)1/(i+1),centroid.getLeaf2NumMap(),centroid.isRooted());	
		
			// For Cauchy sequence convergence check 
		
			if (displayIter > 0) {
				if ((i%displayIter ==0) &&  (i >= displayStart)) {
					System.out.println("Iteration " + i + ". Variance, tree: " + d6o.format(variance(centroid,trees)) + ", " + centroid.getNewick(true));
				}
			}
	
			/* Check for Cauchy sequence to determine convergence. */
			for(int j = cauchyIndex; j >= 0; j--) {
				double dist = calcGeoDist(oldCentroids.get(j), centroid);
				if (dist > epsilon) {
					// don't have a Cauchy sequence.  
					// Move all already checked centroids to the front of the queue.
					for (int k = j; k >= 0; k--) {
						oldCentroids.remove(k);
					}
					break;
				}
			}
			oldCentroids.add(centroid.clone());
			cauchyIndex = oldCentroids.size() - 1;
			
			// Check if we have a Cauchy sequence of right length.
			if (oldCentroids.size() == cauchyLength ) {
				converge = true;
			}
			i++;
		} // end while not converged

		displayCauchy("choose tree randomly",i,cauchyLength,epsilon,converge,centroid,trees,outfile);
		
		if (converge) {
			return centroid;
		}
		else {
			return null;
		}

	}
	
	public static PhyloTree getCentroidViaRandomTwoRuns(PhyloTree[] trees,int numIter, int matchLength, double epsilon,String outfile, int displayIter, int displayStart) {
		DecimalFormat d6o = new DecimalFormat("#0.######");
		Boolean converge = false;
		int numTrees = trees.length;
		int convergenceCounter = 0;  // counts number of times the two runs are within epsilon of each other

		int numNextTree, numNextTree2;
		Random r = new Random();
		PhyloTree centroid = null, centroid2 = null;
		
		// choose a tree at random to start
		centroid = trees[r.nextInt(numTrees)];	
		centroid2 = trees[r.nextInt(numTrees)];			

	
		int i = 1;
		while ((i < numIter) && (!converge)) {
			numNextTree = r.nextInt(numTrees);
			numNextTree2 = r.nextInt(numTrees);

			centroid =  getGeodesic(centroid, trees[numNextTree],null).getTreeAt((double)1/(i+1),centroid.getLeaf2NumMap(),centroid.isRooted());	
			centroid2 =  getGeodesic(centroid2, trees[numNextTree2],null).getTreeAt((double)1/(i+1),centroid.getLeaf2NumMap(),centroid.isRooted());	


			if (displayIter > 0) {
				if ((i%displayIter ==0) &&  (i >= displayStart)) {
					System.out.println("Iteration " + i + ". Run 1 variance, tree: " + d6o.format(variance(centroid,trees)) + ", " + centroid.getNewick(true));
					System.out.println("Iteration " + i + ". Run 2 variance, tree: " + d6o.format(variance(centroid2,trees)) + ", " + centroid2.getNewick(true));
				}
			}
			
			Geodesic geo = getGeodesic(centroid, centroid2, null);
			double dist = geo.getDist();
			// check for convergence
			if (dist <= epsilon) {
				convergenceCounter++;
				if (convergenceCounter == matchLength) {
					converge = true;
					centroid = geo.getTreeAt(0.5, centroid.getLeaf2NumMap(), centroid.isRooted());
				}
			}
			else {
				convergenceCounter = 0;
			}
			i++;
		}
		
		displayTwoRuns("choose tree randomly.",i,matchLength,epsilon,converge,centroid,centroid2,trees,outfile);
		
		if (converge) {
			return centroid;
		}
		else {
			return null;
		}

	}

	
	public static PhyloTree getCentroidViaRandPermTwoRuns(PhyloTree[] trees,int numIter, int matchLength,double epsilon, String outfile, int displayIter, int displayStart) {
		DecimalFormat d6o = new DecimalFormat("#0.######");
		Boolean converge = false;
		int numTrees = trees.length;
		int convergenceCounter = 0;  // counts number of times the two runs are within epsilon of each other

		PhyloTree centroid = null, centroid2 = null;	
		
		PhyloTree [] shuffledTrees1 = permuteTrees(trees);
		PhyloTree [] shuffledTrees2 = permuteTrees(trees);

		
		// find the first centroid candidate
		centroid = shuffledTrees1[0];
		// find the second centroid candidate
		centroid2 = shuffledTrees2[0];
		
		int i = 1;
		while ((i < numIter ) && (!converge)) {
			// Shuffle trees if necessary
			if (i % numTrees == 0) {
				shuffledTrees1 = permuteTrees(trees);
				shuffledTrees2 = permuteTrees(trees);
			}
			
			centroid =  getGeodesic(centroid, shuffledTrees1[i % numTrees],null).getTreeAt((double)1/(i+1),centroid.getLeaf2NumMap(),centroid.isRooted());				
			centroid2 =  getGeodesic(centroid2, shuffledTrees2[i % numTrees],null).getTreeAt((double)1/(i+1),centroid.getLeaf2NumMap(),centroid.isRooted());	

		
			if (displayIter > 0) {
				if ((i%displayIter ==0) &&  (i >= displayStart)) {
					System.out.println("Iteration " + i + ". Run 1 variance, tree: " + d6o.format(variance(centroid,trees)) + ", " + centroid.getNewick(true));
					System.out.println("Iteration " + i + ". Run 2 variance, tree: " + d6o.format(variance(centroid2,trees)) + ", " + centroid2.getNewick(true));
				}
			}
			
			Geodesic geo = getGeodesic(centroid, centroid2, null);
			double dist = geo.getDist();
			// check for convergence
			if (dist <= epsilon) {
				convergenceCounter++;
				if (convergenceCounter == matchLength) {
					converge = true;
					centroid = geo.getTreeAt(0.5, centroid.getLeaf2NumMap(), centroid.isRooted());
				}
			}
			else {
				convergenceCounter = 0;
			}
			i++;
		}
		
		displayTwoRuns("random permutations of the tree set.",i,matchLength,epsilon,converge,centroid,centroid2,trees,outfile);
		
		if (converge) {
			return centroid;
		}
		else {
			return null;
		}
	}
	
	/** 
	 *  Removes elements from the front of the vector so that all remaining elements
	 *   differ by less than epsilon. i.e. last centroid in the vector was the latest one added
	 *   So we can check for Cauchy convergence by checking that the vector returned contains the 
	 *   desired number of elements.
	 * 
	 * 
	 * @param oldCentroids
	 * @param epsilon
	 * @return
	 */
	public static Vector<PhyloTree> checkCauchyConvergence(Vector<PhyloTree> centroids, double epsilon, int cauchyLength) {
		int secondLastCentroidIndex = centroids.size() - 1; 
		int lastIndexToKeep = -1;
		
		for (int i = secondLastCentroidIndex; i >= 0; i-- ) {
			for (int j = i + 1; j < centroids.size(); j++) {
				if (calcGeoDist(centroids.get(i), centroids.get(j)) > epsilon) {
					lastIndexToKeep = i + 1;
					break;
				}
			}
		}
		
		if (lastIndexToKeep > -1) {
			// Remove all centroids with lower index than lastIndexToKeep.
			for (int k = lastIndexToKeep - 1; k >= 0; k--) {
				centroids.remove(k);
			}
		}
		return centroids;
	}
	
	
	
	
	
	/** Reports on the centroid computations.
	 * 
	 * @param i
	 * @param centroid
	 * @param trees
	 */
	public static void report(int i, PhyloTree centroid, PhyloTree[] trees) {
		System.out.println("Centroid candidate " + i + " with score " + sumOfSquareDist(centroid,trees) + " = " + centroid);
	}
	
	/** Outputs and writes to file mean tree and convergence info.
	 *  
	 */
	public static void displayCauchy(String algorithm,long numIter, int cauchyLength, double epsilon, boolean converge,PhyloTree centroid, PhyloTree[] trees, String outfile) {
		DecimalFormat d6o = new DecimalFormat("#0.######");
		double stdDevOfCentroid = stdDev(centroid,trees);
		PrintWriter outputStream;

		
		System.out.println("Tree selection:  " + algorithm);
		System.out.println("Convergence:  Cauchy sequence of length " + cauchyLength + " with epsilon = " + epsilon);

		if (!converge) {
			System.out.println("Did not converge after " + numIter + " iterations.  Displaying last tree found.");
		}
		else {
			System.out.println("Converged after " + numIter + " iterations.  Displaying approximate mean tree.");
		}
		System.out.println("Standard Deviation if tree is mean: " + d6o.format(stdDevOfCentroid) + "  Tree:");
		System.out.println(centroid.getNewick(true));
		
		
		
		//  print centroid to outfile
        if (outfile != null) {
        	try {
        		outputStream = new PrintWriter(new FileWriter(outfile));
        		outputStream.println("Tree selection:  " + algorithm);
        		outputStream.println("Convergence:  Cauchy sequence of length " + cauchyLength + " with epsilon = " + epsilon);
        		
        		if (!converge) {
        			outputStream.println("Did not converge after " + numIter + " iterations.  Displaying last tree found.");
        		}
        		else {
        			outputStream.println("Converged after " + numIter + " iterations.  Displaying approximate mean tree.");
        		}
        		outputStream.println("Standard Deviation if tree is mean: " + d6o.format(stdDevOfCentroid) + "  Tree:");
        		outputStream.println(centroid.getNewick(true));
        		
        		if (outputStream != null) {
        			outputStream.close();
        		}
        	} catch (FileNotFoundException e) {
        		System.out.println("Error opening or writing to " + outfile + ": "+ e.getMessage());
        		System.exit(1);
        	} catch (IOException e) {
        		System.out.println("Error opening or " +
        			"writing to " + outfile + ": " + e.getMessage());
        		System.exit(1);
        	}
        }
	}
	
	/** Displays/writes to file output when convergence checked by two runs.
	 * 
	 * @param algorithm
	 * @param numIter
	 * @param cauchyLength
	 * @param epsilon
	 * @param converge
	 * @param centroid
	 * @param centroid2
	 * @param trees
	 * @param outfile
	 */
	public static void displayTwoRuns(String algorithm,long numIter,int matchLength ,double epsilon,boolean converge, PhyloTree centroid, PhyloTree centroid2,PhyloTree[] trees,String outfile) {
		DecimalFormat d6o = new DecimalFormat("#0.######");
		double stdDevCentroid = stdDev(centroid,trees);
		double stdDevCentroid2 = stdDev(centroid2,trees);
		PrintWriter outputStream;
		
		
		System.out.println("Tree selection:  " + algorithm);
		System.out.println("Convergence:  two runs are within epsilon = " + epsilon + " of each other for " + matchLength + " iterations.");
		
		
		
		if (!converge) {
    		// Print convergence info to console/log file
			System.out.println("Did not converge after " + numIter + " iterations.  Displaying last two trees found.");
    		
			System.out.println("Standard deviation if run 1 tree is mean: " + d6o.format(stdDevCentroid) + "  Run 1 tree:");
			System.out.println(centroid.getNewick(true));
    		
			System.out.println("Standard deviation if run 2 tree is mean: " + d6o.format(stdDevCentroid2) + "  Run 2 tree:");
			System.out.println(centroid2.getNewick(true));
		}
		else {    		
    		// Print above to console/log file
    		System.out.println("Converged after " + numIter + " iterations.  Displaying approximate mean tree.");
    		System.out.println("Standard deviation if tree is mean: " + d6o.format(stdDevCentroid) + "  Tree:");
    		System.out.println(centroid.getNewick(true));
		}
		
		//  print centroid to outfile
        if (outfile != null) {
        	try {
        		outputStream = new PrintWriter(new FileWriter(outfile));
        		outputStream.println("Tree selection:  " + algorithm);
        		outputStream.println("Convergence:  two runs are within epsilon = " + epsilon + " of each other for " + matchLength + " iterations.");
        		
        		
        		
        		if (!converge) {
        			outputStream.println("Did not converge after " + numIter + " iterations.  Displaying last two trees found.");
            		
        			outputStream.println("Standard deviation if run 1 tree is mean: " + d6o.format(stdDevCentroid) + "  Run 1 tree:");
            		outputStream.println(centroid.getNewick(true));
            		
            		outputStream.println("Standard deviation if run 2 tree is mean: " + d6o.format(stdDevCentroid2) + "  Run 2 tree:");
            		outputStream.println(centroid2.getNewick(true));
            		
        		}
        		else {
        			outputStream.println("Converged after " + numIter + " iterations.  Displaying approximate mean tree.");
            		outputStream.println("standard deviation if tree is mean: " + d6o.format(stdDevCentroid) + "  Tree:");
            		outputStream.println(centroid.getNewick(true));
        		}
        		
        		if (outputStream != null) {
        			outputStream.close();
        		}
        	} catch (FileNotFoundException e) {
        		System.out.println("Error opening or writing to " + outfile + ": "+ e.getMessage());
        		System.exit(1);
        	} catch (IOException e) {
        		System.out.println("Error opening or " +
        			"writing to " + outfile + ": " + e.getMessage());
        		System.exit(1);
        	}
        }
	}
	
	
	/**  Calculates an epsilon to use for testing convergence.  
	 * epsilonFactor is the amount by which we divide the square of the variance.
	 * 
	 * @param trees
	 * @return
	 */
	public static double getEpsilon(PhyloTree[] trees, double epsilonFactor) {
		// We'll run through all trees in random order 5 times, compute the variance, and choose the distance between runs based on that.
		int numRounds = 5;	// number of rounds of all trees
		PhyloTree testCentroid = null;
		double epsilon = 0;  // distance that the two runs must be in to converge
		
		int numTrees = trees.length;
		
		
		PhyloTree[] shuffledTrees = permuteTrees(trees);
		
		// find the first centroid candidate
		testCentroid = shuffledTrees[0];
		
		int i = 1;
		while (i < numRounds*numTrees ) {
			testCentroid = getGeodesic(testCentroid, shuffledTrees[i % numTrees], null).getTreeAt((double)1/(i+1),testCentroid.getLeaf2NumMap(), testCentroid.isRooted());
			i++;
			
			if (i % numTrees == 0) {
				shuffledTrees = permuteTrees(trees);
			}
		}
		
		System.out.println("Calculating epsilon...");
		System.out.println("Test mean tree is " + testCentroid.getNewick(true));
		double testStdDev = stdDev(testCentroid,trees);
		System.out.println("Test standard deviation is " + testStdDev);
				
		epsilon = 1.0/epsilonFactor*testStdDev;
			
		System.out.println("epsilon is " + epsilon);
		return epsilon;
	}
	
	/**  Permutes the trees, and returns a new array.
	 * 
	 * @param trees
	 * @return
	 */
	public static PhyloTree[] permuteTrees(PhyloTree[] trees) {
		PhyloTree[] copyTrees = Arrays.copyOf(trees, trees.length);
		List<PhyloTree> treeList = Arrays.asList( copyTrees );
		Collections.shuffle( treeList );
		return treeList.toArray(new PhyloTree[0]);
	}
	
	/** Computes the sum of square distances between the centroid and each tree in trees.
	 * 
	 * @param centroid
	 * @param trees
	 * @return
	 */
	public static double sumOfSquareDist(PhyloTree centroid, PhyloTree[] trees) {
		double sumOfSquares = 0;
		
		for(int i = 0; i < trees.length; i++) {
			sumOfSquares = sumOfSquares + Math.pow(calcGeoDist(centroid, trees[i]),2);
		}
		return sumOfSquares;
	}
	
	/**  Computes the variance of the centroid of the set of trees.
	 *   Uses an unbiased estimator as this is a sample variance.
	 * 
	 * @param centroid
	 * @param trees
	 * @return
	 */
	public static double variance(PhyloTree centroid, PhyloTree[] trees) {
		return 1/( (double) trees.length -1)*sumOfSquareDist(centroid,trees);
	}
	
	/** Computes the (sample) standard deviation of the given centroid as the mean of the trees.
	 *  Standard deviation = sqrt(variance) and is in the same units as the mean.
	 * @param centroid
	 * @param trees
	 * @return
	 */
	public static double stdDev(PhyloTree centroid, PhyloTree[] trees) {
		return Math.sqrt(variance(centroid, trees));
	}
	
	/**  Prints the current state of the algorithm, including the iteration #, centroid estimate 
	 *  (of both runs, if applicable), and the distance between centroid estimates (if two runs).
	 * 
	 * @param iter
	 * @param centroid
	 * @param centroid2
	 */
	public static void printCurrentState(int iter, PhyloTree centroid, PhyloTree centroid2) {
//		DecimalFormat d4o = new DecimalFormat("#0.####");
//		DecimalFormat d6o = new DecimalFormat("#0.####");

			
		
/*		System.out.println("Iteration " + iter);
		System.out.println("run 1 centroid: " + centroid.getNewick(false));
		System.out.println("\t" + centroid);

		// If we have done two runs
		if (centroid2 != null) {
			System.out.println("run 2 centroid: " + centroid2.getNewick(false));
			System.out.println("\t" + centroid2);
			System.out.println("Distance between runs: " + d6o.format(PolyMain.getGeodesic(centroid, centroid2, null).getDist()));
			System.out.println("RF distance between runs: " + d6o.format(Distances.rf(centroid, centroid2)));
			System.out.println("Weighted RF distance between runs: " + d6o.format(Distances.weightedRF(centroid, centroid2)));
		}
		System.out.println();*/
	}
	
	
	public static PhyloTree getPMeanViaRandPermCauchy(PhyloTree[] trees, int p, int numIter, int cauchyLength, double epsilon,String outfile, int displayIter, int displayStart) {
		DecimalFormat d6o = new DecimalFormat("#0.######");
		Boolean converge = false;
		int numTrees = trees.length;
		Vector<PhyloTree> oldCentroids = new Vector<PhyloTree>();	  // stores the last 5 centroids form for checking for cauchy sequence		
			
		// initialize the coefficients for the polynomial we solve each iteration
		// make a double[] containing the coefficients
		// constant, deg 1, deg 2, ...
		// p d(T_k,x_{k-1})^{p-2} k t^{p-1} +t -1 =0
		double[] coeffs= new double[p];
//		coeffs[0] = -1;
//		coeffs[1] = 1;
		for(int j = 2; j < p-1; j++) {
			coeffs[j] = 0;
		}
		
		// initialize the solver
		//BrentSolver solver = new BrentSolver();
		LaguerreSolver solver = new LaguerreSolver();
		
		/*  Permute the trees for this round. */
		PhyloTree [] shuffledTrees = permuteTrees(trees);
			
		// find the first centroid candidate
		PhyloTree centroid = shuffledTrees[0];
		oldCentroids.add(centroid.clone());
		int cauchyIndex = 0;  // last position in the vector that there is a old centroid at
			
		int i = 1;
		while ((i < numIter ) && (!converge)) {
			// Shuffle trees if necessary
			if (i % numTrees == 0) {
				shuffledTrees = permuteTrees(trees);
			}
			
			// make a double[] containing the coefficients
			// constant, deg 1, deg 2, ...
			// p d(T_k,x_{k-1})^{p-2} 1/k t^{p-1} +t -1 =0
			
			double d = getGeodesic(centroid,shuffledTrees[i % numTrees],null).getDist();
			
			// need to cast for when p = 2, and need to add one (the coeffs[1] = 1 that we initialized)
			if (p ==2) {
				coeffs[p-1] = p*Math.pow(d,p)/(double)i + 1;
			}
			else {
				coeffs[p-1] = p*Math.pow(d,p)/(double)i;
			}
			coeffs[0] = -Math.pow(d,2);
			coeffs[1] = Math.pow(d,2);
			// initialize the polynomial
			PolynomialFunction func = new PolynomialFunction(coeffs);
			// arguments are:  max number of iterations (100 suggested), functions, two domain values (min and max)
			double t = solver.solve(100,func,0,1);
			
			
			if (i < 2000) {
				System.out.println("p is " + p + ", p-1 coeff is " + coeffs[p-1] + " and t is " + t);
			}
			centroid = getGeodesic(centroid, shuffledTrees[i%numTrees],null).getTreeAt((1-t),centroid.getLeaf2NumMap(),centroid.isRooted());		

			if (displayIter > 0) {
				if ((i%displayIter ==0) &&  (i >= displayStart)) {
					System.out.println("Iteration " + i + ". Variance, tree: " + d6o.format(variance(centroid,trees)) + ", " + centroid.getNewick(true));
				}
			}
			
			/* Check for Cauchy sequence to determine convergence. */
			
			for(int j = cauchyIndex; j >= 0; j--) {
				double dist = calcGeoDist(oldCentroids.get(j), centroid);
				if (dist > epsilon) {
					// don't have a Cauchy sequence.  Move all already checked centroids to the front of the queue.
					for (int k = j; k >= 0; k--) {
						oldCentroids.remove(k);
					}
					break;
				}
			}
			oldCentroids.add(centroid.clone());
			cauchyIndex = oldCentroids.size() - 1;
			
			// Check if we have a Cauchy sequence of the desired length.
			if (oldCentroids.size() == cauchyLength ) {
				converge = true;
			}
			i++;
		}
		
		displayCauchy("random permutations of the tree set.",i,cauchyLength,epsilon,converge,centroid,trees,outfile);

		if (converge) {
			return centroid;
		}
		else {
			return null;
		}
	}
	
}


