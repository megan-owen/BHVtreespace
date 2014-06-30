package pca;

import polyAlg.PolyMain;
import distanceAlg1.PhyloTree;
import distanceAlg1.Geodesic;
import distances.XORShiftRandom;
import distances.Analysis;
import distances.Presentation;
import static polyAlg.PolyMain.getGeodesic;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class LDA {

	
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
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		Boolean rooted = true;   // default
		String treeFile = null;
		String otherTreeFile = null;
		String thirdTreeFile = null;
		String algorithm = "";
		int numIter = 100;   // default
		String outfile = "output.txt";
		PhyloTree[] otherTrees = null;
		PhyloTree[] thirdTrees = null;
		double epsilon = 0.5;
		String[] means = null;
		
		/* Parse command line arguments */

		if (args.length < 1) {
			System.out.println("Error: Missing input file name"); displayHelp(); throw new Exception();
		}
		treeFile = args[args.length-1];
		for (int i = 0; i < args.length - 1; i++) {
			
			if (!args[i].startsWith("-")) { System.out.println("Invalid command line option"); displayHelp(); throw new Exception(); }
			// specify algorithm
			else if (args[i].equals("-a")) {
				if (i < args.length -2) { algorithm = args[i+1]; i++; }
				else { System.err.println("Error: algorithm not specified"); displayHelp(); throw new Exception(); }
			}
			// epsilon, if desired
			else if (args[i].equals("-e")) {
				if (i < args.length -2) { epsilon = Double.valueOf(args[i+1]); i++; }
				else { System.err.println("Error: value for epsilon not specified"); displayHelp(); throw new Exception(); }
			}
			// other tree file, if desired
			else if (args[i].equals("-f")) {
				if (i < args.length -2) { otherTreeFile = args[i+1]; i++; }
				else { System.err.println("Error: name of second tree file not specified"); displayHelp(); throw new Exception(); }
			}
			// third tree file, if desired
			else if (args[i].equals("-g")) {
				if (i < args.length -2) { thirdTreeFile = args[i+1]; i++; }
				else { System.err.println("Error: name of third tree file not specified"); displayHelp(); throw new Exception(); }
			}
			// number of iterations, if desired
			else if (args[i].equals("-i")) {
				if (i < args.length -2) { numIter = Integer.valueOf(args[i+1]); i++; }
				else { System.err.println("Error: number of iterations not specified"); displayHelp(); throw new Exception(); }
			}
			// means for classification
			else if (args[i].equals("-m")) {
				if (i < args.length -2) { means = args[i+1].split(","); i++; }
				else { System.err.println("Error: means not specified"); displayHelp(); throw new Exception(); }
			}
			// output file
			else if (args[i].equals("-o")) {
				if (i < args.length -2) { outfile = args[i+1]; i++; }
				else { System.err.println("Error: output file not specified");  displayHelp(); throw new Exception(); }
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
		
		
		// Read in trees from file.
		PhyloTree[] trees = PolyMain.readInTreesFromFile(treeFile,rooted);
		
		if (otherTreeFile != null) {
			otherTrees = PolyMain.readInTreesFromFile(otherTreeFile,rooted);
		}
		
		if (thirdTreeFile != null) {
			thirdTrees = PolyMain.readInTreesFromFile(thirdTreeFile, rooted);
		}
		
		// choose random pairs of trees to form geodesic and output them, plus LDA score to file
		if (algorithm.equals("random")) {
			// tree file parameters are:  set1, set2, endpoints
			randomLDA(trees,otherTrees, thirdTrees,numIter,epsilon,outfile);
		}
		else if (algorithm.equals("classify")) {
			if (means != null) {
				if (means.length < 2) {
					System.err.println("Error: 2 means must be specified for classify algorithm\n");
					throw new Exception();
				}
				double mean1 = Double.valueOf(means[0]);
				double mean2 = Double.valueOf(means[1]);
				// otherTrees = endpoints file
				classify(trees,otherTrees,mean1, mean2, epsilon,outfile);
			}
			else {
				System.err.println("Error:  no means specified for classify algorithm\n");
				throw new Exception();
			}
		}
		else {
			System.out.println("Error:  no command specified.\n");
			throw new Exception();
		}
		
		System.out.println("Finished");

	}
	
	

/**  Method used in LDA computations for IPMI 2013 paper.
 * 
 *  Randomly selects two trees from geoEndPts, and projects the trees in set1 and set2 onto the
 *  geodesic formed by the randomly selected trees.
 *  Computes the mean and variance of the indices of the projected trees for each of set1 and set2.
 *  Uses this to compute the 'LDA score', and write the two endpoint trees, the score, 
 *  and the projected means to a file.
 *  Repeats the above for numIter iterations, outputting to a new file outfile_i at the i-th iteration.
 *  
 * @param set1
 * @param set2
 * @param geoEndPts
 * @param numIter
 * @param epsilon
 * @param outfile
 * @throws Exception 
 */
	public static void randomLDA(PhyloTree[] set1, PhyloTree[] set2, PhyloTree[] geoEndPts, int numIter, double epsilon, String outfile) throws Exception {
		XORShiftRandom r = new XORShiftRandom();
		PhyloTree t1, t2;
		Vector<String> leaf2NumMap  = set1[0].getLeaf2NumMap();
		Boolean rooted = set1[0].isRooted();
		
		for(int i = 0; i < numIter; i++) {
			// randomly select a pair of trees from the set of possible end point
			int i1 = r.nextInt(geoEndPts.length);
			int i2 = r.nextInt(geoEndPts.length);
			while (i2 == i1) {
       			i2 = r.nextInt(geoEndPts.length);
       		}

			t1 = geoEndPts[i1];
			t2 = geoEndPts[i2];
		
			Geodesic geo = getGeodesic(t1,t2,null);
		
			// get the indices of the projections of set1 and set2 onto the line specified by the chosen pair
			double[] set1Indices = PCA.getProjIndices(geo,set1,epsilon);
			double[] set2Indices = PCA.getProjIndices(geo,set2,epsilon);
		
			// compute the two means and variances of the indices
			double set1Mean = Analysis.mean(set1Indices);
			double set1Scatter = Analysis.scatter(set1Indices);
		
			double set2Mean = Analysis.mean(set2Indices);
			double set2Scatter = Analysis.scatter(set2Indices);
		
			// compute the LDA score
			double numerator = Math.pow(PolyMain.calcGeoDist(geo.getTreeAt(set1Mean,leaf2NumMap, rooted), geo.getTreeAt(set2Mean,leaf2NumMap, rooted)),2);
			double score = numerator/(Math.pow(set1Scatter,2) + Math.pow(set2Scatter,2));
		
			// write the pair of trees and LDA score to outfile
			String outString = t1.getNewick(true) + "\n" + t2.getNewick(true) + "\nScore: " + score + "\nMean 1: " + set1Mean + "\nMean 2: " + set2Mean;
			Presentation.printStringToFile(outString, outfile + "_" + i);
			
		}
	}
	
	/** Classifies the trees in trees by which mean they are closest too,
	 * on the geodesic between the two trees given in endpoints.
	 * 
	 * @param trees
	 * @param endpoints
	 * @param mean1
	 * @param mean2
	 * @param epsilon
	 * @param outfile
	 * @throws Exception 
	 */
	public static void classify(PhyloTree[] trees, PhyloTree[] endpoints, double mean1, double mean2, double epsilon, String outfile) throws Exception {
		Geodesic geo = getGeodesic(endpoints[0],endpoints[1],null);
		
		// open up outfile for writing
		PrintWriter outputStream = null;
			  
		try {
			outputStream = new PrintWriter(new FileWriter(outfile));
		
			//  for each tree in trees
			for (PhyloTree t: trees) {
		
				// project onto the geodesic
				double index = t.projectToGeoIndex(geo,epsilon);
		
				// figure out which mean is closer
				double abs1 = Math.abs(mean1 - index);
				double abs2 = Math.abs(mean2 - index);
				
				if (abs1 < abs2) {
					outputStream.print("1\n");
				}
				else {
					outputStream.print("2\n");
				}
			}		
			if (outputStream != null) {
       			outputStream.close();
       		}
		} 
		catch (FileNotFoundException e) {
			System.out.println("Error opening or writing to " + outfile + ": "+ e.getMessage());
    		throw new Exception();
		} 
		catch (IOException e) {
			System.out.println("Error opening or writing to " + outfile + ": " + e.getMessage());
			throw new Exception();
		} // end try/catch
	}
	

}
