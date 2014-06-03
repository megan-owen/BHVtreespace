package distances;

import polyAlg.PolyMain;
import distanceAlg1.Geodesic;
import distanceAlg1.PhyloTree;
import distanceAlg1.TreeDistance;
import java.util.*;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.*;

public class Experiments {
	public static int verbose = 0;


	

	
	public static void displayHelp() {
		System.out.println("Command line syntax:");
		System.out.println("distance [options] treefile");
		System.out.println("Optional arguments:");	
		System.out.println("\t -h \t displays this message");
		System.out.println("\t -o <outfile> \t store the output in the file <outfile>");
		System.out.println("\t -u \t unrooted trees (default is rooted trees)");
		System.out.println("\t -v \t verbose output");
	}
	
	
	/**
	 * @param args
	 */
public static void main(String[] args) {
	String treeFile = "";
	
	String outFile = "exp_out_" + getTimeStamp(); // default
	boolean rooted = true;
	

	/*  Parse Arguments */
	if (args.length < 1) {
		System.out.println("Error: Missing input file name"); displayHelp(); System.exit(1);
	}
	treeFile = args[args.length-1];
	for (int i = 0; i < args.length - 1; i++) {		
		if (!args[i].startsWith("-")) {
			System.out.println("Invalid command line option"); displayHelp(); System.exit(1);
		}
		// output file
		else if (args[i].equals("-o")) {
			if (i < args.length -2) {
				outFile = args[i+1]; i++;
			}
			else { System.out.println("Error: Output file not specified"); displayHelp(); System.exit(0); }
		}
		// all other arguments.  Note we can have -vn
		else {
			for (int j = 1; j<args[i].length(); j++) {			
				switch(args[i].charAt(j)) {						
					// display help
					case 'h': Distances.displayHelp(); System.exit(0); break;										
					// unrooted trees?
					case 'u': rooted = false; break;
					// verbose output
					case 'v': verbose = 1; break;	
					default: System.out.println("Illegal command line option.\n"); displayHelp(); System.exit(1); break;
				} // end switch
			} // end for j
		} // end parsing an individual argument
	}  // end for i (looping through arguments)
	/*  End Parse Arguments */


	// Read in trees from file. 
//	PhyloTree[] trees = PolyMain.readInTreesFromFile(treeFile,rooted);
	
	try {
		// File containing trees.  Read one per line.
        BufferedReader inStream = new BufferedReader(new FileReader(treeFile));
        

	// Open output file. 
    	PrintWriter outStream = new PrintWriter(new FileWriter(outFile));
	
    	//  Experiment  
    	
    	PhyloTree speciesTree = (PolyMain.readInTreesFromFile("/Users/meganowen/projects/gene_and_species_trees/experiments/2011-07-22/left128_theta0.01.txt", true))[0];
    	
    	// Print file header
    	outStream.println("RF,weighted RF,geodesic,# common edges,# topologies in geodesic,siDisjointNoLeaves");
    	
    	String line;
    	while ((line = inStream.readLine()) != null) {
    		PhyloTree tree = new PhyloTree(line,rooted);
    		Geodesic geo = PolyMain.getGeodesic(speciesTree,tree,null);

    		PhyloTree disjointSpeciesTree = new PhyloTree(speciesTree);
    		disjointSpeciesTree.setEdges(disjointSpeciesTree.getEdgesIncompatibleWith(tree));
    		
    		PhyloTree disjointGeneTree = new PhyloTree(tree);
    		disjointGeneTree.setEdges(disjointGeneTree.getEdgesIncompatibleWith(speciesTree));
   
    		outStream.println("" + Distances.rf(speciesTree,tree) + "," + Distances.weightedRF(speciesTree,tree) + "," + geo.getDist() + "," + geo.numCommonEdges() + "," + geo.numTopologies() + "," + Distances.siDisjointNoLeaves(disjointSpeciesTree, disjointGeneTree) );
    	}
    	
    	
    	
    	if (outStream != null) {
			outStream.close();
		}
    	inStream.close();
    } catch (FileNotFoundException e) {
    	System.out.println("In Experiments, error opening or writing file: "+ e.getMessage());
    	System.exit(1);
    } catch (IOException e) {
    	System.out.println("In Experiments, error opening or writing file: " + e.getMessage());
    	System.exit(1);
    }

    
    /* Remove common splits and write to output file.
	*PhyloTree[] treesNoCommon = Analysis.removeCommonSplits(trees);
	*Presentation.printTreesToStream(treesNoCommon, outStream,true);
	*/

}



/** Returns the present date and time nicely formatted.
 *  (For use in automatically naming outfiles.)  
 * 
 * @return
 */
public static String getTimeStamp() {
	Date date = new Date();
	
	Format formatter = new SimpleDateFormat("yyMMdd-H'h'mmssSSS");
	return formatter.format(date);
}



	public static void testDec6 (PhyloTree[] trees) {
		int numTrees = trees.length;
		int numRatios, numCommonEdges;
		int numDists = 0;
		double normalizedRF = 0;
		double avgNormRF = 0;
		double weightedRF = 0;
		double avgWeightedRF = 0;
		double geoDist, si2;
		double averageDist = 0;
		double avgSI2 =0;
		double avgNumCommonEdges = 0;
		double avgNumRatios = 0;
		
		
		for (int i = 0; i < numTrees - 1; i++) {
			for (int j = i + 1; j < numTrees; j++) {
				Geodesic geo = PolyMain.getGeodesic(trees[i], trees[j], null);
				numRatios = geo.getRS().size();
				geoDist = geo.getDist();
				si2 = Distances.si2(trees[i], trees[j]);
				normalizedRF = Distances.rf(trees[i], trees[j])/2;
				weightedRF = Distances.weightedRF(trees[i], trees[j]);
				
				numCommonEdges = geo.numCommonEdges();
				
				avgNumRatios = avgNumRatios + numRatios;
				averageDist = averageDist + geoDist;
				avgSI2 = avgSI2 + si2;
				avgNumCommonEdges = avgNumCommonEdges + numCommonEdges;
				avgNormRF = avgNormRF + normalizedRF;
				avgWeightedRF = avgWeightedRF + weightedRF;
				numDists++;
				System.out.println("Tree " + i + " to " + j + ": geo dist = " + geoDist + ", # ratios = " + numRatios + ", si2 = " + si2 + ", # common edges = " + numCommonEdges + ", normalized RF = " + normalizedRF + ", weightedRF = " + weightedRF);
			}
		}
		averageDist = averageDist/numDists;
		avgSI2 = avgSI2/numDists;
		avgNumCommonEdges = avgNumCommonEdges/numDists;
		avgNormRF = avgNormRF/numDists;
		avgWeightedRF = avgWeightedRF/numDists;
		avgNumRatios = avgNumRatios/numDists;
		
		
		System.out.println("Averages: geo dist = " + averageDist + ", # ratios = " + avgNumRatios + ", si2 = " + avgSI2 + ", # common edges = " + avgNumCommonEdges + ", normalized RF = " + avgNormRF + ", weighted RF = " + avgWeightedRF);
		
	}
	
	/*  Testing done Feb. 16 -, 2011. */
	public static void test110216 (PhyloTree[] trees) {
		PhyloTree t1, t2;
		Geodesic geo;
		
		
		for (int i = 0; i < trees.length; i++) {
			for (int j = 0; j <= i; j++) { 
				Vector<PhyloTree> disjointTreePairs = TreeDistance.splitOnCommonEdge(trees[i], trees[j]);
				
				for (int k = 0; k < disjointTreePairs.size() /2 ;k++) {
					//  Compute the SI index.
					t1 = (PhyloTree) disjointTreePairs.get(2*k);
					t2 = (PhyloTree) disjointTreePairs.get(2*k+1);
					geo = PolyMain.getGeodesicNoCommonEdges(t1, t2);
					
					System.out.println("# leaves = " + t1.getLeaf2NumMap().size() + "; # edges = " + t1.getEdges().size() + "; # ratios = " + geo.getRS().size() + "; SI (disjoint, no leaves) = " + Distances.siDisjointNoLeaves(t1,t2));
				}
			}
		}
		
		
		
	}


}

	
