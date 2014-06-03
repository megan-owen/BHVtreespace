/** This file is part of GTP, a program for computing the geodesic distance between phylogenetic trees,
 * and sturmMean, a program for computing the Frechet mean between phylogenetic trees.
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

package polyAlg;


import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Random;
import java.util.Vector; 

import distanceAlg1.EdgeAttribute;
import distanceAlg1.Geodesic;
import distanceAlg1.PhyloTree;
import distanceAlg1.PhyloTreeEdge;
import distanceAlg1.RatioSequence;
import distanceAlg1.Ratio;

//import org.biojavax.bio.phylo.io.nexus.*;
//import org.biojava.bio.seq.io.ParseException;

public class PolyMain {
	// stores pairs of trees with no common edges.  Should be reset at each new distance calculation
	public static Vector<PhyloTree> aTreesNoCommonEdges = new Vector<PhyloTree>();
	public static Vector<PhyloTree> bTreesNoCommonEdges = new Vector<PhyloTree>();
	
//	public static boolean rooted = true;  //holds if the trees are rooted or not.
	public static boolean normalize = false;  // holds if we should normalize the tree split lengths
	
	public static int verbose = 0;
	
	public static String LEAF_CONTRIBUTION_SQUARED_DESCRIPTION = "(Leaf contribution squared = square of the length of the vector" +
	" whose i-th element is the absolute value of the difference between the length of the split ending in leaf i in the first tree" +
	" and the length of the split ending in leaf i in the second tree.)";
	
/** Stores subtrees with no common edges in the global variables aTreesNoCommonEdges (from t1)
 * and in bTreesNoCommonEdges (from t2).
 * 
 * If one of the trees has 0 edges, then all edges in the other tree will be compatible with it.
 * Thus we will not get any more subtree pairs with disjoint leaves, and should return.
 * 
 */
public static void splitOnCommonEdge(PhyloTree t1, PhyloTree t2) {
	int numEdges1 = t1.getEdges().size(); // number of edges in tree 1
	int numEdges2 = t2.getEdges().size(); /// number of edges in tree 2

	if (numEdges1 == 0 || numEdges2 == 0) {
		return;
	}	
	// look for common edges
	Vector<PhyloTreeEdge> commonEdges = PhyloTree.getCommonEdges(t1, t2);	
	
	// if there are no common edges
	// XXX: need to check the following methods don't require the trees to have the same number of edges
	if (commonEdges.size() == 0) {
//			System.out.println("In splitOnCommonEdges, no common edges in " + t1 + " and " + t2);
		aTreesNoCommonEdges.add(t1);
		bTreesNoCommonEdges.add(t2);
		return;
	}
//		System.out.println("At least one common split; edges are " + commonEdges);
	
	// else if there exists a common split: split the trees along the first split in commonEdges
	// and recursively call getDistance for the two new pairs of trees.
	PhyloTreeEdge commonEdge = commonEdges.get(0);
//		System.out.println("Common edges is " + commonEdge);
	
	Vector<PhyloTreeEdge> edgesA1 = new Vector<PhyloTreeEdge>();
	Vector<PhyloTreeEdge> edgesA2 = new Vector<PhyloTreeEdge>();
	Vector<PhyloTreeEdge> edgesB1 = new Vector<PhyloTreeEdge>();
	Vector<PhyloTreeEdge> edgesB2 = new Vector<PhyloTreeEdge>();
	
//		System.out.println("tree 1 is " + t1);
	
	for (PhyloTreeEdge e:  t1.getEdges()) {
//		System.out.println("t1.getEdge(i).getOriginalID() is " + t1.getEdge(i).getOriginalID()  +" \n");
//		edgesA1.add(new PhyloTreeEdge(t1.getEdge(i).getLength(), t1.getEdge(i).getOriginalEdge(), t1.getEdge(i).getOriginalID() ));
		
		
		// tree A is the tree under the common edge (i.e. the common edge is the root)
		if (commonEdge.properlyContains(e)) {
			// add to A tree
			edgesA1.add(e.clone());
		}
		// tree B is the tree above the common edge (i.e. think of all leaves below the common edge as one big leaf)
		else if (!e.sameBipartition(commonEdge)) {
			edgesB1.add(e.clone());
		}
	}
	
	for (PhyloTreeEdge e : t2.getEdges()) {
		// tree A is the tree under the common edge (i.e. the common edge is the root)
		if (commonEdge.properlyContains(e)) {
			// add to A tree
			edgesA2.add(e.clone());
		}
		// tree B is the tree above the common edge (i.e. think of all leaves below the common edge as one big leaf)
		else if (!e.sameBipartition(commonEdge)) {
			edgesB2.add(e.clone());
		}
	}

	
	// make the 4 trees
	PhyloTree tA1 = new PhyloTree(edgesA1, Tools.myVectorCloneString(t1.getLeaf2NumMap()), t1.isRooted());
	PhyloTree tB1 = new PhyloTree(edgesB1, Tools.myVectorCloneString(t1.getLeaf2NumMap()), t1.isRooted());
	PhyloTree tA2 = new PhyloTree(edgesA2, Tools.myVectorCloneString(t1.getLeaf2NumMap()), t1.isRooted());
	PhyloTree tB2 = new PhyloTree(edgesB2, Tools.myVectorCloneString(t1.getLeaf2NumMap()), t1.isRooted());
	
	//System.out.println("Subtrees below common edge are A1 = " + tA1 + " and A2 = " + tA2 + "...");
	//System.out.println("Corresponding B subtrees are B1 = " + tB1 + " and B2 = " + tB2);
	splitOnCommonEdge(tA1, tA2);
	

	splitOnCommonEdge(tB1,tB2);
}	

	
/** Reads in all the phylogenetic trees from the file inFileName.
 *  The trees should be one per line, in the Newick format.
 *  There can also be a ";" at the end of each line.
 *  Exits with error message if there is an error in the Newick of any tree.
 * @param inFileName
 * @return
 */
public static PhyloTree[] readInTreesFromFile (String inFileName, boolean rooted) {
	int numTrees =0;  // count the number of trees read in
	boolean nexus = false;
	Vector<String> stringTrees = new Vector<String>();
	
	BufferedReader inputStream = null;

	// Check to see if the file is in Nexus format.
	try {
		inputStream = new BufferedReader(new FileReader(inFileName));

        String l = inputStream.readLine();
        if (l.equals("#NEXUS")) {
        	nexus = true;
        }
        inputStream.close();
	} catch (FileNotFoundException e) {
        System.out.println("Error opening or reading from " + inFileName + ": " + e.getMessage());
        System.exit(1);
    } catch (IOException e) {
    	System.out.println("Error opening or reading from " + inFileName + ": " + e.getMessage());
    	System.exit(1);
    }
	
    // File is not in NEXUS format.  Assume it is a list of trees in Newick format,
    // possibly separated by blank lines.
    if (nexus == false) {
    try {
        inputStream = new BufferedReader(new FileReader(inFileName));

        String l;
        while ((l = inputStream.readLine()) != null) {
        	// check for blank lines
        	if (!(l.equals("")) && !(l.equals("\n")) ) {
        		if (l.charAt(0) != '#') {  // also ignore lines beginning with # as comment lines
        			// If there is a ( in the line, assume there is a tree starting there
        			// and remove all text before.
        			// Otherwise, ignore the line.
        			int index = l.indexOf("(");
        			if (!(index == -1)) {
        				l = l.substring(index);
        				stringTrees.add(l);
        				numTrees++;
        			}
        		}
        	}
        }
        if (inputStream != null) {
            inputStream.close();
        }
    } catch (FileNotFoundException e) {
        System.out.println("Error opening or reading from " + inFileName + ": " + e.getMessage());
        System.exit(1);
    } catch (IOException e) {
    	System.out.println("Error opening or reading from " + inFileName + ": " + e.getMessage());
    	System.exit(1);
    }
    }
    else {
    	System.out.println("Nexus format not supported.");
    	System.exit(0);
    }
    

    
    // convert the Vector of strings representing trees into an array of PhyloTrees
    PhyloTree[] trees = new PhyloTree[numTrees];
    int treesWithErrors = 0;
    for (int i = 0; i < numTrees; i++) {
    	// check for syntax errors in the string
    	if (PhyloTree.errorInSyntax(stringTrees.get(i))) {
    		treesWithErrors++;
    		System.out.println("Error in tree " + i);
    	}
      	trees[i] = new PhyloTree(stringTrees.get(i),rooted);
       	if (PolyMain.normalize) {
       		trees[i].normalize();
       	}
    }
    if (treesWithErrors > 0) {
    	System.out.println("Exiting: errors in " + treesWithErrors + " trees");
    	System.exit(1);
    }
    
    // verify all trees have the same leaf set
    if (numTrees > 1) {
    	Vector<String> leaf2NumMap = trees[0].getLeaf2NumMap();
    	for (int i = 1; i < numTrees; i++) {
    		if (!leaf2NumMap.equals(trees[i].getLeaf2NumMap())) {
    			System.out.println("Warning:  tree at line " + (i + 1) + " does not have same leaves as first tree in file");
    			System.out.println("Line 1 tree leaf set: " + leaf2NumMap );
    			System.out.println("Line " + (i+1) + " tree leaf set: " + trees[i].getLeaf2NumMap());
    			System.exit(1);
    		}
    	}
    }
    return trees;
}

/** Returns the distance between t1 and t2, accounting for any common edges and leaf edges.
 *  Calls recursive getGeodesic
 *  Does not assume t1 and t2 have the same number of edges.
 *  Pass in null for geoFile to not write to a file.
 * 
 */
public static Geodesic getGeodesic(PhyloTree t1, PhyloTree t2, String geoFile) {
	double leafContributionSquared = 0;
	EdgeAttribute [] t1LeafEdgeAttribs = t1.getLeafEdgeAttribs();
	EdgeAttribute [] t2LeafEdgeAttribs = t2.getLeafEdgeAttribs();
	Geodesic geo = new Geodesic(new RatioSequence(),t1.getLeafEdgeAttribs(),t2.getLeafEdgeAttribs());
	
	String verboseOutput = "";
	
	// get the leaf contributions
	for(int i = 0; i < t1.getLeaf2NumMap().size(); i++ ) {
		if ( !(t1.getLeaf2NumMap().get(i).equals(t2.getLeaf2NumMap().get(i)) ) ) {
			System.out.println("Error getting geodesic: trees do not have the same sets of leaves");
			System.out.println("Starting tree leaves: " + t1.getLeaf2NumMap());
			System.out.println("Target tree leaves: " + t2.getLeaf2NumMap());
			
			System.out.println("Starting tree: " + t1.getNewick(true));
			System.out.println("Target tree: " + t2.getNewick(true));
			
			System.exit(1);
		}
//		System.out.println("leaf: " + t1.getLeaf2NumMap().get(i) + " | " + t1LeafEdgeLengths[i] + " - " + t2LeafEdgeLengths[i] + "| = " + (t1LeafEdgeLengths[i] - t2LeafEdgeLengths[i]) );

//		leafContributionSquared = leafContributionSquared+ Math.pow(Math.abs(t1LeafEdgeLengths[i] - t2LeafEdgeLengths[i]), 2);
		leafContributionSquared = leafContributionSquared+ Math.pow(EdgeAttribute.difference(t1LeafEdgeAttribs[i],t2LeafEdgeAttribs[i]).norm(), 2);
	}
	geo.setLeafContributionSquared(leafContributionSquared);
	
	if (PolyMain.verbose == 1) {
		System.out.println("Starting tree: " + t1.getNewick(true));
		verboseOutput = verboseOutput + "Starting tree: " + t1.getNewick(true) + "\n";
		System.out.println("Target tree: " + t2.getNewick(true));
		verboseOutput = verboseOutput + "Target tree: " + t2.getNewick(true) + "\n";
		
		System.out.println("\nStarting tree edges:");
		verboseOutput = verboseOutput + "\nStarting tree edges:\n";
		verboseOutput = verboseOutput + PhyloTreeEdge.printEdgesVerbose(t1.getEdges(), t1.getLeaf2NumMap(), true);
		
		System.out.println("\nTarget tree edges:");
		verboseOutput = verboseOutput + "\nTarget tree edges:\n"; 
		verboseOutput = verboseOutput + PhyloTreeEdge.printEdgesVerbose(t2.getEdges(), t2.getLeaf2NumMap(), true);
		
		// leaf contributions
		System.out.println("\nLeaf contribution squared " + Tools.truncate(leafContributionSquared,6));
		verboseOutput = verboseOutput + "\nLeaf contribution squared " + Tools.truncate(leafContributionSquared,6) + "\n";
		System.out.println(LEAF_CONTRIBUTION_SQUARED_DESCRIPTION);
		verboseOutput = verboseOutput + LEAF_CONTRIBUTION_SQUARED_DESCRIPTION + "\n";
	}
	
	aTreesNoCommonEdges = new Vector<PhyloTree>();
	bTreesNoCommonEdges = new Vector<PhyloTree>();
	
	// get the pairs of trees with no common edges put into aTreesNoCommonEdges and bTreesNoCommonEdges
	//  aTreesNoCommonEdges.get(i) goes with bTreesNoCommonEdges.get(i)
	splitOnCommonEdge(t1,t2);
	
	//set the common edges
	Vector<PhyloTreeEdge> commonEdges = PhyloTree.getCommonEdges(t1,t2);
	geo.setCommonEdges(commonEdges);
	
	// set the t1 and t2 attributes of the common edges
	Vector<PhyloTreeEdge> eCommonEdges = Tools.myVectorClonePhyloTreeEdge(commonEdges);
	Vector<PhyloTreeEdge> fCommonEdges = Tools.myVectorClonePhyloTreeEdge(commonEdges);
	for (int i = 0; i < eCommonEdges.size(); i++) {
		eCommonEdges.get(i).setAttribute(t1.getAttribOfSplit(eCommonEdges.get(i)));
		fCommonEdges.get(i).setAttribute(t2.getAttribOfSplit(fCommonEdges.get(i)));
	}
	geo.seteCommonEdges(eCommonEdges);
	geo.setfCommonEdges(fCommonEdges);
	
	if (verbose == 1) {
		System.out.println("\nCommon edges are:  (Length = abs. value of difference in length between the two trees)");
		verboseOutput = verboseOutput +"\nCommon edges are:  (Length = abs. value of difference in length between the two trees)\n";
		
		verboseOutput = verboseOutput + PhyloTreeEdge.printEdgesVerbose(commonEdges, t1.getLeaf2NumMap(), true);
		double commonEdgeContributionSquared = 0;
		for (int i = 0; i < commonEdges.size(); i++) {
				commonEdgeContributionSquared = commonEdgeContributionSquared + Math.pow(commonEdges.get(i).getNorm(),2);
		}
		System.out.println("\nCommon edges contribution squared: " + Tools.truncate(commonEdgeContributionSquared, 6));
		verboseOutput = verboseOutput +"\nCommon edges contribution squared: " + Tools.truncate(commonEdgeContributionSquared, 6) + "\n";
		System.out.println("(sum of squares of above differences in length)");
		verboseOutput = verboseOutput +"(sum of squares of above differences in length)\n";
		System.out.println("=============================================================================================================================");
		verboseOutput = verboseOutput +"=============================================================================================================================\n";
		
		System.out.println("\nNow finding the geodesic between the following subtrees, which have no edges in common:");
		verboseOutput = verboseOutput +"\nNow finding the geodesic between the following subtrees, which have no edges in common:\n";
	}
		
	// find the geodesic between each pair of subtrees found by removing the common edges
	for(int i = 0; i < aTreesNoCommonEdges.size(); i++) {
		PhyloTree subTreeA = aTreesNoCommonEdges.get(i);
		PhyloTree subTreeB = bTreesNoCommonEdges.get(i);
		
		if (verbose == 1) {
			//  Commented out Nov. 8, 2012, as use parent leafset for the subtrees now
//			System.out.println("Leaves or subtree representatives in subtrees:");
//			verboseOutput = verboseOutput +"Leaves or subtree representatives in subtrees:\n";
//			for (int j = 0; j <subTreeA.getLeaf2NumMap().size(); j++) {
//				System.out.println("" + subTreeA.getLeaf2NumMap().get(j));
//				verboseOutput = verboseOutput + subTreeA.getLeaf2NumMap().get(j) + "\n";
//			}
			
			System.out.println("\nStarting subtree edges:");
			verboseOutput = verboseOutput + "\nStarting subtree edges:\n";
			verboseOutput = verboseOutput + PhyloTreeEdge.printEdgesVerbose(subTreeA.getEdges(), subTreeA.getLeaf2NumMap(), false);
			
			System.out.println("\nTarget subtree edges:");
			verboseOutput = verboseOutput + "\nTarget subtree edges:\n";
			verboseOutput = verboseOutput + PhyloTreeEdge.printEdgesVerbose(subTreeB.getEdges(),subTreeB.getLeaf2NumMap(), false);
		}
		
		Geodesic newGeo = getGeodesicNoCommonEdges(subTreeA, subTreeB);
		
		if (verbose == 1) {
			System.out.println("\nGeodesic distance between above subtrees, ignoring edges ending in leaves: " + Tools.truncate(newGeo.getRS().getNonDesRSWithMinDist().getDistance(), 6) );
			verboseOutput = verboseOutput + "\nGeodesic distance between above subtrees, ignoring edges ending in leaves: " + Tools.truncate(newGeo.getRS().getNonDesRSWithMinDist().getDistance(), 6) + "\n";
			
			System.out.println("Ratio sequence corresponding to the geodesic:\nCombinatorial type: " + newGeo.getRS().getNonDesRSWithMinDist().toStringCombType() );
			verboseOutput = verboseOutput + "Ratio sequence corresponding to the geodesic:\nCombinatorial type: " + newGeo.getRS().getNonDesRSWithMinDist().toStringCombType() + "\n";
		    
			System.out.println(newGeo.getRS().getNonDesRSWithMinDist().toStringVerbose(subTreeA.getLeaf2NumMap()) );
			verboseOutput = verboseOutput + newGeo.getRS().getNonDesRSWithMinDist().toStringVerbose(subTreeA.getLeaf2NumMap()) + "\n";
			
			System.out.println("------------------------------------------------------------------------------------------------------------");
			verboseOutput = verboseOutput + "------------------------------------------------------------------------------------------------------------\n";
		}
		
		geo.setRS(RatioSequence.interleave(geo.getRS(), newGeo.getRS()));
	}

	if (verbose == 1) {
		 System.out.println("\nGeodesic distance between start and target tree is " + Tools.truncate(geo.getDist(),6));
		 verboseOutput = verboseOutput + "\nGeodesic distance between start and target tree is " + Tools.truncate(geo.getDist(),6) + "\n";
	
	//	 verboseOutput = verboseOutput + "\n\nGeodesic         : " + geo;
	//	 verboseOutput = verboseOutput + "\n\nGeodesic reversed: " + geo.reverse();
	}
	
	// write verbose output to geofile, if in verbose mode
	if (verbose == 1 && geoFile != null) {
        PrintWriter outputStream = null;
        
        try {
        	outputStream = new PrintWriter(new FileWriter(geoFile));
        	
        	outputStream.println(verboseOutput);
        
    		if (outputStream != null) {
                outputStream.close();
            }
        } catch (FileNotFoundException e) {
            System.out.println("Error opening or writing to " + geoFile + ": "+ e.getMessage());
            System.exit(1);
        } catch (IOException e) {
        	System.out.println("Error opening or writing to " + geoFile + ": " + e.getMessage());
        	System.exit(1);
        }
	}
	
	// write alternative output to geofile if selected -w switch
	if (verbose ==2 && geoFile != null) {
		PrintWriter outputStream = null;
        
        try {
        	outputStream = new PrintWriter(new FileWriter(geoFile));
        	
        	outputStream.println("Starting tree: " + t1.getNewick(true));
        	outputStream.println("Target tree: " + t2.getNewick(true));
        	outputStream.println();
        	
        	outputStream.println("Boundary trees are: ");
        	Vector<PhyloTree> boundaryTrees = Geodesic.getBoundaryTrees(t1,t2);
        	for (PhyloTree t: boundaryTrees) {
        		outputStream.println("" + t.getNewick(true));
        	}
        	
        	outputStream.println();
        	
        	outputStream.println("Topology sequence of geodesic: ");
        	outputStream.println("" + t1.getNewick(false));
        	for (PhyloTree t: geo.getMidOrthantTrees(t1.getLeaf2NumMap(), t1.isRooted())) {
        		outputStream.println("" + t.getNewick(false));
        	}
        	outputStream.println("" + t2.getNewick(false));
        
    		if (outputStream != null) {
                outputStream.close();
            }
        } catch (FileNotFoundException e) {
            System.out.println("Error opening or writing to " + geoFile + ": "+ e.getMessage());
            System.exit(1);
        } catch (IOException e) {
        	System.out.println("Error opening or writing to " + geoFile + ": " + e.getMessage());
        	System.exit(1);
        }
	}
	return geo;
}

/** Returns the geodesic between t1 and t2, which are assumed to have no common edges.
 *  Does not assume t1 and t2 have the same number of edges.
 *  Does not take into account the leaf edges.
 *  Uses polynomial algorithm.
 *  XXX: how to deal with multifurcating trees
 * 
 *  Returns:  a Geodesic with just the ratio sequence set 
 */
public static Geodesic getGeodesicNoCommonEdges(PhyloTree t1, PhyloTree t2 ) {
	int numEdges1 = t1.getEdges().size(); // number of edges in tree 1
	int numEdges2 = t2.getEdges().size(); // number of edges in tree 2
	RatioSequence rs = new RatioSequence();
	int [] aVertices, bVertices;
	Vector<Ratio> queue = new Vector<Ratio>();
	Ratio ratio;
	int[][] cover;
	
	if (numEdges1 == 0 && numEdges2 == 0) {
		return new Geodesic(new RatioSequence());
	}
		
	// double-check no common edges
	Vector<PhyloTreeEdge> commonEdges = PhyloTree.getCommonEdges(t1, t2);
	if (commonEdges.size() != 0) {
		System.out.println("Exiting: tried to compute geodesic between subtrees that should not have common edges, but do!  t1 = " + t1 + " and t2 = " + t2);
		System.exit(1);
	}
	
	// double-check that both trees have splits.  Otherwise didn't remove a common edge.
	if (numEdges1 ==0 || numEdges2 == 0) {
		System.out.println("Exiting: tried to compute geodesic between subtrees that should not have common/compatible edges, but do!  t1 = " + t1 + " and t2 = " + t2);
		System.exit(1);
	}
	
	// if we can't split the ratio because it has too few edges in either the numerator or denominator
	if ((numEdges1 ==1) || (numEdges2 ==1)) {
		rs.add( new Ratio(t1.getEdges(), t2.getEdges()) );
		return new Geodesic(rs);
	}
	
	// initialize BipartiteGraph
	boolean[][] incidenceMatrix = Tools.getIncidenceMatrix(t1.getEdges(), t2.getEdges());

	
	BipartiteGraph bg = new BipartiteGraph(incidenceMatrix, t1.getIntEdgeAttribNorms(), t2.getIntEdgeAttribNorms());
	
	queue.add(new Ratio(t1.getEdges(), t2.getEdges() ));
	
	while(queue.size() >0) {
		ratio = queue.remove(0);
		
		aVertices = new int[ratio.getEEdges().size()];
		bVertices = new int[ratio.getFEdges().size()];
		
		// convert the ratio to what we pass to vertex cover
		for (int i = 0; i < ratio.getEEdges().size(); i++) {
			aVertices[i] = t1.getEdges().indexOf(ratio.getEEdges().get(i));
		}
		for (int i = 0; i < ratio.getFEdges().size(); i++) {
			bVertices[i] = t2.getEdges().indexOf(ratio.getFEdges().get(i));
		}
		
		// get the cover
		cover = bg.vertex_cover(aVertices, bVertices);
		
		// check if cover is trivial
		if ( (cover[0][0] == 0) || (cover[0][0] == aVertices.length) ){
			// add ratio to geodesic
			rs.add(ratio);

			
		}
		else {  // cover not trivial
			// make two new ratios
			Ratio r1 = new Ratio();
			Ratio r2 = new Ratio();
			
			int j = 0;  // for index in cover array
			
			// split the ratio based on the cover
			for (int i = 0; i < aVertices.length; i++) {
				if ( (j < cover[2].length) && (aVertices[i] == cover[2][j]) ) {
					r1.addEEdge( t1.getEdge(aVertices[i]) );
					j++;
				}
				else { // the split is not in the cover, and hence dropped first
					r2.addEEdge( t1.getEdge(aVertices[i]) );
				}
			}
			
			j = 0;   // reset index
			// split the ratio based on the cover
			for (int i = 0; i < bVertices.length; i++) {
				if ( (j < cover[3].length) && (bVertices[i] == cover[3][j]) ) {	
					r2.addFEdge( t2.getEdge(bVertices[i]) );
					j++;
				}
				else { // the split is not in the cover, and hence dropped first
					r1.addFEdge( t2.getEdge(bVertices[i]) );
				}
			}
			
			// add ratios to the queue
			queue.add(0, r2);
			queue.add(0, r1);
		}
	}
	
	return new Geodesic(rs);
}



/** Computes all the inter-tree distances between the trees in trees using algorithm, and returns them in a 
 * matrix.  Prints the average time on the screen.  If doubleCheck is true, computes each distance both ways,
 * and displays a message if they differ.
 * @param trees
 * @param algorithm
 * @return
 */
public static Geodesic[][] getAllInterTreeGeodesics(PhyloTree[] trees, boolean doubleCheck) {
	Date startTime;
	Date endTime;
	int numTrees = trees.length;
	long[][] compTimes = new long[numTrees][numTrees];
	long avgCompTime = 0;
	double[][] dists = new double[numTrees][numTrees];
	Geodesic[][] geos = new Geodesic[numTrees][numTrees];
	
	for (int i = 0; i < numTrees; i++) {
		for (int j = i+1 ; j < numTrees;j++) {
			startTime = new Date();
			geos[i][j] = getGeodesic(trees[i], trees[j], "geo_" + i + "_" + j);
			dists[i][j] = geos[i][j].getDist();
			endTime = new Date();
			compTimes[i][j] = endTime.getTime() - startTime.getTime();
			// sum up all the times, then divide by number of computations
			avgCompTime = avgCompTime + compTimes[i][j];
			//	System.out.println("computed in time " + compTimes[i][j] + " geos[" + i + "][" + j + "].getRS().getMinAscRSDistance() is " + geos[i][j].getRS().getMinAscRSDistance() + "; commonEdges is " +geos[i][j].getCommonEdges() + "; and leafContributionSquared is " + geos[i][j].getLeafContributionSquared());
		}
	}
	
	// compute average
	avgCompTime = avgCompTime/ (numTrees * (numTrees - 1)/2);
	System.out.println("Average dist. computation was " + avgCompTime + " ms for " + numTrees * (numTrees - 1)/2 + " trees.");
			
	// we want to doublecheck
	if (doubleCheck) {
		avgCompTime = 0;
		for (int i = 0; i < numTrees; i++) {
			for (int j = i+1 ; j < numTrees;j++) {
				startTime = new Date();
				geos[j][i] = getGeodesic(trees[j], trees[i], "geo_" + j + "_" + i);
				dists[j][i] = geos[j][i].getDist();
				endTime = new Date();
				compTimes[j][i] = endTime.getTime() - startTime.getTime();
				// sum up all the times, then divide by number of computations
				avgCompTime = avgCompTime + compTimes[j][i];

				if (Tools.truncate(geos[i][j].getDist(), 10) != Tools.truncate(geos[j][i].getDist(),10) ) {
					System.out.println("*** Distances don't match for trees " + i + " and " + j + "***");
					System.out.println("Dist " + i + " -> " + j + " is " + geos[i][j].getDist() + " but dist " + j + " -> " + i+ " is " + geos[j][i].getDist());
					System.out.println("RS " + i + " -> " + j + "           : " + geos[i][j]);
					System.out.println("geos[" + i + "][" + j + "].getRS().getAscRSWithMinDist().getDistance() is " + geos[i][j].getRS().getNonDesRSWithMinDist().getDistance() + "; commonEdges is " +geos[i][j].getCommonEdges() + "; and leafContributionSquared is " + geos[i][j].getLeafContributionSquared());

					
					System.out.println("RS " + j + " -> " + i + " (reversed): " + geos[j][i].reverse());
					System.out.println("geos[" + j + "][" + i + "].getRS().getAscRSWithMinDist().getDistance() is " + geos[j][i].getRS().getNonDesRSWithMinDist().getDistance() + "; commonEdges is " +geos[j][i].getCommonEdges() + "; and leafContributionSquared is " + geos[j][i].getLeafContributionSquared());

				}
			}
		}
//		 compute average
		avgCompTime = avgCompTime/ (numTrees * (numTrees - 1)/2);
		System.out.println("In doubleCheck, average dist. computation was " + avgCompTime + " ms for " + numTrees * (numTrees - 1)/2 + " trees.");
	}
	else {
		// XXX: to avoid null pointers in outputting.  Fix this!
		for (int i = 0; i < numTrees; i++) {
			for (int j = i+1; j < numTrees; j++) {
				geos[j][i] = geos[i][j];
			}
		}
	}
	
	return geos;
}

	
/** Open file fileName, reads in the trees, and outputs the distances computed by the polynomial distance algorithm.
 *  Assumes first line of file is number of trees, and then one tree per line.
 *
 */
public static void computeAllInterTreeGeodesicsFromFile(String inFileName, String outFileName, boolean doubleCheck, boolean rooted){

    
    PhyloTree[] trees = readInTreesFromFile(inFileName,rooted);
    int numTrees = trees.length;
    if (verbose >= 1 ) {
    	System.out.println("" + numTrees + " trees read in from " + inFileName);
    }
    Geodesic[][] geos = getAllInterTreeGeodesics(trees, doubleCheck);
    
    // print distances to file
        PrintWriter outputStream = null;
  
        // Outputs the distances in a column, with the first two columns being the trees numbers and the third
    // number the geodesic distance between those trees
        try {
        	outputStream = new PrintWriter(new FileWriter(outFileName));
 
    		for (int i = 0; i < numTrees -1 ; i++) {
    			for (int j = i + 1; j< numTrees; j++) {
/*    				if (verbose >0) {
    					System.out.println("Geo " + i + " -> " + j + " is " + geos[i][j]);
    					System.out.println(geos[i][j].toStringVerbose(trees[i], trees[j]));
    				}*/
				outputStream.println(i + "\t" + j + "\t" + Tools.roundSigDigits(geos[i][j].getDist(), 6));
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
	

/** Open file fileName, reads in the trees, and outputs the distances computed by the polynomial distance algorithm.
 *  Assumes first line of file is number of trees, and then one tree per line.
 *
 */
public static void computeGeodesicsLargeFile(String inFileName, String outFileName, boolean rooted){

    
    PhyloTree[] trees = readInTreesFromFile(inFileName,rooted);
    int numTrees = trees.length;
    if (verbose >= 1 ) {
    	System.out.println("" + numTrees + " trees read in from " + inFileName);
    }
    
    // print distances to file
    PrintWriter outputStream = null;
  
    // Outputs the distances in a column, with the first two columns being the trees numbers and the third
    // number the geodesic distance between those trees
    try {
        outputStream = new PrintWriter(new FileWriter(outFileName));
 
    	for (int i = 0; i < numTrees -1 ; i++) {
    		for (int j = i + 1; j< numTrees; j++) {
    			double geoDist = calcGeoDist(trees[i],trees[j]);
				outputStream.println(i + "\t" + j + "\t" + Tools.roundSigDigits(geoDist, 6));
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

/** Help message (ie. which arguments can be used, etc.)
 * 
 */
public static void displayHelp() {
	System.out.println("Command line syntax:");
	System.out.println("gtp [options] treefile");
	System.out.println("Optional arguments:");	
	System.out.println("\t -d \t double check results, by computing each distance with the target tree as the starting tree and vice versa; default is false");
	System.out.println("\t -h || --help \t displays this message");
	System.out.println("\t -n \t normalize (vector of the lengths of all edges has length 1)");
	System.out.println("\t -o <outfile> \t store the output in the file <outfile>");
	System.out.println("\t -u \t unrooted trees (default is rooted trees)");
	System.out.println("\t -v || --verbose \t verbose output");
}
	
	
/**
 * @param args
 */
public static void main(String[] args) {
	String treeFile = "";
	String outFile = "output.txt"; // default
	boolean doubleCheck = false;
	boolean minLabelling = false;   // used for relabelling one tree to minimize the geodesic distance
	boolean rooted = true;
	boolean largeFile = true;
	
	if (args.length < 1) {
		displayHelp();
		System.exit(0);
	}
	treeFile = args[args.length-1];
	for (int i = 0; i < args.length - 1; i++) {
		
		if (!args[i].startsWith("-")) {
			System.out.println("Invalid command line option");
			displayHelp();
			System.exit(0);
		}
			
		if (args[i].equals("--verbose")) {
			verbose = 1;
		}
		else if (args[i].equals("--help")) {
			displayHelp();
			System.exit(0);
		}

		// output file
		else if (args[i].equals("-o")) {
			if (i < args.length -2) {
				outFile = args[i+1];
				i++;
			}
			else {
				displayHelp();
				System.exit(0);
			}
		}
			
		// all other arguments.  Note we can have -vn
		else {
			for (int j = 1; j<args[i].length(); j++) {
				
				switch(args[i].charAt(j)) {						
				// doublecheck distances
				case 'd':
					doubleCheck = true;
					break;
					
				// display help
				case 'h':
					displayHelp();
					System.exit(0);
					break;
					
				// signify the tree file is large
				case 'l':
					largeFile = true;
					break;
					
				// relabel one tree to minimize the geodesic distance
				case 'm':
					minLabelling = true;
					break;	
					
					
				// normalize trees?
				case 'n':
					normalize = true;
					break;
					
					// unrooted trees?
				case 'u':
					rooted = false;
					break;
					
				// verbose output
				case 'v':
					verbose = 1;
					break;
				
				// output the geodesic information using Newick trees and giving the boundary trees
				case 'w':
					verbose =2;
					break;
						
				default:
					System.out.println("Illegal command line option.\n");
					displayHelp();
					System.exit(0);
					break;
				} // end switch
			} // end for j
		} // end parsing an individual argument
	}  // end for i (looping through arguments)
	
	if (minLabelling) {
		PhyloTree [] trees = readInTreesFromFile(treeFile,rooted);
		getMinLabelling(trees[0],trees[1], outFile);
		System.exit(0);
	}
	
	if (largeFile) {
		computeGeodesicsLargeFile(treeFile, outFile, rooted);
		System.exit(0);
	}
	
	computeAllInterTreeGeodesicsFromFile(treeFile, outFile, doubleCheck, rooted);
	
	System.exit(0);

}

/** Tries to find the labelling of the second tree to minimize the geodesic distance to the first tree.
 *  XXX:  swapLeaves has a bug right now
 * 
 * @param treeFile
 * @param outFile
 */
public static PhyloTree getMinLabelling(PhyloTree tree1, PhyloTree tree2, String outFile) {
	int numIter = 100;
	double potentialGeo;
	PhyloTree potentialTree;
	double controlParam = 1;  // used in computing the probability.  Need to figure out how to update.
	Random r = new Random();
	
	PhyloTree currentTree = tree2.clone();
	double currentGeo = calcGeoDist(tree1, currentTree);
	
	double minGeo = currentGeo;
	PhyloTree minLabelledTree = currentTree.clone();
	
	// iterate the specified number of times to find the best labelling possible in that time.
	for (int i= 0; i < numIter; i++) {
		System.out.println("current geo: " + currentGeo + "; tree: " + currentTree);
		// from the currentTree, generate a potential tree
		potentialTree = currentTree.clone();
		potentialTree.swapleaves( r.nextInt(currentTree.numLeaves()) , r.nextInt(currentTree.numLeaves()) );
		
		
		potentialGeo = PolyMain.getGeodesic(tree1, potentialTree, null).getDist();
		
		// check if we have a shorter geo
		if (potentialGeo <= currentGeo) {
			currentTree = potentialTree.clone();
			currentGeo = potentialGeo;
		}
		else {
			// change the current tree to the potential with probability exp(- (potentialGeo - currentGeo)/controlParm )
			double prob = Math.exp(-(potentialGeo - currentGeo)/controlParam );
			if (Math.random() < prob) {
				currentTree = potentialTree.clone();
				currentGeo = potentialGeo;
			}
		}
		
		// set the minLabelledTree, if relevant
		if (currentGeo < minGeo) {
			minLabelledTree = currentTree.clone();
			minGeo = currentGeo;
		}
		
		// update controlParam if desired.
	}
	
	System.out.println("Min labelled tree with dist " + minGeo + " is: " + minLabelledTree);
	return minLabelledTree;
}


public static double calcGeoDist(PhyloTree t1, PhyloTree t2) {
	return getGeodesic(t1, t2, null).getDist();
}

}

