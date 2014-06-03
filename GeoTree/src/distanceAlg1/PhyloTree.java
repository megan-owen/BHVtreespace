/** This file is part of GeodeMAPS and GTP, programs for computing the geodesic distance between phylogenetic trees,
 *  and sturmMean, a program for computing the Frechet mean between phylogenetic trees.

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


package distanceAlg1;

import java.util.*;

import polyAlg.*;
import static polyAlg.PolyMain.getGeodesic;
import static polyAlg.PolyMain.calcGeoDist;

/* XXX: (general problems)
1) we need a method to convert trees without the same leaf2NumMap to having the same one, if possible.

*/

public class PhyloTree {
	private Vector<PhyloTreeEdge> edges;	 // each element of edges represents an interior split of the tree.
	private Vector<String> leaf2NumMap;   // the position the leaf name occurs in this vector corresponds to its coordinate in an split vector.
							// ie. if leaf A is in position 5, then, a 1 in an split vector at position 5 means that edges contains leaf A
							// standardized to always be in the natural order for strings
	private EdgeAttribute[] leafEdgeAttribs;
	private String newick;  // the newick format of the tree, if available
	private Boolean rooted;  // if the tree is rooted
	
	// constructor
	// assume these trees are already rooted according to convention; 
	// this shouldn't be a problem, because we root any tree read in from a string according to convention
	public PhyloTree(Vector<PhyloTreeEdge> edges, Vector<String> leaf2NumMap, EdgeAttribute[] leafEdgeLengths, Boolean rooted) {
		this.edges = edges;
		this.leaf2NumMap = leaf2NumMap;
		this.leafEdgeAttribs = leafEdgeLengths;
		this.rooted = rooted;
	}
	
	// constructor
	// assume these trees are already rooted according to convention; 
	// this shouldn't be a problem, because we root any tree read in from a string according to convention
	// XXX:  Should we be checking that none of the edges have length 0?
	public PhyloTree(Vector<PhyloTreeEdge> edges, Vector<String> leaf2NumMap, Boolean rooted) {
		this.edges = edges;
		this.leaf2NumMap = leaf2NumMap;
		this.rooted = rooted;
	}
	
	// copy constructor
	public PhyloTree(PhyloTree t) {
		this.edges = Tools.myVectorClonePhyloTreeEdge(t.edges);
		this.leaf2NumMap = Tools.myVectorCloneString(t.leaf2NumMap);
		if (t.leafEdgeAttribs != null) {
			this.leafEdgeAttribs = Arrays.copyOf(t.leafEdgeAttribs, t.leafEdgeAttribs.length);
		}
		if (t.newick != null) {
			this.newick = new String(t.newick);
		}
		this.rooted = t.rooted;
	}

	
	private void setLeaf2NumMapFromNewick() {
		// go through the string and pull out all the leaf labels:  any string between '(' and ':' or ',' and ':'
		int i =0;
		while (i < newick.length()) {
			// however, the first character might be the beginning of a leaf label
			if ((newick.charAt(i) == '(' || newick.charAt(i) == ',')&&(newick.charAt(i+1) != '(')) {
				leaf2NumMap.add(newick.substring(i+1, newick.substring(i).indexOf(":")+i));
				i = nextIndex(newick, i, ",)");
			}
			else {
				i++;
			}
		}
		// sort the elements of leaf2NumMap
		Collections.sort(leaf2NumMap);
	}
	
	
	/** Constructor - turns t into a PhyloTree, by storing the split data in edges and creating a conversion chart between the leaf names
	 * and their position in split vectors.  Basically each name is assigned a different number from 1 to # of leaves.
	 * Example of Newick Standard:  (B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);
	 * If tree is really unrooted, as passed in as argument and stored in TreeDistance class variable, 
	 * then reroot the tree so that the last leaf in leaf2NumMap is the root (and then everything automatically becomes one leaf smaller)
	 * Handles multifurcating trees.
	 * Convention:  Do not add any interior edge with length 0.
	 * If the Newick syntax may not be correct,  call ErrorInSyntax first to get a better error message (indeed, this constructor
	 * may still parse incorrect syntax without returning an error (particularly, the situation where a bracket pair does not have its own comma))
	 * Assumes that the string starts with (
	 * TODO:  insert check for same vertex (so can't have two of the same vertices)
	 * TODO:  Do we want to keep the length of any root vertex?
	 * @param t representation of a tree in Newick standard.
	 */
	public PhyloTree(String t, boolean rooted) {
		int leafNum = 0;
		this.rooted = rooted;

		LinkedList<PhyloTreeEdge> queue = new LinkedList<PhyloTreeEdge>();   // queue for keeping track of edges we are still adding leaves to
		edges = new Vector<PhyloTreeEdge>();
		leaf2NumMap = new Vector<String>();
		int endOfLeafEdgeLength = 0;
				
		newick = t;
			
		//pull off ';' if at end
		if (t.charAt(t.length()-1) == ';') {
			t = t.substring(0, t.length() -1);
		}
		
		// pull off the first and last brackets (and a root length, between the last bracket and ;, if there is one.
		t = t.substring(t.indexOf('(') + 1);
		t = t.substring(0, t.lastIndexOf(')') );
			

try{  // for stringIndexOutOfBoundsException	
		this.setLeaf2NumMapFromNewick();  
		
		leafEdgeAttribs = new EdgeAttribute[leaf2NumMap.size()];
		
		// now go back, and associate leaves with edges.
		int i = 0;
		while (i < t.length() && i > -1) {
			switch(t.charAt(i)) {
			case '(': 
				queue.addFirst(new PhyloTreeEdge());
				i++;
				break;
			
			case ')':
				// Extract this split's attribute, which is between the : following the ), and either a , or a ).
				EdgeAttribute attrib = new EdgeAttribute(t.substring(i+2, nextIndex(t, i, ",)") ) );
				// Only add this edge if norm is > 0
				if (attrib.norm() > 0) {
					queue.peek().setAttribute(attrib ); 
					edges.add(queue.poll());
				}
				// else remove the edge from the list
				else {
					queue.poll();
				}
				// increment i, so it is the next index after the length
				i = nextIndex(t, i, ",)");  // while loop will end if this is -1, and hence we are at the end of the string
				break;
				
			case ',':
				i++;
				break;
			
				// this character is the beginning of a leaf name
			default:
				// get following leaf label, which ends with a :
				leafNum = leaf2NumMap.indexOf(t.substring(i, t.substring(i).indexOf(":")+i));
				
				// the following two lines get the leaf split length for this leaf, and store it
				if ((t.substring(i+1).indexOf(',') > -1) && (t.substring(i+1).indexOf(')') > -1) ) { 
					endOfLeafEdgeLength = Math.min( t.substring(i+1).indexOf(','), t.substring(i+1).indexOf(')') ) + i + 1;		
				}
				else if (t.substring(i+1).indexOf(')') > -1) {
					endOfLeafEdgeLength = t.substring(i+1).indexOf(')') + i+1;
				}
				else if (t.substring(i+1).indexOf(',') > -1) {
					endOfLeafEdgeLength = t.substring(i+1).indexOf(',') + i+1;
				}
				else {
					// we have removed the end bracket, so we are not finding that
					endOfLeafEdgeLength = t.length();
				}

				leafEdgeAttribs[leafNum] = new EdgeAttribute( t.substring( (t.substring(i+1).indexOf(':') + i + 2),endOfLeafEdgeLength  ) );
				
				// we now want to add this leaf to all the edges in our queue.
				for (PhyloTreeEdge e : queue) {
					e.addOne(leafNum);
				}

				i = nextIndex(t, i, ",)");  // while loop will end if this is -1, and hence we are at the end of the string
				break;
				
			}  //switch
				
		}	// while
			
	} 	//try
	catch(StringIndexOutOfBoundsException e) {
		System.err.println("Error reading in tree:  invalid Newick string: (" + t + ");");
		System.exit(1);
	}
	
	// if tree is really unrooted, reroot so that the last leaf in leaf2NumMap is the root
	if (!rooted) {
		for (PhyloTreeEdge e : edges) {
			int lastLeaf = leaf2NumMap.size()-1;
			if (e.contains(lastLeaf)) {
				e.complement(lastLeaf +1);
			}
		}
	}
	
	for (int k = 0; k < edges.size(); k++) {
		edges.get(k).setOriginalEdge(edges.get(k).asSplit());
		edges.get(k).setOriginalID(k);
	}
	
	} // constructor
			

	
		
	// Getters and Setters
	public Vector<PhyloTreeEdge> getEdges() {
		return edges;
	}

	/** Returns the split at position i in the split vector
	 * 
	 * @param i
	 * @return
	 */
	public PhyloTreeEdge getEdge(int i) {
		return edges.get(i);
	}
	
	
	/**  Checks the Newick syntax.  
	 *   Exits with message if problem.
	 *   TODO:  nice error if missing leaf edge
	 *   
	 */
	public static boolean errorInSyntax(String t) {
		boolean error = false;
		int numCommas = 0;
		int numRightParens = 0;  // to be able to say where the problem is
		
		// strip everything after the last ), including ;
		t = t.substring(0,t.lastIndexOf(')')+1);
		
		// count the number of commas
		for (int i = 0; i < t.length(); i++) {
			if (t.charAt(i) == ',') {
				numCommas++;
			}
		}
				
		int[] commaCount = new int[numCommas];
		// fill with 0's.
		Arrays.fill(commaCount, 0);
		
		// Sanity check brackets by adding 1 for each left bracket encountered
		// and subtracting 1 for each right bracket encountered.
		// Bracket "total" should always be >= 0.
				
		// At the same time, check that there is a comma between each pair 
		// of brackets.
		int index = -1;  // add one for each open bracket and subtract one for each closed bracket
		numCommas = 0;  // reset comma count
		for (int i =0; i < t.length(); i++) {
			if (t.charAt(i) == '(') {
				index++;
			}
			else if (t.charAt(i) == ',') {
				commaCount[index]++;
			}
			else if (t.charAt(i) == ')') {
				numRightParens++;
				// check the number of commas for this bracket
				if (commaCount[index] < 1) {
					System.out.println("Error parsing Newick string: missing comma for right parenthesis number " + numRightParens);
					error = true;
					break;
				}
				// reset the comma count for this index, so it can be used again
				commaCount[index] = 0;
				index--;
				// Check that each ) is followed by : or that ) is the last character
				if (i != t.length() - 1) {
					if (!((t.charAt(i+1) == ':') || (t.charAt(i+1) == ';'))) {  
						System.out.println("Error parsing Newick string:  missing ':' after right parenthesis number " + numRightParens);
						error=true;
						break;
					}
				}
				// Check that we still have enough (
				if ((index < 0) && (i != t.length()-1)) {
					System.out.println("Error parsing Newick string:  not enough left parentheses by right parenthesis number " + numRightParens);
					error = true;
					break;
				}
			}
		} // end for
		// index count might be off if there is already an error (since we broke out of the loop)
		// so don't do this test in that case (to avoid misleading error messages)
		if ((!error) && (index != -1)) {
			System.out.println("Error parsing Newick string: uneven number of parentheses");
			error = true;
		}
		return error;
	}
	
	public void setEdges(Vector<PhyloTreeEdge> edges) {
		this.edges = edges;
	}


	public Vector<String> getLeaf2NumMap() {
		return leaf2NumMap;
	}

	public EdgeAttribute getAttribOfSplit(Bipartition edge) {
		Iterator<PhyloTreeEdge> edgesIter = edges.iterator();
		while (edgesIter.hasNext()){
			PhyloTreeEdge e = (PhyloTreeEdge) edgesIter.next();
			if (e.sameBipartition(edge)) {
				return e.getAttribute();
			}
		}
		return null;
	}
	
	/** Returns a vector containing the splits (as Bipartitions) corresponding to the edges of this tree.
	 * 
	 * @return
	 */
	public Vector<Bipartition> getSplits() {
		Vector<Bipartition> splits = new Vector<Bipartition>();
		
		for(int i = 0; i < edges.size(); i++) {
			splits.add(getEdge(i).asSplit() );
		}
		
		return splits;
	}
	
	
	public Boolean isRooted() {
		return rooted;
	}
	
	/** Returns true if this tree is binary.
	 * 
	 * @return
	 */
	public Boolean isBinary() {
		if (isRooted()) {
			if (numEdges() < numLeaves() - 2) {
				return false;
			}
		}
		else {
			if (numEdges() < numLeaves() - 3) {
				return false;
			}
		}
		return true;
	}
	
	
	public Boolean hasSameTopology(PhyloTree t) {
		Vector<Bipartition> splitSet1 = this.getSplits();
		Vector<Bipartition> splitSet2 = t.getSplits();
		
		for (int i = 0; i < this.leaf2NumMap.size(); i++) {
			if (!(this.leaf2NumMap.get(i).equals(t.leaf2NumMap.get(i)))) {
				return false;
			}
		}
		
		return splitSet1.containsAll(splitSet2) && splitSet2.containsAll(splitSet1);
	}
	
	/* * Normalizes so that the vector of all edges (both internal and ones ending in leaves) has length 1.
	 * 
	 */
	public void normalize() {
		double vecLength = getDistanceFromOrigin();
		
		// divide by the length of the split length vector to normalize
		for (int i =0; i < leafEdgeAttribs.length; i++) {
			leafEdgeAttribs[i].scaleBy(1.0/vecLength);
		}
		for (int i = 0; i < edges.size(); i++) {
			edges.get(i).getAttribute().scaleBy(1.0/vecLength);
		}
		// reset getNewick
		newick = null;
	}

	public int numLeaves() {
		return leaf2NumMap.size();
	}
	
	public void setLeaf2NumMap(Vector<String> leaf2NumMap) {
		this.leaf2NumMap = leaf2NumMap;
	}
	

	/**  Returns the distance from this tree to the origin of tree space (i.e. includes leaf edges).
	 * 
	 * @return
	 */
	public double getDistanceFromOrigin() {
		double dist = 0;
		for (int i = 0; i< edges.size(); i++) {
			dist = dist + Math.pow(edges.get(i).getNorm(), 2);
		}
		
		for (int i = 0; i < leafEdgeAttribs.length; i++) {
			dist = dist + Math.pow(leafEdgeAttribs[i].norm(),2);
		}
					
		return Math.sqrt(dist);
	}
	
	/**  Returns the distance from this tree to the origin of tree space \ leaf space.  (i.e. doesn't include leaves)
	 * 
	 * @return
	 */
	public double getDistanceFromOriginNoLeaves() {
		double dist = 0;
		for (int i = 0; i< edges.size(); i++) {
			dist = dist + Math.pow(edges.get(i).getNorm(), 2);
		}
					
		return Math.sqrt(dist);
	}
	
	
	
	/** Returns a vector containing edges with the same partition in t1 and t2.  
	 *  If t1 and t2 do not have the same leaf2NumMap, then returns null.
	 *  If t1 and t2 have no common edges, then returns a vector of size 0. 
	 *  We set the length of the returned edges to be the 
	 *  attribute of the split in t1 minus the attribute of the split in t2.
	 *  Edges from one tree that are compatible with all edges of the other tree are returned as common edges. 
	 *  By convention, no tree contains 0 length common edges.
	 * @param t1
	 * @param t2
	 * @return
	 */
	public static Vector<PhyloTreeEdge> getCommonEdges(PhyloTree t1, PhyloTree t2) {	
		Vector<PhyloTreeEdge> commonEdges = new Vector<PhyloTreeEdge>();
		
		// if the two trees do not have the same leaf2NumMap
		if (!(t1.getLeaf2NumMap().equals(t2.getLeaf2NumMap()))){
			System.out.println("Error: the two trees do not have the same leaves!");
			System.out.println("First tree's leaves are " + t1.getLeaf2NumMap() );
			System.out.println("Second tree's leaves are " + t2.getLeaf2NumMap() );
			System.exit(1);
		}
		
		for (PhyloTreeEdge e1 : t1.edges) {
			if (t2.getSplits().contains(e1.asSplit() ) ){
				// then we have found the same split in each tree
				EdgeAttribute commonAttrib = EdgeAttribute.difference(e1.getAttribute(), t2.getAttribOfSplit(e1.asSplit()));
				commonEdges.add(new PhyloTreeEdge(e1.asSplit(), commonAttrib,e1.getOriginalID() ));
			}
			// otherwise check if the split is compatible with all splits in t2
			else if (e1.isCompatibleWith(t2.getSplits()) ) {
				EdgeAttribute commonAttrib = EdgeAttribute.difference(e1.getAttribute(), null);
				commonEdges.add(new PhyloTreeEdge(e1.asSplit(),commonAttrib,e1.getOriginalID() ));
			}
		}
		// check for splits in t2 that are compatible with all splits in t1	
		for (PhyloTreeEdge e2 : t2.getEdges()) {
			if (e2.isCompatibleWith(t1.getSplits()) && !(t1.getSplits().contains(e2.asSplit()))) {
				EdgeAttribute commonAttrib = EdgeAttribute.difference(null,e2.getAttribute());
				commonEdges.add(new PhyloTreeEdge(e2.asSplit(),commonAttrib,e2.getOriginalID() ));
			}
		}		
		return commonEdges;
	}
	
	/**  Returns the edges that are not in common with edges in t, excluding edges of length 0.
	 * 
	 * @param t
	 * @return
	 */
	public Vector<PhyloTreeEdge> getEdgesIncompatibleWith(PhyloTree t) {
		Vector<PhyloTreeEdge> incompEdges = new Vector<PhyloTreeEdge>();
		
		// if the two trees do not have the same leaf2NumMap
		if (!(this.getLeaf2NumMap().equals(t.getLeaf2NumMap()))){
			System.out.println("Error: the two trees do not have the same leaves!");
			System.out.println("First tree's leaves are " + this.getLeaf2NumMap() );
			System.out.println("Second tree's leaves are " + t.getLeaf2NumMap() );
			System.exit(1);
		}		
		
		for (PhyloTreeEdge e : edges) {
			if (! (e.isCompatibleWith(t.getSplits()) )) { 
				incompEdges.add(e.clone());
			}
		}
		return incompEdges;
	}
	
	
	
	public void permuteLeaves() {
		int numLeaves = leaf2NumMap.size();
		ArrayList<Integer> permutation = new ArrayList<Integer>();
		
		for (int i = 0; i < numLeaves; i++) {
			permutation.add(i);
		}
		
		Collections.shuffle(permutation);
		
		this.permuteLeaves(permutation.toArray(new Integer[numLeaves]));
	}
	
	/** Permutes the leaves of this tree, as specified by the permutation array.
	 * Specicially leaf i is changed to be in position permutation(i).
	 * @param permutation
	 */
	public void permuteLeaves(Integer[] permutation) {
		if (permutation.length != leaf2NumMap.size()) {
			System.err.println("Error: size of permutation map does not match the number of leaves");
			System.exit(1);
		}
		int numEdges = this.numEdges();
		
		Vector<PhyloTreeEdge> permutedV = new Vector<PhyloTreeEdge>();
		// set all the original edges and lengths
		for (PhyloTreeEdge e : this.edges) {
			permutedV.add(new PhyloTreeEdge(e.getAttribute(), e.getOriginalID() ));
		}
		
		for (int j = 0; j < numEdges; j++ ) {
			for (int i = 0; i < permutation.length; i++) {
				//  changing leaf i to be in position permutation(i)
				if(this.getEdge(j).getPartition().get(permutation[i])) {
					permutedV.get(j).getPartition().set(i);
				}
			}
		}
		
		// Permute leaf edge attributes.
		EdgeAttribute[] permutedAttribs = new EdgeAttribute[this.leaf2NumMap.size()];
		for (int i = 0; i < permutation.length; i++ ) {
			permutedAttribs[i] = this.leafEdgeAttribs[permutation[i]]; 
		}
		
		// if tree is really unrooted, reroot so that the last leaf in leaf2NumMap is the root
		if (!rooted) {
			for (PhyloTreeEdge e : permutedV) {
				int lastLeaf = this.leaf2NumMap.size()-1;
				if (e.contains(lastLeaf)) {
					e.complement(lastLeaf +1);
				}
			}
		}
		
//		System.out.println("permutedV is " + permutedV);
		edges = TreeDistance.myVectorClonePhyloTreeEdge(permutedV);
		
		leafEdgeAttribs = permutedAttribs;
	}
	
	// XXX:  need to swap the edge lengths too
	public void swapleaves(int l1,int l2) {
		Vector<PhyloTreeEdge> v = this.edges;
		int numEdges = v.size();
		
		int l1swap;   // saves the value of position l1 when swapping leaves
		
		for (int i = 0; i < numEdges; i++) {
		//  changing leaf l1 to be in position l2
			if(v.get(i).getPartition().get(l1)) {
				l1swap = 1;
			}
			else {
				l1swap = 0;
			}
			if (v.get(i).getPartition().get(l2) ) {
				v.get(i).getPartition().set(l1);
			}
			else {
				v.get(i).getPartition().clear(l1);
			}
			if (l1swap == 1) {
				v.get(i).getPartition().set(l2);
			}
			else {
				v.get(i).getPartition().clear(l2);
			}
		}
	}
	
	
	/** Deletes the split corresponding to bipartition e from the tree.
	 * 
	 * @param e
	 */
	public boolean removeSplit(Bipartition e){
		boolean removed = false;
		
		for(int i = 0; i < edges.size(); i++) {
			if (edges.get(i).sameBipartition(e)) {
				edges.remove(i);
				removed = true;
				break;
			}
		}
		return removed;
	}
	
	/** Remove the edges with bipartitions corresponding to those in Vector<Bipartitions>
	 * 
	 */
	public void removeSplits(Vector<Bipartition> splits) {
		for (int i = 0; i < splits.size(); i++) {
			this.removeSplit(splits.get(i));
		}
	}

	/**  Removes the edges corresponding to the ones in v.
	 * 
	 * @param v
	 */
	public void removeEdgesIndicatedByOnes(Bipartition v) {
		int numRemovedSoFar = 0;
		for (int i = 0; i < leaf2NumMap.size() -2; i++) {
			if (v.contains(i)) {
				edges.remove(i - numRemovedSoFar);
				numRemovedSoFar++;
			}
		}
	}
	
	/** Adds the split e to the tree
	 * 
	 * @param e
	 */
	public void addEdge(PhyloTreeEdge e) {
		edges.add(e);
	}
	
	/** Returns the smallest index in a string t after position i, at which one of the characters in s is located.
	 * 
	 * @param i
	 * @param s
	 * @return t.length() if no character in s is located after position i in t; otherwise return the min index
	 */
	private static int nextIndex(String t, int i, String s) {
		int minIndex = t.length();
		int tempIndex = -1;
		
		for (int j = 0; j < s.length(); j++) {
			tempIndex = t.substring(i+1).indexOf(s.charAt(j));
			if ((tempIndex != -1)  && (tempIndex + i + 1 < minIndex)) {
				minIndex = tempIndex + i +1;
			}
		}
		return minIndex;
	}
	
	
	/** Returns a vector containing representations of 0-1 vectors which
	 * have a 1 in coordinate i if split i of tree t crosses the split of this tree, which this 0-1 
	 * vector represents.  The entries of this vector correspond to the edges of this tree.
	 * Doesn't assume the trees have the same number of edges.
	 * XXX maybe should move exiting if trees don't have same leaf2numMap to this method??
	 * @param t
	 * @return
	 */
	public Vector<Bipartition> getCrossingsWith(PhyloTree t) {
		// if the two trees do not have the same leaf2NumMap
		if (!(this.getLeaf2NumMap().equals(t.getLeaf2NumMap()))){
			return null;
		}
		
/*		Vector<Bipartition> edges = new Vector<Bipartition>();
		Iterator<PhyloTreeEdge> thisIter = this.edges.iterator();
		int i = 0;  // keeps track of which vector coordinate the split 
		// that we are comparing with will correspond to
		
		while (thisIter.hasNext()) {
			PhyloTreeEdge thisTreeEdge = thisIter.next();
			Bipartition newEdge = new Bipartition();
			
			Iterator<PhyloTreeEdge> tIter = t.edges.iterator();
			i = 0;  // reset i
			while (tIter.hasNext()) {
				PhyloTreeEdge e = (PhyloTreeEdge) tIter.next();
				if (e.getEdge() != 0 && e.crosses(thisTreeEdge)) {
					newEdge.addOne(i);
				}
				i++;
			}
			edges.add(newEdge);
		}
		*/
//		System.out.println("in getCrossings, this tree is " + this);
		Vector<Bipartition> edges = new Vector<Bipartition>();
		// for each vector in m 
//		System.out.println("this.edges.size() is " + this.edges.size());
		for (int i = 0; i < this.edges.size(); i++) {
			
			PhyloTreeEdge thisTreeEdge = this.edges.get(i);;
			Bipartition newSplit = new Bipartition();
			
			// add the appropriate 1's
			for (int j = 0; j < t.edges.size(); j++) {
				PhyloTreeEdge e = t.edges.get(j);
				// e should not equal 0.
				if (e.crosses(thisTreeEdge)) {
					newSplit.addOne(j);
				}
			}
			edges.add((Bipartition) newSplit.clone());
//			System.out.println("Added split " + newEdge.clone());
		}
//		System.out.println("in getCrossingsWith returning " + edges);
		return edges;
	}
	
	public String toString() {
		return "Leaves: " + leaf2NumMap + "; edges: " + edges + "; leaf edges: " + Arrays.toString(leafEdgeAttribs); 
	}


	// TODO:  not actually overriding the object clone method.  Also, clone should not call constructors.
	public PhyloTree clone() {
		return new PhyloTree(TreeDistance.myVectorClonePhyloTreeEdge(edges), TreeDistance.myVectorCloneString(leaf2NumMap), leafEdgeAttribs.clone(), rooted);
	}
	
	/* Returns true of the two trees are equal, 
	 * namely if they have the same leaf-label set, the same
	 * leaf edge attributes, and the same set of internal edges (including attributes),
	 * where order doesn't matter.
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
  	@Override public boolean equals(Object t) {
		if (t == null) {
			return false;
		}
		if (this == t) {
			return true;
		}
		
		if (!(t instanceof PhyloTree)) {
			return false;
		}
		
		return leaf2NumMap.equals( ((PhyloTree) t).leaf2NumMap) && 
			Arrays.equals(leafEdgeAttribs, ((PhyloTree)t).leafEdgeAttribs ) && 
				edges.containsAll( ((PhyloTree) t).edges) &&
					((PhyloTree) t).edges.containsAll(edges);
	}

	public EdgeAttribute[] getLeafEdgeAttribs() {
		return leafEdgeAttribs;
	}
	
	public EdgeAttribute[] getCopyLeafEdgeAttribs() {
		if (leafEdgeAttribs == null) {
			return null;
		}
		
		EdgeAttribute[] copy = new EdgeAttribute[leafEdgeAttribs.length];
		
		for (int i = 0; i < leafEdgeAttribs.length; i++) {
			copy[i] = leafEdgeAttribs[i].clone();
		}
		
		return copy; 
	}
	
	/** Returns vector with the weights of the edges in the Vector of interior edges.
	 * 
	 * @return
	 */
	public double[] getIntEdgeAttribNorms() {
		double[] norms = new double[edges.size()];
			
		for (int i = 0; i < edges.size(); i++) {
			norms[i] = edges.get(i).getAttribute().norm();
		}
		return norms;
	}

	public void setLeafEdgeAttribs(EdgeAttribute[] leafEdgeAttribs) {
		this.leafEdgeAttribs = leafEdgeAttribs;
	}

	/** Returns the tree in Newick format, either with or without branch lengths.
	 *  If the tree was constructed using branch lengths, then 
	 * @param branchLengths
	 * @return
	 */
	public String getNewick(Boolean branchLengths) {
		String newNewick = "";
		Vector<String> strPieces = new Vector<String>();	// stores the pieces of the Newick format constructed so far
		Vector<PhyloTreeEdge> corrEdges = new Vector<PhyloTreeEdge>();	// stores the top split for each piece of Newick constructed so far
		
		if ((newick != null) && branchLengths) {
			return newick;
		}
		
		// otherwise we have to figure out the Newick representation for this tree
		// (since even if we have a Newick representation with branch lengths, but want it
		// without, we want to be sure all 0 length branches are excluded)
		if (edges.size() == 0) {
			// we have the star tree
			newNewick = "(";
			for (int i = 0;i < leaf2NumMap.size(); i++) {
				if (i >0 ) {
					newNewick = newNewick + ",";
				}
				newNewick = newNewick + leaf2NumMap.get(i);
				if (branchLengths) {
					newNewick = newNewick + ":" + leafEdgeAttribs[i];
				}
			}
			newNewick = newNewick  + ");";
			return newNewick;
		}
		
		Vector<PhyloTreeEdge> edgesLeft = TreeDistance.myVectorClonePhyloTreeEdge(edges);
		
		while (edgesLeft.size() > 0) {
			
			// pick the next minimal split to add
			PhyloTreeEdge minEdge = edgesLeft.get(0);
		
			for (int i = 1; i < edgesLeft.size(); i++) {
				if (minEdge.contains(edgesLeft.get(i))) {
					// minEdge was not the min elements
					minEdge = edgesLeft.get(i);
				}
			}
			// remove minEdge from edgesLeft
			edgesLeft.remove(minEdge);
					
			corrEdges.add(0,minEdge.clone());
			strPieces.add(0,"");
			
			// now we have a min element.
			// Start its Newick string.
			String str1 = "(";
			// Find out if it contains one of the min elements we have already processed.
			// Since we are allowing degenerate trees, there could be an arbitrary number of such edges.
			
			int k = 1;
			while (k < corrEdges.size()) {
				if (minEdge.contains(corrEdges.get(k))) {
					// add it to the string
					str1 = str1 + strPieces.get(k) + ",";
					
					// remove each leaf in this split from minEdge
					minEdge.getPartition().andNot(corrEdges.get(k).getPartition());
					
					// remove this split and its corresponding string from the vectors
					strPieces.remove(k);
					corrEdges.remove(k);
				}
				else {
					k++;
				}
//				System.out.println("in here");
			}
				
			// add all the elements still in minEdge (These are leaves that weren't already added as part of 
			// a min split contained by minEdge.)
			if (!minEdge.getPartition().isEmpty()) {
				for (int i = 0; i < minEdge.getPartition().length(); i++) {
					if (minEdge.getPartition().get(i)) {
						str1 = str1 + leaf2NumMap.get(i);
						if (branchLengths) {
//							str1 = str1 + ":" + d6o.format(leafEdgeAttribs[i]);
							str1 = str1 + ":" + leafEdgeAttribs[i];
						}
						str1 = str1 +  ",";
					}
				}
			}
			// remove the last , and add the bracket and minEdge length
			str1 = str1.substring(0, str1.length()-1) + ")";
			if (branchLengths) {
//				str1 = str1 + ":" + d6o.format(minEdge.getLength());
				str1 = str1 + ":" + minEdge.getAttribute();
			}
					
			
			// store str1 
			strPieces.set(0,new String(str1));
		}	
				
		// now we need to combine all edges in corrEdges and all remaining leaves
		BitSet allLeaves = new BitSet();
		allLeaves.set(0,leaf2NumMap.size(),true);
		
		String newickString = "(";
		// add all the string pieces we've accumulated
		for (int i = 0; i < corrEdges.size(); i++) {
			newickString = newickString + strPieces.get(i) + ",";
			allLeaves.andNot(corrEdges.get(i).getPartition());
		}
		// add all remaining leaves
		if (!allLeaves.isEmpty()) {
			for (int i = 0; i < allLeaves.length(); i++) {
				if (allLeaves.get(i)) {
					newickString = newickString + leaf2NumMap.get(i);
					if (branchLengths) {
						newickString = newickString + ":" + leafEdgeAttribs[i];
//						newickString = newickString + ":" + d6o.format(leafEdgeAttribs[i]);
					}
					newickString = newickString + ",";
				}
			}
		}
		
		// remove the last , 
		newickString = newickString.substring(0, newickString.length()-1) + ");";
		
		return newickString;
	}

	public void setNewick(String newick) {
		this.newick = newick;
	}
	
	/** Returns number of edges in this tree, including any 0 length ones.
	 * 
	 * @return
	 */
	public int numEdges() {
		return edges.size();
	}
	
	
	
	/**  Projects this tree onto the geodesic between trees t1 and t2;
	 * i.e. returns the index (between 0 and 1) of the tree on the geodesic
	 *  that is closest (= has min geo distance) to this tree.
	 * Algorithm ends when distance between left and right is < epsilon
	 * @param geo
	 * @return
	 */
	public double projectToGeoIndex(Geodesic geo, double epsilon) {
		PhyloTree t1 = geo.getTreeAt(0,leaf2NumMap,rooted);
		PhyloTree t2 = geo.getTreeAt(1,leaf2NumMap,rooted);
		
		// test that the projection is not one of the ends of the geodesic
		// check left end
		if (calcGeoDist(this,t1) < calcGeoDist(this,geo.getTreeAt(epsilon,leaf2NumMap,rooted))) {
			return 0;
		}
		else if (getGeodesic(this,t2,null).getDist() < calcGeoDist(this,geo.getTreeAt(1-epsilon, leaf2NumMap, rooted))) {
			return 1;
		}
		
		
		final double PHI = (1 + Math.sqrt(5))/2;
		double a = 0;
		double b = 1;
		double c1 = 2 - PHI;
		double c2 = PHI - 1;
		
		double distToc1 = calcGeoDist(this, geo.getTreeAt(c1, leaf2NumMap, rooted));
		double distToc2 = calcGeoDist(this, geo.getTreeAt(c2, leaf2NumMap, rooted ));
		while ( (b -a) > epsilon) {
			if (distToc1 < distToc2) {
				// c2 becomes b
				b = c2;
				// c1 becomes c2
				c2 = c1;
				distToc2 = distToc1;
				// compute new c1
				c1 = (PHI - 1)*a + (2 - PHI)*b;
				distToc1 = calcGeoDist(this, geo.getTreeAt(c1, leaf2NumMap, rooted));
			}
			else if (distToc2 < distToc1) {
				// c1 becomes a
				a = c1;
				// c2 becomes c1
				c1 = c2;
				distToc1 = distToc2;
				// compute new c2
				c2 = (2 - PHI)*a + (PHI-1)*b;
				distToc2 = calcGeoDist(this, geo.getTreeAt(c2, leaf2NumMap, rooted));	
			}
			else {
				// distToc2 = distToc1
				// since the distance function is stictly convex, this implies that the 
				// minimum is in between c1 and c2
				// c1 becomes a
				a = c1;
				// c2 becomes b
				b = c2;
				// compute new c1
				c1 = (PHI - 1)*a + (2 - PHI)*b;
				distToc1 = calcGeoDist(this, geo.getTreeAt(c1, leaf2NumMap, rooted));
				// compute new c2
				c2 = (2 - PHI)*a + (PHI-1)*b;
				distToc2 = calcGeoDist(this, geo.getTreeAt(c2, leaf2NumMap, rooted));
			}
		}
		if (distToc1 < distToc2) {
			return (a + c2)/2;
		}
		else {
			// distToc1 > distToc2
			return (c1 + b)/2;
		}
	}
	
	/**  Projects this tree onto the geodesic between trees t1 and t2;
	 * i.e. returns the tree on the geodesic that is closest (= has min geo distance) to this tree.
	 * Algorithm ends when distance between left and right is < epsilon
	 * @param geo
	 * @return
	 */
	public PhyloTree projectToGeo(Geodesic geo, double epsilon) {
		return geo.getTreeAt(this.projectToGeoIndex(geo,epsilon), leaf2NumMap, rooted);
	}
	
	
	
	/** Resamples with replacement n trees from an input sample of trees.
	 *  
	 * @return
	 */
	public static PhyloTree[] resample(PhyloTree[] trees, int n) {
		if (trees.length == 0) {
			return trees;
		}
		
		PhyloTree[] newSample = new PhyloTree[n];
		Random r = new Random();		
		for(int i = 0; i < n; i++) {
			newSample[i] = trees[r.nextInt(trees.length)].clone() ;	
		}
		return newSample;
	}

	
	/**  Get the direction vector from this tree to tree t.
	 *   This is the direction of the geodesic leaving from this tree.
	 *   TODO:  This is not normalized.
	 *   TODO:  For now, this tree must be binary.
	 * @param t
	 * @return
	 */
	public double[] getDirectionTo(PhyloTree t) {
		int numCoords = numEdges() + numLeaves();
		int dimAttrib = this.getDimAttribs();
		double[] coords = new double[numCoords*dimAttrib];
		
		// Check this tree is binary
		// If not, return error.
		if (!this.isBinary()) {
			System.out.println("Error: cannot compute the direction to a non-binary tree " + getNewick(true));
			System.exit(1);
		}
		
		// if tree t doesn't have the same topology, get the boundary tree
		if (!this.hasSameTopology(t)) {
			// if t is on the boundary of this tree's orthant, then it is actually ok
			if (!this.getSplits().containsAll(t.getSplits())) {
				// otherwise get the boundary tree
				t = Geodesic.getBoundaryTrees(this,t).get(0);
			}
		}
		
		// Coordinates will be in the order:  all internal edges, then all leaf edges
		// (In the same order as in this tree)
		int i = 0;
		// internal edges
		for (PhyloTreeEdge edge: getEdges()) {
			// if edge is in t, add in the difference
			if (t.getSplits().contains(edge.asSplit())) {
				EdgeAttribute attrib = EdgeAttribute.difference(t.getAttribOfSplit(edge),edge.getAttribute());
				for (int j = 0; j < dimAttrib; j++) {
					coords[i] = attrib.get(j);
					i++;
				}
			}
			// otherwise, if edge is not in t (because it's multi-furcating),
			// then add the negative of the this tree's edge attribute
			else {
				for (int j = 0; j < dimAttrib; j++) {
					coords[i] = -edge.getAttribute().get(j);
					i++;
				}
			}
		}
		// leaf edges
		for (int k = 0; k< numLeaves(); k++ ) {
			EdgeAttribute attrib = EdgeAttribute.difference(t.getLeafEdgeAttribs()[k],getLeafEdgeAttribs()[k]);
			for (int j = 0; j < dimAttrib; j++) {
				coords[i] = attrib.get(j);
				i++;
			}
		}
		
		return coords;
	}
	
	/**  Get the dimension of an edge attribute in this trees.
	 * 
	 * @return
	 */
	public int getDimAttribs() {
		return getLeafEdgeAttribs()[0].size();
	}

	/** Returns the coordinates of tree t in the log map centred at this tree.
	 *  Note:  If need to speed this up, could check if t has the same topology, 
	 *  and then just return getDirectionTo(t)
	 * 
	 * @param t
	 * @return
	 */
	public double[] getLogMap(PhyloTree t) {
		// get the geodesic distance between the trees
		double geoDist = calcGeoDist(this,t);
			
		// get the direction of the geodesic
		double[] directionCoords = getDirectionTo(t);
			
		return Tools.scaleBy(Tools.normalize(directionCoords),geoDist);	
	}
}
