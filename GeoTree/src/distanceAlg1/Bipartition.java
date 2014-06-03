/** This file is part of GeodeMAPS and GTP, programs for computing the geodesic distance between phylogenetic trees,
 *  and sturmMean, a program for computing the Frechet mean between phylogenetic trees.
    Copyright (C) 2008, 2009  Megan Owen, Scott Provan

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

// TODO:  create sameAs method and fix getCommonEdges
import java.util.*;

public class Bipartition implements Cloneable {
	
	protected BitSet partition;

	public Bipartition(){
		partition = new BitSet();
	}
	
	public Bipartition(BitSet edge) {
		this.partition = edge;
	}
	
	public Bipartition(String s) {
		partition = new BitSet();
		
		for (int i = 0; i < s.length(); i++)
			if (s.charAt(i) == '1') {
				partition.set(i);
			}
			else if (s.charAt(i) != '0') {
				System.err.println("Error creating bipartition: input string " + s + " should only contain 0s and 1s");
				System.exit(1);
		}
	}
	
	public BitSet getPartition() {
		return partition;
	}

	public boolean isEmpty() {
/*		if (split.isEmpty()) {
			System.out.println("in isZero, split is " + this.edge + " and will return " + split.isEmpty());
		}*/
		return partition.isEmpty();
	}
	
	public void setPartition(BitSet edge) {
		this.partition = edge;
	}

	/** Adds a 1 to the split, or changes an element from one block to the other in the partition.
	 * Adds the 1 at coordinate one, from the right.
	 * Precondition:  this split contains a 0 at position one
	 * 
	 * @param one
	 */
	public void addOne(int one) {
//		System.out.println("split is " + split);
		partition.set(one);
//		System.out.println("split with flipped bit in " + one + " is " + split);
	}

	/** Removes a 1 from the split.  In other words, transfers the element one to the other (zero)
	 * block in the partition.
	 * Precondition:  this split contains a 1 at position one 
	 * @param leafNum
	 */
	public void removeOne(int one) {
		partition.clear(one);
	}

	/** = e1 and e2 are disjoint edges
	 * Precondition:  the leaf2NumMap for the trees containing e1 and e2 are the same.
	 * @param e1
	 * @param e2
	 * @return
	 */
	public boolean disjointFrom(Bipartition e) {
//		return (this.edge.and(e.edge).compareTo(BigInteger.ZERO) == 0 );
//		System.out.println("" + this.edge + " and " + e.edge + " are disjoint: " + !this.edge.intersects(e.edge));
		return !this.partition.intersects(e.partition);
//		return (this.edge & e.edge) == 0;
	}

	/** = this split contains e
	 * (in particular, returns true if the two splits are equal)
	 * Precondition:  the leaf2NumMap for the tree containing this split and the tree containing e are the same.
	 * @param e
	 * @return
	 */
	public boolean contains(Bipartition e) {
		
		BitSet edgeClone = (BitSet) e.partition.clone();
		edgeClone.and(this.partition);
		return edgeClone.equals(e.partition);

	}

	/** = this split contains a 1 in the bit corresponding to the i-th column
	 * 
	 * @param i
	 * @return
	 */
	public boolean contains(int i) {
//		System.out.println("" + this.edge + " contains split " + i + ": " + this.edge.get(i));
		return this.partition.get(i);
		// return this.edge.testBit(i);
//		return this.contains( new Bipartition((long) Math.pow(2,i)) );
	}
		
	/** = this split properly contains e  (ie. contains e, but is not equal to e)
	 * 
	 * @param e
	 * @return
	 */
	public boolean properlyContains(Bipartition e) {
		return this.contains(e) && !e.contains(this);
	}
	
	/** = this split crosses e
	 * (in particular, will return false if this split is equal to e)
	 * 
	 * @param e
	 * @return
	 */
	public boolean crosses(Bipartition e) {
		return !(disjointFrom(e) || this.contains(e) || e.contains(this));
	}

	/** Returns true if e contains ones exactly where this edges doesn't.
	 * numLeaves gives the total number of ones each should contain.
	 * @param e
	 * @param numLeaves
	 * @return
	 */
/*	public boolean isComplementOf(Bipartition e, int numLeaves) {
//		System.out.println("Sum of edges " + e.edge + " and " + this.edge + " is " + (e.edge + this.edge) + " and want " + ((long)Math.pow(2,numLeaves) - 1));
		if ( (e.edge + this.edge) == ((long)Math.pow(2,numLeaves) - 1) )  {
			return true;
		}
		else {
			return false;
		}
	}*/
	
	/** Changes e to be its complement.
	 *  XXX:  I think this is a slow method; make faster
	 * @param e
	 * @param numLeaves
	 * @return
	 */
	public void complement (int numLeaves) {
//		this.edge = (long)Math.pow(2,numLeaves) - 1 - this.edge;
		for (int i = 0; i < numLeaves; i++) {
			this.partition.flip(i);
		}
	}
	
	/** Converts split into a 0-1 vector.  XXX: leading zeros not shown.
	 *  XXX:  maybe fix???
	 */
	public String toString() {
//		return Long.toBinaryString(split);
//		return "" + this.edge.toByteArray();
		return this.partition.toString();
	}

	public static String toStringVerbose(BitSet edge, Vector<String> leaf2NumMap) {
		String toDisplay = "";
		for (int i = 0; i < edge.length(); i++) {
			if (edge.get(i)) {
				toDisplay = toDisplay + leaf2NumMap.get(i) + ",";
			}
		}
		// remove the last ,
		return toDisplay.substring(0, toDisplay.length()-1);
	}
	
	
	public static String toStringReroot(BitSet edge, Vector<String> leaf2NumMap, String newRoot) {
		// convert newRoot into index
		int rootIndex = leaf2NumMap.indexOf(newRoot);
		
		if (rootIndex == -1) {
			System.err.println("Warning: specified root is not a leaf name, not re-rooting");
			return toStringVerbose(edge,leaf2NumMap);
		}
		
		String toDisplay = "";
		
		if (edge.get(rootIndex) == true) {
			// then we need to reroot the edge, so print the leaf if there is a 0 in that position
			// (vs. printing if there is a 1 normally)
			for (int i = 0; i < edge.length(); i++) {
				if (!(edge.get(i))) {
					toDisplay = toDisplay + leaf2NumMap.get(i) + ",";
				}
			}
			// also, print the root, which is the last leaf in leaf2NumMap
			toDisplay = toDisplay + leaf2NumMap.lastElement() + ",";
		}
		else {
			// display as usual
			for (int i = 0; i < edge.length(); i++) {
				if (edge.get(i)) {
					toDisplay = toDisplay + leaf2NumMap.get(i) + ",";
				}
			}
		}
		// remove the last ,
		return toDisplay.substring(0, toDisplay.length()-1);
	}

	
	/**
	 * @param e
	 * @return
	 */
	
	@Override public boolean equals(Object e) {
		if (e == null) {
			return false;
		}
		if (this == e) {
			return true;
		}
		
		if (!(e instanceof Bipartition) || e instanceof PhyloTreeEdge) {
			return false;
		}
		return partition.equals(((Bipartition) e).partition);
	}
	
	//TODO:  currently this method does NOT override the object one, which has header "public Object clone()"
	// Do we care?  It seems like there are problems with making clones, and one could just use a copy constructor.
	public Bipartition clone() {
		return new Bipartition((BitSet) partition.clone());
	}	
	

	
	/** Returns true if this edge is compatible with the given splits,
	 * that is, if it does not cross any of the given splits.
	 * (in particular if this edge is equal to any of the given splits,
	 * and compatible with the rest then this returns true)
	 * 
	 * @param splits
	 * @return
	 */
	public boolean isCompatibleWith(Vector<Bipartition> splits) {
		boolean compatible =true;
		for (Bipartition s : splits) {
			if (this.crosses(s)) {
				compatible = false;
			}
		}
		return compatible;
	}
	
}
