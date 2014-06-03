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

import polyAlg.PolyMain;

public class Geodesic {
	private RatioSequence rs;
	private Vector<PhyloTreeEdge> commonEdges;   // stores the common edges where their length is the difference in lengths of that edge between the two trees
	private double leafContributionSquared = 0;
	private Vector<PhyloTreeEdge> eCommonEdges;
	
	public static final double TOLERANCE = 0.0000000001;  // compares double values up to TOLERANCE (currently 10 decimal places)

	@Override public boolean equals(Object g) {
		if (g == null) {
			return false;
		}
		if (this == g) {
			return true;
		}
		
		if (!(g instanceof Geodesic)) {
			return false;
		}
		
		// check common edges 
		return (this.commonEdges.containsAll(((Geodesic) g).commonEdges)
				&& ((Geodesic) g).commonEdges.containsAll(this.commonEdges))
				// check e commonEdges
				&& this.eCommonEdges.containsAll(((Geodesic) g).eCommonEdges)
				&& ((Geodesic) g).eCommonEdges.containsAll(this.eCommonEdges)
				// check f commonEdges
				&& this.fCommonEdges.containsAll(((Geodesic) g).fCommonEdges)
				&& ((Geodesic) g).fCommonEdges.containsAll(this.fCommonEdges)
				// check leafAttrib arrays
				&& Arrays.equals(this.eLeafAttribs,((Geodesic) g).eLeafAttribs)
				&& Arrays.equals(this.fLeafAttribs,((Geodesic) g).fLeafAttribs)
				// check leafContributionSquared
				&& Math.abs(this.leafContributionSquared - ((Geodesic) g).leafContributionSquared) < TOLERANCE
				// check ratio sequence
				&& this.rs.equals( ( (Geodesic) g).rs);
	}
	
	public Vector<PhyloTreeEdge> geteCommonEdges() {
		return eCommonEdges;
	}

	public void seteCommonEdges(Vector<PhyloTreeEdge> eCommonEdges) {
		this.eCommonEdges = eCommonEdges;
	}

	public Vector<PhyloTreeEdge> getfCommonEdges() {
		return fCommonEdges;
	}

	public void setfCommonEdges(Vector<PhyloTreeEdge> fCommonEdges) {
		this.fCommonEdges = fCommonEdges;
	}

	public EdgeAttribute[] geteLeafAttribs() {
		return eLeafAttribs;
	}

	public void seteLeafAttribs(EdgeAttribute[] eLeafAttribs) {
		this.eLeafAttribs = eLeafAttribs;
	}

	public EdgeAttribute[] getfLeafAttribs() {
		return fLeafAttribs;
	}

	public void setfLeafAttribs(EdgeAttribute[] fLeafAttribs) {
		this.fLeafAttribs = fLeafAttribs;
	}

	private Vector<PhyloTreeEdge> fCommonEdges;
	private EdgeAttribute[] eLeafAttribs;
	private EdgeAttribute[] fLeafAttribs;
	
	
	
	// constructors
	public Geodesic(RatioSequence rs) {
		this.rs = rs;
		//		this.rs = rs.getAscRSWithMinDist();
		commonEdges = new Vector<PhyloTreeEdge>();
//		dist = this.rs.getDistance();
	}
	
	public Geodesic(RatioSequence rs, EdgeAttribute[] eLeaves, EdgeAttribute[] fLeaves) {
		this.rs = rs;
		this.eLeafAttribs = eLeaves;
		this.fLeafAttribs = fLeaves;
		commonEdges = new Vector<PhyloTreeEdge>();
	}
	
	public Geodesic(RatioSequence rs, Vector<PhyloTreeEdge> cEdges) {
		this.rs = rs;
		commonEdges = cEdges;
	}
	
	public Geodesic(RatioSequence rs, Vector<PhyloTreeEdge> cEdges, double leafContributionSquared) {
		this.rs = rs;
		commonEdges = cEdges;
		this.leafContributionSquared = leafContributionSquared;
	}
	
	
	
	public RatioSequence getRS() {
		return rs;
	}
	
	public void setRS(RatioSequence rs) {
		this.rs = rs;
	}
	
	public double getDist() {
		double commonEdgeDistSquared = 0;
		for(int i = 0; i < commonEdges.size(); i++) {
			commonEdgeDistSquared = commonEdgeDistSquared + Math.pow(commonEdges.get(i).getNorm(),2);
		}
/*		if (TreeDistance.verbose > 0) {
			System.out.println("rs.getAscRSWithMinDist().getDistance() is " + rs.getAscRSWithMinDist().getDistance() + "; commonEdgeDistSquared is " +commonEdgeDistSquared + "; and leafContributionSquared is " + leafContributionSquared);
		} */
		return Math.sqrt(Math.pow(rs.getNonDesRSWithMinDist().getDistance(), 2) + commonEdgeDistSquared + leafContributionSquared);
	}
	
	
	/**  Displays the geodesic in user-friendly form.
	 * 
	 * @return
	 */
	public String toStringVerboseOld(PhyloTree t1, PhyloTree t2) {
		Vector<PhyloTreeEdge> commonEdges = this.getCommonEdges();
		Boolean cEdge = false;
		
		// display T1 only splits
		String toDisplay = "\nSplits only in T1:\n";
		for (int i = 0; i < t1.getEdges().size(); i++) {
			cEdge = false;
			for (int j = 0; j < commonEdges.size(); j++) {
				if (commonEdges.get(j).sameBipartition(t1.getEdge(i))) {
					cEdge = true;
					break;
				}
			}
			if (!cEdge) {
				toDisplay = toDisplay + t1.getEdge(i).toStringVerbose(t1.getLeaf2NumMap()) + "\n" + t1.getEdge(i) + "\n";
			}
		}
		
		//display T2 only splits
		toDisplay = toDisplay + "\n\nSplits only in T2:\n";
		for (int i = 0; i < t2.getEdges().size(); i++) {
			if (! commonEdges.contains(t2.getEdge(i)) ) {
				// if this is not a common split, display
//				toDisplay = toDisplay + t2.getEdge(i).toStringVerbose(t2.getLeaf2NumMap()) + "\n";
				toDisplay = toDisplay + t2.getEdge(i).toStringVerbose(t2.getLeaf2NumMap()) + "\n" + t2.getEdge(i) + "\n";
			}
		}
		
		// display common splits
		toDisplay = toDisplay + "\n\nCommon splits:\n";
		for (int i = 0; i < commonEdges.size(); i++) {
//			toDisplay = toDisplay + commonEdges.get(i).toStringVerbose(t1.getLeaf2NumMap()) + "\n";
			toDisplay = toDisplay + commonEdges.get(i).toStringVerbose(t1.getLeaf2NumMap()) + "\n" + commonEdges.get(i) + "\n";
		}
		
		return toDisplay;
	}
	
	public Geodesic clone() { 
		return new Geodesic(rs.clone(), TreeDistance.myVectorClonePhyloTreeEdge(commonEdges));
	}
	
	public String toString() {
		return "" + getDist() + "; " + rs.getNonDesRSWithMinDist();
//		return "" + getDist() + "; " + rs;
	}

	public Vector<PhyloTreeEdge> getCommonEdges() {
		return commonEdges;
	}

	public void setCommonEdges(Vector<PhyloTreeEdge> commonEdges) {
		this.commonEdges = commonEdges;
	}
	
	public int numCommonEdges() {
		return commonEdges.size();
	}
	
	/** Returns the number of orthants/topologies that the geodesic passes through (not including boundaries between orthants).
	 *  = # ratios in the strictly ascending ratio sequences + 1
	 * @return
	 */
	public int numTopologies() {
		return rs.getAscRSWithMinDist().size() + 1;
	}
	
	/** Returns the geodesic with the ratio sequence (and ratios) reversed.
	 * 
	 * @return
	 */
	public Geodesic reverse() {
		return new Geodesic(rs.reverse(), commonEdges, leafContributionSquared);
	}

	public double getLeafContributionSquared() {
		return leafContributionSquared;
	}

	public void setLeafContributionSquared(double leafContributionSquared) {
		this.leafContributionSquared = leafContributionSquared;
	}
	
	
	/** Returns the tree at the given position (as a number between 0 and 1).
	 * 
	 * @param position
	 * @return
	 */
	public PhyloTree getTreeAt(double position, Vector<String> leaf2NumMap, Boolean isRooted) {
		int lowerRatioIndex = -2;	// the index of the ratio containing all f edges in the tree we want
									// i.e. the index of the ratio with time < position, but such that the next ratio has time >= position 
									// if position is in the starting orthant, we don't want any f edges.
		int higherRatioIndex = -2;  // the index of the ratio containing all e edges in the tree we want
									// i.e. the index of the ratio with time >= position, but such that the next ratio has time >= position
									// if position is in the target orthant, we don't want any e edges
		Vector<PhyloTreeEdge> eEdges = new Vector<PhyloTreeEdge>();
		Vector<PhyloTreeEdge> fEdges =  new Vector<PhyloTreeEdge>();
		
		
		PhyloTree tree = new PhyloTree(this.getCommonEdges(position),leaf2NumMap, isRooted);

		// set the leaf lengths
		EdgeAttribute[] newLeafEdgeAttribs = new EdgeAttribute[eLeafAttribs.length];
		for (int i =0; i < newLeafEdgeAttribs.length; i++) {
			//newLeafEdgeAttribs[i] = ((1-position)*t1.getLeafEdgeAttribs()[i] + position*t2.getLeafEdgeAttribs()[i]);
			newLeafEdgeAttribs[i] = EdgeAttribute.weightedPairAverage(eLeafAttribs[i], fLeafAttribs[i], position);

			//System.out.println("t1 leaf length: " +  t1.getLeafEdgeLengths()[i] + "; t2 leaf length: " +  t2.getLeafEdgeLengths()[i] + "; new length: " + ((1-position)*t1.getLeafEdgeLengths()[i] + position*t2.getLeafEdgeLengths()[i]));

		}
		tree.setLeafEdgeAttribs(newLeafEdgeAttribs);

		if (this.getRS().size() == 0) {
			// then we are done, because the two trees are in the same orthant
			return tree;
		}
		// figure out what orthant the new tree is in
		// first check if the new tree is in the starting orthant
		if (this.getRS().getRatio(0).getTime() > position ) {
			// new tree is in the interior of the starting orthant
			lowerRatioIndex = -1;
			//System.out.println("in starting orthant: setting lower ratio index to -1.  First ratio time is " + geo.getRS().getRatio(0).getTime() + " and position is " + position);

			higherRatioIndex = 0;
			//System.out.println("in starting orthant: setting higher ratio index to 0.  First ratio time is " + geo.getRS().getRatio(0).getTime() + " and position is " + position);

		}
		// if the new tree is in the last orthant
		else if (this.getRS().getRatio(this.getRS().size()-1).getTime() < position) {
			lowerRatioIndex = this.getRS().size()-1;
			higherRatioIndex = this.getRS().size();
			//System.out.println("in target orthant: setting lower ratio index to " + (geo.getRS().size()-1) + ".  Final ratio time is " + geo.getRS().getRatio(geo.getRS().size()-1).getTime() + " and position is " + position);
			//System.out.println("in target orthant: setting higehr ratio index to " + (geo.getRS().size()) + ".  Final ratio time is " + geo.getRS().getRatio(geo.getRS().size()-1).getTime() + " and position is " + position);

		}
		// the new tree is in an intermediate orthant
		else {
			for (int i = 0; i < this.getRS().size(); i++) {
				// note:  want < instead of <= so we are in an orthant and not still on the boundary,
				// if we have a string of equalities
				double ratioTime = this.getRS().getRatio(i).getTime();
				//System.out.println("ratio: " + geo.getRS().getRatio(i));
				//System.out.println("time: " + ratioTime);

				if ((lowerRatioIndex == -2) && (ratioTime >= position)) {
					lowerRatioIndex = i - 1;
					//System.out.println("setting lower ratio index to " + (i-1) + ".  Time is " + ratioTime + " and position is " + position);
				}
				if ((higherRatioIndex == -2) && (lowerRatioIndex != -2) && (ratioTime > position )) {
					higherRatioIndex = i;
					//System.out.println("setting higher ratio index to " + i + ".  Time is " + ratioTime + " and position is " + position);t hot
				}
			}
		}
	// if we didn't set the higherRatioIndex, then we are on the boundary with the target orthant.
		// we want all no e edges, so set higherRatioIndex to 
		if (higherRatioIndex == -2) {
			higherRatioIndex = this.getRS().size();
		}

		// add the edges for all f edges in ratios indexed <= lowerRatioIndex
		for (int i = 0; i <= lowerRatioIndex; i++) {
			fEdges = this.getRS().getRatio(i).getFEdges();

			for (PhyloTreeEdge f : fEdges) {
				
				//double newLength = ((position*geo.getRS().getRatio(i).getFLength() - (1-position)*geo.getRS().getRatio(i).getELength())/geo.getRS().getRatio(i).getFLength())*fEdges.get(j).getLength();
				EdgeAttribute newAttrib = f.getAttribute().clone();
				newAttrib.scaleBy( (position*this.getRS().getRatio(i).getFLength() - (1-position)*this.getRS().getRatio(i).getELength())/this.getRS().getRatio(i).getFLength() );
				// don't have to clone newAttrib because a new object is created each iteration of this loop
				// add the edge only if the length is > 0
				if (newAttrib.norm() > TOLERANCE) {
					tree.addEdge(new PhyloTreeEdge(f.asSplit(), newAttrib,f.getOriginalID() ));
				}
			}
	}

		// to the new tree, add the e edges in the ratios indexed >= higherRatioIndex
		for (int i = higherRatioIndex; i < this.getRS().size(); i++) {
			eEdges = this.getRS().getRatio(i).getEEdges();

			for (PhyloTreeEdge e : eEdges) {
				
				//double newLength = ((1-position)*geo.getRS().getRatio(i).getELength() - position*geo.getRS().getRatio(i).getFLength())/geo.getRS().getRatio(i).getELength()*eEdges.get(j).getLength();
				EdgeAttribute newAttrib = e.getAttribute().clone();
				newAttrib.scaleBy( ((1-position)*this.getRS().getRatio(i).getELength() - position*this.getRS().getRatio(i).getFLength())/this.getRS().getRatio(i).getELength() );
				// add the edge only if it has positive length
				if (newAttrib.norm() > TOLERANCE) {
					tree.addEdge(new PhyloTreeEdge(e.asSplit(), newAttrib,e.getOriginalID() ));
				}
			}
		}
		return tree;
	}
	

	/** Returns all common edges in the geodesic, where  
	 *  we set the length of the returned edges to be (1-position)*e_Attrib_length + position*f_Attrib_length.
	 *  xxx:  assume the common edges are in the same order in all three vectors
	 * @param t1
	 * @param t2
	 * @return
	 */
 	public Vector<PhyloTreeEdge> getCommonEdges(double position) {	
		EdgeAttribute commonEdgeAttribute;
		Bipartition commonSplit;
		
		Vector<PhyloTreeEdge> commonEdgesToReturn = new Vector<PhyloTreeEdge>();
		
		if (position < 0 || position > 1) {
			System.err.println("Error:  position " + position + " must be between 0 and 1");
			System.exit(1);
		}
				
		for (int i = 0; i < commonEdges.size(); i++) {
			commonSplit = commonEdges.get(i).asSplit();
			commonEdgeAttribute = EdgeAttribute.weightedPairAverage(eCommonEdges.get(i).getAttribute(),fCommonEdges.get(i).getAttribute(),position);
			commonEdgesToReturn.add(new PhyloTreeEdge(commonSplit,commonEdgeAttribute.clone(),-1));
		}
		
		
//		// end error checking
//		
//		for (PhyloTreeEdge e : t1.getEdges() ) {
//			if (t2.getSplits().contains(e.asSplit() ) ){
//				// then we have found the same split in each tree
//				commonSplit = e.asSplit();
//				commonEdgeAttribute = EdgeAttribute.weightedPairAverage(e.getAttribute(),t2.getAttribOfSplit(commonSplit),position);
//				commonEdges.add(new PhyloTreeEdge(commonSplit.clone(),commonEdgeAttribute.clone(),-1));
//			}
//			// otherwise check if the split is compatible with all splits in t2
//			else if (e.isCompatibleWith(t2.getSplits()) ) {
//				commonEdgeAttribute = EdgeAttribute.weightedPairAverage(e.getAttribute(),null, position);
//				commonEdges.add(new PhyloTreeEdge(e.asSplit(),commonEdgeAttribute.clone(),-1));
//			}
//		}
//		// check for splits in t2 that are compatible with all splits in t1	
//		for (PhyloTreeEdge e : t2.getEdges()) {
//			if (e.isCompatibleWith(t1.getSplits()) && !(t1.getSplits().contains(e.asSplit()))) {
//				commonEdgeAttribute = EdgeAttribute.weightedPairAverage(null,e.getAttribute(),position);
//				commonEdges.add(new PhyloTreeEdge(e.asSplit(),commonEdgeAttribute.clone(),-1));
//			}
//		}
			
		return commonEdgesToReturn;	
	}

	/**  Returns a vector of the trees in the intersection of the geodesic with all boundary orthants.
	 * 
	 * @return
	 */
	public static Vector<PhyloTree> getBoundaryTrees(PhyloTree t1, PhyloTree t2) {
		Vector<PhyloTree> boundaryTrees = new Vector<PhyloTree>();
		
		Geodesic geo = PolyMain.getGeodesic(t1, t2, null);
		RatioSequence minAscRS = geo.getRS().getAscRSWithMinDist();
		for (Ratio ratio: minAscRS) {
			boundaryTrees.add(geo.getTreeAt(ratio.getTime(), t1.getLeaf2NumMap(), t1.isRooted() ));
		}
		return boundaryTrees;
	}

	/** Returns an ArrayList of the trees that are in the middle of each geodesic 
	 * segment in each othant.  Does not include the starting or ending tree/orthant.
	 *
	 * @return
	 * 
	 */
	public ArrayList<PhyloTree> getMidOrthantTrees(Vector<String> leaf2NumMap, Boolean isRooted) {
		ArrayList<PhyloTree> midTrees = new ArrayList<PhyloTree>();
		RatioSequence minAscRS = this.getRS().getAscRSWithMinDist();
		
		for (int i = 0; i < minAscRS.size() - 1; i++) {
			double timeEnter = minAscRS.get(i).getTime();  // time the geodesic enters the orthant
			double timeLeave = minAscRS.get(i + 1).getTime();	// time the geodesic leaves the orthant
			midTrees.add(this.getTreeAt(timeEnter + (timeLeave - timeEnter)/2 , leaf2NumMap, isRooted ));
		}
		return midTrees;
	}
}
