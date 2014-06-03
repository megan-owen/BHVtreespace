package polyAlg;

import static org.junit.Assert.*;

import java.util.*;

import distanceAlg1.*;
import static polyAlg.PolyMain.calcGeoDist;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class PolyMainTest {
	private static PhyloTree t1;
	private static PhyloTree t1_one_diff_edge;
	private static PhyloTree t2;
	private static PhyloTree t3;
	
	private static PhyloTree v1;
	private static PhyloTree v2;
	private static PhyloTree v3;

	
	private static PhyloTree s1;
	private static PhyloTree s2;
	
	
	@Before
	public void setUp() throws Exception {
		// reset global variables holding subtrees
		PolyMain.aTreesNoCommonEdges = new Vector<PhyloTree>();
		PolyMain.bTreesNoCommonEdges = new Vector<PhyloTree>();
		
		t1 = new PhyloTree("(F:1,(((A:1,B:1):1,C:1):1,(D:1,E:1):1):1);",true);	// rooted; 
								// splits: AB|CDEF0, ABC|DEF0, DE|ABCF0, ABCDE|F0
		t1_one_diff_edge = new PhyloTree("(F:1,(((A:1,B:1):5,C:1):1,(D:1,E:1):1):1);",true);
		
		t2 = new PhyloTree("((((B:1,(A:1,C:1):1):1,D:1):1,F:1):1,E:1);",true); // rooted;
								// splits: AC|BDEF0, ABC|DEF0, ABCD|EF0, ABCDF|E0
								// split ABC|DEF0 in common with t1
		
		t3 = new PhyloTree("(((A:1,B:1):1,C:1):1,((E:1,F:1):1,D:1):1);", true); // rooted;
								// splits:  AB|CDEF0, ABC|DEF0, EF|ABCD0, DEF|ABC0
								// 2 splits AB|CDEF0, ABC|DEF0 in common with t2
		
		
		v1 = new PhyloTree("((a:1,b:1):1,c:1,(d:1,e:1):2);",false); 	// unrooted; splits: ab|cde, de|abc   (1 in common with v2)
		v2 = new PhyloTree("((a:1,b:1):2,d:1,(c:1,e:1):1.5);",false);	// unrooted; splits: ab|cde, ce|abd
//		v4 = new PhyloTree("(((a:1,d:1):0.3971780669,b:1):0.0661963445,c:1,e:1);",false);  //caused problem with v2
//		v4 = new PhyloTree("(((a:1,d:1):0.4,b:1):0.07,c:1,e:1);",false);   // unrooted; splits: ad|bce, abd|ce

		v3 = new PhyloTree("((a:1,d:1):3,b:1,(c:1,e:1):0.5);",false);  // (-3, -0.5)

		
		
		s1 = new PhyloTree("(((b:1,c:1):1,a:1):1,d:1);",true);   //rooted; splits bc|ad0, abc|d0
		s2 = new PhyloTree("((c:1,d:1):1,(a:1,b:1):1);",true);	 // rooted; splits cd|ab0, ab|cd0
		
		
	}

	@After
	public void tearDown() throws Exception {
		t1 = null;
		t1_one_diff_edge = null;
		t2 = null;
		t3 = null;
		
		v1 = null;
		v2 = null;
		v3 = null;
		
		s1 = null;
		s2 = null;
	}

	@Test
	public void testGetGeodesic() {
		assertEquals("Test 1 failed",4,calcGeoDist(t1,t1_one_diff_edge),0.0000001);
		
		
		// Test 2: rooted, binary, no common edges (s1 and s2)
		// Construct the answer
		// There is only 1 ratio in the ratio sequence;
		Vector<PhyloTreeEdge> eEdges = new Vector<PhyloTreeEdge>();
		eEdges.add(new PhyloTreeEdge(new Bipartition("0110"), new EdgeAttribute("1"),-1));	// split bc|ad0
		eEdges.add(new PhyloTreeEdge(new Bipartition("1110"), new EdgeAttribute("1"),-1));  // split abc|d0
		
		Vector<PhyloTreeEdge> fEdges = new Vector<PhyloTreeEdge>();
		fEdges.add(new PhyloTreeEdge(new Bipartition("1100"), new EdgeAttribute("1"),-1));   // split ab|cd0
		fEdges.add(new PhyloTreeEdge(new Bipartition("0011"), new EdgeAttribute("1"),-1));	 // split cd|ab0

		Ratio ratio = new Ratio(eEdges,fEdges);
		RatioSequence answerRS = new RatioSequence();
		answerRS.add(ratio);
		
		Geodesic geo = PolyMain.getGeodesic(s1,s2,null);
		
		// check each element of the geodesic:
		// commonEdges
		assertEquals("Test 2 (commonEdges) failed; ", new Vector<PhyloTreeEdge>(), geo.getCommonEdges());
		
		// eCommonEdges
		assertEquals("Test 2 (eCommonEdges) failed; ", new Vector<PhyloTreeEdge>(), geo.geteCommonEdges());

		// fCommonEdges
		assertEquals("Test 2 (fCommonEdges) failed; ", new Vector<PhyloTreeEdge>(), geo.getfCommonEdges());

		// eLeafAttribs
		EdgeAttribute[] answer_leafAttribs = {new EdgeAttribute("1"),new EdgeAttribute("1"),new EdgeAttribute("1"),new EdgeAttribute("1")};
		assertArrayEquals("Test 2 (eLeafAttribs) failed; ", answer_leafAttribs, geo.geteLeafAttribs());
		
		// fLeafAttribs
		assertArrayEquals("Test 2 (fLeafAttribs) failed; ", answer_leafAttribs, geo.getfLeafAttribs());
		
		// leafContributionSquared
		assertEquals("Test 2 (leafContributionSquared) failed; ", 0, geo.getLeafContributionSquared(),0);

		// rs
		assertEquals("Test 2 (rs) failed; ", answerRS, geo.getRS());
		
		PolyMain.getGeodesic(v2,v1,null);
		
		PolyMain.getGeodesic(v2,v3,null);
	}

	
	@Test
	public void testSplitOnCommonEdge() {
		// reset global variables holding subtrees
		PolyMain.aTreesNoCommonEdges = new Vector<PhyloTree>();
		PolyMain.bTreesNoCommonEdges = new Vector<PhyloTree>();
		
		// Test 1:  1 common edge, rooted
		PolyMain.splitOnCommonEdge(t1, t2);
		PhyloTree tA1 = new PhyloTree("((A:1,B:1):1,C:1,D:1,E:1,F:1);",true);
		tA1.setLeafEdgeAttribs(null);
		PhyloTree tA2 = new PhyloTree("((A:1,C:1):1,B:1,D:1,E:1,F:1);",true);
		tA2.setLeafEdgeAttribs(null);
		PhyloTree tB1 = new PhyloTree("(F:1,(A:1,B:1,C:1,(D:1,E:1):1):1);", true);
		tB1.setLeafEdgeAttribs(null);
		PhyloTree tB2 = new PhyloTree("(((B:1,A:1,C:1,D:1):1,F:1):1,E:1);",true);
		tB2.setLeafEdgeAttribs(null);
		assertEquals("Test 1 (tA1) failed; ", tA1, PolyMain.aTreesNoCommonEdges.get(0));
		assertEquals("Test 1 (tA2) failed; ", tA2, PolyMain.bTreesNoCommonEdges.get(0));
		assertEquals("Test 1 (tB1) failed; ", tB1, PolyMain.aTreesNoCommonEdges.get(1));
		assertEquals("Test 1 (tB2) failed; ", tB2, PolyMain.bTreesNoCommonEdges.get(1));
		
		// reset global variables holding subtrees
		PolyMain.aTreesNoCommonEdges = new Vector<PhyloTree>();
		PolyMain.bTreesNoCommonEdges = new Vector<PhyloTree>();
		
		// Test 2: all edges in common, rooted
		PolyMain.splitOnCommonEdge(t1, t1_one_diff_edge);
		assertEquals("Test 2 (A) failed; ", new Vector<PhyloTree>(), PolyMain.aTreesNoCommonEdges);
		assertEquals("Test 2 (B) failed; ", new Vector<PhyloTree>(), PolyMain.bTreesNoCommonEdges);
		
		// reset global variables holding subtrees
		PolyMain.aTreesNoCommonEdges = new Vector<PhyloTree>();
		PolyMain.bTreesNoCommonEdges = new Vector<PhyloTree>();
		
		// Test 3: 2 common edges, rooted
		PolyMain.splitOnCommonEdge(t1,t3);
		Vector<PhyloTree> answer = new Vector<PhyloTree>();
		answer.add(new PhyloTree("(F:1,(A:1,B:1,C:1,(D:1,E:1):1):1);",true));
		answer.get(0).setLeafEdgeAttribs(null);
		assertEquals("Test 3 (A) failed; ", answer, PolyMain.aTreesNoCommonEdges);
		
		answer = new Vector<PhyloTree>();
		answer.add(new PhyloTree("(A:1,B:1,C:1,((E:1,F:1):1,D:1):1);", true));
		answer.get(0).setLeafEdgeAttribs(null);
		assertEquals("Test 3 (B) failed; ", answer, PolyMain.bTreesNoCommonEdges);
		
		// reset global variables holding subtrees
		PolyMain.aTreesNoCommonEdges = new Vector<PhyloTree>();
		PolyMain.bTreesNoCommonEdges = new Vector<PhyloTree>();
		
		// Test 4: 1 common edge, unrooted
		PolyMain.splitOnCommonEdge(v1,v2);
		answer = new Vector<PhyloTree>();
		answer.add(new PhyloTree("(a:1,b:1,c:1,(d:1,e:1):2);",false));
		answer.get(0).setLeafEdgeAttribs(null);
		assertEquals("Test 4 (A) failed; ", answer, PolyMain.aTreesNoCommonEdges);
		
		answer = new Vector<PhyloTree>();
		answer.add(new PhyloTree("(a:1,b:1,d:1,(c:1,e:1):1.5);", false));
		answer.get(0).setLeafEdgeAttribs(null);
		assertEquals("Test 3 (B) failed; ", answer, PolyMain.bTreesNoCommonEdges);
	}
	
}