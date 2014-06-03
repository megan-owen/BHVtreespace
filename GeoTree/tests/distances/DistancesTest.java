package distances;
import distanceAlg1.PhyloTree;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.After;
import org.junit.Test;

public class DistancesTest {
	private static PhyloTree t1;
	private static PhyloTree t2;
	private static PhyloTree t3;
	private static PhyloTree t4;
	private static PhyloTree t5;
	private static PhyloTree multi;
	private static PhyloTree multi_2;
	private static PhyloTree star1;
	private static PhyloTree star2;
	
	
	@Before
	public void setUp() {
		t1 = new PhyloTree("((C:1,(A:1,B:1):1):2,((D:1,E:0.5):3,F:1):1);", true);  // rooted
		t2 = new PhyloTree("((B:1,((A:1,F:1):1,D:1):2):1,(C:1,E:1):1);", true);  // rooted
		t3 = new PhyloTree("((((C:1,B:1):1,A:1):1,E:1):1,(D:1,F:1):1);", true);  // rooted; split ABC in common with t1
		multi = new PhyloTree("(A:1,B:1,(C:1,D:1):1,(E:1,F:1):1);",true); //rooted, multifurcating
		multi_2 = new PhyloTree("((A:1,C:1):1,(B:1,D:1):1,(E:1,F:1):1);",true); //rooted, multifurcating, disjoint with t1
		
		t4 = new PhyloTree("(A:0.083,B:0.033,(C:0.138,(D:0.090,(E:0.086,((F:0.817,G:0.539):0.068,H:1.134):0.187):0.077):0.078):0.055);", false);  //unrooted, 99th tree in Rokas file nt_8_stripped.tre
		t5 = new PhyloTree("(A:0.047,B:0.016,(C:0.051,(D:0.064,(E:0.051,(F:0.244,(G:0.216,H:0.805):0.105):0.252):0.036):0.034):0.021);", false);  //unrooted, 100th tree in Rokas file nt_8_stripped.tre
	
		star1 = new PhyloTree("(A:1,B:1,C:1,D:1,E:2,F:1);", true);
		star2 = new PhyloTree("(C:2,A:1,B:1,D:1,E:1,F:1);",true);
	}
	
	@After
	public void tearDown() {
		t1 = null;
		t2 = null;
		t3 = null;
		t4 = null;
		t5 = null;
		
		multi = null;
		multi_2 = null;
		star1 = null;
		star2 = null;
	}
	
	
	@Test
	public void testGetLeafSpaceL2Dist() {
		/*  assertEquals(String, expected, actual, allowed error) */
		assertEquals("getLeafSpaceL2Dist(t1,t2)",0.5, Distances.getLeafSpaceL2Dist(t1,t2),0.00000001);
		// sqrt( (0.083 - 0.047)^2 + (0.033-0.016)^2 + (0.138 - 0.051)^2 + (0.09 - 0.064)^2 + (0.086 - 0.051)^2 + (0.817 - 0.244)^2 + (0.539 - 0.216)^2 + (1.134 - 0.805)^2 )
		assertEquals("getLeafSpaceL2Dist(t3,t4)",0.74293606723593652, Distances.getLeafSpaceL2Dist(t4, t5),0.00000001);
	}
	
	@Test
	public void testConePath() {
		// sqrt( (sqrt(15) + sqrt(7) )^2 + 0.5^2) 
		assertEquals("conePath(t1,t2)", 6.5378820371676328346, Distances.conePath(t1,t2),0.00000001);
		
		// sqrt(0.74293606723593652^2 + (0.187 - 0.252)^2 + (0.077 - 0.036)^2 + (0.078 - 0.034)^2+  (0.055 - 0.021)^2 + (0.068 + 0.105)^2)
		assertEquals("conePath(t3,t4)",0.768687843015615793, Distances.conePath(t4,t5),0.00000001 );
	}

	@Test
	public void testConePathDisjointNoLeaves() {
		assertEquals("Test 1 (disjoint bifurcating trees) failed",Math.sqrt(15.0) + Math.sqrt(7.0), Distances.conePathDisjointNoLeaves(t1, t2), 0.00000001 );
		assertEquals("Test 2 (disjoint bifurcating/multifurcating trees) failed",Math.sqrt(15.0) + Math.sqrt(3.0), Distances.conePathDisjointNoLeaves(t1, multi_2), 0.00000001 );
		
		assertEquals("Test 3 (bifurcating and star tree) failed", Math.sqrt(15.0), Distances.conePathDisjointNoLeaves(t1,star1),0.00000001);
		assertEquals("Test 3 (2 star trees) failed", 0, Distances.conePathDisjointNoLeaves(star1,star2),0.00000001);
		
	}
	

	@Test
	public void testSi() {
		assertEquals("si(t1,t2)", 1.0, Distances.si(t1, t2),0.00000001 );
	}
	
	@Test
	public void testRf() {
		assertEquals("Test 1 (binary trees, rooted, no common edges) failed;",8, Distances.rf(t1,t2) );
		assertEquals("Test 2 (binary trees, rooted, 1 common edges) failed;",6, Distances.rf(t1,t3) );
		assertEquals("Test 3 (1 binary tree, 1 multifurcating tree, rooted, no common edges) failed;",6, Distances.rf(t1,multi) );
		assertEquals("Test 4 (binary trees, unrooted, 4 common edges) failed;", 2, Distances.rf(t4, t5));
		
		PhyloTree t3_zero = new PhyloTree("((((C:1,B:1):1,A:1):0,E:1):1,(D:1,F:1):1);", true);  // rooted; split ABC in common with t1
		assertEquals("Test 5 (binary trees, rooted, common edge in second trees is 0) failed;", 7, Distances.rf(t1,t3_zero));
		
		PhyloTree t1_zero = new PhyloTree("((C:1,(A:1,B:1):1):0,((D:1,E:0.5):3,F:1):1);", true);  // rooted; split ABC in common with t3
		assertEquals("Test 6 (binary trees, rooted, common edge in first trees is 0) failed;", 7, Distances.rf(t1_zero,t3));


	}
	
	@Test
	public void testWeightedRF() {
		assertEquals("Test 1 (binary trees, rooted, no common edges) failed;",12.5, Distances.weightedRF(t1,t2),0.000000001 );
		assertEquals("Test 2 (binary trees, rooted, 1 common edges) failed;",9.5, Distances.weightedRF(t1,t3),0.0000000001 );
		assertEquals("Test 3 (1 binary tree, 1 multifurcating tree, rooted, no common edges) failed;",9.5, Distances.weightedRF(t1,multi), 0.000000001 );
		assertEquals("Test 4 (binary trees, unrooted, 4 common edges) failed;", 1.783, Distances.weightedRF(t4, t5), 0.000000001);
	}

}
