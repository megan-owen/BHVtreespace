package distances;

import static org.junit.Assert.*;

import org.junit.Test;
import org.junit.Before;
import org.junit.After;
import distanceAlg1.*;

import java.util.*;


public class AnalysisTest {
	private static PhyloTree t1;
	private static PhyloTree t2;
	private static PhyloTree t3;
	private static PhyloTree multi;
	
	private static PhyloTree t1_noCommon;
	private static PhyloTree t3_noCommon;
	
	private static PhyloTree multi_len3;
	private static PhyloTree t1_len3;
	
	private static PhyloTree[] treesNoCommonEdges;
	private static PhyloTree[] treesCommonEdges;
	private static PhyloTree[] treesCommonRemoved;
	private static PhyloTree[] treesEmpty;
	private static PhyloTree[] treesLen3;



//	@BeforeClass
//	public static void setUpBeforeClass() throws Exception {
	
	/* Set up before each test. (To avoid one test interfering with another.) */
 @Before
	public void setUp() {
		t1 = new PhyloTree("((C:1,(A:1,B:1):1):2,((D:1,E:0.4):3,F:1):1);", true);  // rooted
		t2 = new PhyloTree("((B:1,((A:1,F:1):1,D:1):2):1,(C:1,E:1):1);", true);  // rooted
		t3 = new PhyloTree("((((C:1,B:1):1,A:1):1,E:1):1,(D:1,F:1):1);", true);  // rooted; split ABC in common with t1
		
		t1_noCommon = new PhyloTree("(C:1,(A:1,B:1):1,((D:1,E:0.4):3,F:1):1);", true);  // rooted;  t1 with split ABC removed
		t3_noCommon = new PhyloTree("(((C:1,B:1):1,A:1,E:1):1,(D:1,F:1):1);", true);  // rooted; t3 with split ABC removed 
		
		multi = new PhyloTree("(A:1,B:1,(C:1,D:1):1,(E:1,F:1):1);",true); //rooted, multifurcating
		
		multi_len3 = new PhyloTree("(A:[4 2 9],(B:[-2.4 1 0],(C:[1 1 1.4],D:[-1 2 3],(E:[5.6666 1 2],F:[2 2 2]):[22 33 44]):[-1 -1 -1.44]):[100 1 1.2]);",true);
		t1_len3 = new PhyloTree("((C:[1 1 1],(A:[1 1 1],B:[1 1 1]):[1 1 1]):[2 2 2],((D:[1 1 1],E:[0.4 0.4 0.4]):[1 1 1],F:[1 1 1]):[1 1 1]);", true);
		
		/*  Initialize treesNoCommonEdges - set of trees with no edges in common, including multifurcating tree  */
		treesNoCommonEdges = new PhyloTree[3];
		treesNoCommonEdges[0] = t1;
		treesNoCommonEdges[1] = t2;
		treesNoCommonEdges[2] = multi;
		
		/* Initialize treesCommonEdges - two trees with 1 common split ABC |DEF */
		treesCommonEdges = new PhyloTree[2];
		treesCommonEdges[0] = t1;
		treesCommonEdges[1] = t3;
		
		/* Initialize treesCommonRemoved - the same two trees as in treesCommonEdges, 
		 * but with the common split ABC |DEF removed */
		treesCommonRemoved = new PhyloTree[2];
		treesCommonRemoved[0] = t1_noCommon;
		treesCommonRemoved[1] = t3_noCommon;
		
		treesEmpty = new PhyloTree[15];
		
		treesLen3 = new PhyloTree[2];
		treesLen3[0] = t1_len3;
		treesLen3[1] = multi_len3;
		
		
		
	}
	
	/* Tear down after each test. */
	@After
	public void tearDown() {
		t1 = null;
		t2 = null;
		t3 = null;
		t1_noCommon = null;
		t3_noCommon = null;
		
		multi_len3 = null;
		
		
		multi = null;
		treesNoCommonEdges = null;
		treesCommonEdges = null;
		treesCommonRemoved = null;
		treesLen3 = null;
		
		treesEmpty = null;
	}

	@Test
	public void testGetCommonSplits() {
		assertNull("Test 1 (tree array empty) failed", Analysis.getCommonSplits(treesEmpty));
		
		assertEquals("Test 2 (no common edges) failed", 0, Analysis.getCommonSplits(treesNoCommonEdges).size());
		
		Bipartition split = new Bipartition("111000");		
		
		assertEquals("Test 3a (1 common split) failed", 1, Analysis.getCommonSplits(treesCommonEdges).size());
		assertEquals("Test 3b (1 common split) failed", split, Analysis.getCommonSplits(treesCommonEdges).get(0) );		
	}
	
	@Test
	public void testRemoveCommonSplits() {
		assertArrayEquals("Test 1 (tree array empty) failed", treesEmpty, Analysis.removeCommonSplits(treesEmpty));
		
		assertArrayEquals("Test 2 (no common edges) failed", treesNoCommonEdges, Analysis.removeCommonSplits(treesNoCommonEdges));
		
		assertArrayEquals("Test 3 (1 common split) failed;", treesCommonRemoved, Analysis.removeCommonSplits(treesCommonEdges) );
	}
	
	@Test
	public void testGetStarTree() {
		assertNull("Test 1 (empty trees[]) failed", Analysis.getStarTree(treesEmpty));
		

		assertEquals("Test 2 (no common edges, bi & multi, length 1) failed; ",new PhyloTree("(A:1,B:1,C:1,D:1,E:0.8,F:1);",true),Analysis.getStarTree(treesNoCommonEdges));
	
		assertEquals("Test 3 (common edges, bi, length 1) failed; ",new PhyloTree("(A:1,B:1,C:1,D:1,E:0.7,F:1);",true),Analysis.getStarTree(treesCommonEdges));
		
		assertEquals("Test 4 (bi & multi, length 3) failed; ", new PhyloTree("(A:[2.5 1.5 5],B:[-0.7 1 0.5],C:[1 1 1.2],D:[0 1.5 2],E:[3.0333 0.7 1.2],F:[1.5 1.5 1.5]);",true),Analysis.getStarTree(treesLen3) );
	}

	@Test
	public void testMean() {
		// 1 positive number
		double[] numsTest1 = {4.5};
		assertEquals("Test 1 failed; ", 4.5, Analysis.mean(numsTest1),0.0000000001);
		
		// 1 negative number
		double[] numsTest2 = {-0.0000006};		
		assertEquals("Test 2 failed; ", -0.0000006, Analysis.mean(numsTest2),0.0000000001);

		// multiple numbers both positive and negative
		double[] numsTest3 = {0, 3, -9.1, 57, 0.000982, - 1.33};
		assertEquals("Test 3 failed; ", 8.26183033333, Analysis.mean(numsTest3),0.0000000001);
	}
	
	@Test
	public void testVar() {
		// 1 positive number
		double[] numsTest1 = {4.5};
		assertEquals("Test 1 failed; ", 0, Analysis.var(numsTest1),0.0000000001);
				
		// 1 negative number
		double[] numsTest2 = {-0.0000006};		
		assertEquals("Test 2 failed; ", 0, Analysis.var(numsTest2),0.0000000001);

		// 2 positive numbers
		double[] numsTest3 = {0, 48.3};
		assertEquals("Test 3 failed; ", 583.2225, Analysis.var(numsTest3),0.0000000001);
		
		// multiple numbers both positive and negative
		double[] numsTest4 = {0, 3, -9.1, 57, 0.000982, - 1.33};
		assertEquals("Test 4 failed; ", 488.838643037267, Analysis.var(numsTest4),0.0000000001);
	}
}
