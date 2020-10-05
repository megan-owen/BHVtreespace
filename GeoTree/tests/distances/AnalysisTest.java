package distances;

import org.junit.jupiter.api.*;

import static org.junit.Assert.*;

import distanceAlg1.*;
import distances.Analysis;

import java.io.ByteArrayOutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.*;


public class AnalysisTest {
	private static PhyloTree t1;
	private static PhyloTree t2;
	private static PhyloTree t3;
	private static PhyloTree multi;
	
	private static PhyloTree t1_different_internal_lengths;
	private static PhyloTree t1_different_leaf_lengths;
	private static PhyloTree t1_multi;
	private static PhyloTree t1_multi_e;
	private static PhyloTree t1_one_off;
	private static PhyloTree t1_one_off2;
	
	private static PhyloTree t1_noCommon;
	private static PhyloTree t3_noCommon;
	
	private static PhyloTree multi_len3;
	private static PhyloTree t1_len3;
	
	private static PhyloTree[] treesNoCommonEdges;
	private static PhyloTree[] treesCommonEdges;
	private static PhyloTree[] treesCommonRemoved;
	private static PhyloTree[] treesEmpty;
	private static PhyloTree[] treesLen3;
	
	private static PhyloTree k1;
	private static PhyloTree k2;



//	@BeforeClass
//	public static void setUpBeforeClass() throws Exception {
	
	/* Set up before each test. (To avoid one test interfering with another.) */
 @BeforeEach
	public void setUp() {
		t1 = new PhyloTree("((C:1,(A:1,B:1):1):2,((D:1,E:0.4):3,F:1):1);", true);  // rooted: AB|CDEF0, ABC|DEF0, DE|ABCF0, EDF|ABC0
		t2 = new PhyloTree("((B:1,((A:1,F:1):1,D:1):2):1,(C:1,E:1):1);", true);  // rooted: AF|BCDE0, ADF|BCE0, ABDF|CE0, CE|ABDF0
		t3 = new PhyloTree("((((C:1,B:1):1,A:1):1,E:1):1,(D:1,F:1):1);", true);  // rooted; split ABC in common with t1
		
		t1_different_internal_lengths = new PhyloTree("((C:1,(A:1,B:1):1):7,((D:1,E:0.4):3,F:1):1);", true);  // rooted
		t1_different_leaf_lengths = new PhyloTree("((C:1,(A:9,B:1):1):2,((D:1,E:0.4):3,F:1):1);", true);  // rooted	
		
		t1_multi = new PhyloTree("((C:1,A:1,B:1):2,((D:10,E:0.4):3,F:1):1);", true);  // rooted, splits:  ABC|DEF0, DE|ABCF0, DEF|ABC0
		  // same splits as t1, except AB|CDEFO
		  // no change in edge lengths, except leaf edge D is 10 instead of 1
		
		t1_multi_e = new PhyloTree("((C:1,A:1,B:1):2,((D:10,E:0.5):3,F:1):1);", true);  // same as t1_multi but E edge is 0.5 (as in PhyloTreeTest) instead of 0.4
		
		t1_one_off = new PhyloTree("(((C:1,A:1):1,B:1):2,((D:1,E:0.4):3,F:1):1);", true);  // rooted, same splits as t1, except has 
		 // split AC|BDEF0 instead of AB|CDEF0

		t1_one_off2 = new PhyloTree("((C:1,(A:1,B:1):1):2,((D:1,F:1):2,E:1):1);",true); // rooted, same splits as t1 except has
		// split DF|ABCD instead of DE|ABCF
		// also edge E has length 1

		
		
		t1_noCommon = new PhyloTree("(C:1,(A:1,B:1):1,((D:1,E:0.4):3,F:1):1);", true);  // rooted;  t1 with split ABC removed
		t3_noCommon = new PhyloTree("(((C:1,B:1):1,A:1,E:1):1,(D:1,F:1):1);", true);  // rooted; t3 with split ABC removed 
		
		multi = new PhyloTree("(A:1,B:1,(C:1,D:1):1,(E:1,F:1):1);",true); //rooted, multifurcating: CD|ABEF, EF|ABCD
		
		multi_len3 = new PhyloTree("(A:[4 2 9],(B:[-2.4 1 0],(C:[1 1 1.4],D:[-1 2 3],(E:[5.6666 1 2],F:[2 2 2]):[22 33 44]):[-1 -1 -1.44]):[100 1 1.2]);",true);
		t1_len3 = new PhyloTree("((C:[1 1 1],(A:[1 1 1],B:[1 1 1]):[1 1 1]):[2 2 2],((D:[1 1 1],E:[0.4 0.4 0.4]):[1 1 1],F:[1 1 1]):[1 1 1]);", true);
		
		k1 = new PhyloTree("((a:[1 1],b:[1 1]):[1 1],(c:[1 1],d:[1 1]):[1 1],e:[1 1]);",false);  // unrooted, splits: ab|cde, cd|abe
		k2 = new PhyloTree("((a:[-1 0],b:[1 1]):[2 -3],(c:[1 1],e:[1 1]):[1 -1],d:[1 1])",false); 		// unrooted, splits: ab|cde, ce|abd  (1 in common with k1)

		
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
	@AfterEach
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
	
	@Test
	public void testGetLogMapWithOrder() {
		StringWriter stringWriter = new StringWriter();
		PrintWriter printWriter = new PrintWriter(stringWriter);
		
		// Test 1: rooted, attrib length =  1, trees in same orthant, only one internal edge changes
		String answer = "0.0 5.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \n";
		answer = answer + "A,B\nA,B,C\nD,E\nD,E,F\n";
		Analysis.getLogMapWithOrder(t1_different_internal_lengths,t1,printWriter);
		assertEquals("Test 1 failed; ", answer, stringWriter.toString());
		
		stringWriter = new StringWriter();
		printWriter = new PrintWriter(stringWriter);
		// Test 2: rooted, attrib length = 1, trees in same orthant, only leaf edges change
		String answer2 = "0.0 0.0 0.0 0.0 8.0 0.0 0.0 0.0 0.0 0.0 \n";
		answer2 = answer2 + "A,B\nA,B,C\nD,E\nD,E,F\n";
		Analysis.getLogMapWithOrder(t1_different_leaf_lengths, t1, printWriter);
		assertEquals("Test 2 failed; ", answer2, stringWriter.toString());
		
		stringWriter = new StringWriter();
		printWriter = new PrintWriter(stringWriter);
		// Test 3: rooted, attrib length = 1, trees in same orthant, but second tree is on boundary
		String answer3 = "-1.0 0.0 0.0 0.0 0.0 0.0 0.0 9.0 0.0 0.0 \n";
		answer3 = answer3 + "A,B\nA,B,C\nD,E\nD,E,F\n";
		Analysis.getLogMapWithOrder(t1_multi, t1,  printWriter);
		assertEquals("Test 3 failed; ", answer3, stringWriter.toString());
		
		
		stringWriter = new StringWriter();
		printWriter = new PrintWriter(stringWriter);
		// Test 4: rooted, attrib length = 1, trees in different orthant, geodesic leaves through codim 1 boundary
		String answer4 = "-2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \n";
		answer4 = answer4 + "A,B\nA,B,C\nD,E\nD,E,F\n";
		Analysis.getLogMapWithOrder(t1_one_off, t1,  printWriter);
		assertEquals("Test 4 failed; ", answer4, stringWriter.toString());

		stringWriter = new StringWriter();
		printWriter = new PrintWriter(stringWriter);
		// Test 5: rooted, attrib length = 1, trees in different orthants, geodesic passes through origin
		String answer5 = "-1.6831300510639735 -3.366260102127947 -5.04939015319192 -1.6831300510639735 0.0 0.0 0.0 0.0 0.6 0.0 \n";
		answer5 = answer5 + "A,B\nA,B,C\nD,E\nD,E,F\n";
		Analysis.getLogMapWithOrder(t2, t1, printWriter);
		assertEquals("Test 5 failed; ", answer5, stringWriter.toString());
		
		stringWriter = new StringWriter();
		printWriter = new PrintWriter(stringWriter);
		// Test 6:  unrooted, attrib length > 1, second tree in adjacent orthant
		//String answer6 = "1.0 -4.0 -2.0 -2.0 -2.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \n";
		String answer6 = "0.9999999999999999 -3.9999999999999996 -1.9999999999999998 ";
		answer6 = answer6 + "-1.9999999999999998 -1.9999999999999998 -0.9999999999999999 ";
		answer6 = answer6 + "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 \n";
		answer6 = answer6 + "a,b\nc,d\n";
		Analysis.getLogMapWithOrder(k2, k1, printWriter);
		assertEquals("Test 6 failed; ", answer6, stringWriter.toString());
		
		stringWriter = new StringWriter();
		printWriter = new PrintWriter(stringWriter);
		// Test 7: rooted, attrib length = 1, trees in different orthants, base tree is non-binary, 
		// geodesic passes through two orthants
		String answer7 = "0.0 -5.0 0.0 0.9999999999999999 0.0 0.0 0.0 -9.000000000000002 ";
		answer7 = answer7 + "0.5000000000000001 0.0 \n";
		answer7 = answer7 + "A,B,C\nD,E\nD,E,F\nA,B\n";
		Analysis.getLogMapWithOrder(t1_one_off2, t1_multi_e, printWriter);
		assertEquals("Test 7 failed; ", answer7, stringWriter.toString());

	}
	
}
