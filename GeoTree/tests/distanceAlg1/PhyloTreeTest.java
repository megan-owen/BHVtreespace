
package distanceAlg1;

import static org.junit.Assert.*;

import java.util.Vector;

import org.junit.Test;
import org.junit.Before;
import org.junit.After;

import distanceAlg1.PhyloTree;
import java.util.*;
import polyAlg.PolyMain;

public class PhyloTreeTest {
	private static PhyloTree t1;
	private static PhyloTree t1_equals;
	private static PhyloTree t1_different_internal_lengths;
	private static PhyloTree t1_different_leaf_lengths;
	private static PhyloTree t2;
	private static PhyloTree t3;
	private static PhyloTree multi;
	private static PhyloTree multi_unrooted;
	private static PhyloTree multi_len3;
	private static PhyloTree t1_multi;
	private static PhyloTree t1_one_off;
	
	private static PhyloTree s1;
	private static PhyloTree s2;
	private static PhyloTree s3;
	
	private static PhyloTree s1_equals;
	private static PhyloTree s1_different_leaves;
	
	private static PhyloTree k1;
	private static PhyloTree k2;

	
	
	
//	@BeforeClass
//	public static void setUpBeforeClass() throws Exception {
	@Before
	public void setUp() {
		t1 = new PhyloTree("((C:1,(A:1,B:1):1):2,((D:1,E:0.5):3,F:1):1);", true);  // rooted, splits:  AB|CDEF0, ABC|DEF0, DE|ABCF0, DEF|ABC0
		t2 = new PhyloTree("((B:1,((A:1,F:1):1,D:1):2):1,(C:1,E:1):1);", true);  // rooted, no edges in common with t1, multi
																				// splits:  AF|BCDE0, ADF|BCE0, ABDF|CE0, CE|ABDF0
																				// geodesic to t1 is the cone path
																				// all leaf edges same as t1, except E is 1 instead of 0.5
		t3 = new PhyloTree("(((A:1,C:1):1,B:1):1,(E:1,(D:1,F:1):1):1);", true);  // rooted, splits: AC|BDEF0, ABC|DEF0, DF|ABCE0, DEF|ABC0
														// splits ABC|DEF0, DEF|ABC0 in common with t1
														// no splits in common with t2
		
		t1_equals = new PhyloTree("((C:1,(B:1,A:1):1):2,(F:1,(D:1,E:0.5):3):1);", true);  // rooted
		t1_different_internal_lengths = new PhyloTree("((C:1,(A:1,B:1):1):7,((D:1,E:0.5):3,F:1):1);", true);  // rooted
		t1_different_leaf_lengths = new PhyloTree("((C:1,(A:9,B:1):1):2,((D:1,E:0.5):3,F:1):1);", true);  // rooted		
		t1_one_off = new PhyloTree("(((C:1,A:1):1,B:1):2,((D:1,E:0.5):3,F:1):1);", true);  // rooted, same splits as t1, except has 
																						 // split AC|BDEF0 instead of AB|CDEF0

		
		
		t1_multi = new PhyloTree("((C:1,A:1,B:1):2,((D:10,E:0.5):3,F:1):1);", true);  // rooted, splits:  ABC|DEF0, DE|ABCF0, DEF|ABC0
																					  // same splits as t1, except AB|CDEFO
																					  // no change in edge lengths, except leaf edge D is 10 instead of 1
		multi = new PhyloTree("(A:1,B:1,(C:1,D:1):1,(E:1,F:1):1);",true); //rooted, multifurcating, splits: CD|ABEF0, EF|ABCD0
		
		multi_unrooted = new PhyloTree("(Alligator:9.45,(Bobcat:0.0034,Camel:2.3,Dromedary:11):0.789,(Elephant:0.00023,Fox:7.2):7);", false);
	
		multi_len3 = new PhyloTree("(aa:[4 2 9],(bb:[-2.4 1 0],(cc:[1 1 1.4],dd:[-1 2 3],(ee:[5.6666 1 2],ff:[2 2 2]):[22 33 44]):[-1 -1 -1.44]):[100 1 1.2]);",true);	
	
		s1 = new PhyloTree("((a:1,b:1):1,c:1,(d:1,e:1):2);",false);   // splits: ab|cde, de|abc
		s2 = new PhyloTree("((a:1,b:1):2,d:1,(c:1,e:1):1.5);",false);   // splits:  ab|cde, ce|abd
		s3 = new PhyloTree("((a:1,d:1):3,b:1,(c:1,e:1):0.5);",false);  // splits: ad|bce, ce|abd
		
		s1_different_leaves = new PhyloTree("((A:1,B:1):1,C:1,(D:1,E:1):2);",false);
		s1_equals = new PhyloTree("((a:1,b:1):1,c:1,(d:1,e:1):2);",false);
		
		k1 = new PhyloTree("((a:[1 1],b:[1 1]):[1 1],(c:[1 1],d:[1 1]):[1 1],e:[1 1]);",false);  // unrooted, splits: ab|cde, cd|abe
		k2 = new PhyloTree("((a:[-1 0],b:[1 1]):[2 -3],(c:[1 1],e:[1 1]):[1 -1],d:[1 1])",false); 		// unrooted, splits: ab|cde, ce|abd  (1 in common with k1)

	}
	
	@After
	public void tearDown() {
		t1 = null;
		t2 = null;
		t3 = null;
		
		t1_equals = null;
		t1_different_internal_lengths = null;
		t1_different_leaf_lengths = null;
		t1_multi = null;
		t1_one_off = null;
		
		multi = null;
		multi_unrooted = null;
		multi_len3 = null;
		
		s1 = null;
		s2 = null;
		s3 = null;
		
		s1_equals = null;
		s1_different_leaves = null;
		
		k1 = null;
		k2 = null;
	}
	
	@Test
	public void testConstructor() { 
		assertEquals("Test 1 (dendropy Newick format) failed;", t1, new PhyloTree("[&R] ((C:1,(A:1,B:1):1):2,((D:1,E:0.5):3,F:1):1):0.0;", true) );
	
		assertEquals("Test 2 (different representations) failed;", t1, new PhyloTree("((C:1,(B:1,A:1):1):2,(F:1,(D:1,E:0.5):3):1);", true) );
	
		assertEquals("Test 3 (multi, unrooted, longer names) failed;", multi_unrooted, new PhyloTree("(Camel:2.3,Bobcat:0.0034,Dromedary:11,(Alligator:9.45,(Fox:7.2,Elephant:0.00023):7):0.789);",false));
	
		assertEquals("Test 4 (multi,rooted, length > 1) failed;", multi_len3, new PhyloTree("((bb:[-2.4 1 0],(dd:[-1 2 3],(ee:[5.6666 1 2],ff:[2 2 2]):[22 33 44],cc:[1 1 1.4]):[-1 -1 -1.44]):[100 1 1.2],aa:[4 2 9]);",true));
	
		// Tree is t1, with AB split having length 0.
		PhyloTree answer = new PhyloTree("((C:1,A:1,B:1):2,((D:1,E:0.5):3,F:1):1);",true);
		assertEquals("Test 5 (0 length edge) failed; ",answer, new PhyloTree("((C:1,(A:1,B:1):0):2,((D:1,E:0.5):3,F:1):1);",true));
	}
	
	@Test
	public void testGetNewick() {
		// Test 1:  tree defined by string, no edge lengths
		assertEquals("Test 1 failed;", "(((D,E),F),((A,B),C));", t1.getNewick(false));
		
		// Test 2: multi, unrooted, tree defined by string
		assertEquals("Test 2 failed;", "(((Bobcat,Camel,Dromedary),Alligator),Elephant,Fox);", multi_unrooted.getNewick(false));

		// Test 3: vector on edge, tree defined by string
		assertEquals("Test 3 failed;", "((((ee,ff),cc,dd),bb),aa);", multi_len3.getNewick(false));
	}
	
	@Test
	public void testErrorInSyntax() {
		// Extra pair of outside brackets
		System.out.println("ErrorInSyntax test 1");
		assertEquals("Test 1 failed; ", true, PhyloTree.errorInSyntax("((a:1,b:1,c:1,d:1):1)"));
		// Missing final bracket
		System.out.println("ErrorInSyntax test 2");
		assertEquals("Test 2 failed; ", true, PhyloTree.errorInSyntax("(a:1,b:1,(c:1,d:1):1"));
		// Extra final bracket
		System.out.println("ErrorInSyntax test 3");
		assertEquals("Test 3 failed; ", true, PhyloTree.errorInSyntax("(a:1,b:1,(c:1,d:1):1):1)"));
		// Missing outside brackets
		System.out.println("ErrorInSyntax test 4");
		assertEquals("Test 4 failed; ", true, PhyloTree.errorInSyntax("(a:1,b:1):1,(c:1,d:1):1"));
		// Extra pair of middle brackets
		System.out.println("ErrorInSyntax test 5");
		assertEquals("Test 5 failed; ", true, PhyloTree.errorInSyntax("(a:1,(b:1,c:1):1,d:1,((e:1,f:1):1):1)"));
		// Extra pair of middle brackets
		System.out.println("ErrorInSyntax test 6");
		assertEquals("Test 6 failed; ", true, PhyloTree.errorInSyntax("(((((L8:1,L7:1):1,(L10:1,L9:1):1):1,L6:1):1,(((L1:1,L2:1):1,L3:1):1,(L4:1,L5:1):1):1):1,((((R3:1,R1:1):1,R2:1):1,(((R7:1,((R10:1,R9:1):1,R8:1):1):1,R6:1):1,(R4:1,R5:1):1):1):1):1)"));
		// Missing : after )
		System.out.println("ErrorInSyntax test 7");
		assertEquals("Test 7 failed; ", true, PhyloTree.errorInSyntax("((a:1,b:1),(c:1,d:1))"));
		// Correct
		assertEquals("Test 8 failed; ", false, PhyloTree.errorInSyntax("((a:1,b:1):1,(c:1,d:1):1);"));
		assertEquals("Test 9 failed; ", false, PhyloTree.errorInSyntax("(d:1,(e:1,(a:1,b:1):1):1,c:1);"));
	}

	@Test
	public void testNumEdges() {
		assertEquals(4,t1.numEdges());
		assertEquals(4,t2.numEdges());
		assertEquals(2,multi.numEdges());
	}
	
	@Test
	public void testEquals() {
		assertEquals("Test 1 (tree equals itself) failed;",t1, t1);
		assertEquals("Test 2 (equals trees) failed;", t1, t1_equals);
		assertFalse("Test 3 (different topologies) failed;", t1.equals(t2));
		assertFalse("Test 4 (different internal split length) failed;", t1.equals(t1_different_internal_lengths));
		assertFalse("Test 5 (different leaf split lengths) failed;", t1.equals(t1_different_leaf_lengths));
	}
	
	@Test 
	public void testGetSplits() {
		Vector<Bipartition> splitsOf_t1 = new Vector<Bipartition>();
	
		splitsOf_t1.add(new Bipartition("110000") );  // split AB
		splitsOf_t1.add(new Bipartition("111000") );  // split ABC
		splitsOf_t1.add(new Bipartition("000110") );  // split DE
		splitsOf_t1.add(new Bipartition("000111") );  // split DEF
		
		assertTrue("Test 1 (get splits from tree t1) failed;", splitsOf_t1.containsAll(t1.getSplits()) && t1.getSplits().containsAll(splitsOf_t1) );
	}

	
	@Test
	public void testCopyConstructor() {
		assertEquals("Test 1 (copy tree created from Newick) failed;", t1, new PhyloTree(t1));
	}
	
	@Test
	public void testGetCommonEdges() {
		// Test 1:  no common edges; binary, rooted trees
		assertEquals("Test 1 (no common edges) failed;", new Vector<PhyloTreeEdge>(), PhyloTree.getCommonEdges(t1, t2));
		
		
		// Test 2:  1 edge in binary tree compatible with all edges in multi tree; rooted
		Vector<PhyloTreeEdge> answer = new Vector<PhyloTreeEdge>();
		answer.add(new PhyloTreeEdge(new Bipartition("110000"), new EdgeAttribute("1"), -1));
		assertEquals("Test 2 failed;", answer, PhyloTree.getCommonEdges(t1, multi));
	
		// Test 3: same tree, rooted, binary
		answer = TreeDistance.myVectorClonePhyloTreeEdge(t1.getEdges());
		answer.get(0).setAttribute(new EdgeAttribute("0"));
		answer.get(1).setAttribute(new EdgeAttribute("0"));
		answer.get(2).setAttribute(new EdgeAttribute("0"));
		answer.get(3).setAttribute(new EdgeAttribute("0"));
		assertEquals("Test 3 (same tree) failed;", answer, PhyloTree.getCommonEdges(t1,t1_equals) );
		
		// Test 4:  same tree, different internal edge lengths
		answer.get(1).setAttribute(new EdgeAttribute("-5"));
		assertEquals("Test 4 failed;", answer, PhyloTree.getCommonEdges(t1,t1_different_internal_lengths));
		
		// Test 5: binary trees with 2 edges in common, rooted
		// common splits:  ABC|DEF0, DEF|ABC0
		answer = new Vector<PhyloTreeEdge>();
		// split ABC|DEF0
		answer.add(new PhyloTreeEdge(new Bipartition("111000"),new EdgeAttribute("1"), -1));
		// split DEF|ABC0
		answer.add(new PhyloTreeEdge(new Bipartition("000111"),new EdgeAttribute("0"),-1));
		Vector<PhyloTreeEdge> commonEdges = PhyloTree.getCommonEdges(t1,t3);
		Boolean test = answer.containsAll(commonEdges) && commonEdges.containsAll(answer);
		assertTrue("Test 5 failed; ", test );
		
		// Test 6: binary and multi tree with no edges in common, rooted
		assertEquals("Test 6 failed; ", new Vector<PhyloTreeEdge>(), PhyloTree.getCommonEdges(multi, t2));
		
		// Test 7: binary, unrooted, 1 common edge
		answer = new Vector<PhyloTreeEdge>();
		answer.add(new PhyloTreeEdge(new Bipartition("11000"),new EdgeAttribute("-1"),-1));
		assertEquals("Test 7 failed; ", answer, PhyloTree.getCommonEdges(s1,s2));
		
		// Test 8: binary, unrooted, 1 common edge
		answer = new Vector<PhyloTreeEdge>();
		answer.add(new PhyloTreeEdge(new Bipartition("11010"),new EdgeAttribute("-1"),-1));
		assertEquals("Test 8 failed; ", answer, PhyloTree.getCommonEdges(s3,s2));
	}
	
	@Test
	public void testGetDistanceFromOrigin() {
		assertEquals("Test 1 (bifurcating, rooted, length 1) failed; ",Math.sqrt(1 + 1 + 1 + 1 + 4 + 1 + 0.25 + 9 + 1 + 1),t1.getDistanceFromOrigin() ,0.00000001);
		assertEquals("Test 2 (multifurcating, rooted, length 1) failed; ", Math.sqrt(8), multi.getDistanceFromOrigin(),0.00000001);
		
		assertEquals("Test 3 (multifurcating, unrooted, length 1) failed; ", Math.sqrt(Math.pow(9.45,2) + Math.pow(0.0034,2) +Math.pow(2.3,2) + 121 + Math.pow(0.789,2) + Math.pow(0.00023,2) + Math.pow(7.2,2) + 49),multi_unrooted.getDistanceFromOrigin(),0.00000001);
	
		assertEquals("Test 4 (multi, rooted, length >1) failed: ", Math.sqrt( Math.pow(-2.4,2) + 15 + Math.pow(5.6666,2) + 5  + 12 + Math.pow(22,2) + Math.pow(33,2) + Math.pow(44,2) + 2 + Math.pow(1.4,2) + 2 + Math.pow(-1.44,2) + 10001 + Math.pow(1.2,2) + 20 +81), multi_len3.getDistanceFromOrigin(),0.00000001);
	}
	
	@Test
	public void testGetCopyLeafEdgeAttribs()  {
		assertArrayEquals("Test 1 (bifurcating, rooted, length 1) failed; ", t1.getLeafEdgeAttribs(), t1.getCopyLeafEdgeAttribs());
		assertArrayEquals("Test 2 (multifurcating, unrooted, length 1) failed; ", multi_unrooted.getLeafEdgeAttribs(), multi_unrooted.getCopyLeafEdgeAttribs());
		assertArrayEquals("Test 2 (multifurcating, rooted, length >1) failed; ", multi_len3.getLeafEdgeAttribs(), multi_len3.getCopyLeafEdgeAttribs());

	
	}
	
	@Test
	public void testProjectToGeo() {
		// test it projects correctly onto endpoints
		PhyloTree t1p = new PhyloTree("((a:1,b:1):1,(c:1,d:1):1);", true);
		PhyloTree t2p = new PhyloTree("(((a:1,c:1):1,b:1):2,d:1);",true);
		Geodesic geo = PolyMain.getGeodesic(t1p, t2p,null);
		double epsilon = 0.0000000000001;
		
		// project tree in same orthant as t1p, but far away
		// answer is t1p
		PhyloTree x1 = new PhyloTree("((a:1,b:1):10,(c:1,d:1):5);",true);
		assertEquals("Test 1 (projection is t1) failed; ", t1p, x1.projectToGeo(geo,epsilon) );
		
		
		// project tree in orthant next to t1p
		// answer is t1p
		PhyloTree x2 = new PhyloTree("(((c:1,d:1):5,b:1):2,a:1)",true);
		assertEquals("Test 2 failed; ", t1p, x2.projectToGeo(geo,epsilon) );

		// project tree in same orthant as t2p, but far away
		// answer is t2p
		PhyloTree x3 = new PhyloTree("(((a:1,c:1):7,b:1):8,d:1)",true);
		assertEquals("Test 3 failed; ", t2p, x3.projectToGeo(geo,epsilon) );

		// project tree in same orthant as t2p, but close to geodesic (so projects non-trivially)
		// answer is (((a,c):7/13,b):17/13,d)
		PhyloTree x4 = new PhyloTree("(((a:1,c:1):1,b:1):1,d:1)", true);
		PhyloTree answer = new PhyloTree("(((a:1,c:1):0.53846153846,b:1):1.30769230769,d:1)",true);
		PhyloTree projection = x4.projectToGeo(geo,epsilon);
		assertEquals("Test 4 failed. Computed projection is " + projection.getNewick(true), answer, projection);
		System.out.println("distance is " + PolyMain.calcGeoDist(x4,x4.projectToGeo(geo,epsilon)));


		// project tree in middle orthant
		// answer is (((a,b):1/13,c):5/13,d)
		PhyloTree x5 = new PhyloTree("(((a:1,b:1):1,c:1):1,d:1)",true);
		answer = new PhyloTree("(((a:1,b:1):0.07692307692,c:1):0.38461538461,d:1);",true);
		projection = x5.projectToGeo(geo,epsilon);
		assertEquals("Test 5 failed. Computed projection is " + projection.getNewick(true), answer, projection );

		// project tree in flap attached to boundary between middle and last orthant (along geodesic)
		// answer is ((a,b,c):0.5,d)
		PhyloTree x6 = new PhyloTree("(((b:1,c:1):1,a:1):1,d:1);",true);
		answer = new PhyloTree("((a:1,b:1,c:1):0.5,d:1)",true);
		projection = x6.projectToGeo(geo,epsilon);
		assertEquals("Test 6 failed. Computed projection is " + projection.getNewick(true), answer, projection );
		
		// the two trees are the same
		geo = PolyMain.getGeodesic(t1p,t1p,null);
		projection = x6.projectToGeo(geo, epsilon);
		assertEquals("Test 7 failed. Computed projection is " + projection.getNewick(true), t1p, projection);

	}

	@Test
	public void testResample() {
		PhyloTree[] trees = new PhyloTree[0];
		// Test if tree array has length 0
		assertArrayEquals("Test 1 (length 0 tree array) failed; ", new PhyloTree[0], PhyloTree.resample(trees, 10));
		
		trees = new PhyloTree[] {t1};
		PhyloTree[] answer = new PhyloTree[] {t1};
		
		// Test if tree array has length 1
		assertArrayEquals("Test 2a (length 1 tree array, n = 1) failed; ", answer, PhyloTree.resample(trees,1));
		
		answer = new PhyloTree[] {t1, t1, t1};
		assertArrayEquals("Test 2b (length 1 tree array, n > 1) failed; ", answer, PhyloTree.resample(trees,3));
		
		// Test if tree array has length > 1
		trees = new PhyloTree[] {t1, t2, multi};
		System.out.println("Check visually: " + Arrays.toString(PhyloTree.resample(trees,6)));
		
	}
	
	@Test
	public void testPermuteLeaves() {
		
		// Test 1:  rooted, 3 leaves
		PhyloTree t_l3 = new PhyloTree("(A:1,(B:1,C:2):1);", true); 	// rooted
		// swap leaves A and B
		Integer[] permutation = new Integer[] {new Integer(1), new Integer(0), new Integer(2)};
		PhyloTree t_l3_permuted = new PhyloTree("(B:1,(A:1,C:2):1);", true); 	// rooted
		
		t_l3.permuteLeaves(permutation);
		assertEquals("Test 1 (rooted, 3 leaves) failed; ", t_l3_permuted, t_l3);
		
		// Test 2:  unrooted, 4 leaves
		PhyloTree t_l4_unrooted = new PhyloTree("(A:1,B:3,(C:1,D:4):2);", false);		// unrooted
		// swap: A -> C, C -> D, D -> A  (leaf i goes to permutation(i))
		permutation = new Integer[] {new Integer(2), new Integer(1), new Integer(3), new Integer(0)};
		PhyloTree t_l4_unrooted_permuted = new PhyloTree("(D:1,B:3,(A:1,C:4):2);", false);		// unrooted
		
		t_l4_unrooted.permuteLeaves(permutation);
		assertEquals("Test 2 (unrooted, 4 leaves) failed; ", t_l4_unrooted_permuted, t_l4_unrooted);

	}
	
	@Test
	public void testGetEdgesNotInCommonWith() {
		// Test 1:  rooted, no edges in common
		assertEquals("Test 1 (rooted, no edges in common) failed;", t1.getEdges(), t1.getEdgesIncompatibleWith(t2));
		
		// Test 2:  rooted, all edges in common
		assertEquals("Test 2 (rooted, all edges in common) failed;", new Vector<PhyloTreeEdge>(),t1.getEdgesIncompatibleWith(t1_different_internal_lengths));

		// Test 3: rooted, 1 edges compatible with multifurcating tree
		Vector<PhyloTreeEdge> answer = TreeDistance.myVectorClonePhyloTreeEdge(t1.getEdges());
		answer.remove(0);
		assertEquals("Test 3 (rooted, edge compatible with multifurcating tree) failed;", answer, t1.getEdgesIncompatibleWith(multi));
	}

	@Test
	public void testHasSameTopology() {
		// Test 1: same tree
		assertTrue("Test 1 failed; ",t1.hasSameTopology(t1));
		
		// Test 2: same tree but defined separately
		assertTrue("Test 2 failed; ", t1.hasSameTopology(t1_equals));
		
		// Test 3:  same tree, different internal edges
		assertTrue("Test 3 failed; ", t1.hasSameTopology(t1_different_internal_lengths));
		
		// Test 4: same tree, different leaf lengths
		assertTrue("Test 4 failed; ", t1.hasSameTopology(t1_different_leaf_lengths));
		
		// Test 5: same tree, but different leaves
		assertTrue("Test 5 failed; ", !s1.hasSameTopology(s1_different_leaves));
		
		// Test 6: same tree, different definition, unrooted
		assertTrue("Test 6 failed; ", s1.hasSameTopology(s1_equals));
		
		// Test 7: different trees
		assertTrue("Test 7 failed; ", !t1.hasSameTopology(t2));

	}
	
	@Test
	public void testIsBinary() {
		// Test 1: binary, rooted
		assertTrue("Test 1 failed; ", t1.isBinary());
		
		// Test 2: multi, rooted
		assertTrue("Test 2 failed; ", !multi.isBinary());
		
		// Test 3: binary, unrooted
		assertTrue("Test 3 failed; ", s1.isBinary());
		
		// Test 4: multi, unrooted
		assertTrue("Test 4 failed; ", !multi_unrooted.isBinary());

		
	}

	@Test
	public void testGetDimAttribs() {
		// Test 1: rooted, length 1 attribs
		assertEquals("Test 1 failed; ", 1, t1.getDimAttribs());
		
		// Test2: unrooted, length 1 attribs
		assertEquals("Test 2 failed; ", 1, s1.getDimAttribs());
		
		// Test 3: rooted, length > 1 attribs
		assertEquals("Test 3 failed; ", 3, multi_len3.getDimAttribs());
	}

	@Test
	public void testGetDirectionTo() {
		// Test 1: rooted, attrib length =  1, trees in same orthant, only one internal edge changes
		double[] answer = {0,5,0,0,0,0,0,0,0,0};
		for (int i = 0; i < answer.length; i++) {
			assertEquals("Test 1-" + i + " failed; ", answer[i], t1.getDirectionTo(t1_different_internal_lengths)[i],0);
		}
		
		// Test 2: rooted, attrib length = 1, trees in same orthant, only leaf edges change
		double[] answer2 = {0,0,0,0,8,0,0,0,0,0};
		for (int i = 0; i < answer2.length; i++) {
			assertEquals("Test 2-" + i + " failed; ", answer2[i], t1.getDirectionTo(t1_different_leaf_lengths)[i],0);
		}
		
		// Test 3: rooted, attrib length = 1, trees in same orthant, but second tree is on boundary
		double[] answer3 = {-1,0,0,0,0,0,0,9,0,0};
		for (int i = 0; i < answer3.length; i++) {
			assertEquals("Test 3-" + i + " failed; ", answer3[i], t1.getDirectionTo(t1_multi)[i],0);
		}
		
		// Test 4: rooted, attrib length = 1, trees in different orthant, geodesic leaves through codim 1 boundary
		double[] answer4 = {-1,0,0,0,0,0,0,0,0,0};
		for (int i = 0; i < answer4.length; i++) {
			assertEquals("Test 4-" + i + " failed; ", answer4[i], t1.getDirectionTo(t1_one_off)[i],0);
		}
		
		// Test 5: rooted, attrib length = 1, trees in different orthants, geodesic passes through origin
		double[] answer5 = {-1,-2,-3,-1,0,0,0,0,0.29706557712,0};
		for (int i = 0; i < answer5.length; i++) {
			assertEquals("Test 5-" + i + " failed; ", answer5[i], t1.getDirectionTo(t2)[i],0.00000000001);
		}
		
		// Test 6:  unrooted, attrib length > 1, second tree in adjacent orthant
		double[] answer6 = {0.5,-2,-1,-1,-1,-0.5,0,0,0,0,0,0,0,0};
		for (int i = 0; i < answer6.length; i++) {
			assertEquals("Test 6-" + i + " failed; ", answer6[i], k1.getDirectionTo(k2)[i],0);
		}
		
	}
	@Test
	public void testGetLogMap() {
		// Test 1: rooted, attrib length =  1, trees in same orthant, only one internal edge changes
		double[] answer = {0,5,0,0,0,0,0,0,0,0};
		for (int i = 0; i < answer.length; i++) {
			assertEquals("Test 1-" + i + " failed; ", answer[i], t1.getLogMap(t1_different_internal_lengths)[i],0.0000000001);
		}
		
		// Test 2: rooted, attrib length = 1, trees in same orthant, only leaf edges change
		double[] answer2 = {0,0,0,0,8,0,0,0,0,0};
		for (int i = 0; i < answer2.length; i++) {
			assertEquals("Test 2-" + i + " failed; ", answer2[i], t1.getLogMap(t1_different_leaf_lengths)[i],0.0000000001);
		}
		
		// Test 3: rooted, attrib length = 1, trees in same orthant, but second tree is on boundary
		double[] answer3 = {-1,0,0,0,0,0,0,9,0,0};
		for (int i = 0; i < answer3.length; i++) {
			assertEquals("Test 3-" + i + " failed; ", answer3[i], t1.getLogMap(t1_multi)[i],0.0000000001);
		}
		
		// Test 4: rooted, attrib length = 1, trees in different orthant, geodesic leaves through codim 1 boundary
		double[] answer4 = {-2,0,0,0,0,0,0,0,0,0};
		for (int i = 0; i < answer4.length; i++) {
			assertEquals("Test 4-" + i + " failed; ", answer4[i], t1.getLogMap(t1_one_off)[i],0.0000000001);
		}
		
		// Test 5: rooted, attrib length = 1, trees in different orthants, geodesic passes through origin
		double[] answer5 = {-1.68313005106,-3.36626010213,-5.04939015319,-1.68313005106,0,0,0,0,0.5,0};
		for (int i = 0; i < answer5.length; i++) {
			assertEquals("Test 5-" + i + " failed; ", answer5[i], t1.getLogMap(t2)[i],0.00000000001);
		}
		
		// Test 6:  unrooted, attrib length > 1, second tree in adjacent orthant
		double[] answer6 = {1,-4,-2,-2,-2,-1,0,0,0,0,0,0,0,0};
		for (int i = 0; i < answer6.length; i++) {
			assertEquals("Test 6-" + i + " failed; ", answer6[i], k1.getLogMap(k2)[i],0.0000000001);
		}
		
	}
}