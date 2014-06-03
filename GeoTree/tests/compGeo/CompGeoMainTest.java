package compGeo;


import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;


import distanceAlg1.PhyloTree;

import java.util.*;

public class CompGeoMainTest {

	private static PhyloTree t3_1;
	private static PhyloTree t3_2;
	private static PhyloTree t3_3;
	private static PhyloTree[] threeLeafTrees;
	private static Vector<PhyloTree> threeLeafTreesVect;
	
	private static PhyloTree t4_1;
	private static PhyloTree t4_2;
	private static PhyloTree t4_3;
	private static PhyloTree[] fourLeafTrees;
	private static Vector<PhyloTree> fourLeafTreesVect;
	
	
	
	
	@Before
	public void setUp() throws Exception {
		t3_1 = new PhyloTree("((a:1,b:1):1,c:1);",true);  
		t3_2 = new PhyloTree("((a:1,c:1):1,b:1);",true);  
		t3_3 = new PhyloTree("((b:1,c:1):1,a:1);",true);  
		threeLeafTrees = new PhyloTree[] {t3_1, t3_2, t3_3};
		threeLeafTreesVect = new Vector<PhyloTree>();
		threeLeafTreesVect.add(t3_1);
		threeLeafTreesVect.add(t3_2);
		threeLeafTreesVect.add(t3_3);
		
		t4_1 = new PhyloTree("((a:1,b:1):2,(c:1,d:1):1);",true);  // rooted; splits ab|cd0, cd|ab0
		t4_2 = new PhyloTree("(((a:1,c:1):1,b:1):2,d:1);",true);	// rooted; splits ac|bd0, abc|d0
		t4_3 = new PhyloTree("(((a:1,b:1):2,c:1):2,d:1);",true);  // rooted; splits ab|cd0, abc|d0
		
		fourLeafTrees = new PhyloTree[] {t4_1,t4_2,t4_3};
		fourLeafTreesVect = new Vector<PhyloTree>();
		fourLeafTreesVect.add(t4_1);
		fourLeafTreesVect.add(t4_2);
		fourLeafTreesVect.add(t4_3);
	}

	@After
	public void tearDown() throws Exception {
		t3_1 = null;
		t3_2 = null;
		t3_3 = null;
		threeLeafTrees = null;
		threeLeafTreesVect = null;
		
		t4_1 = null;
		t4_2 = null;
		t4_3 = null;
		
		fourLeafTrees = null;
		fourLeafTreesVect = null;
	}
	
	@Test
	public void testRestrictedWeightedComb() {
//		assertEquals("Test 1 failed", new PhyloTree(triangleMean.approxEquals(mean,0.005));

	}
	
	@Test
	public void testGetBoundaryTrees() {
		Vector<PhyloTree> answer = new Vector<PhyloTree>();
		answer.add(new PhyloTree("(a:1,b:1,c:1);",true));
		assertEquals("Test 1 (t3, input is array) failed",answer,CompGeoMain.getBoundaryTrees(threeLeafTrees));
		assertEquals("Test 1b (t3, input is vector) failed",answer,CompGeoMain.getBoundaryTrees(threeLeafTreesVect));

		
		answer.removeAllElements();
		answer.add(new PhyloTree("((a:1,b:1):1,c:1,d:1);",true));
		answer.add(new PhyloTree("((a:1,b:1,c:1):1,d:1);",true));
		answer.add(new PhyloTree("((a:1,b:1):2,c:1,d:1);",true));
		answer.add(new PhyloTree("((a:1,b:1,c:1):2,d:1);",true));
		
		// need to test this way, since we don't care what order they are in the vector
		Vector<PhyloTree> boundaryTrees = CompGeoMain.getBoundaryTrees(fourLeafTrees);
		Boolean test = answer.containsAll(boundaryTrees) && boundaryTrees.containsAll(answer);
		assertTrue("Test 2 (t4, input is array) failed",test);
		
		boundaryTrees = CompGeoMain.getBoundaryTrees(fourLeafTreesVect);
		test = answer.containsAll(boundaryTrees) && boundaryTrees.containsAll(answer);
		assertTrue("Test 2b (t4, input is vector) failed",test);
		

		
	}
	
	@Test
	public void testGetExtremalTrees() {
		Vector<PhyloTree> answer = new Vector<PhyloTree>();
		answer.add(new PhyloTree("((a:1,b:1):1,c:1,d:1);",true));
		answer.add(new PhyloTree("((a:1,b:1,c:1):1,d:1);",true));
		answer.add(new PhyloTree("((a:1,b:1):2,c:1,d:1);",true));
		answer.add(new PhyloTree("((a:1,b:1,c:1):2,d:1);",true));
		answer.add(t4_1);
		answer.add(t4_2);
		answer.add(t4_3);
		
		Vector<PhyloTree> extremalTrees = CompGeoMain.getExtremalTrees(fourLeafTreesVect,5);
		Boolean test = answer.containsAll(extremalTrees) && extremalTrees.containsAll(answer);
		System.out.println("Boundary trees are: " + extremalTrees);
		assertTrue("Test 1 (t4, input is vector) failed",test);
		
		
		
		
		
		
	}

	
	
		
}
