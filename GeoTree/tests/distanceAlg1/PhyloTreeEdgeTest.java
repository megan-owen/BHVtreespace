package distanceAlg1;

import static org.junit.Assert.*;

//import org.junit.BeforeClass;
import org.junit.jupiter.api.*;


public class PhyloTreeEdgeTest {
	private static PhyloTree t1;
	private static PhyloTree t1_different_length;
	private static PhyloTree t2;
	private static PhyloTree t3;

//	@BeforeClass
//	public static void setUpBeforeClass() throws Exception {

	// Set up before each test.
	@BeforeEach
	public void setUp() {
		t1 = new PhyloTree("((A:1,B:1):1,C:1);", true);	// rooted
		t1_different_length = new PhyloTree("((A:1,B:1):8,C:1);", true);	// rooted
		t2 = new PhyloTree("((A:1,C:1):1,B:1);", true);	// rooted
		t3 = new PhyloTree("(C:1,(A:1,B:1):1);", true);	// rooted
	}
	
	@AfterEach
	public void tearDown() {
		t1 = null;
		t1_different_length = null;
		t2 = null;
		t3 = null;
		
	}

	@Test
	public void testEqualsBipartition() {
		assertEquals("Test 1 (equal Bipartitions) failed", t1.getEdge(0).asSplit(), t3.getEdge(0).asSplit());
		
		
		assertFalse("Test 2 (not equal Bipartitions) failed:", t1.getEdge(0).asSplit().equals(t2.getEdge(0).asSplit()) );
		
		assertFalse("Test 3a (comparing Bipartition with PhyloTreeEdge) failed", t1.getEdge(0).asSplit().equals(t1.getEdge(0)));
		assertFalse("Test 4 (comparing Bipartition with random object) failed", t1.getEdge(0).asSplit().equals(new EdgeAttribute()));
	}
	
	@Test
	public void testEqualsPhyloTreeEdge() {
		assertEquals("Test 0 (PhyloTreeEdge equals itself) failed", t1.getEdge(0), t1.getEdge(0));
		
		assertEquals("Test 1a (equal PhyloTreeEdges) failed", t1.getEdge(0), t3.getEdge(0));
		assertEquals("Test 1b (equal PhyloTreeEdges - opposite order) failed", t3.getEdge(0), t1.getEdge(0));
		
		assertFalse("Test 2 (not equal PhyloTreeEdges) failed:", t1.getEdge(0).equals(t2.getEdge(0)) );
		
		assertFalse("Test 3 (equal splits but different lengths) failed:", t1.getEdge(0).equals(t1_different_length.getEdge(0)) );
		
		assertFalse("Test 4 (comparing PhyloTreeEdge with Bipartition) failed", t1.getEdge(0).equals(t1.getEdge(0).asSplit()));
	}
	

}