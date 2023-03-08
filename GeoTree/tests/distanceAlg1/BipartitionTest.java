package distanceAlg1;

import static org.junit.Assert.*;

import java.util.Vector;

import org.junit.jupiter.api.*;

public class BipartitionTest {
	private static PhyloTree t1;
	
	// Set up before each test.
	@BeforeEach
	public void setUp() {
		t1 = new PhyloTree("((C:1,(A:1,B:1):1):2,((D:1,E:0.5):3,F:1):1);", true);  // rooted, splits:  AB|CDEF0, ABC|DEF0, DE|ABCF0, DEF|ABC0
	}
	
	@AfterEach
	public void tearDown() {
		t1 = null;
	}
	
	@Test
	public void testToStringVerbose() {
		Vector<Bipartition> bipartitionsVec = t1.getSplits();
		assertEquals("Test 1 failed", bipartitionsVec.get(0).toStringVerbose(t1.getLeaf2NumMap()), "A,B");
		assertEquals("Test 2 failed", bipartitionsVec.get(1).toStringVerbose(t1.getLeaf2NumMap()), "A,B,C");
		assertEquals("Test 3 failed", bipartitionsVec.get(2).toStringVerbose(t1.getLeaf2NumMap()), "D,E");
		assertEquals("Test 4 failed", bipartitionsVec.get(3).toStringVerbose(t1.getLeaf2NumMap()), "D,E,F");

	}
}
