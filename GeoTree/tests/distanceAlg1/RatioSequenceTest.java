package distanceAlg1;

import static org.junit.Assert.*;

import org.junit.jupiter.api.*;

import polyAlg.*;

public class RatioSequenceTest {
	private static RatioSequence rs_s1s2;
	
	private static PhyloTree s1;
	private static PhyloTree s2;
	
	private static PhyloTree v1;
	private static PhyloTree v2;
	
	private static RatioSequence rs_v1v2;
	
	@BeforeEach
	public void setUp() throws Exception {
		s1 = new PhyloTree("(((b:1,c:1):1,a:1):1,d:1);",true);   //rooted; splits bc|ad0, abc|d0
		s2 = new PhyloTree("((c:1,d:1):1,(a:1,b:1):1);",true);	 // rooted; splits cd|ab0, ab|cd0
		
		rs_s1s2 = PolyMain.getGeodesic(s1,s2,null).getRS();  // rs with no common edges
		
		v1 = new PhyloTree("((a:1,b:1):1,c:1,(d:1,e:1):2);",false); 	// unrooted; splits: ab|cde, de|abc   (1 in common with s2)
		v2 = new PhyloTree("((a:1,b:1):2,d:1,(c:1,e:1):1.5);",false);	// unrooted; splits: ab|cde, ce|abd
		
		rs_v1v2 = PolyMain.getGeodesic(v1,v2,null).getRS();  // rs with common edges
		
	}

	@AfterEach
	public void tearDown() throws Exception {
		s1 = null;
		s2 = null;	
		
		rs_s1s2 = null;
		
		v1 = null;
		v2 = null;
		
		rs_v1v2 = null;
	}

	@Test
	public void testEquals() {
		// Test 1:  ratio sequence with itself  
		assertEquals("Test 1 failed; ", rs_s1s2, rs_s1s2);
		
		// Test 2:  ratio sequence between the same two trees (rooted, no common edges) 
		assertEquals("Test 2 failed; ", rs_s1s2, PolyMain.getGeodesic(s1,s2,null).getRS() );
		
		// Test 3: ratio sequence between the same two trees (unrooted, one common edge)
		assertEquals("Test 3 failed; ", rs_v1v2, PolyMain.getGeodesic(v1,v2,null).getRS() );

	}

}
