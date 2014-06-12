package distanceAlg1;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class EdgeAttributeTest {
	private static EdgeAttribute a1;
	private static EdgeAttribute a1_same;
	private static EdgeAttribute a2;
	
	private static EdgeAttribute b1;
	private static EdgeAttribute b1_same;
	private static EdgeAttribute b2;
	
	private static EdgeAttribute c1;

	
	// Set up before each test.
	@Before
	public void setUp() {
		
		a1 = new EdgeAttribute( new double[] {1.0});
		a1_same = new EdgeAttribute(new double[] {1.0});
		a2 = new EdgeAttribute(new double[] {4.0});
		
		b1 = new EdgeAttribute(new double[] {1.0, 2.0, 3.0});
		b1_same = new EdgeAttribute(new double[] {1.0, 2.0, 3.0});
		b2 = new EdgeAttribute(new double[] {2.0, 2.0, 2.0});
		
		c1 = new EdgeAttribute(new double[] {0.945,0,-7.129463692});
	}
	
	@After
	public void tearDown() {
		a1 = null;
		a1_same = null;
		a2 = null;
		
		b1 = null;
		b1_same = null;
		b2 = null;
		
		c1 = null;
	}
	
	@Test
	public void testConstructor() {
		assertEquals("Test 1 (length 1, no vector) failed", a1, new EdgeAttribute("1.0"));
		
		assertEquals("Test 2 (length 1, vector) failed", a2, new EdgeAttribute("[4.0]"));
		
		assertEquals("Test 3 (length >1) failed", b1, new EdgeAttribute("[1.0 2.0 3.0]"));

		assertEquals("Test 4 (length >1, random doubles) failed", c1, new EdgeAttribute("[0.945 0.0 -7.129463692]"));

		// TODO:  how to test that code exits on bad string
//		EdgeAttribute e = new EdgeAttribute("hello");
	}
	
	@Test
	public void testToString() {
		assertEquals("Test 1 (null) failed", "", (new EdgeAttribute()).toString() );
		assertEquals("Test 2 (double array, length 1) failed","1", a1.toString());
		
		assertEquals("Test 3 (double array, length > 1) failed","[1 2 3]", b1.toString());
	}
	
	@Test
	public void testNorm() {
		assertEquals("Test 1 (null) failed", 0.0, (new EdgeAttribute()).norm(),0.000001 );
		assertEquals("Test 2 (double array, length 1) failed", 1, a1.norm(),0.000001);
		assertEquals("Test 3 (double array, length > 1) failed", Math.sqrt(14),b1.norm(),0.000001);
	}
	
	@Test
	public void testEquals() {
		assertEquals("Test 1 (itself) failed", a1, a1);
		assertEquals("Test 2a (equal attributes - double array, length 1) failed", a1, a1_same );
		assertEquals("Test 2b (equal attributes, opposite order - double array, length 1) failed", a1_same,a1);

		assertFalse("Test 3 (unequal attributes, double array, length 1) failed", a1.equals(a2 ) );
		
		assertEquals("Test 4a (equal attributes - double array, length > 1) failed", b1, b1_same);
		assertEquals("Test 4b (equal attributes, opposite order - double array, length >1) failed", b1_same, b1);
		assertFalse("Test 5 (unequals attributes - double array, length >1) failed", b1.equals(b2));	
	
		assertEquals("Test 6 (equal attributes, one is a function) failed; ", new EdgeAttribute(new double[] {0.8}), new EdgeAttribute(new double[] {2.4/3}));
	}
	
	@Test
	public void testDifference() {	
		
		assertEquals("Test 1 (same attribute - double array, length 1) failed",new EdgeAttribute("0.0"),EdgeAttribute.difference(a1, a1) );
		assertEquals("Test 2 (equal attributes - double array, length 1) failed", new EdgeAttribute("0.0"), EdgeAttribute.difference(a1,a1_same));
		assertEquals("Test 3a (unequal attributes - double array, length 1) failed", new EdgeAttribute("-3.0"), EdgeAttribute.difference(a1, a2));
		assertEquals("Test 3b (unequal attributes, opposite order - double array, length 1) failed", new EdgeAttribute("3.0"), EdgeAttribute.difference(a2, a1));

		EdgeAttribute zero_l3 = new EdgeAttribute(new double[] {0.0,0.0,0.0});
		
		EdgeAttribute diff_l3 = new EdgeAttribute(new double[] {-1.0, 0.0, 1.0});
		
		EdgeAttribute diff_opp_l3 = new EdgeAttribute(new double[] {1.0, 0.0, -1.0} );
		
		assertEquals("Test 4 (equal attributes - double array, length >1) failed", zero_l3, EdgeAttribute.difference(b1,b1_same));
		assertEquals("Test 3a (unequal attributes - double array, length >1) failed", diff_l3, EdgeAttribute.difference(b1, b2));
		assertEquals("Test 3b (unequal attributes, opposite order - double array, length 1) failed", diff_opp_l3, EdgeAttribute.difference(b2, b1));

	}
	
	@Test
	public void testSum() {	
		
		assertEquals("Test 1 (same attribute - double array, length 1) failed",new EdgeAttribute(new double[] {2.0}) ,EdgeAttribute.add(a1, a1) );
		assertEquals("Test 2 (unequal attributes - double array, length 1) failed", new EdgeAttribute(new double[] {5.0}), EdgeAttribute.add(a1, a2));
		
		assertEquals("Test 3 (double array, length >1) failed", new EdgeAttribute(new double[] {2.0,4.0,6.0}), EdgeAttribute.add(b1,b1_same));
		assertEquals("Test 3b (different double array, length >1) failed", (new EdgeAttribute(new double[] {2.945,2.0,-5.129463692})), EdgeAttribute.add(b2, c1));
	}
	
	@Test
	public void testProduct() {
		
		assertEquals("Test 1 (same attribute - double array, length 1) failed",new EdgeAttribute(new double[] {1.0}),EdgeAttribute.product(a1, a1_same));
		assertEquals("Test 2 (unequal attribute - double array, length 1) failed",new EdgeAttribute(new double[] {4.0}),EdgeAttribute.product(a1, a2));
	
		assertEquals("Test 3 (same attribute - double array, length >1) failed", new EdgeAttribute(new double[] {1.0,4.0,9.0}), EdgeAttribute.product(b1,b1_same));
		assertEquals("Test 4 (unequal attribute - double array, length >1) failed", new EdgeAttribute(new double[] {2.0,4.0,6.0}), EdgeAttribute.product(b1,b2));	
	}
	
	@Test
	public void testSumOfVector() {
		assertEquals("Test 1 (a1 singular vector) failed", 1, a1.sumOfAttributeVector(),0);
		assertEquals("Test 2 (a2 singular vector) failed", 4, a2.sumOfAttributeVector(),0);
		assertEquals("Test 3 (b1 vectored) failed", 6, b1.sumOfAttributeVector(),0);
		assertEquals("Test 4 (b2 vectored) failed", 6, b2.sumOfAttributeVector(),0);	
	}
	
	
	@Test
	public void testWeightedPairAverage() {
		
		assertEquals("Test 1 (same attribute - double array, length 1) failed; ", a1, EdgeAttribute.weightedPairAverage(a1,a1,0.4));
		assertEquals("Test 2 (same attribute, different instances - double array, length 1) failed; ", a1, EdgeAttribute.weightedPairAverage(a1_same,a1,0.9));
		assertEquals("Test 3a (different attributes - double array, length 1) failed; ", a1, EdgeAttribute.weightedPairAverage(a1,a2,0.0));
		assertEquals("Test 3b (different attributes - double array, length 1) failed; ", new EdgeAttribute(new double[]{2.2}), EdgeAttribute.weightedPairAverage(a1,a2,0.4));
		assertEquals("Test 3c (different attributes - double array, length 1) failed; ", new EdgeAttribute(new double[]{2.5}), EdgeAttribute.weightedPairAverage(a1,a2,0.5));
		assertEquals("Test 3d (different attributes - double array, length 1) failed; ", a2, EdgeAttribute.weightedPairAverage(a1,a2,1));
		assertEquals("Test 3e (different attributes - double array, length 1) failed; ", new EdgeAttribute(new double[]{4.0}), EdgeAttribute.weightedPairAverage(a2,a1,0.0));
		assertEquals("Test 3f (different attributes - double array, length 1) failed; ", new EdgeAttribute(new double[]{1.9}), EdgeAttribute.weightedPairAverage(a2,a1,0.7));
		assertEquals("Test 3g (different attributes - double array, length 1) failed; ", new EdgeAttribute(new double[]{2.5}), EdgeAttribute.weightedPairAverage(a2,a1,0.5));
		assertEquals("Test 3h (different attributes - double array, length 1) failed; ", a1, EdgeAttribute.weightedPairAverage(a2,a1,1));

		assertEquals("Test 4 (same attribute - double array, length >1) failed; ", b1, EdgeAttribute.weightedPairAverage(b1,b1,0.99));
		assertEquals("Test 5 (same attribute, different instances - double array, length >1) failed; ", b1, EdgeAttribute.weightedPairAverage(b1,b1_same,0.3));
		assertEquals("Test 6a (different attributes - double array, length >1) failed; ", b1, EdgeAttribute.weightedPairAverage(b1,c1,0.0));
		assertEquals("Test 6b (different attributes - double array, length >1) failed; ", new EdgeAttribute(new double[]{0.978,1.2,-1.0517854768}), EdgeAttribute.weightedPairAverage(b1,c1,0.4));
		assertEquals("Test 6c (different attributes - double array, length >1) failed; ", new EdgeAttribute(new double[]{0.9725,1,-2.064731846}), EdgeAttribute.weightedPairAverage(b1,c1,0.5));
		assertEquals("Test 6d (different attributes - double array, length >1) failed; ", c1, EdgeAttribute.weightedPairAverage(b1,c1,1));
		assertEquals("Test 6e (different attributes - double array, length >1) failed; ", new EdgeAttribute(new double[]{0.945,0.0,-7.129463692}), EdgeAttribute.weightedPairAverage(c1,b1,0.0));
		assertEquals("Test 6f (different attributes - double array, length >1) failed; ", new EdgeAttribute(new double[]{0.9835,1.4,-0.0388391076}), EdgeAttribute.weightedPairAverage(c1,b1,0.7));
		assertEquals("Test 6g (different attributes - double array, length >1) failed; ", new EdgeAttribute(new double[]{0.9725,1,-2.064731846}), EdgeAttribute.weightedPairAverage(c1,b1,0.5));
		assertEquals("Test 6h (different attributes - double array, length >1) failed; ", b1, EdgeAttribute.weightedPairAverage(c1,b1,1));
	}
	
	@Test
	public void testGetZeroAttribute() {
		assertEquals("Test 1 (length 1) failed; ", new EdgeAttribute(new double[] {0.0}), EdgeAttribute.zeroAttribute(1));
		assertEquals("Test 2 (length >1) failed; ", new EdgeAttribute(new double[] {0.0,0.0,0.0}), EdgeAttribute.zeroAttribute(3));

	}

	
	@Test
	public void testGet() {
		// Test 1: length 1
		assertEquals("Test 1 failed; ", 1.0, a1.get(0), 0);
		
		// Test 2: length >1, first position
		assertEquals("Test 2 failed; ", 0.945, c1.get(0), 0);
		
		// Test 3: length >1, middle position
		assertEquals("Test 3 failed; ", 0, c1.get(1), 0);
				
		// Test 4: length >1, last position
		assertEquals("Test 4 failed; ", -7.129463692, c1.get(2), 0);
	}
}