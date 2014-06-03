package polyAlg;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class ToolsTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testNormalize() {
		// Test 1 - all 0 coords
		double [] coords1 = {0.0, 0.0, 0.0};
		for (int i = 0; i < coords1.length; i++) {
			assertEquals("Test 1-" + i + " failed; ", coords1[i], Tools.normalize(coords1)[i],0);
		}
		
		// Test 2 - all positive coords
		double [] coords2 = {3,4};
		double [] coords2_answer = {3.0/5, 4.0/5};
		for (int i = 0; i < coords2.length; i++) {
			assertEquals("Test 2-" + i + " failed; ", coords2_answer[i], Tools.normalize(coords2)[i],0.0000000000001);
		}
		
		// Test 3 - mixed coords
		double [] coords3 = {1, -1, 0};
		double [] coords3_answer = {0.707106781186548,-0.707106781186548,0};
		for (int i = 0; i < coords3.length; i++) {
			assertEquals("Test 3-" + i + " failed; ", coords3_answer[i], Tools.normalize(coords3)[i],0.0000000001);
		}
	}

	@Test
	public void testScaleBy() {
		// Test 1 - all 0 coords
		double [] coords1 = {0.0, 0.0, 0.0};
		for (int i = 0; i < coords1.length; i++) {
			assertEquals("Test 1-" + i + " failed; ", coords1[i], Tools.scaleBy(coords1,5)[i],0);
		}
		
		// Test 2 - mixed coords, positive scale
		double [] coords2 = {1.1, -1, 0};
		double [] coords2_answer = {7.7,-7,0};
		for (int i = 0; i < coords2.length; i++) {
			assertEquals("Test 2-" + i + " failed; ", coords2_answer[i], Tools.scaleBy(coords2,7)[i],0.0000000001);
		}
		
		// Test 3 - mixed coords, negative scale
		double [] coords3 = {1.1, -1, 0};
		double [] coords3_answer = {-3.3,3,0};
		for (int i = 0; i < coords3.length; i++) {
			assertEquals("Test 3-" + i + " failed; ", coords3_answer[i], Tools.scaleBy(coords3,-3)[i],0.0000000001);
		}
		
		// Test 4 - mixed coords, 0 scale
		double [] coords4 = {1.1, -1, 0};
		double [] coords4_answer = {0,0,0};
		for (int i = 0; i < coords4.length; i++) {
			assertEquals("Test 4-" + i + " failed; ", coords4_answer[i], Tools.scaleBy(coords4,0)[i],0.0000000001);
		}
		
	}

	@Test
	public void testRoundSigDigits() {
		// Test 1: no decimals, 1 sig. digit, round down
		assertEquals("Test 1 failed; ", 1000,Tools.roundSigDigits(1222,1),0.000000001);
		
		// Test 2: decimals, > 1 sig. digit, round up
		assertEquals("Test 2 failed; ", 0.0783,Tools.roundSigDigits(0.07829976,3),0.0000000001);
	}
}
