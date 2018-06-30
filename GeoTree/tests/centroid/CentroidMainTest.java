package centroid;

import static org.junit.Assert.*;

import org.junit.jupiter.api.*;

import distanceAlg1.*;
import distances.Analysis;
import static polyAlg.PolyMain.calcGeoDist;


public class CentroidMainTest {
	
	private static PhyloTree v1;
	private static PhyloTree v2;
	private static PhyloTree v3;
	private static PhyloTree[] triangleTrees;
	private static PhyloTree triangleMean;
	

	private static PhyloTree w1;
	private static PhyloTree w2;
	private static PhyloTree w3;
	private static PhyloTree w4;
	private static PhyloTree[] squareTrees;
	private static PhyloTree squareMean;
	
	@BeforeEach
	public void setUp() throws Exception {
	
		// example of triangle in three quadrants, where mean is Euclidean mean
		// This is equivalent to the triangle in R^2 with vertices (1,2), (-3, -0.5), and (2, -1.5)
		v1 = new PhyloTree("((a:1,b:1):1,c:1,(d:1,e:1):2);",false);  // (1,2)
		v2 = new PhyloTree("((a:1,b:1):2,d:1,(c:1,e:1):1.5);",false);  //(2, -1.5)
		v3 = new PhyloTree("((a:1,d:1):3,b:1,(c:1,e:1):0.5);",false);  // (-3, -0.5)
		triangleTrees = new PhyloTree[] {v1, v2, v3};
		triangleMean = new PhyloTree("(a:1,b:1,c:1,d:1,e:1);",false);
		
		//example of a (weighted) quadrilateral in three orthants that meet at a common axis
		// Equivalent to the weighted quadrilateral 2*(2,3), 1*(1,1), 1*(6,4), 1/2*(5, 0.5).
		// Translated to tree space so the shared axis is a x = 4.
		// Weighting will be done by including the appropriate number of trees in the set passed to the 
		// mean algorithm.
		w1 = new PhyloTree("((a:1,b:1):3,c:1,(d:1,e:1):1);",false);  // (1,1)
		w2 = new PhyloTree("((a:1,b:1):2,c:1,(d:1,e:1):3);",false);	 // (2,3)
		w3 = new PhyloTree("((a:1,c:1):1,b:1,(d:1,e:1):0.5);",false);	// (5,0.5)
		w4 = new PhyloTree("((b:1,c:1):2,a:1,(d:1,e:1):4);",false);		// (6,4)
		squareTrees = new PhyloTree[] {w1, w1, w2, w2, w2, w2, w3, w4, w4};  // include each tree twice the amount of its weighting
		squareMean = new PhyloTree("((a:1,b:1):1,c:1,(d:1,e:1):2.5);",false);	// (3,2.5)
		
	}

	@AfterEach
	public void tearDown() throws Exception {
		
		triangleTrees = null;
		v1 = null;
		v2 = null;
		v3 = null;
		triangleMean = null;
		
		squareTrees = null;
		w1 = null;
		w2 = null;
		w3 = null;
		w4 = null;
		squareMean = null;
		
		System.out.println();
	}


	
	@Test
	public void testCentroidViaRandPermCauchy() {
		// Test the mean of three vertices of a triangle is the origin. 
		// (Exact mean can be computed by Euclidean methods).
		PhyloTree mean = CentroidMain.getCentroidViaRandPermCauchy(triangleTrees,5000000,5,CentroidMain.getEpsilon(triangleTrees,50000),null,0,0);
		assertEquals("Test 1 failed.  Calculated mean was " + mean, 0, calcGeoDist(mean, triangleMean), 0.001);
	
		// Test the mean of a weighted quadrilateral in three quadrants sharing an axis.
		// Exact mean can be computed by Euclidean methods, by putting the two flaps that don't contain the mean into the same space.
		mean = CentroidMain.getCentroidViaRandPermCauchy(squareTrees,5000000,5,CentroidMain.getEpsilon(squareTrees,50000),null,0,0);
		assertEquals("Test 2 failed.  Calculated mean was " + mean, 0, calcGeoDist(mean, squareMean), 0.001);
	}
	
	@Test
	public void testCentroidViaRandomCauchy() {
		// Test the mean of three vertices of a triangle is the origin. 
		// (Exact mean can be computed by Euclidean methods).
//		PhyloTree mean = CentroidMain.getCentroidViaRandomCauchy(triangleTrees,5000000,25,CentroidMain.getEpsilon(triangleTrees,1000000),null,5000);
		PhyloTree mean = CentroidMain.getCentroidViaRandomCauchy(triangleTrees,10000000,25,0.000001,null,0,0);
		System.out.println("Average distance to origin: " + Analysis.avgDistanceFromOrigin(triangleTrees));
		assertEquals("Test 1 failed.  Calculated mean tree is " + mean, 0, calcGeoDist(triangleMean, mean), 0.005); 
		
		// Test the mean of a weighted quadrilateral in three quadrants sharing an axis.
		// Exact mean can be computed by Euclidean methods, by putting the two flaps that don't contain the mean into the same space.
		mean = CentroidMain.getCentroidViaRandomCauchy(squareTrees,10000000,25,0.000001,null,0,0);
		assertEquals("Test 2 failed,", 0, calcGeoDist(mean, squareMean), 0.005);
	}
	
	@Test
	public void testCentroidViaRandPermTwoRuns() {
		// Test the mean of three vertices of a triangle is the origin. 
		// (Exact mean can be computed by Euclidean methods).
		PhyloTree mean = CentroidMain.getCentroidViaRandPermTwoRuns(triangleTrees,10000000,25,0.000001,null,0,0);
		assertEquals("Test 1 failed.  Calculated mean was " + mean,0, calcGeoDist(triangleMean, mean) ,0.005);
	
		// Test the mean of a weighted quadrilateral in three quadrants sharing an axis.
		// Exact mean can be computed by Euclidean methods, by putting the two flaps that don't contain the mean into the same space.
		mean = CentroidMain.getCentroidViaRandPermTwoRuns(squareTrees,10000000,25,0.000001,null,0,0);
		assertEquals("Test 2 failed.  Calculated mean was " + mean, 0, calcGeoDist(squareMean, mean), 0.005);
	}
	
	@Test
	public void testCentroidViaRandomTwoRuns() {
		// Test the mean of three vertices of a triangle is the origin. 
		// (Exact mean can be computed by Euclidean methods).
		PhyloTree mean = CentroidMain.getCentroidViaRandomTwoRuns(triangleTrees,80000000,20,0.0001,null,0,0);
		assertEquals("Test 1 failed.  Calculated mean was " +mean, 0, calcGeoDist(triangleMean, mean),0.005);
	
		// Test the mean of a weighted quadrilateral in three quadrants sharing an axis.
		// Exact mean can be computed by Euclidean methods, by putting the two flaps that don't contain the mean into the same space.
		mean = CentroidMain.getCentroidViaRandomTwoRuns(squareTrees,80000000,20,0.0001,null,0,0);
		assertEquals("Test 2 failed.  Calculated mean was " + mean, 0, calcGeoDist(squareMean, mean), 0.005);
	}
}