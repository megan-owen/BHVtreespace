package pca;

import polyAlg.PolyMain;
import distanceAlg1.PhyloTree;
import distanceAlg1.Geodesic;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class PCATest {

	@Before
	public void setUp() throws Exception {
		// Trees from triangle example in the SturmMean manual
		PhyloTree r1 = new PhyloTree("((a:1,b:1):1,(d:1,e:1):2);",true);
		PhyloTree r2 = new PhyloTree("(((a:1,b:1):2,d:1):1.5,e:1);",true);
		PhyloTree r3 = new PhyloTree("(((a:1,d:1):3,b:1):0.5,e:1);",true);
		
		PhyloTree[] triangleTrees = {r1,r2,r3};
		
		// Trees from weighted quadrilateral in SturmMean manual
		PhyloTree q1 = new PhyloTree("((a:1,b:1):3,c:1,(d:1,e:1):1);",false);
		PhyloTree q2 = new PhyloTree("((a:1,b:1):3,c:1,(d:1,e:1):1);",false);
		PhyloTree q3 = new PhyloTree("((a:1,b:1):2,c:1,(d:1,e:1):3);",false);
		PhyloTree q4 = new PhyloTree("((a:1,b:1):2,c:1,(d:1,e:1):3);",false);
		PhyloTree q5 = new PhyloTree("((a:1,b:1):2,c:1,(d:1,e:1):3);",false);
		PhyloTree q6 = new PhyloTree("((a:1,b:1):2,c:1,(d:1,e:1):3);",false);
		PhyloTree q7 = new PhyloTree("((a:1,c:1):1,b:1,(d:1,e:1):0.5);",false);
		PhyloTree q8 = new PhyloTree("((b:1,c:1):2,a:1,(d:1,e:1):4);",false);
		PhyloTree q9 = new PhyloTree("((b:1,c:1):2,a:1,(d:1,e:1):4);",false);
		
		PhyloTree[] quadTrees = {q1, q2, q3, q4, q5, q6,q7, q8, q9};
		
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testScoreProjsToGeodesic() {
		// Test 1:  Trees are from Example 2 in SturmMean manual.
		// Geodesic endpoints are the two trees in the same orthant.
		PhyloTree e1 = new PhyloTree("((a:1,b:1):3,c:1,(d:1,e:1):1);", false);
		PhyloTree e2 = new PhyloTree("((a:1,b:1):2,c:1,(d:1,e:1):3);", false);
		
		PhyloTree[] trees = new PhyloTree[2];
		trees[0] = new PhyloTree("((b:1,c:1):2,a:1,(d:1,e:1):4);", false);
		trees[1] = new PhyloTree("((a:1,c:1):1,b:1,(d:1,e:1):0.5)", false);
		
		Geodesic geo = PolyMain.getGeodesic(e1,e2,null);
		
		// trees[0] projects to e2, contributing a distance of sqrt(17)
		// trees[1] projects to ((a:1,b:1):12/5,c:1,(d:1,e:1):11/5), 
		//          contributing a distance of sqrt((17/5)^2 + (17/10)^2)
		assertEquals("Test 1 failed.", 31.45, PCA.scoreProjsToGeodesic(geo,trees,0.01),0.00001);
		
	}
	
	@Test
	public void testProjToAllGeodesics() {
		// Test 1:  Trees from triangle example (example 1) in SturmMean manual.
//		PCA.projToAllGeodesics(triangleTrees,triangleTrees)
		
	}

}
