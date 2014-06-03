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

}
