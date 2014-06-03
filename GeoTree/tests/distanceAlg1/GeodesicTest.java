package distanceAlg1;

import static org.junit.Assert.*;

import java.util.Vector;

import distanceAlg1.PhyloTree;
import polyAlg.PolyMain;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class GeodesicTest {
	private static PhyloTree u1;
	private static PhyloTree u2;
	private static PhyloTree u3;

	private static PhyloTree t1;
	private static PhyloTree t1_different_edge_lengths;
	private static PhyloTree t2;
	private static PhyloTree multi;
	
	private static PhyloTree v1;
	private static PhyloTree v2;
	
	private static PhyloTree s1;
	private static PhyloTree s2;
	private static PhyloTree s_multi;
	
	private static PhyloTree t4_1;
	private static PhyloTree t4_2;
	

	@Before
	public void setUp() throws Exception {
		u1 = new PhyloTree("((C:1,(A:1,B:1):1):2,((D:1,E:0.5):3,F:1):1);", true);  // rooted; splits AB|CDEF0, ABC|DEF0, DE|ABCF0, DEF|ABC0
		u2 = new PhyloTree("((B:1,((A:1,F:1):1,D:1):2):1,(C:1,E:1):1);", true);  // rooted; splits AF|BCDE0, ADF|BCE0, CE|ABDF0, ABDF|CE0 
		u3 = new PhyloTree("((((C:1,B:1):1,A:1):1,E:1):1,(D:1,F:1):1);", true);  // rooted; split ABC in common with t1
	
		t1 = new PhyloTree("(a:1,(b:1,c:1):0.3,(e:1,(d:10,f:0.22):1):1);",false);
		t1_different_edge_lengths = new PhyloTree("((b:1.04,c:1):0.39,a:5,((f:10,d:0.22):1,e:1):1);",false);
		t2 = new PhyloTree("((a:1,b:2.44):1,(c:1,d:1):1,(e:1,f:1):1);",false);
		
		// unrooted trees.  Caused problems in getTreeAt
		v1 = new PhyloTree("((a:1,b:1):1,c:1,(d:1,e:1):2);",false); 	// unrooted; splits: ab|cde, de|abc   (1 in common with s2)
		v2 = new PhyloTree("((a:1,b:1):2,d:1,(c:1,e:1):1.5);",false);	// unrooted; splits: ab|cde, ce|abd
		
		multi = new PhyloTree("((a:1,b:1,c:1):3,f:1,(d:1,e:1):1);",false);
		
		// same topology as t1
		s1 = new PhyloTree("(a:[1 1 1],(b:[1 1 1 ],c:[1 1 1]):[0.3 -3.7778 -4],(e:[1 1 1],(d:[10 -0.6 2],f:[0.22 0 0]):[1 1 1]):[2 2 2]);",false); 
		// same topology as t2
		s2 = new PhyloTree("((a:[1 1 1],b:[2.44 2.44 2.44]):[1 -1 1],(c:[1 1 1],d:[1 0 -1]):[1 1 1],(e:[1 1 1],f:[1 1 1]):[2 2 2]);",false);
		// same topology as multi
		s_multi = new PhyloTree("((a:[1 1 1],b:[1 1 1],c:[1 1 1]):[3 3 3],f:[1 1 1],(d:[1 1 1],e:[1 1 1]):[1 1 1]);",false);
		
		t4_1 = new PhyloTree("((a:1,b:1):2,(c:1,d:1):1);",true);  // rooted; splits ab|cd0, cd|ab0
		t4_2 = new PhyloTree("(((a:1,c:1):1,b:1):2,d:1);",true);	// rooted; splits ac|bd0, abc|d0
	
	}

	@After
	public void tearDown() throws Exception {
		t1 = null;
		t2 = null;
		
		u1 = null;
		u2 = null;
		u3 = null;
		
		v1 = null;
		v2 = null;
		
		t1_different_edge_lengths = null;
		multi = null;
		
		s1 = null;
		s2 = null;
		s_multi = null;
	}

	@Test
	public void testNumTopologies() {
		assertEquals("Test 1 (no common edges, both bifurcating, cone path) failed;", 2,PolyMain.getGeodesic(u1,u2,null).numTopologies());
		assertEquals("Test 2 (common edges, both bifurcating, not cone path) failed;", 3,PolyMain.getGeodesic(u1,u3,null).numTopologies());
	}

	@Test
	public void testGetCommonEdges() {
		// same tree, different edges, edge attrib vector length 1
		
		Geodesic geo = PolyMain.getGeodesic(t1, t1_different_edge_lengths,null);
		assertEquals("Test 1a failed; ",t1.getEdges(),geo.getCommonEdges(0));
		
		Vector<PhyloTreeEdge> answer = TreeDistance.myVectorClonePhyloTreeEdge(t1.getEdges());
		answer.get(0).setAttribute(new EdgeAttribute("0.345"));
		
		assertEquals("Test 1b failed; ",answer,geo.getCommonEdges(0.5));
		assertEquals("Test 1c failed; ",t1_different_edge_lengths.getEdges(),geo.getCommonEdges(1));
	
		// different trees, no common edges, length 1
		geo = PolyMain.getGeodesic(t1,t2,null);
		assertEquals("Test 2 failed;", new Vector<PhyloTreeEdge>(), geo.getCommonEdges(0.7));
	
		// different trees (one multi), 1 common edge (different lengths), 1 edge compatible with all edges in the other tree, length 1
		answer = new Vector<PhyloTreeEdge>();
		answer.add(new PhyloTreeEdge(new Bipartition("011000"), new EdgeAttribute("0.3"),-1));
		answer.add(new PhyloTreeEdge(new Bipartition("111000"),new EdgeAttribute("1"),-1));
		
		geo = PolyMain.getGeodesic(t1,multi,null);
		assertEquals("Test 3a failed;", answer,geo.getCommonEdges(0));
		
		answer.get(1).setAttribute(new EdgeAttribute("1.6"));
		answer.get(0).setAttribute(new EdgeAttribute("0.21"));
		assertEquals("Test 3b failed;", answer,geo.getCommonEdges(0.3));
		
		answer.get(1).setAttribute(new EdgeAttribute("2.2"));
		answer.get(0).setAttribute(new EdgeAttribute("0.12"));
		assertEquals("Test 3c failed;", answer,geo.getCommonEdges(0.6));
		
		answer.get(1).setAttribute(new EdgeAttribute("3"));
		answer.get(0).setAttribute(new EdgeAttribute("0"));
		assertEquals("Test 3d failed;", answer,geo.getCommonEdges(1));
		
		// same trees, length >1
		geo = PolyMain.getGeodesic(s1,s1,null);
		assertEquals("Test 4 failed;", s1.getEdges(), geo.getCommonEdges(0.97));
		
		// different trees, no common edges, length >1
		geo = PolyMain.getGeodesic(s1,s2,null);
		assertEquals("Test 5 failed;", new Vector<PhyloTreeEdge>(), geo.getCommonEdges(0.7));
	
		// different trees (one multi), 1 edge compatible with edges in the other tree, length >1
		answer = new Vector<PhyloTreeEdge>();
		answer.add(new PhyloTreeEdge(new Bipartition("110000"), new EdgeAttribute("[0 0 0]"),-1));
		
		geo = PolyMain.getGeodesic(s_multi,s2,null);
		assertEquals("Test 6a failed;", answer,geo.getCommonEdges(0));
		
		answer.get(0).setAttribute(new EdgeAttribute("[0.3 -0.3 0.3]"));
		assertEquals("Test 6b failed;", answer,geo.getCommonEdges(0.3));
		
		answer.get(0).setAttribute(new EdgeAttribute("[0.6 -0.6 0.6]"));
		assertEquals("Test 6c failed;", answer,geo.getCommonEdges(0.6));
		
		answer.get(0).setAttribute(new EdgeAttribute("[1 -1 1]"));
		assertEquals("Test 6d failed;", answer,geo.getCommonEdges(1));

		
	}
	
	@Test
	public void testGetBoundaryTrees() {
		Vector<PhyloTree> answer = new Vector<PhyloTree>();
		
		assertEquals("Test 1 (same orthants) failed; ", answer, Geodesic.getBoundaryTrees(t1,t1_different_edge_lengths));
		
		answer.add(new PhyloTree("((a:1,b:1):1,c:1,d:1);",true));
		answer.add(new PhyloTree("((a:1,b:1,c:1):1,d:1);",true));
		
		assertEquals("Test 2 (different orthants) failed;", answer, Geodesic.getBoundaryTrees(t4_1,t4_2));
	}
	
    @Test
    public void testGetTreeAt() {
    	// Test 1:  end points of the geodesic are the same (attribute vector length is 1), rooted
    	Geodesic geo = PolyMain.getGeodesic(u1,u1, null);
    	assertEquals("Test 1 failed; ", u1, geo.getTreeAt(0.73,u1.getLeaf2NumMap(), u1.isRooted()) );
    	
    	// Test 2: different geodesic endpoints, position is 0 (attribute vector length is 1), rooted
    	geo = PolyMain.getGeodesic(u1, u2, null);
    	assertEquals("Test 2 failed; ", u1, geo.getTreeAt(0, u1.getLeaf2NumMap(), u1.isRooted() ) );
    	
    	// Test 3: different geodesic endpoints, position is 1 (attribute vector length is 1), rooted
    	// same geo as Test 2
    	assertEquals("Test 3 failed; ", u2, geo.getTreeAt(1, u1.getLeaf2NumMap(), u1.isRooted() ) );
    	
    	// Test 4: different geodesic endpoints, position is random (attribute vector length is 1), rooted
    	// geo is cone path
    	PhyloTree answer = new PhyloTree("((C:1,(A:1,B:1):0.49506098467):0.99012196935,((D:1,E:0.65):1.48518295403,F:1):0.49506098467);", true);
    	assertEquals("Test 4 failed; ", answer, geo.getTreeAt(0.3, u1.getLeaf2NumMap(), u1.isRooted() ) );
    	
    	//Test 5: different geodesic endpoints, position is the cone point (attribute vector length is 1), rooted
    	answer = new PhyloTree("(A:1,B:1,C:1,D:1,E:0.79706557712,F:1);",true);
    	assertEquals("Test 5 failed; ", answer, geo.getTreeAt(0.59413115425,u1.getLeaf2NumMap(),u1.isRooted()) );
    	
    	// Test 6: unrooted trees, with one split in common (caused problems in code)
    	System.out.println("************* Start of Test 6 ***********");
    	answer = new PhyloTree("((a:1,b:1):1.5,c:1,(d:1,e:1):0.25);",false);
    	geo = PolyMain.getGeodesic(v1,v2,null);
    	assertEquals("Test 6 failed; ", answer, geo.getTreeAt(0.5, v1.getLeaf2NumMap(), v1.isRooted() ) );
    }

    @Test
    public void testEquals() {
    	// Test 1: basic test
    	assertEquals("Test 1 failed; ", PolyMain.getGeodesic(t1,t2, null), PolyMain.getGeodesic(t1,t2,null));
    }

}