package distances;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import distanceAlg1.PhyloTree;

import java.io.*;

public class PresentationTest {
	private static StringWriter buffer;
	private static PrintWriter outStream;
	
	private static PhyloTree t1;
	private static PhyloTree t2;
	
	private static PhyloTree[] trees;
	private static PhyloTree[] treesWithNull;

	@Before
	public void setUp() throws Exception {
		buffer = new StringWriter();    // will wrap in a PrintWriter
		outStream = new PrintWriter(buffer);  // wraps the StringWriter
		
		t1 = new PhyloTree("(((C:1,(A:1,B:1):1):2,((D:1,E:0.5):3,F:1):1);", true);  // rooted
		t2 = new PhyloTree("(((B:1,((A:1,F:1):1,D:1):2):1,(C:1,E:1):1);", true);  // rooted
		
		trees = new PhyloTree[2];
		trees[0] = t1;
		trees[1] = t2;
		
		treesWithNull = new PhyloTree[3];
		treesWithNull[0] = t1;
		treesWithNull[2] = t2;
	}

	@After
	public void tearDown() throws Exception {
		buffer = null;
		outStream = null;
		
		t1 = null;
		t2 = null;
		trees = null;
		treesWithNull = null;
	}

	@Test
	public void testPrintTreesToStream_notNullTrees() {
		Presentation.printTreesToStream(trees,outStream,true);
		assertEquals("Test 1 (2 non-null trees) failed;", "" + t1.getNewick(true) + "\n" + t2.getNewick(true) + "\n", buffer.toString() );		
	}
	
	@Test
	public void testPrintTreesToStream_nullTree() {
		Presentation.printTreesToStream(treesWithNull,outStream,true);
		assertEquals("Test 1 (one tree null) failed;", "" + t1.getNewick(true) + "\nnull\n" + t2.getNewick(true) + "\n", buffer.toString() );		
	}

}
