package distances;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import distanceAlg1.PhyloTree;

public class Presentation {

	
	/** Prints the given trees to a file, one per line.
	 * 
	 * @param trees
	 * @param file
	 * @throws Exception 
	 */
	public static void printTreesToFile(PhyloTree[] trees, String file) throws Exception {
		int numTrees = trees.length;
		
		PrintWriter outputStream = null;
		  
        try {
        	outputStream = new PrintWriter(new FileWriter(file));
			
        	for (int i = 0; i < numTrees; i++) {
        		if (trees[i] != null) {
        			outputStream.println(trees[i].getNewick(true));
        		}
        		else {
        			outputStream.println("null");
        		}
        	}
	
			if (outputStream != null) {
				outputStream.close();
			}
        } catch (FileNotFoundException e) {
        	System.out.println("Error opening or writing to " + file + ": "+ e.getMessage());
        	throw new Exception();
        } catch (IOException e) {
        	System.out.println("Error opening or writing to " + file + ": " + e.getMessage());
        	throw new Exception();
        }
	}
	
	/** Prints the trees, one per line in Newick format, to outStream.
	 *  Pass in true to have trees written with branch lengths.
	 * 
	 * @param trees
	 * @param outStream
	 */
	public static void printTreesToStream(PhyloTree[] trees, PrintWriter outStream, Boolean branchLengths) {
		for (int i = 0; i < trees.length; i++) {
    		if (trees[i] != null) {
    			outStream.println(trees[i].getNewick(branchLengths));
    		}
    		else {
    			outStream.println("null");
    		}
    	}
	}
	
	/** Prints the string s to the file outfile.
	 * 
	 * @param s
	 * @param outfile
	 * @throws Exception 
	 */
	public static void printStringToFile(String s, String outfile) throws Exception {
		PrintWriter outputStream = null;
		  
        try {
        	outputStream = new PrintWriter(new FileWriter(outfile));
			
        	outputStream.println(s);
	
			if (outputStream != null) {
				outputStream.close();
			}
        } catch (FileNotFoundException e) {
        	System.out.println("Error opening or writing to " + outfile + ": "+ e.getMessage());
        	throw new Exception();
        } catch (IOException e) {
        	System.out.println("Error opening or writing to " + outfile + ": " + e.getMessage());
        	throw new Exception();
        }
		
	}
}
