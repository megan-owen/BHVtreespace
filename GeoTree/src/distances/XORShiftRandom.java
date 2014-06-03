package distances;

import java.util.Random;

/** Code from Javamex:  http://www.javamex.com/tutorials/random_numbers/java_util_random_subclassing.shtml
 * Implements the XORShift Random Number Generator (RNG) of George Marsaglia. (Journal of Statistical Software, Vol. 8, Jul. 2003.)
 * 
 * Notes: - The implementation above is not thread-safe, whereas the base java.lang.Random is. If we wanted
 *  to fix this, we'd declare seed as an AtomicLong and use an atomic update: see the source code to java.lang.Random
 *   for an example, plus this site's tutorial on the Java atomic classes.
    * Our implementation wouldn't serialise its state (current seed), unlike java.lang.Random.
    * We might want to override nextLong(), since we can generate 64 bits at a time. This method would then not be able
    *  to generate the value zero, of course, but as it stands, there'll always be one value that nextLong() can't generate,
    *   since the generator has a period of 264-1, and there are 264 possible longs. 
 *
 */


public class XORShiftRandom extends Random {
	  private long seed = System.nanoTime();
	  private static final long serialVersionUID = 42L;  // to get rid of warning that serialVersionUID is not declared

	  public XORShiftRandom() {
	  }
	  protected int next(int nbits) {
	    // N.B. Not thread-safe!
	    long x = this.seed;
	    x ^= (x << 21);
	    x ^= (x >>> 35);
	    x ^= (x << 4);
	    this.seed = x;
	    x &= ((1L << nbits) -1);
	    return (int) x;
	  }
	}
