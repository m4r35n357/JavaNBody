/**
 * 
 */
package uk.me.doitto.orbits;

/**
 * @author ian
 *
 */
public class Particle {
	
	public double mass;
	
	public double qX;
	
	public double qY;
	
	public double qZ;

	public double pX;
	
	public double pY;
	
	public double pZ;

	/**
	 * 
	 */
	public Particle (double qX, double qY, double qZ, double pX ,double pY ,double pZ, double mass) {
		this.qX = qX;
		this.qY = qY;
		this.qZ = qZ;
		this.pX = pX;
		this.pY = pY;
		this.pZ = pZ;
		this.mass = mass;
	}

}
